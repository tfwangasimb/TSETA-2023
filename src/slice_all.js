// @ts-check

const os = require("os");
const fs = require("fs");
const child_process = require("child_process");
const Path = require("path");

const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide } = require("./blastn_util.js");
const { Dataset } = require("./dataset.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");

const WRITE_VERBOSE_LOG = process.env.WRITE_VERBOSE_LOG;

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}

const argv_start_chr            = parseSafeIntegerLike(argv["--chr"], 1);
const argv_end_chr              = parseSafeIntegerLike(argv["--chr-end"], Infinity);
const argv_slice_size           = parseSafeIntegerLike(argv["--s"],   10000);
const argv_add_size             = parseSafeIntegerLike(argv["--a"],   10000);
const argv_init_delta           = parseSafeIntegerLike(argv["--d"],   argv_slice_size * 3);// local diff
const argv_init_max_delta       = parseSafeIntegerLike(argv["--x"],   argv_init_delta * 1);// telomere diff
const argv_minimum_align_length = parseSafeIntegerLike(argv["--m"],   1000);

/**
 * @param {any} value 
 */
function parseSafeIntegerLike(value, default_value) {
	let num = Number(value);
	if (Number.isFinite(num) && !Number.isNaN(num) && Number.isSafeInteger(num)) {
		return num;
	}
	else {
		return default_value;
	}
}

const VERBOSE = !!argv["--verbose"];

if (VERBOSE) {
	console.log({
		argv_slice_size, argv_add_size, argv_init_delta, argv_init_max_delta, argv_minimum_align_length,
	});
}

const genome_info_list = dataset.loadGenomeInfoList();

if (!fs.existsSync(`${dataset.tmp_path}/ma_util_blastn`)) {
	fs.mkdirSync(`${dataset.tmp_path}/ma_util_blastn`);
}

if (fs.realpathSync(process.argv[1]) == __filename) {
	main();
}

async function main() {
	program_log(`${dataset.name}.log.txt`, "start");

	if (!fs.existsSync(`${dataset.tmp_path}/seq_frag`)) {
		fs.mkdirSync(`${dataset.tmp_path}/seq_frag`);
	}
	if (WRITE_VERBOSE_LOG) {
		if (!fs.existsSync(`${dataset.tmp_path}/meta`)) {
			fs.mkdirSync(`${dataset.tmp_path}/meta`);
		}
	}

	const start_chr = Math.max(1, argv_start_chr);
	const end_chr = Math.min(genome_info_list[0].chr_list.length, argv_end_chr);

	if (VERBOSE) {
		console.log({
			start_chr,
			end_chr,
		});
	}

	for (let nChr = start_chr; nChr <= end_chr; ++nChr) {
		if (VERBOSE) {
			console.log("start ch", nChr);
		}
		try {
			fs.writeFileSync(`${dataset.tmp_path}/multi_coord_ch${nChr}.txt`, "");//clear
			fs.writeFileSync(`${dataset.tmp_path}/frag_length_ch${nChr}.txt`, "");//clear
			
			const count = await multi_coord(nChr);

			if (VERBOSE) {
				console.log("number of fragments", count);
			}
		}
		catch (ex) {
			fs.writeFileSync(`${dataset.tmp_path}/ma_util.error.txt`, "ch" + nChr + "\n" + ex.stack + "\n", { flag: "a" });
			console.error(ex);
		}
		console.log(`finished ch${nChr}`);
	}

	console.log("next step:", "execute MAFFT");
	console.log("command:", `node ${__dirname}/run_mafft.js -dataset ${argv_dataset_path} --mafft-algorithm localpair --mafft-maxiterate 1000 --mafft-thread 20`);

	program_log(`${dataset.name}.log.txt`, "exit");
}

/**
 * @param {number} nChr
 */
async function multi_coord(nChr) {
	const chr_idx = nChr - 1;
	
	const chr_info_list = genome_info_list.map(genome_info => genome_info.chr_list[chr_idx]);

	const chr_name_list = chr_info_list.map(chr_info => chr_info.chr);
	const chr_fastaPath_list = chr_info_list.map(chr_info => chr_info.path);

	const chr_seq_list = chr_fastaPath_list.map((filepath, i) => readFasta(filepath)[chr_name_list[i]]);
	
	if (chr_info_list.some(a => a.length <= argv_slice_size || a.length <= argv_init_delta)) {
		if (VERBOSE) {
			console.log(`chr size less then ${argv_slice_size}`);
		}

		const fragId = 1;
		const search_start = 1;

		const extract_length_list = chr_info_list.map(a => a.length);

		const search_end = Math.max(...extract_length_list);

		{
			const coord_text_1 = fragId + "\tstart\t" + [search_start, ...chr_info_list.map(_ => 1)].join("\t|\t");
			const coord_text_2 = fragId + "\t  end\t" + [search_end,   ...chr_info_list.map(chrInfo => chrInfo.length)].join("\t|\t");

			if (VERBOSE) {
				console.log(coord_text_1);
				console.log(coord_text_2);
			}
			
			fs.writeFileSync(`${dataset.tmp_path}/multi_coord_ch${nChr}.txt`, coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
		}

		const output_coord_list = chr_info_list.map(chrInfo => [1, chrInfo.length]);
		const extract_seq_list = chr_seq_list;
		{
			const fa_seq_name_list = chr_name_list.map((chrName, i) => {
				const [start, end] = output_coord_list[i];// charAt
				return `${chrName} ${start + 1}-${end} ${extract_seq_list[i].length}`;
			});

			const output_fasta_file_name = `${dataset.tmp_path}/seq_frag/ch${nChr}_${fragId}.fa`;
			
			if (VERBOSE) {
				console.log(output_fasta_file_name);
			}

			saveFasta(output_fasta_file_name, fa_seq_name_list.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_list[i];
				return fasta;
			}, {}));
		}

		{
			const min_len = Math.min(...extract_length_list);
			const max_len = Math.max(...extract_length_list);
			fs.writeFileSync(`${dataset.tmp_path}/frag_length_ch${nChr}.txt`, [
				fragId,
				...extract_length_list,
				...extract_length_list.map(len => (len / max_len).toFixed(2)),
				min_len, max_len,
				(min_len / max_len).toFixed(2),
			].join("\t") + "\n", { flag: "a" });
		}

		return;
	}

	const pos_start_list = dataset.genomeNameList.map(a => 1);
	const pos_end_list = dataset.genomeNameList.map(a => -1);
	const find_next_start_list = dataset.genomeNameList.map(a => 0);

	let max_delta = argv_init_max_delta;

	let search_start = 1;
	let search_end = search_start + argv_slice_size;
	
	pos_start_list[0] = search_start;
	pos_end_list[0] = search_end;

	let fragId = 1;

	const seq_from_ref1 = dataset.genomeNameList.map(gName => gName == dataset.ref);

	let _global_search_align, _local_overlap_align;
	let pos_search_start_list = dataset.genomeNameList.map(a => 1);

	const has_centromere = dataset.centromere && dataset.centromere[nChr];

	/**
	 * after: after slice centromere, once("after", add centeromere tag), remove centro_status
	 * @type {"before"|"in"|"after"}
	 */
	let centro_status = has_centromere ? "before" : null;

	const [centro_start, centro_end] = has_centromere ? dataset.centromere[nChr] : [];

	for (; search_start <= chr_info_list[0].length && search_end <= chr_info_list[0].length; ++fragId) {
		try {
			/**
			 * @param {BlastnCoord} coord_1
			 * @param {BlastnCoord} coord_2
			 */
			function check_q_start_end(coord_1, coord_2) {
				return coord_1.qstart <= coord_2.qend && coord_1.qend >= coord_2.qstart;
			}
			
			/**
			 * @param {BlastnCoord} row
			 * @param {any[]} params
			 */
			function next_start_filter(row, params) {
				return (
					row.qstart < row.qend &&
					row.sstart < row.send
					// //row.sstart >= next_start
					// row.send >= params[0]
				);
			}

			for (let i = 1; i < pos_search_start_list.length; ++i) {
				pos_search_start_list[i] = Math.max(pos_start_list[i], find_next_start_list[i]);
			}
			//skip ref, pos_search_start[0] => 1
			if (pos_search_start_list.some((pos_search_start, i) => pos_search_start >= chr_info_list[i].length)) {
				//translocation
				++fragId;
				break;
			}

			let pos_search_end_list = pos_search_start_list.map((pos_search_start, i) => Math.min(pos_search_start + max_delta, chr_info_list[i].length))
			
			// remove ref1-ref1 repeat and duplicate
			let r1_coords_task = blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[0], search_start, search_end, search_start, search_end, next_start_filter, [], nChr + "_" + fragId);

			let coords_task_list = [];
			for (let i = 1; i < pos_search_start_list.length; ++i) {
				let task = blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[i], search_start, search_end, pos_search_start_list[i], pos_search_end_list[i], next_start_filter, [], nChr + "_" + fragId);
				coords_task_list.push(task);
			}
			let result_coords_list = await Promise.all([r1_coords_task, ...coords_task_list]);
			_global_search_align = result_coords_list;
			
			/**
			 * @param {BlastnCoord[][]} coords_list - [r1[], r2[], s1[], s2[], s3[], s4[]]
			 */
			function check_overlap_ref(coords_list) {
				let _coords_list = coords_list.map(coords => coords.filter(r2c => r2c.align >= argv_minimum_align_length));
				
				if (_coords_list[1] && _coords_list[1].length) {
					let all_match_group = _coords_list[1].map(r2_coord => {
						let _ret = [
							_coords_list[0].sort((a, b) => (Math.abs(a.qstart - r2_coord.qstart) + Math.abs(a.qend - r2_coord.qend)) - (Math.abs(b.qstart - r2_coord.qstart) + Math.abs(b.qend - r2_coord.qend)))[0],
							r2_coord,
							..._coords_list.slice(2).map(_coords => _coords.sort((a, b) => (Math.abs(a.qstart - r2_coord.qstart) + Math.abs(a.qend - r2_coord.qend)) - (Math.abs(b.qstart - r2_coord.qstart) + Math.abs(b.qend - r2_coord.qend)))[0]),
						];
						return _ret;
					}).sort((a, b) => b[1].send - a[1].send);

					return {
						/** returns [r1, r2, s1, s2, s3, s4] */
						best_match: all_match_group[0],
						all_match_groups: all_match_group,
					};
				}
				return {
					best_match: [],
					all_match_groups: [],
				};
			}

			let match_results = check_overlap_ref(result_coords_list);
			_local_overlap_align = match_results.best_match;

			if (match_results.best_match.length && match_results.best_match.every((match, i) => match)) {
				if (match_results.best_match.some((match, i) => (match.sstart - pos_start_list[i]) >= max_delta)) {

					if (VERBOSE) {
						console.log(fragId, "out of range:", "max_delta", max_delta);
					}
					// if (VERBOSE) {
					// 	let oor = {
					// 		r2D: Math.max(0, (r2_coord.sstart - r2_start) - max_delta),
					// 		s1D: Math.max(0, (s1_coord.sstart - s1_start) - max_delta),
					// 		s2D: Math.max(0, (s2_coord.sstart - s2_start) - max_delta),
					// 		s3D: Math.max(0, (s3_coord.sstart - s3_start) - max_delta),
					// 		s4D: Math.max(0, (s4_coord.sstart - s4_start) - max_delta)
					// 	};
					// 	console.table(oor);
					// }

					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
					max_delta += argv_add_size;
					continue;
				}
				for (let i = 1; i < chr_name_list.length; ++i) {
					find_next_start_list[i] = match_results.best_match[i].send + 1;
				}

				//check qend
				let _min_qend = Math.min(...match_results.best_match.map(match => match.qend));

				if (!match_results.best_match.every(match => check_q_start_end(match_results.best_match[1], match)) ||
					match_results.best_match.some(match => match_results.best_match[1].qstart == _min_qend)
				) {
					// if (VERBOSE) {
					// 	console.log(fragId, "out of q start end:", {
					// 		search_start, search_end,
					// 		r2_coord, s1_coord, s2_coord, s3_coord, s4_coord
					// 	});
					// }
					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
					max_delta += argv_add_size;
					if (VERBOSE) {
						console.log("max_delta", max_delta);
					}
					debugger;

					continue;
				}

				let max_iterate = 10;
				while ((--max_iterate) > 0) {
					/**
					 * @param {BlastnCoord} row
					 * @param {any[]} params
					 */
					function q_end_filter(row, params) {
						return (
							row.qstart < row.qend &&
							row.sstart < row.send &&
							//row.sstart >= next_start
							row.qend == params[0]
						);
					}

					const _coords_tasks = match_results.best_match.map((match, i) => {
						if (match.qend != _min_qend) {
							return blastn_coord(chr_fastaPath_list[0], chr_fastaPath_list[i], search_start, _min_qend, pos_search_start_list[i], match.send, q_end_filter, [_min_qend], nChr + "_" + fragId);
						}
						else {
							return null;
						}
					});
					let _coords_list = await Promise.all(_coords_tasks);
					
 					if (_coords_list.some(_coords => _coords != null && _coords.length <= 0)) {
						if (VERBOSE) {
							console.log({ _min_qend });
						}
						--_min_qend;
						continue;
					}

					_coords_list.forEach((_coords, i) => {
						if (_coords) {
							const n = _coords.sort((a, b) => b.send - a.send)[0];
							const p = match_results.best_match[i];
							match_results.best_match[i] = n || p;
						}
					});
					
					match_results.best_match.forEach((coord, i) => {
						seq_from_ref1[i] = coord != null;
					});//last identical loaction

					//check qend
					const qend_list = match_results.best_match.slice(1).map(match => match.qend);
					const min_qend = Math.min(...qend_list);
					if (!qend_list.every(a => min_qend == a)) {
						if (VERBOSE) {
							console.log(fragId, "gap or mis");
						}
						continue;
					}
					
					if (match_results.best_match.every((match, i) => {
						return match.strand <= 0 || (pos_start_list[i] - 1) >= match.send;
					})) {
						if (VERBOSE) {
							console.log(fragId, "inv:", {
								search_start, search_end,
								best_match: match_results.best_match,
								pos_start_list,
							});
						}
						continue;
					}

					const a_list = match_results.best_match.map((coord, i) => chr_seq_list[i][coord.send - 1]);
					const a1_s = [...a_list[1]];
					if (a_list.slice(2).every(a => a1_s.every((v, i) => v == a[i]))) {
						pos_start_list[0] = search_start;
						pos_end_list[0] = min_qend;
						break;
					}
					else {
						if (VERBOSE) {
							console.log({ _min_qend, ...a_list });
						}
						--_min_qend;
					}
				}
				if (max_iterate <= 0) {
					search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);
					max_delta += argv_add_size;
					if (VERBOSE) {
						console.log("max_delta", max_delta);
					}
					continue;
				}
				
				match_results.best_match[0].send = pos_end_list[0];
				const output_coord_list = match_results.best_match.map((match, i) => [pos_start_list[i] - 1, match.send]);// charAt
				const extract_seq_list = output_coord_list.map(([start, end], i) => chr_seq_list[i].slice(start, end));
				const extract_length_list = extract_seq_list.map(seq => seq.length);

				if (centro_status) {
					const { [0]: ref1_start_pos, [1]: ref1_end_pos } = output_coord_list[0];

					const over_centromere = ref1_start_pos <= centro_start && ref1_end_pos >= centro_end;
					if (centro_status == "before" && over_centromere) {
						if (nChr != 6) {
							console.log({
								ref1_start_pos, centro_start, ref1_end_pos, centro_end
							});
							throw new Error(nChr + "over_centromere");
						}
						centro_status = "in";
						console.log("v3", `over_centromere -> centro_status = "in"`);
					}

					if (centro_status == "before") {
						const is_centromere = ref1_end_pos >= centro_start;
						if (is_centromere) {
							console.log("v3", "before centromere");
							console.log({
								ref1_start_pos, ref1_end_pos,
								centro_start, centro_end,
							});
							centro_status = "in";
						}
					}
					else if (centro_status == "in") {
						console.log("v3", "in centromere");
						console.log({
							ref1_start_pos, ref1_end_pos,
							centro_start, centro_end,
						});

						const leave_centromere = ref1_end_pos >= centro_end;
						if (leave_centromere) {
							centro_status = "after";
							console.log("v3", "after centromere");
						}
						else {
							continue;//skip centromere
						}
					}
					// else if (centro_status == "after") {
					// 	// nothing
					// }
				}

				const min_len = Math.min(...extract_length_list);
				const max_len = Math.max(...extract_length_list);

				{
					const coord_text_1 = fragId + "\tstart\t" + [search_start, ...match_results.best_match.map((match, i) => pos_start_list[i])].join("\t|\t") + "\t" + (centro_status == "after" ? "centromere" : "");
					const coord_text_2 = fragId + "\t  end\t" + [search_end,   ...match_results.best_match.map((match, i) => match.send)].join("\t|\t") + "\t" + (centro_status == "after" ? "centromere" : "");
					if (VERBOSE) {
						console.log(coord_text_1);
						console.log(coord_text_2);
					}
					
					fs.writeFileSync(`${dataset.tmp_path}/multi_coord_ch${nChr}.txt`, coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
				}

				fs.writeFileSync(`${dataset.tmp_path}/frag_length_ch${nChr}.txt`, [
					fragId,
					...extract_length_list,
					...extract_length_list.map(len => (len / max_len).toFixed(2)),
					min_len, max_len,
					(min_len / max_len).toFixed(2),
					"centromere_" + centro_status == "after",
					"remove_" + centro_status == "after",
				].join("\t") + "\n", { flag: "a" });

				if (WRITE_VERBOSE_LOG) {
					fs.writeFileSync(`${dataset.tmp_path}/meta/meta_ch${nChr}_${fragId}.json`, JSON.stringify({
						id: fragId,
						best_match: match_results.best_match,
						match_group: match_results.all_match_groups,
						coord: output_coord_list,
						length: extract_length_list,
						centromere: centro_status == "after",
						remove: centro_status == "after",
					}, null, "\t"));
				}
				
				{
					let cc = extract_seq_list.map(a => a.slice(-2));
					if (cc.some(a => cc[0] != a)) {
						if (VERBOSE) {
							console.log("diff align len fragId=", fragId);
						}
						debugger;
					}
				}
				
				{
					const fa_seq_name_list = chr_name_list.map((chrName, i) => {
						const [start, end] = output_coord_list[i];// charAt
						return `${chrName} ${start + 1}-${end} ${extract_seq_list[i].length}`;
					});

					const output_fasta_file_name = `${dataset.tmp_path}/seq_frag/ch${nChr}_${fragId}.fa`;
					
					if (VERBOSE) {
						console.log(output_fasta_file_name);
					}

					saveFasta(output_fasta_file_name, fa_seq_name_list.reduce((fasta, seqName, i) => {
						fasta[seqName] = extract_seq_list[i];
						return fasta;
					}, {}));
				}
				
				if (centro_status == "after") {
					centro_status = null;
				}

				//search next range
				search_start = Math.min(pos_end_list[0] + 1, chr_info_list[0].length);
				search_end = Math.min(search_start + argv_slice_size, chr_info_list[0].length);//
				//
				match_results.best_match.forEach((match, i) => {
					pos_start_list[i] = match.send + 1;
				});

				if (search_start == search_end || search_end >= chr_info_list[0].length) {
					++fragId;
					break;
				}

				max_delta = argv_init_delta;
			}
			else {
				//search next range
				if (search_end >= chr_info_list[0].length) {
					++fragId;
					break;
				}

				let _seq_from_ref1 = result_coords_list.map(coords => coords.length > 0);
				if (_seq_from_ref1.every(a => a)) {
					// every has align results, but all query loc no overlap
				}
				else {//save 1:3, 2:2, 3:1
					_seq_from_ref1.forEach((value, i) => {
						seq_from_ref1[i] = value;
					});//copy value
				}
				search_end = Math.min(search_end + argv_slice_size, chr_info_list[0].length);//inc search range
				max_delta += argv_add_size;
				
				if (VERBOSE) {
					console.log("not found");
				}

				if (VERBOSE) {
					console.table({
						search_start, search_end,
						pos_start_list,
						find_next_start_list,
						max_delta,
					});

					console.table(match_results.best_match.map(match => match != null));
				}
			}
		}
		catch (ex) {
			console.error(fragId, ex.stack);
			throw ex;
		}
	}//for fragId
	
	if (VERBOSE) {
		console.log({
			max_delta,
			search_start, search_end,
		});
	}

	{
		pos_start_list[0] = pos_end_list[0] + 1;
		let output_coord_list = pos_start_list.map((start, i) => {
			return [start - 1, chr_info_list[i].length];// charAt
		});
		const extract_seq_list = output_coord_list.map(([start, end], i) => chr_seq_list[i].slice(start, end));
		const fa_seq_name_list = chr_name_list.map((chrName, i) => {
			const [start, end] = output_coord_list[i];// charAt
			return `${chrName} ${start + 1}-${end} ${extract_seq_list[i].length}`;
		});
		const extract_length_list = extract_seq_list.map(seq => seq.length);
		const extract_seq_map = fa_seq_name_list.reduce((fasta, seqName, i) => {
			fasta[seqName] = extract_seq_list[i];
			return fasta;
		}, {});

		{
			const coord_text_1 = fragId + "\tstart\t" + [search_start, ...pos_start_list].join("\t|\t");
			const coord_text_2 = fragId + "\t  end\t" + [search_end,   ...chr_info_list.map(chrInfo => chrInfo.length)].join("\t|\t");

			if (VERBOSE) {
				console.log(coord_text_1);
				console.log(coord_text_2);
			}
			
			fs.writeFileSync(`${dataset.tmp_path}/multi_coord_ch${nChr}.txt`, coord_text_1 + "\n" + coord_text_2 + "\n", { flag: "a" });
		}

		if (VERBOSE) {
			console.log("seq_from_ref1", seq_from_ref1);
				
			console.log({
				_global_search_align: _global_search_align.map(a => a.length),
				local_overlap_align: Object.keys(_local_overlap_align).map(key => _local_overlap_align[key] != null),
			});
		}

		if (extract_length_list.every(len => len <= argv_init_max_delta) ||
			seq_from_ref1.every(a => a) ||
			seq_from_ref1.slice(1).every(a => !a)
		) {//seq_from_ref1[0] -> true
			const output_fasta_file_name = `${dataset.tmp_path}/seq_frag/ch${nChr}_${fragId}.fa`;

			if (VERBOSE) {
				console.log(output_fasta_file_name);
			}
			
			saveFasta(output_fasta_file_name, extract_seq_map);
		}
		else {//check translocation
			let like_r1_seq = fa_seq_name_list.filter((seq, idx) => {
				return seq_from_ref1[idx];
			});
			let like_r2_seq = fa_seq_name_list.filter((seq, idx) => {
				return !seq_from_ref1[idx];
			});

			const output_r1_fasta_file_name = `${dataset.tmp_path}/seq_frag/ch${nChr}_${fragId}_ref1.fa`;
			const output_r2_fasta_file_name = `${dataset.tmp_path}/seq_frag/ch${nChr}_${fragId}_ref2.fa`;

			if (VERBOSE) {
				console.log(output_r1_fasta_file_name);
				console.log(output_r2_fasta_file_name);
			}
			saveFasta(output_r1_fasta_file_name, like_r1_seq.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_map[seqName];
				return fasta;
			}, {}));
			saveFasta(output_r2_fasta_file_name, like_r2_seq.reduce((fasta, seqName, i) => {
				fasta[seqName] = extract_seq_map[seqName];
				return fasta;
			}, {}));
		}
		if (WRITE_VERBOSE_LOG) {
			fs.writeFileSync(`${dataset.tmp_path}/meta/meta_ch${nChr}_${fragId}.json`, JSON.stringify({
				id: fragId,
				best_match: [],
				match_group: [],
				coord: output_coord_list,
				length: extract_length_list,
			}, null, "\t"));
		}
		{
			const min_len = Math.min(...extract_length_list);
			const max_len = Math.max(...extract_length_list);
			fs.writeFileSync(`${dataset.tmp_path}/frag_length_ch${nChr}.txt`, [
				fragId,
				...extract_length_list,
				...extract_length_list.map(len => (len / max_len).toFixed(2)),
				min_len, max_len,
				(min_len / max_len).toFixed(2),
			].join("\t") + "\n", { flag: "a" });
		}
	}
	return fragId;
}
