//@ts-check

if (typeof String.prototype.matchAll != "function") {
	const matchAll = require("string.prototype.matchall");
	matchAll.shim();
}

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");

const { readFasta, saveFasta, joinFastaSeq, multialign_to_chrPos_posMap, chrPos_to_multialign_posMap } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);
const argv_output_summary_only = !!argv["--output-summary-only"];

// /** output adp_map.unit8 file */
// const argv_output_map_only = !!argv["--output-map-only"];

const argv_adp_count_window = Number(argv["--adp-count-window"] | 0);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}

const VERBOSE = !!argv["--verbose"];

const genome_info_list = dataset.loadGenomeInfoList();
dataset.load_GC_content();
dataset.load_rDNA_info();

const output_prefix = (dataset.name);


class ChrSummary {
	constructor() {
		/** @type {number|string} */
		this.nChr = "";

		this.snp = 0;
		this.indel = 0;
		this.adp = 0;// ADP + indel
		
		this.snp_clear = 0;
		this.indel_clear = 0;
		this.adp_clear = 0;// ADP + indel
	}

	static get header_list() {
		const aa = new ChrSummary();
		aa.nChr = "ch";
		aa.snp = "SNP";
		aa.indel = "InDel";
		aa.adp = "ADP";
		aa.snp_clear = "SNP (all - rDNA - telomere)";
		aa.indel_clear = "InDel (all - rDNA - telomere)";
		aa.adp_clear = "ADP (all - rDNA - telomere)";
		return aa;
	}
}


const split_adp_chr = false;

/**
 * @param {number} nChr
 * @param {number} min_pos
 * @param {number} max_pos
 * @param {number[][]} ref_pos_map_list
 */
function make_ignore_region(nChr, min_pos, max_pos, ref_pos_map_list) {
	/**
	 * @type {{ start: number; end: number; tag: string; }[]}
	 */
	const range_list = [];

	try {
		{
			let left = min_pos;
			let right = max_pos;

			if (dataset.all_telomere[0][nChr].length) {
				dataset.all_telomere.slice(0, ref_pos_map_list.length).forEach((telomere, genome_idx) => {
					if (ref_pos_map_list[genome_idx] != null) {
						extends_telomere(telomere[nChr], genome_idx);
					}
				});
			}
			if (dataset.telomere[nChr]) {
				extends_telomere(dataset.telomere[nChr], 0);
			}
			range_list.push({
				start: 1,
				end: left,
				tag: "left-telomere",
			});
			range_list.push({
				start: right,
				end: max_pos,
				tag: "right-telomere",
			});
			
			function extends_telomere(telomere, genome_idx) {
				if (!(telomere && telomere[0] && telomere[1])) {
					return false;
				}
				try {
					if (telomere[0][0] != telomere[0][1]) {
						const l_telo_end = ref_pos_map_list[genome_idx][telomere[0][1] - 1 + 1];
						// range_list.push([1, l_telo_end]);
						left = Math.max(left, l_telo_end);
					}
					if (telomere[1][0] != telomere[1][1]) {
						const r_telo_start = ref_pos_map_list[genome_idx][telomere[1][0] - 1];
						// range_list.push([r_telo_start, seq_list[genome_idx].length]);
						right = Math.min(right, r_telo_start);
					}
				}
				catch (ex) {
					console.error(ex);
				}
			}
		}

		if (dataset.rDNA && dataset.rDNA.region && (
			nChr == dataset.rDNA.nChr ||
			nChr == dataset.rDNA.chr
		)) {
			const ref1 = dataset.rDNA.region.map(p => ref_pos_map_list[0][p - 1]);
			const ref2 = dataset.rDNA.region_ref2.map(p => ref_pos_map_list[1][p - 1]);

			range_list.push({
				start: Math.min(ref1[0], ref2[0]),
				end: Math.max(ref1[1], ref2[1]),
				tag: "rDNA",
			});

			// range_list.push([dataset.rDNA_info.region_start, dataset.rDNA_info.region_end]);
		}
	}
	catch (ex) {
	}

	return range_list;
}

main();

// function unit_test() {
// 	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
// 		try {
// 			if (dataset.telomere == null || dataset.telomere[nChr] == null) {
// 				throw { nChr };
// 			}
// 			const [[x11, x12], [x21, x22]] = dataset.telomere[nChr];
// 			if (!is_telomere(nChr, (x12 + x11) / 2)) {
// 				console.log(nChr, "L");
// 				return false;
// 			}
// 			if (!is_telomere(nChr, (x22 + x21) / 2)) {
// 				console.log(nChr, "R");
// 				return false;
// 			}
// 			if (is_telomere(nChr, (x21 + x12) / 2)) {
// 				console.log(nChr, "C");
// 				return false;
// 			}
// 		}
// 		catch (ex) {
// 			console.log(nChr, JSON.stringify(dataset.telomere[nChr], null, "\t"));

// 			for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
// 				console.log(`typeof dataset.telomere[${nChr}] => ${typeof dataset.telomere[nChr]}`);
// 			}

// 			throw ex;
// 		}
// 	}
// 	return true;
// }

function main() {
	if (dataset.mode != "SNP") {
		console.warn("switch SNP mode");
		dataset.mode = "SNP";
	}
	// if (dataset.telomere && Object.values(dataset.telomere).length) {
	// 	if (!unit_test()) {
	// 		throw "unit_test";
	// 	}
	// }

	if (!argv_output_summary_only) {
		main_output_table();
	}
	output_viewer();
}

function main_output_table() {
	const output_adp_map_path = `${dataset.output_path}/adp_map.uint8`;
	//const stream = fs.createWriteStream(output_adp_map_path);
	
	if (!split_adp_chr) {
		fs.writeFileSync(`${dataset.output_path}/snp_map.txt`, [
			"ch",
			"length(bp)",
			"n",
		].join("\t") + "\n");

		fs.writeFileSync(`${dataset.output_path}/adp_map.txt`, [
			"ch",
			"length(bp)",
			"n",
		].join("\t") + "\n");
	}

	const chr_summary_list = genome_info_list[0].chr_list.map(_ => new ChrSummary());
	const total_summary = new ChrSummary();

	const final_summary = {
		chr: chr_summary_list,
		total: total_summary,
	};

	// genome_telomere_total_length(genome_idx)

	/**
	 * sss_pos_len_rows = sss_pos_len_rows.concat(...)
	 */
	let sss_pos_len_rows = [];

	// /** @type {{ [ref_target:string]: { snp:number, indel:number, adp:number } }[]} */
	// const xCmp_list = [];
	
	/** @type {{ snp:number, indel:number, adp:number }[][][]} */
	const xCmp_list = [];

	genome_info_list[0].chr_list.forEach((chrInfo, chrIdx) => {
		const nChr = chrIdx + 1;
		const input_fa = readFasta(`${dataset.output_path}/mafft_ch${nChr}.fa`);
		
		// const chr_name_list = genome_info_list.map(genomeInfo => genomeInfo.chr_list[chrIdx].chr);
		// const seq_list = chr_name_list.map(chrName => input_fa[chrName]);
		
		const chr_name_list = Object.keys(input_fa);
		const seq_list = Object.values(input_fa);

		let check_adp = 0;
		for (let i = 0; i < seq_list[0].length; ++i) {
			if (seq_list[0][i] != seq_list[1][i] && (seq_list[1][i] != "N" && seq_list[0][i] != "N")) {
				++check_adp;
			}
		}
		console.log(nChr, " check_adp:", check_adp);
		
		const map_to_seq_1 = multialign_to_chrPos_posMap(seq_list[0]);
		const map_to_seq_2 = multialign_to_chrPos_posMap(seq_list[1]);

		const ref1_map_to_pos = chrPos_to_multialign_posMap(seq_list[0]);
		const ref2_map_to_pos = chrPos_to_multialign_posMap(seq_list[1]);

		const ignore_region = make_ignore_region(nChr, 0, seq_list[0].length, [ref1_map_to_pos, ref2_map_to_pos]);
		// if (nChr == 4) {
		// 	console.table(ignore_region);
		// 	console.log(is_telomere_or_rDNA(nChr, [
		// 		ref1_map_to_pos[19379],
		// 		ref2_map_to_pos[798]
		// 	]));
		// }

		/**
		 * @param {number} nChr
		 * @param {number[]} pos_array
		 */
		function is_telomere_or_rDNA(nChr, pos_array) {
			return ignore_region.some(region => {
				return pos_array.some(pos => {
					return region.start <= pos && pos <= region.end;
				});
			});
			// nChr, pos_list
		}

		// console.log(chr_name_list);
		// console.log(Object.keys(input_fa));
		// console.log(seq_list.map(a => a.length));

		/** @type {Uint8Array[]} */
		const ui8_ab = seq_list.slice(1).map(seq => new Uint8Array(seq.length));

		const snp_rows = [];
		const indel_rows = [];
		const _adp_list = [];

		const snp_clear_rows = [];
		const indel_clear_rows = [];
		const _adp_clear_list = [];

		let pos_list = chr_name_list.map(_ => 1);
		for (let tseta_pos = 0; tseta_pos < seq_list[0].length; ++tseta_pos) {
			const ref1 = seq_list[0][tseta_pos];
			const columns = [];
			let has_adp, is_indel;

			columns.push(pos_list[0], ref1);
			
			const current_pos = seq_list.map((target_seq, target_seqIdx) => {
				let pos;

				if (seq_list[target_seqIdx][tseta_pos] == "N") {
					// skip
				}
				else if (seq_list[target_seqIdx][tseta_pos] == "-") {
					is_indel = true;
					pos = pos_list[target_seqIdx] - 1;
				}
				else {
					pos = pos_list[target_seqIdx];
					++pos_list[target_seqIdx];
				}

				if (target_seqIdx >= 1) {// not ref
					const target = seq_list[target_seqIdx][tseta_pos];
					if (target != ref1) {
					// if (target != ref1 && !is_telomere_or_rDNA(nChr, tseta_pos)) {
						columns.push(pos, target);
						has_adp = true;
						
						// if ((pos_list[target_seqIdx] - 1) > 0) {
						// 	console.log(target_seqIdx, pos_list[target_seqIdx] - 1);
						// }
						ui8_ab[target_seqIdx - 1][pos_list[target_seqIdx] - 1] = 1;
					}
					else {
						columns.push(null, null);
					}
				}

				return pos;
			});

			// const sub_list = seq_list.slice(1);
			// sub_list.map((sub_seq, sub_seqIdx) => {
			// 	const pos_1 = current_pos[sub_seqIdx + 1];
				
			// 	sub_list.map((sub2_seq, sub2_seqIdx) => {
			// 		const x_key = [sub_seqIdx, sub2_seqIdx].join("_");
					
			// 		if (sub_seq[i] != sub2_seq[1]) {
			// 			xCmp[x_key] = xCmp[x_key] || {
			// 				snp: 0,
			// 				indel: 0,
			// 				adp: 0,
			// 			};

			// 			if (sub_seq[i] == "-" || sub2_seq[i] == "-") {
			// 				xCmp[x_key].indel += 1;
			// 			}
			// 			else {
			// 				xCmp[x_key].snp += 1;
			// 			}
			// 			xCmp[x_key].adp += 1;
			// 		}
			// 	});
			// });

			if (has_adp) {
				if (is_indel) {
					indel_rows.push(columns);
					if (!is_telomere_or_rDNA(nChr, [
						ref1_map_to_pos[pos_list[0]],
						ref2_map_to_pos[pos_list[1]],
					])) {
						indel_clear_rows.push(columns);
					}
				}
				else {//is snp
					snp_rows.push(columns);
					if (!is_telomere_or_rDNA(nChr, [
						ref1_map_to_pos[pos_list[0]],
						ref2_map_to_pos[pos_list[1]],
					])) {
						snp_clear_rows.push(columns);
					}
				}
				_adp_list.push(columns);
				if (!is_telomere_or_rDNA(nChr, [
					ref1_map_to_pos[pos_list[0]],
					ref2_map_to_pos[pos_list[1]],
				])) {
					_adp_clear_list.push(columns);
				}
			}
		}//for each bp

		xCmp_list[chrIdx] = [];

		for (let i = 0; i < seq_list.length; ++i) {
			const ref_seq = [...seq_list[i]];

			xCmp_list[chrIdx][i] = [];

			for (let j = 0; j < seq_list.length; ++j) {
				if (i != j) {
					// const k = i + "_" + j;
					
					const adps = ref_seq.map((v, pos) => {
						if (seq_list[j][pos] == "N" ||
							v == "N"
						) {
							return null;
						}
						else if (seq_list[j][pos] != v) {
							return [pos, v != "-" && seq_list[j][pos] != "-"];
						}
						else {
							return null;
						}
					}).filter(a => a);

					xCmp_list[chrIdx][i][j] = {
						snp: 0,
						indel: 0,
						adp: 0,
					};
					xCmp_list[chrIdx][i][j].adp = adps.length;
					xCmp_list[chrIdx][i][j].snp = adps.filter(aa => aa[1]).length;
					xCmp_list[chrIdx][i][j].indel = adps.length - xCmp_list[chrIdx][i][j].snp;
				}
			}
		}
		// xCmp_list.push(xCmp);
		
		// if (argv_output_map_only) {
		// 	ui8_ab.forEach(arr => {
		// 		//const output_path = `${genomeInfo.name}_adp_map_${nChr}.uint8`;
		// 		stream.write(arr);
		// 	});
		// }
		// else {
			try {
				final_summary.total.snp += snp_rows.length;
				final_summary.chr[chrIdx].snp = snp_rows.length;
				final_summary.total.snp_clear += snp_clear_rows.length;
				final_summary.chr[chrIdx].snp_clear = snp_clear_rows.length;
				output_table("snp", snp_rows);
			}
			catch (ex) {
				console.error(ex);
			}
			
			try {
				final_summary.total.indel += indel_rows.length;
				final_summary.chr[chrIdx].indel = indel_rows.length;
				final_summary.total.indel_clear += indel_clear_rows.length;
				final_summary.chr[chrIdx].indel_clear = indel_clear_rows.length;
				output_table("indel", indel_rows);
			}
			catch (ex) {
				console.error(ex);
			}
			
			try {
				final_summary.total.adp += _adp_list.length;
				final_summary.chr[chrIdx].adp = _adp_list.length;
				final_summary.total.adp_clear += _adp_clear_list.length;
				final_summary.chr[chrIdx].adp_clear = _adp_clear_list.length;
				output_table("adp", _adp_list);
			}
			catch (ex) {
				console.error(ex);
			}
		//}

		if (argv_adp_count_window != null && argv_adp_count_window > 0) {
			let count_list = [];

			_adp_list.forEach(([pos, ...a]) => {
				let wId = Math.trunc(Number(pos) / argv_adp_count_window);
				count_list[wId] = (count_list[wId] | 0) + 1;
			});

			const output_path = `${dataset.output_path}/adp_per_${argv_adp_count_window}_ch${nChr}.txt`;
			const ws = fs.createWriteStream(output_path);
			ws.write(["chr", "start", "end", "count"].join("\t") + "\n");
			count_list.forEach((count, wId) => {
				const start = wId * argv_adp_count_window + 1;
				const end = start - 1 + argv_adp_count_window;
				let row = [nChr, start, end, count].join("\t") + "\n";
				ws.write(row);
			});
			ws.end();
			console.log("output table:", output_path);
		}

		/**
		 * 
		 * @param {"snp"|"indel"|"adp"} type
		 * @param {(string|number)[][]} rows
		 */
		function output_table(type, rows) {
			if (VERBOSE) {
				console.log(nChr, type + ":", rows.length);
			}
	
			const output_path = `${dataset.output_path}/${type}_ch${nChr}.txt`;
			
			const output_header = ["ref_pos", "ref"];
			seq_list.slice(1).map((_, targetIdx) => output_header.push("target" + (1 + targetIdx) + "_pos", "target" + (1 + targetIdx)));
			
			const ws = fs.createWriteStream(output_path);
			ws.write(output_header.join("\t") + "\n");
			rows.forEach(cols => {
				let row = cols.map(v => v != null ? v : "").join("\t") + "\n";
				ws.write(row);
			});
			ws.end();
			console.log("output table:", output_path);

			//fs.writeFileSync(output_path, output_header.join("\t") + "\n");
			
			// rows.forEach(cols => {
			// 	let row = cols.map(v => v != null ? v : "").join("\t") + "\n";
			// 	fs.writeFile(output_path, row, { flag: "a" }, function (err) {
			// 		if (err) {
			// 			console.error(err);
			// 		}
			// 	});
			// });
			
			// const output_table = rows.map(cols => cols.map(v => v != null ? v : "").join("\t")).join("\n");
			// //console.log(output_table.length);
			// fs.writeFile(output_path, output_table, { flag: "a" }, function (err) {
			// 	if (err) {
			// 		console.error(err);
			// 	}
			// 	console.log("output:", output_path);
			// });
		}

		// 20200713
		sss_pos_len_rows = make_sss_region(seq_list, map_to_seq_1, map_to_seq_2, nChr, final_summary, chrIdx, sss_pos_len_rows, is_telomere_or_rDNA);
	});// chr_list.forEach()
	
	//stream.end();
	
	// 20200713
	xCmp_table(xCmp_list);// xCmp

	{
		const file_path = `${dataset.output_path}/${output_prefix}_final_summary.json`;
		fs.writeFile(file_path, JSON.stringify(final_summary, null, "\t"), make_writeFile_callback(file_path));
	}
	{
		final_summary.chr.forEach((v, i) => {
			v.nChr = romanize(Number(i) + 1);
		});
		final_summary.total.nChr = "total";
		final_summary.chr.push(final_summary.total);

		let text = "";
		text += Object.assign(
			Object.keys(final_summary.chr[0]),
			Object.values(ChrSummary.header_list)
		).join("\t") + "\n";// ["Chromosome", "SNP", "INDEL", "ADP"].join("\t") + "\n";
		// text += final_summary.chr.map((v, i) => [romanize(Number(i) + 1), v.snp, v.indel, v.adp].join("\t")).join("\n") + "\n";
		text += final_summary.chr.map((v, i) => Object.values(v).join("\t")).join("\n") + "\n";
		// text += ["total", final_summary.total.snp, final_summary.total.indel, final_summary.total.adp].join("\t") + "\n";

		const file_path = `${dataset.output_path}/${output_prefix}_final_summary.txt`;
		fs.writeFile(file_path, text, make_writeFile_callback(file_path));
	}

	{
		const text = sss_pos_len_rows.join("\n");
		const file_path = `${dataset.output_path}/${output_prefix}_sss_pos_len.txt`;
		fs.writeFile(file_path, text, make_writeFile_callback(file_path));
	}
}

/**
 * @param {string[]} seq_list
 * @param {number[]} map_to_seq_1
 * @param {number[]} map_to_seq_2
 * @param {number} nChr
 * @param {{ chr: ChrSummary[]; total: ChrSummary; }} final_summary
 * @param {number} chrIdx
 * @param {*} sss_pos_len_rows
 * @param {*} is_telomere_or_rDNA
 */
function make_sss_region(seq_list, map_to_seq_1, map_to_seq_2, nChr, final_summary, chrIdx, sss_pos_len_rows, is_telomere_or_rDNA) {
	if (dataset.progeny_list.length >= 1) {
		// && seq_list.slice(2).filter(ss => ss[i] != "-").length == 2
		// const output_map = {};
		const min_indel_len = 1;

		const qqq = [...seq_list[0]].map((a, i) => a != "-" && a != "N" && seq_list[1][i] == "-" ? "1" : "0").join("");
		/** @type {RegExpMatchArray[]} */
		// 20200714 // const qqq_sss = [...qqq.matchAll(/1{10,}/g)];
		const qqq_sss = [...qqq.matchAll(new RegExp(String.raw`1{${min_indel_len},}`, "g"))];
		const qqq_sss_t_len = qqq_sss.reduce((t, v) => t + v[0].length, 0);

		const ccc = [...seq_list[0]].map((a, i) => a == "-" && seq_list[1][i] != "-" && seq_list[1][i] != "N" ? "1" : "0").join("");
		/** @type {RegExpMatchArray[]} */
		// 20200714 // const ccc_sss = [...ccc.matchAll(/1{10,}/g)];
		const ccc_sss = [...ccc.matchAll(new RegExp(String.raw`1{${min_indel_len},}`, "g"))];
		const ccc_sss_t_len = ccc_sss.reduce((t, v) => t + v[0].length, 0);

		// output_map.Chromosome = nChr;
		// output_map.QC_SNP = ref_snp.length;
		// output_map.QC_InDel = ref_adp_indel.length;
		// // output_map.QC_ADP = ref_adp.length;
		if (true && "SNP") {
			function aaaa() {
				return [...seq_list[0]].map((q, i) => [
					map_to_seq_1[i],
					q,
					seq_list[1][i],
					map_to_seq_2[i],
					i
				]);
			}

			function bbbb() {
				return aaaa().map(([ref1_pos, q, c, ref2_pos, i]) => [
					ref1_pos,
					q, c,
					ref2_pos,
					i,
					!is_telomere_or_rDNA(nChr, [Number(i) + 1]) && q != c && (q != "-" && c != "-") && (q != "N" && c != "N")
				]);
			}

			function cccc(bbb) {
				return bbb.map(([ref1_pos, q, c, ref2_pos, i, isSNP]) => isSNP ? "1" : "0").join("");
			}

			/** @returns {RegExpMatchArray[]} */
			function mmmm(bbbb) {
				return [...cccc(bbbb).matchAll(/1+/g)];
			}

			const rrrr = (function () {
				return mmmm(bbbb()).map(a => a[0].length).sort((a, b) => a - b);
			})();

			const max_n_indel = Infinity; //50 + 
			const indel_map = {};
			rrrr.forEach(n_indel => {
				const n_indel_ID = Math.min(n_indel, max_n_indel);
				indel_map[n_indel_ID] = (indel_map[n_indel_ID] | 0) + 1;
			});
			// fs.writeFileSync(`${dataset.output_path}/ch${nChr}_snp_map.json`, JSON.stringify(rrr));
			console.log({
				"snp:": Object.keys(indel_map).reduce((prev, len) => prev + Number(len) * indel_map[len], 0),
			});

			const out_text = Object.keys(indel_map).map(n_indel => [nChr, n_indel, indel_map[n_indel]]).filter(([nChr, n_indel, n]) => n > 0).map(row => row.join("\t")).join("\n");
			if (split_adp_chr) {
				fs.writeFileSync(`${dataset.output_path}/ch${nChr}_snp_map.txt`, `ch\tlength(bp)\tn\n${out_text}`);
			}
			else {
				fs.writeFileSync(`${dataset.output_path}/snp_map.txt`, out_text + "\n", { flag: "a" });
			}
		}
		if (true && "InDel") {
			function aaaa() {
				return [...seq_list[0]].map((q, i) => [
					map_to_seq_1[i],
					q,
					seq_list[1][i],
					map_to_seq_2[i],
					i
				]);
			}

			function bbbb() {
				return aaaa().map(([ref1_pos, q, c, ref2_pos, i]) => [
					ref1_pos,
					q, c,
					ref2_pos,
					i,
					!is_telomere_or_rDNA(nChr, [Number(i) + 1]) && q != c && (q != "N" && c != "N"),
					!is_telomere_or_rDNA(nChr, [Number(i) + 1]) && q != c && ((q == "-" && c != "-") || (q != "-" && c == "-")) && (q != "N" && c != "N") //q is del or c is del
				]);
			}

			function cccc(bbb) {
				return bbb.map(([ref1_pos, q, c, ref2_pos, i, isADP, isIndel]) => isIndel ? "1" : "0").join("");
			}

			/** @returns {RegExpMatchArray[]} */
			function mmmm(bbbb) {
				return [...cccc(bbbb).matchAll(/1+/g)];
			}

			const rrrr = (function () {
				return mmmm(bbbb()).map(a => a[0].length).sort((a, b) => a - b);
			})();

			const max_n_indel = Infinity; //50 + 
			const indel_map = {};
			rrrr.forEach(n_indel => {
				const n_indel_ID = Math.min(n_indel, max_n_indel);
				indel_map[n_indel_ID] = (indel_map[n_indel_ID] | 0) + 1;
			});
			// fs.writeFileSync(`${dataset.output_path}/ch${nChr}_adp_map.json`, JSON.stringify(rrr));
			console.log({
				"adp:": Object.keys(indel_map).reduce((prev, len) => prev + Number(len) * indel_map[len], 0),
			});

			const out_text = Object.keys(indel_map).map(n_indel => [nChr, n_indel, indel_map[n_indel]]).filter(([nChr, n_indel, n]) => n > 0).map(row => row.join("\t")).join("\n");
			if (split_adp_chr) {
				fs.writeFileSync(`${dataset.output_path}/ch${nChr}_adp_map.txt`, `ch\tlength(bp)\tn\n${out_text}`);
			}
			else {
				fs.writeFileSync(`${dataset.output_path}/adp_map.txt`, out_text + "\n", { flag: "a" });
			}
		}

		// output_map[`${dataset.parental_list[0]} specific`] = qqq_sss_t_len;
		// output_map[`${dataset.progeny_list[0]} specific`] = ccc_sss_t_len;
		final_summary.chr[chrIdx][`${dataset.parental_list[0]} specific`] = qqq_sss_t_len;
		final_summary.chr[chrIdx][`${dataset.progeny_list[0]} specific`] = ccc_sss_t_len;

		// const s_tbl_qc_bp_len = Object.keys(output_map).map(k => `${output_map[k]}`).join("\t");
		// fs.writeFile(`${dataset.output_path}/${argv_output_prefix}_sss_total.txt`, s_tbl_qc_bp_len, (err) => err ? console.error(err) : void 0);
		// const s_tbl_qc_sss = [
		// 	qqq_sss.map(m => [dataset.parental_list[0], nChr, data.pos_ref1_uint32array[m.index], m[0].length].join("\t")).join("\n"),
		// 	ccc_sss.map(m => [dataset.parental_list[1], nChr, data.pos_ref1_uint32array[m.index], m[0].length].join("\t")).join("\n")
		// ].join("\n");
		// fs.writeFile(`${dataset.output_path}/${argv_output_prefix}_sss_pos_len.txt`, s_tbl_qc_sss, (err) => err ? console.error(err) : void 0);
		sss_pos_len_rows = sss_pos_len_rows.concat(
			qqq_sss.map(m => [dataset.parental_list[0], nChr, map_to_seq_1[m.index], m[0].length].join("\t")),
			ccc_sss.map(m => [dataset.progeny_list[0], nChr, map_to_seq_2[m.index], m[0].length].join("\t"))
		);
	}
	return sss_pos_len_rows;
}

/**
 * @param {{ snp:number, indel:number, adp:number }[][][]} xCmp_list [chrIdx][target_1][target_2]
 */
function xCmp_table(xCmp_list) {
	if (dataset.progeny_list.length > 1) {
		const total_xCmp = genome_info_list.map(_ => {
			return genome_info_list.map(_ => {
				return {
					snp: 0,
					indel: 0,
					adp: 0,
				};
			});
		});

		//
		genome_info_list[0].chr_list.forEach((chrInfo, chrIdx) => {
			const nChr = chrIdx + 1;

			output_xCmp(`Ch${romanize(nChr)}`, xCmp_list[chrIdx], `${dataset.output_path}/${output_prefix}_xCmp_ch${nChr}.txt`);

			for (let i = 0; i < genome_info_list.length; ++i) {
				for (let j = 0; j < genome_info_list.length; ++j) {
					if (i != j) {
						total_xCmp[i][j].snp += xCmp_list[chrIdx][i][j].snp;
						total_xCmp[i][j].indel += xCmp_list[chrIdx][i][j].indel;
						total_xCmp[i][j].adp += xCmp_list[chrIdx][i][j].adp;
					}
					else { // if (i == j) {
						total_xCmp[i][j].snp = null; // clear
						total_xCmp[i][j].indel = null; // clear
						total_xCmp[i][j].adp = null; // clear
					}
				}
			}
		});
		output_xCmp("total", total_xCmp, `${dataset.output_path}/${output_prefix}_xCmp_total.txt`);
		/**
		 * @param {string} tag
		 * @param {xCmp_list[0]} xCmp
		 * @param {string} output_path
		 */
		function output_xCmp(tag, xCmp, output_path) {
			let ttt = `${tag}\t${output_path}\n\n`;

			ttt += "snp\n";
			ttt += "\t" + genome_info_list.map(a => a.name).join("\t") + "\n";
			ttt += xCmp.map((a, j) => genome_info_list[j].name + "\t" + a.map(b => b.snp).join("\t")).join("\n") + "\n\n";

			ttt += "indel\n";
			ttt += "\t" + genome_info_list.map(a => a.name).join("\t") + "\n";
			ttt += xCmp.map((a, j) => genome_info_list[j].name + "\t" + a.map(b => b.indel).join("\t")).join("\n") + "\n\n";

			ttt += "adp\n";
			ttt += "\t" + genome_info_list.map(a => a.name).join("\t") + "\n";
			ttt += xCmp.map((a, j) => genome_info_list[j].name + "\t" + a.map(b => b.adp).join("\t")).join("\n") + "\n\n";

			fs.writeFile(output_path, ttt, make_writeFile_callback(output_path));
		}

		// {
		// 	fs.writeFile(`${dataset.output_path}/${output_prefix}_xCmp_list.json`, JSON.stringify(xCmp_list, null, "\t"), make_writeFile_callback("xCmp_list.json"));
		// }
	}
}

// /**
//  * @param {number} nChr
//  * @param {number[]} ref_pos_list
//  */
// function is_telomere_or_rDNA(nChr, ref_pos_list) {
// 	return ref_pos_list.some(ref_pos => {
// 		return is_telomere(nChr, ref_pos) || is_rDNA(nChr, ref_pos);
// 	});
// }

// /**
//  * @param {number} nChr
//  * @param {number} ref1_pos
//  */
// function is_rDNA(nChr, ref1_pos) {
// 	// return false;

// 	const [
// 		start, end
// 	] = dataset.rDNA.region;

// 	if (start == end) {
// 		return false;
// 	}
	
// 	if (nChr == dataset.rDNA.nChr) {
// 		if (ref1_pos >= start &&
// 			ref1_pos <= end
// 		) {
// 			return true;
// 		}
// 	}

// 	// if (dataset.rDNA_info && nChr == dataset.rDNA_info.chr) {
// 	// 	if (pos >= dataset.rDNA_info.alignment_start && pos <= dataset.rDNA_info.alignment_end) {
// 	// 		return true;
// 	// 	}
// 	// }
// }

// /**
//  * @param {number} nChr
//  * @param {number} ref1_pos
//  */
// function is_telomere(nChr, ref1_pos) {
// 	const telomere = dataset.telomere[nChr];
// 	if (telomere) {
// 		const [
// 			[
// 				left_start, left_end,
// 			], 
// 			[
// 				right_start, right_end,
// 			]
// 		] = telomere;

// 		// console.log({
// 		// 	nChr,
// 		// 	ref1_pos,
// 		// 	telomere
// 		// });

// 		return (
// 			ref1_pos <= left_end || ref1_pos >= right_start
// 		);
// 	}
// 	else {
// 		return false;
// 	}
// }

/**
 * @param {string} file_path
 */
function make_writeFile_callback(file_path) {
	const foo = {
		/**
		 * @param {Error} err
		 */
		[file_path]: function (err) {
			if (err) {
				console.error(err);
			}
		}
	};
	return foo[file_path];
}

function output_viewer() {
	let template_html = fs.readFileSync(`${__dirname}/template_viewer.html`).toString();
	let analyser_js = fs.readFileSync(`${__dirname}/analyser.js`).toString();
	let web_ui_js = fs.readFileSync(`${__dirname}/web_ui.js`).toString();

	dataset.genome_info_list = genome_info_list;
	
	{
		dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
			const nChr = chrIdx + 1;
			return `mafft_ch${nChr}.fa`;
		});
		const output_path = `${dataset.output_path}/debug_${output_prefix}_snp.html`;
		
		output_html(template_html, output_path, false);
		if (VERBOSE) {
			console.log("output debug viewer", output_path);
		}
	}

	//single html
	{
		let all_seq = load_all_seq();
		dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
			return all_seq[chrIdx];
		});
		const output_path = `${dataset.output_path}/${output_prefix}_snp.html`;

		output_html(template_html, output_path, true);
	
		console.log("output viewer:", output_path);
	}
	
	function output_html(input_html, output_path, inline_script) {
		let output_html = input_html.replace(`<script id="dataset.json" type="application/json"></script>`, `<script id="dataset.json" type="application/json">${JSON.stringify(dataset)}</script>`);
		//output_html = output_html.replace(`<script id="all_seq.json" type="application/json"></script>`, `<script id="all_seq.json" type="application/json">${JSON.stringify(all_seq)}</script>`);

		if (inline_script) {
			output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script>${analyser_js}</script>`);
			output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script>${web_ui_js}</script>`);
		}
		else {
			output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script src="../src/analyser.js"></script>`);
			output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script src="../src/web_ui.js"></script>`);
		}

		fs.writeFileSync(output_path, output_html);
	}
}

function load_all_seq() {
	let all_seq = [];
	genome_info_list[0].chr_list.forEach(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const seq = `${dataset.output_path}/mafft_ch${nChr}.fa`;
		
		const input_fasta = readFasta(seq);

		all_seq[chrIdx] = input_fasta;
	});
	return all_seq;
}

/**
 * @param {number[][]} chr_telo
 */
	function chr_telomere_total_length(chr_telo) {
	/**
	 * @param {number[]} tel
	 */
	function get_tel_len(tel) {
		try {
			return tel[1] - tel[0] + 1;
		}
		catch (ex) {
			console.log(tel);
			return 0;
		}
	}
	try {
		const [left, right] = chr_telo;
		const l = get_tel_len(left);
		const r = get_tel_len(right);
		return l + r;
	}
	catch (ex) {
		console.log(chr_telo);
		return 0;
	}
}

/**
 * @param {number} genome_idx
 * @see {@link table_telomere_total_length}
 */
function genome_telomere_total_length(genome_idx) {
	const len = genome_info_list[genome_idx].chr_list.map((_, chr_idx) => {
		const nChr = chr_idx + 1;
		return chr_telomere_total_length(dataset.all_telomere[genome_idx][nChr]);
	}).reduce((aa, v) => aa + v, 0);
	return len;
}

function table_telomere_total_length() {
	return genome_info_list.map((genome_info, genome_idx) => [genome_info.name, genome_telomere_total_length(genome_idx)].join("\t")).join("\n");
}

/**
 * @see {@link https://stackoverflow.com/a/9083076|stackoverflow}
 * @see {@link http://blog.stevenlevithan.com/archives/javascript-roman-numeral-converter}
 * @param {number} num
 * @returns {string}
 */
function romanize(num) {
	if (isNaN(num))
		return String(num);
	var digits = String(+num).split(""),
		key = ["","C","CC","CCC","CD","D","DC","DCC","DCCC","CM",
			"","X","XX","XXX","XL","L","LX","LXX","LXXX","XC",
			"","I","II","III","IV","V","VI","VII","VIII","IX"],
		roman = "",
		i = 3;
	while (i--)
		roman = (key[+digits.pop() + (i * 10)] || "") + roman;
	return Array(+digits.join("") + 1).join("M") + roman;
}

/**
 * @see {@link http://blog.stevenlevithan.com/archives/javascript-roman-numeral-converter}
 * @param {string} str
 * @returns {number}
 */
function deromanize (str) {
	var	str = str.toUpperCase(),
		validator = /^M*(?:D?C{0,3}|C[MD])(?:L?X{0,3}|X[CL])(?:V?I{0,3}|I[XV])$/,
		token = /[MDLV]|C[MD]?|X[CL]?|I[XV]?/g,
		key = {M:1000,CM:900,D:500,CD:400,C:100,XC:90,L:50,XL:40,X:10,IX:9,V:5,IV:4,I:1},
		num = 0, m;
	if (!(str && validator.test(str)))
		return NaN;
	while (m = token.exec(str))
		num += key[m[0]];
	return num;
}

