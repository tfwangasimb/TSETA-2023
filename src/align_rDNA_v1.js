//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { BlastnCoord, execAsync, exec_blastn, exec_blastn_Ex, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { run_mafft } = require("./run_mafft.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { validation_chr } = require("./validation_seq.js");
const { Dataset, RibosomalDNA_Data } = require("./dataset.js");

const argv = argv_parse(process.argv);

const DEBUG = !!argv["--debug"];
const VERBOSE = process.argv.indexOf("--verbose") >= 0;
const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}

const genome_info_list = dataset.loadGenomeInfoList();

main();

async function main() {
	program_log(`${dataset.name}.log.txt`, "start");

	if (argv["-chr"] != null && argv["-rDNA"] != null) {
		let nChr = Number(argv["-chr"]);
		let fasta_filepath = String(argv["-rDNA"] || "");

		if (!Number.isSafeInteger(nChr)) {
			console.error("error", "rDNA", "number of chromosome:", nChr);
		}

		if (!fs.existsSync(fasta_filepath)) {
			console.error("error", "rDNA", "fasta file path:", fasta_filepath);
			return;
		}

		dataset.rDNA = new RibosomalDNA_Data();
		dataset.rDNA.nChr = nChr;
		dataset.rDNA.sequence = fasta_filepath;
		fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));
	}

	try {
		await re_align_rdna();
	
		validation_chr(Number(dataset.rDNA.nChr), dataset.output_path, true);
	}
	catch (ex) {
		throw ex;
	}
	finally {
		fs.writeFileSync(argv_dataset_path, JSON.stringify(dataset, null, "\t"));
		console.log("save:", argv_dataset_path);
	}

	console.log("next step:", "analysis");
	console.log("snp calling command:", `node ${__dirname}/snp_summary.js -dataset ${argv_dataset_path}`);
	console.log("tetrad analysis command:", `node ${__dirname}/tetrad_summary.js -dataset ${argv_dataset_path} -min-co 5000`);

	program_log(`${dataset.name}.log.txt`, "exit");
}

async function re_align_rdna() {
	const nChr = Number(dataset.rDNA.nChr);
	const chrIdx = nChr - 1;
	let rDNA_filePath = dataset.rDNA.sequence;

	let mafft_fa = readFasta(`${dataset.tmp_path}/mafft_ch${nChr}.fa`);

	/** chrPos_to_multialign_posMap */
	let chr_to_ma_posmap = {};
	/** multialign_to_chrPos_posMap */
	let ma_to_chr_posmap = {};

	/** @type {{ [seqName: string]: rDNA_Data }} */
	let rdna_data = {};

	let chr_list = genome_info_list.map(genome_info => genome_info.chr_list[chrIdx].chr);
	
	let mafft_fa_keys = Object.keys(mafft_fa);//200416
	chr_list.forEach((genomeChrSeqName, genomeIndex) => {
		if (genomeChrSeqName != mafft_fa_keys[genomeIndex]) {
			console.warn(`mafft:${mafft_fa_keys[genomeIndex]} -> info:${genomeChrSeqName}`)
			
			mafft_fa[genomeChrSeqName] = mafft_fa[mafft_fa_keys[genomeIndex]];
			delete mafft_fa[mafft_fa_keys[genomeIndex]];
		}
	});

	let promise = chr_list.map(async function (genomeChrSeqName, genomeIndex) {
		const mafft_seq = mafft_fa[genomeChrSeqName];

		chr_to_ma_posmap[genomeChrSeqName] = _chrPos_to_multialign_posMap(mafft_seq);
		ma_to_chr_posmap[genomeChrSeqName] = _multialign_to_chrPos_posMap(mafft_seq);
		
		rdna_data[genomeChrSeqName] = await find_rDNA_use_blastn(rDNA_filePath, genomeIndex, nChr);
	});
	await Promise.all(promise);
	
	/** @type {{ [seqName: string]: string }} */
	let final_mafft_seq = {};

	function pos_to_ma(seqName) {
		const rd = rdna_data[seqName];

		let posmap = chr_to_ma_posmap[seqName];

		let ma_start = posmap[rd.start - 1];
		let ma_end = posmap[rd.end - 1];

		const ma_seq = mafft_fa[seqName];

		while (ma_seq[ma_start - 1] == "-") {
			ma_start = ma_start - 1;
			// console.log(seqName, ma_start, ma_seq[ma_start]);
		}

		while (ma_seq[ma_end + 1] == "-") {
			ma_end = ma_end + 1;
			// console.log(seqName, ma_end, ma_seq[ma_end]);
		}
		// console.log(seqName, ma_start, ma_seq.slice(ma_start - 1, ma_start + 2));
		// console.log(seqName, ma_end, ma_seq.slice(ma_end - 1, ma_end + 2));

		return {
			ma_start,
			ma_end,
		};
	}

	let ma_rdna_range_list = chr_list.map(function (seqName, seqIndex) {
		let { ma_start, ma_end } = pos_to_ma(seqName);
		return {
			ma_start, ma_end,
		};
	});
	let min_ma_start = Math.min(...ma_rdna_range_list.map(a => a.ma_start));
	let max_ma_end = Math.max(...ma_rdna_range_list.map(a => a.ma_end)) + 1;//???

	let ma_start_delta = ma_rdna_range_list.map(a => a.ma_start - min_ma_start);
	
	let ref1_rDNA_strand = rdna_data[chr_list[0]].strand;
	/** @type {string[][]} */
	let rDNA_repeatSeqList_genomeList = [];
	// /** @type {string[]} */
	// let rDNA_last_IGS_genomeList = [];

	/** @type {{start:number, end:number}[][]} */
	let ex_ranges = chr_list.map(_ => []);

	// maybe rDNA
	let ex_range = [];
	let rdna_seq_list = chr_list.map(function (seqName, seqIndex) {//raw rDNA seq
		let posmap = ma_to_chr_posmap[seqName];
		const rd = rdna_data[seqName];
		
		let start = posmap[min_ma_start - 1];
		let end = posmap[max_ma_end - 1];

		if (ref1_rDNA_strand == 1) {
			throw new Error("if (ref1_rDNA_strand == 1) {");
		}
		else if (ref1_rDNA_strand == -1) {
			let prev_rep_end_pos = 0;
			
			//repeats: (IGS + rDNA) + (IGS + rDNA) + ...
			rDNA_repeatSeqList_genomeList[seqIndex] = rd.info.repeats.map((_range, _range_idx) => {
				let range = [..._range].sort((a, b) => a - b);

				if (!prev_rep_end_pos) {
					prev_rep_end_pos = start + 1;
				}

				let _rep_start = range[0];
				let rep_end = Math.min(range[1] + 1, end);
				
				if (prev_rep_end_pos > _rep_start) {
					throw new Error("????");
				}

				ex_ranges[seqIndex].push({
					start: prev_rep_end_pos,
					end: rep_end,
				});//IGS + rDNA
				let seq = rd.chrSeq.slice(prev_rep_end_pos, rep_end);//IGS + rDNA

				prev_rep_end_pos = rep_end;

				return seq;
			});
			if (prev_rep_end_pos < end) {//IGS
				let igs_seq = rd.chrSeq.slice(prev_rep_end_pos, end);
				rDNA_repeatSeqList_genomeList[seqIndex].push(igs_seq);
				//rDNA_last_IGS_genomeList[seqIndex] = igs_seq;
				console.log("IGS", seqName, seqIndex, igs_seq.length);
			}
			// else {
			// 	rDNA_last_IGS_genomeList[seqIndex] = null;//padding null
			// }
		}

		ex_range[seqIndex] = [start + 1, end];
		let mpos_seq = rd.chrSeq.slice(start + 1, end);
		
		return mpos_seq;
	});

	let max_num_repeats = Math.max(...rDNA_repeatSeqList_genomeList.map(repeats => repeats.length));
	//let max_num_repeats = Math.max(...rDNA_repeatSeqList_genomeList.map(repeats => repeats.length)) + 1;// + last_IGS
	let multiAlign_groups = [...Array(max_num_repeats)].map((_, repIdx) => {
		return rDNA_repeatSeqList_genomeList.map((repeats_list, seqIdx) => {
			return repeats_list[repIdx];
		});
	});
	let no_align_group_list = [];
	///multiAlign_groups.push(rDNA_last_IGS_genomeList);
	let multiAlign_groups_fileNameList = multiAlign_groups.map((group, grpIdx) => {
		let grp_fa = group.reduce((grp_fa, seq, seqIdx) => {
			if (seq) {
				const seq_name = chr_list[seqIdx];
				let out_name;

				if (!rdna_data[seq_name].info.repeats[grpIdx]) {
					// console.error(seq_name, grpIdx, rdna_data[seq_name].info.repeats.length, seq.length);
					out_name = `${seq_name} IGS-${grpIdx} len=${seq.length}`;
					console.info(`rDNA-${grpIdx}-IGS`, grpIdx, out_name);
				}
				else {
					const { [0]: start, [1]: end } = rdna_data[seq_name].info.repeats[grpIdx];
					out_name = `${seq_name} ${start}-${end}=${Math.abs(end - start)} len=${seq.length}`;
					console.info("rDNA", grpIdx, out_name);
				}

				grp_fa[out_name] = seq;
			}
			return grp_fa;
		}, {});

		if (group.filter(a => a).length >= 2) {
			//console.info("grpIdx", (grpIdx + 1));

			let fileName = `rDNA_repeat_${(grpIdx + 1)}.fa`;
			saveFasta(`${dataset.tmp_path}/${fileName}`, grp_fa);
			return fileName;
		}
		else {
			no_align_group_list.push(grp_fa);
			if (VERBOSE) {
				console.info("1 seq, grpIdx:", (grpIdx + 1));
			}
			return null;//skip
		}
	});
	
	if (VERBOSE) {
		const _tmp = Object.assign({}, rdna_data, {
			ma_rdna_range_list,
			min_ma_start,
			max_ma_end,
			ma_start_delta,
			ref1_rDNA_strand,
			ex_ranges,
		});
		fs.writeFileSync(`${dataset.output_path}/_rDNA_debug_data.json`, JSON.stringify(_tmp, function (key, value) {
			if (key == "chrSeq") {
				return undefined;
			}
			else {
				return value;
			}
		}, "\t"));
		console.log("output debug data:", `${dataset.output_path}/_rDNA_debug_data.json`);
	}

	const reAlign = true;
	await multiAlign_groups_fileNameList.filter(a => a).reduce(async function (acc, file_name, grpIdx) {
		await acc;
		let output_path = `${dataset.tmp_path}/mafft_${file_name}`;
		try {
			return await run_mafft(`${dataset.tmp_path}/${file_name}`, output_path, dataset.mafft.algorithm, nChr, `rDNA_${(grpIdx + 1)}`, reAlign);
		}
		catch (ex) {
			console.error(ex);
			return await run_mafft(`${dataset.tmp_path}/${file_name}`, `${dataset.tmp_path}/mafft_${file_name}`, dataset.mafft.default_algorithm, nChr, `rDNA_${(grpIdx + 1)}`, reAlign);
		}
	}, Promise.resolve(true));

	//let max_total_length = Math.max(...rdna_seq_list.map(a => a.length));//raw max length

	/** @type {string[]} */
	let newSeq_genomeList = Array(chr_list.length).fill("");
	multiAlign_groups_fileNameList.filter(a => a).forEach(function (new_fileName) {
		let fa = readFasta(`${dataset.tmp_path}/mafft_${new_fileName}`);

		let mmm = Math.max(...chr_list.map(seqName => fa[seqName] ? fa[seqName].length : 0));
		
		chr_list.forEach(function (seqName, seqIndex) {
			if (fa[seqName]) {
				if (ref1_rDNA_strand == 1) {
					throw new Error("WIP");
				}
				else if (ref1_rDNA_strand == -1) {
					//console.log(`rdna_seq_list[${seqIndex}].length => ${rdna_seq_list[seqIndex].length}`);// re-align seq length
					newSeq_genomeList[seqIndex] += fa[seqName].padStart(mmm, "-");
				}
			}
			else {
				newSeq_genomeList[seqIndex] += "".padStart(mmm, "-");//fill -
			}
		});
	});
	no_align_group_list.forEach(no_align_rdna => {
		Object.keys(no_align_rdna).forEach(sName => {
			let seq = no_align_rdna[sName];
			if (seq) {
				chr_list.forEach((seqName, seqIndex) => {
					console.log({
						"no_align_group_list": "no_align_group_list",
						sName,
						seqName,
						"seq.length": seq.length,
					});
					if (sName.indexOf(seqName) >= 0) {//202000721
						newSeq_genomeList[seqIndex] += seq;
					}
					else {
						newSeq_genomeList[seqIndex] += "".padEnd(seq.length, "-");
					}
				});
			}
		});
	});

	let max_total_length = 0;//init
	chr_list.forEach(function (seqName, seqIndex) {
		const ma_seq = mafft_fa[seqName];
		//const rd = rdna_data[seqName];

		//console.log("rdna_seq_list[" + seqIndex + "].length", rdna_seq_list[seqIndex].length);// raw seq length
		//let new_seq = rdna_seq_list[seqIndex].padEnd(max_total_length, "-");// raw seq

		let new_seq = newSeq_genomeList[seqIndex];//re-align seq
		max_total_length = Math.max(max_total_length, new_seq.length);

		final_mafft_seq[seqName] = replace_seq_range(ma_seq, min_ma_start, max_ma_end, new_seq);

		if (VERBOSE) {
			console.log(seqName, new_seq.slice(-3));
		}
	});
	
	// console.log("rDNA multi align range", {
	// 	start: min_ma_start,
	// 	end: min_ma_start + max_total_length,
	// 	before_end: max_ma_end,
	// });

	let rDNA_info = {
		chr: nChr,
		
		//alignment index to position
		alignment_start: min_ma_start + 1,
		alignment_end: min_ma_start + max_total_length + 1,
		alignment_range_list: ma_rdna_range_list.map(a => [a.ma_start + 1, a.ma_end + 1]),

		alignment_delta_list: ma_start_delta,

		data: chr_list.map(function (seqName, seqIndex) {
			let info = rdna_data[seqName].info;

			//rebuild pos map
			{
				let posmap = __chrPos_to_multialign_posMap(final_mafft_seq[seqName]);

				info.alignment_repeats = [];
				for (let i = 0; i < info.repeats.length; ++i) {
					let start = info.repeats[i][0] - 1;
					let end = info.repeats[i][1] - 1;
					let ma_start = posmap[start - 1] + 1;//pos(bp)
					let ma_end = posmap[end - 1] + 1;//pos(bp)

					info.alignment_repeats[i] = [];
					info.alignment_repeats[i][0] = ma_start;
					info.alignment_repeats[i][1] = ma_end;
				}

				let ex_s = posmap[ex_range[seqIndex][0]];
				let ex_e = posmap[ex_range[seqIndex][1]];

				info.alignment_ex = [ex_s, ex_e];
			}

			return info;
		}),
	};
	fs.writeFileSync(`${dataset.output_path}/rDNA_info.json`, JSON.stringify(rDNA_info));

	if (DEBUG) {
		console.log("debug no output file");
	}
	else {
		const output_path = `${dataset.output_path}/mafft_ch${nChr}.fa`;

		console.log("output:", output_path);

		saveFasta(output_path, final_mafft_seq);
		
		const output_tmp_path = `${dataset.tmp_path}/re-align_ch${nChr}.fa`;
		saveFasta(output_tmp_path, final_mafft_seq);
	}
}

class rDNA_Data {
	constructor() {
		this.strand = 0;
		this.start = 0;
		this.end = 0;
		this.chrSeq = "";
		
		/** @type {{ range: number[], repeats: number[][], alignment_repeats: number[][] }} */
		this.info = null;
	}
}

/**
 * @param {string} rDNA_filePath
 * @param {number} genomeIndex
 * @param {number} nChr
 * @returns {Promise<rDNA_Data>}
 */
async function find_rDNA_use_blastn(rDNA_filePath, genomeIndex, nChr) {
	const subject_genome_name = dataset.genomeNameList[genomeIndex];
	const subject_chrInfo = genome_info_list[genomeIndex].chr_list;
	const subject_chr_name = subject_chrInfo[nChr - 1].chr;

	const subject_fa_filename = subject_chrInfo[nChr - 1].path;

	const result_text = await exec_blastn_Ex(rDNA_filePath, subject_fa_filename, undefined, undefined, undefined, undefined, "-evalue 1e-5");
	const _table = parseBlastnResults(result_text);
	
	const max_len = Math.max(..._table.map(a => a.slen));
	
	/**
	 * TODO: Check Ribosomal DNA structure, IGS 18S ITS 5.8S ITS 26S IGS
	 * CBS1-1 min:7364, max:7835, 7835 / 7364 = 0.94
	 * range include IGS
	 */
	const group = _table.filter(a => (a.slen / max_len) >= 0.9).sort((a, b) => a.s_min - b.s_max);
	if (!group.every(a => a.strand == group[0].strand)) {
		console.warn("found reverse rDNA");
	}
	
	const near_bp = 1;
	const region_list = _table.filter(aln => {
		const [ts, te] = [aln.sstart, aln.send].sort((a, b) => a - b);
		return group.some(rep => {
			const [rs, re] = [rep.sstart, rep.send].sort((a, b) => a - b);
			if (ts <= (re + near_bp) &&
				te >= (rs - near_bp)
			) {
				return true;
			}
		});
	})
	
	if (!DEBUG) {
		fs.writeFileSync(
			`${dataset.tmp_path}/blastn_rDNA_${subject_chr_name}.txt`,
			group.map(a => a.toArray().join("\t")).join("\n")
		);
	}

	let min_sstart = Math.min(...region_list.map(a => a.s_min));
	let max_send = Math.max(...region_list.map(a => a.s_max)) + 1;

	if (VERBOSE) {
		console.log({
			subject: subject_genome_name,
			min_sstart, max_send,
			len: max_send - min_sstart,
			"_table.length": _table.length,
			"group.length": group.length,
		});
	}

	if (!(max_send > min_sstart)) {
		throw new Error("if (!(max_send > min_sstart)) {");
	}
	else {
		let raw_fa = subject_chrInfo[nChr - 1].loadSeq();

		let info = {// save rDNA repeat position
			range: [min_sstart, max_send],
			repeats: group.map(a => [a.sstart, a.send]),
			alignment_repeats: null,
		};
		return {
			strand: group[0].strand,
			chrSeq: raw_fa,
			//seq: getSeq(raw_fa, min_sstart, max_send),
			start: min_sstart, // rDNA region start
			end: max_send,     // rDNA region end
			info,
		};
	}
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function _chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 0;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index;
			++nPos;
		}
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function _multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 0;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		posmap[index] = nPos;

		if (element != "-") {
			++nPos;
		}
	}
	return posmap;
}



/**
 * pos: 1 ~ length
 * @param {string} ma_seq
 * @returns {number[]}
 */
function __chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index;
			++nPos;
		}
	}
	return posmap;
}

/**
 * pos: 1 ~ length
 * @param {string} ma_seq
 * @returns {number[]}
 */
function __multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		if (element == "-") {
			posmap[index] = Math.max(1, nPos - 1);
		}
		else {
			posmap[index] = nPos;
			++nPos;
		}
	}
	return posmap;
}

/**
 * @param {string} ma_seq
 * @param {number} ma_start
 * @param {number} ma_end
 * @param {string} replace_seq
 */
function replace_seq_range(ma_seq, ma_start, ma_end, replace_seq) {
	return getSeq(ma_seq, 1, ma_start) + replace_seq + getSeq(ma_seq, ma_end, ma_seq.length);
}

/**
 * @param {string} seq
 * @param {number} start - bp
 * @param {number} end - bp
 */
function getSeq(seq, start, end) {
	return seq.slice(start - 1, end);
}



