//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { BlastnCoord, execAsync, exec_blastn, exec_blastn_Ex, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { run_mafft } = require("./run_mafft.js");
const { readFasta, saveFasta, chrPos_to_multialign_posMap } = require("./fasta_util.js");
const { validation_chr } = require("./validation_seq.js");
const { Dataset, RibosomalDNA_Data } = require("./dataset.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");
const { join_chr_frag } = require("./join_chr_frag.js");

const argv = argv_parse(process.argv);

const DEBUG = !!argv["--debug"];
const VERBOSE = process.argv.indexOf("--verbose") >= 0;
const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const genome_info_list = dataset.loadGenomeInfoList();

main();

// const region_list = [
// 	{
// 		nChr: 1,
// 		frag_range: [1, 3],
// 	}
// ]

async function main() {
	program_log(`${dataset.name}.log.txt`, "start");

	await re_align(1, 1, 5);
}

/**
 * @param {number} nChr
 * @param {number} fragId_start
 * @param {number} fragId_end
 */
async function re_align(nChr, fragId_start, fragId_end) {
	const all_chr_frag_list = loadFragIdList(dataset);

	const chrIdx = nChr - 1;
	
	let chr_list = genome_info_list.map(genome_info => genome_info.chr_list[chrIdx].chr);
	
	/** @type {number[]} */
	const raw_frag_idx_list = [];
	const raw_frags = all_chr_frag_list[nChr].filter((coord, coord_idx) => {
		if (Number(coord.id) >= fragId_start && Number(coord.id) <= fragId_end) {
			raw_frag_idx_list.push(coord_idx);
			return true;
		}
	});

	// load frags

	/** @type {{ [chr: string]: string }} */
	const concat_raw_seq = {};
	//
	raw_frags.forEach((coord, frag_idx) => {
		const fragId = coord.id;
		const fasta_filename = `ch${nChr}_${fragId}.fa`;
		const input_path = `${dataset.tmp_path}/seq_frag/${fasta_filename}`;
		const fa = readFasta(input_path);
		
		// const g = ["r1", "r2", "s1", "s2", "s3", "s4"];
		// const start_list = g.map(k => raw_frags[0].start[k]);
		// const end_list = g.map(k => raw_frags[raw_frags.length - 1].end[k]);

		chr_list.forEach(function (genomeChrSeqName, genomeIndex) {
			if (!concat_raw_seq[genomeChrSeqName]) {
				concat_raw_seq[genomeChrSeqName] = "";
				// concat_raw_seq_aa[genomeChrSeqName] = [];
			}
			concat_raw_seq[genomeChrSeqName] += fa[genomeChrSeqName];
			// concat_raw_seq_aa[genomeChrSeqName] = concat_raw_seq_aa[genomeChrSeqName].concat(...)
		});
	});
	const frag_id = `${raw_frags[0].id}-${raw_frags[raw_frags.length - 1].id}`;
	const merged_raw_fa_name = `ch${nChr}_${frag_id}.fa`;
	const merged_raw_fa_path = `${dataset.tmp_path}/merged_fa/${merged_raw_fa_name}`;

	saveFasta(merged_raw_fa_path, concat_raw_seq);

	// run mafft

	const merged_fa_output_path = `${dataset.tmp_path}/mafft_seq_frag/mafft_${merged_raw_fa_name}`;
	
	const reAlign = true;
	await run_mafft(
		merged_raw_fa_path,
		merged_fa_output_path,
		dataset.mafft.algorithm,
		nChr,
		frag_id,
		reAlign
	);

	// c f
	
	/**
	 * @type {MyCoord[]}
	 */
	const new_frags = [
	];

	new_frags[0] = new MyCoord();
	new_frags[0].id = frag_id;

	// modify list

	const new_frag_list = all_chr_frag_list[nChr].slice(0);// clone

	if (VERBOSE) {
		console.log("o frags", all_chr_frag_list[nChr].length);
	}

	// modify
	new_frag_list.splice(raw_frag_idx_list[0], raw_frag_idx_list.length, ...new_frags);

	if (VERBOSE) {
		console.log("n frags", new_frag_list.length);

		console.log("r frags", new_frags.length);
	}
	
	const output_tmp_path = `${dataset.output_path}/mafft_ch${nChr}.fa`;
	await fs.promises.copyFile(output_tmp_path, `${dataset.output_path}/${new Date().getTime()}_mafft_ch${nChr}.fa`);

	join_chr_frag(dataset, nChr, new_frag_list, output_tmp_path, { padding: true, chr_list: chr_list });

	validation_chr(nChr, dataset.output_path, true);
}

// 20200727
class rDNA_Data {
	constructor() {
		/**
		 * nChr -> [1,n]
		 */
		this.nChr = 0;

		this.strand = 0;

		/**
		 * @typedef seq_range
		 * @type {[number, number] | number[]}
		 */

		/**
		 * region
		 * @type {seq_range}
		 */
		this.region = [0, 0];

		/**
		 * alignment region
		 * @type {seq_range}
		 */
		this.alignment_region = [0, 0];

		/**
		 * repeat
		 * @type {seq_range[]}
		 */
		this.repeats = [];
		
		/**
		 * alignment repeat
		 * @type {seq_range[]}
		 */
		this.alignment_repeats = [];
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
	
	// const plus_score = group.filter(a => a.strand > 0).reduce((t, v) => t + v.strand * v.score, 0);
	// const minus_score = group.filter(a => a.strand < 0).reduce((t, v) => t + v.strand * v.score, 0);
	
	const strand = (function detect_seq_strand() {
		const strand_score = group.reduce((t, v) => t + (v.strand * v.score), 0);
		return strand_score > 0 ? 1 : (strand_score < 0 ? -1 : 0);
	})();

	if (strand == 0) {
		console.warn("Unknow rDNA strand");
	}

	// if (!group.every(a => a.strand == group[0].strand)) {
	// 	console.warn("found reverse rDNA");
	// }
	
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
		// let raw_fa = subject_chrInfo[nChr - 1].loadSeq();
		// let info = {// save rDNA repeat position
		// 	region: [min_sstart, max_send],
		// 	repeats: group.map(a => [a.sstart, a.send]),
		// 	alignment_repeats: null,
		// };
		// return {
		// 	strand: strand,
		// 	// chrSeq: raw_fa,
		// 	//seq: getSeq(raw_fa, min_sstart, max_send),
		// 	// start: min_sstart, // rDNA region start
		// 	// end: max_send,     // rDNA region end
		// 	info,
		// };

		const data = new rDNA_Data();
		data.region = [min_sstart, max_send];
		data.repeats = group.map(a => [a.sstart, a.send]);
		data.strand = strand;
		data.nChr = nChr;
		return data;
	}
}

