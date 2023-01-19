
// @ts-check

if (typeof String.prototype.matchAll != "function") {
	const matchAll = require("string.prototype.matchall");
	matchAll.shim();
}

const fs = require("fs");
const child_process = require("child_process");

const { argv_parse } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser");
const { Dataset } = require("./dataset.js");
const { Crossover } = require("./crossover_util.js");
const { readFasta, multialign_to_chrPos_posMap, saveFasta } = require("./fasta_util.js");
const { SegRow } = require("./SegFile.js");

const { FinalTableRow } = require("./final_table_format.js");

const { initAnalyser, loadCmpData } = require("./analyser.js");


const { parseGFF, GFF_ROW } = require("./gff.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_output_chr = Number(argv["-chr"]) | 0;
const argv_seq = String(argv["--seq"] || "");
const argv_co_list = String(argv["--co-list"] || "");
const argv_nco_list = String(argv["--nco-list"]);
const argv_output_prefix = String(argv["--output-prefix"] || "");

const argv_snp_count_window = Number(argv["--snp-count-window"]) | 0;
const argv_point_filter = argv["--point_filter"] == "true" || argv["--point_filter"];


const merge_RNA_TPM_repeat_methyl = false;


if (!argv_dataset_path || !argv_seq || !argv_output_chr || !argv_output_prefix || !argv_co_list) {
	console.log("Usage:", "node tetrad_chr_summary.js -dataset <dataset.json> -chr <nChr> --seq <multi_align.fa> --co-list <chrN_co_list.json> --output-prefix <output file prefix name>");
	console.log({
		argv_dataset_path,
		argv_seq,
		argv_output_chr,
		argv_output_prefix,
	})
	throw new Error("agv");
}

const co_list = ((co_list_filename) => {
	const text = fs.readFileSync(co_list_filename).toString();
	const list = JSON.parse(text);
	/**
	 * @param {{ before: string; after: string; }} co
	 */
	list.forEach(co => {
		/**
		 * @param {any} n
		 */
		co.before = co.before.split(",").map(n => Number(n));
		/**
		 * @param {any} n
		 */
		co.after = co.after.split(",").map(n => Number(n));
	});
	return list;
})(argv_co_list);

const nco_list = ((co_list_filename) => {
	const text = fs.readFileSync(co_list_filename).toString();
	const list = JSON.parse(text);
	// list.forEach(co => {
	// 	co.before = co.before.split(",").map(n => Number(n));
	// 	co.after = co.after.split(",").map(n => Number(n));
	// });
	return list;
})(argv_nco_list);

const input_fasta = readFasta(argv_seq);
const dataset = Dataset.loadFromFile(argv_dataset_path);
dataset.load_GC_content();
dataset.load_rDNA_info();

const genome_info_list = dataset.loadGenomeInfoList();

function main() {
	console.log({
		ref: dataset.ref,
		chr: argv_output_chr
	});
	const nChr = argv_output_chr;

	/** @type {string[]} */
	let seq_list = [];
	let seq_id_list = Object.keys(input_fasta);
	seq_id_list.forEach((id, i) => seq_list[i] = input_fasta[id]);

	// // const gff_sIdx = 0;// QM6a
	// // const gff = load_gff(gff_sIdx, `${dataset.output_path}/QM6a.gff3`, nChr);

	// // const gff_sIdx = 1;// CBS
	// // const gff = load_gff(gff_sIdx, `${dataset.output_path}/Trichoderma_reesei_CBS.gff3`, nChr);

	// const gff_files = [
	// 	// `${dataset.output_path}/QM6a.gff3`,
	// 	// `${dataset.output_path}/Trichoderma_reesei_CBS.gff3`,
	// ];
	// const gff = gff_files.length > 0 ? dataset.parental_list.reduce((obj, gff_ref, gff_idx) => {
	// 	if (gff_files[gff_idx]) {
	// 		obj[gff_ref] = load_gff(gff_idx, gff_files[gff_idx], nChr);
	// 	}
	// 	return obj;
	// }, {}) : null;

	// /**
	//  * @type {{ [refId:string]: RepeatSegment[][] }}
	//  */
	// const repeat_segment = {};
	// if (merge_RNA_TPM_repeat_methyl) {
	// 	dataset.parental_list.forEach((ref_name, ref_idx) => {
	// 		if (!repeat_segment[ref_name]) {
	// 			repeat_segment[ref_name] = [];
	// 		}
	// 		const rows = LoadRepeatSegment(ref_name, nChr);
	// 		// console.log({
	// 		// 	refId,
	// 		// 	nChr: nChr,
	// 		// 	"rows.length": rows.length,
	// 		// });
	// 		repeat_segment[ref_name].push(rows);
	// 	});
	// 	// console.log(
	// 	// 	typeof repeat_segment,
	// 	// 	Object.keys(repeat_segment),
	// 	// 	typeof repeat_segment[Object.keys(repeat_segment)[0]],
	// 	// 	repeat_segment[Object.keys(repeat_segment)[0]].length,
	// 	// 	);
	// }

	// const methratio = merge_RNA_TPM_repeat_methyl ? {
	// 	[dataset.parental_list[0]]: [
	// 		LoadMethylRatio(0, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_Q_wt.methratio.float32`),
	// 		LoadMethylRatio(0, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_Q_d4.methratio.float32`),
	// 		LoadMethylRatio(0, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_Q_d8.methratio.float32`),
	// 	],
	// 	[dataset.parental_list[1]]: [
	// 		LoadMethylRatio(1, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_C_wt.methratio.float32`),
	// 		LoadMethylRatio(1, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_C_d4.methratio.float32`),
	// 		LoadMethylRatio(1, nChr, `${dataset.output_path}/20200902_wt_methyl/20200902_C_d8.methratio.float32`),
	// 	],
	// } : null;

	// /**
	//  * @param {number[]} _list
	//  */
	// function average_stdev(_list) {
	// 	const list = _list.filter(a => a != null && Number.isFinite(a) && !Number.isNaN(a));
	// 	if (list.length <= 0) {
	// 		return "";
	// 	}
	// 	if (list.length == 1) {
	// 		return list[0].toFixed(2);
	// 	}

	// 	const average = list.reduce((acc, v) => acc + v, 0) / list.length;
	// 	const stdev = (list.reduce((acc, v) => acc + ((average - v) ** 2), 0) / (list.length - 1)) ** 0.5;
	// 	// return {
	// 	// 	average,
	// 	// 	stdev,
	// 	// };
	// 	return `"${average.toFixed(2)}\n${stdev.toFixed(2)}"`;
	// }

	// const _RNA_sample_name_list = merge_RNA_TPM_repeat_methyl ? [
	// 	["QM6a-1","QM6a-2","QM6a-3",],
	// 	["CBS_1-1-1","CBS_1-1-2","CBS_1-1-3",],
	// 	["WT-D1-1","WT-D1-2","WT-D1-3",],
	// 	["WT-D2-1","WT-D2-2","WT-D2-3",],
	// 	["WT-D3-1","WT-D3-2","WT-D3-3",],
	// 	["WT-D4-1","WT-D4-2","WT-D4-3",],
	// 	["WT-D5-1","WT-D5-2","WT-D5-3",],
	// 	["WT-D6-1","WT-D6-2","WT-D6-3",],
	// 	["WT-D7-1","WT-D7-2","WT-D7-3",],
	// 	["WT-D8-1","WT-D8-2","WT-D8-3",],
	// 	["D1"],
	// 	["D2"],
	// 	["D3"],
	// 	["D4"],
	// 	["D5"],

	// 	// // rid1
	// 	// ["rid-Q-1","rid-Q-2","rid-Q-3",],
	// 	// ["rid1_1-1-1","rid1_1-1-2","rid1_1-1-3",],
	// 	// ["rid-D1-1","rid-D1-2","rid-D1-3",],
	// 	// ["rid-D2-1","rid-D2-2","rid-D2-3",],
	// 	// ["rid-D3-1","rid-D3-2","rid-D3-3",],
	// 	// ["rid-D4-1","rid-D4-2","rid-D4-3",],
	// 	// ["rid-D5-1","rid-D5-2","rid-D5-3",],
	// 	// ["rid-D6-1","rid-D6-2","rid-D6-3",],
	// 	// ["rid-D7-1","rid-D7-2","rid-D7-3",],
	// 	// ["rid-D8-1","rid-D8-2","rid-D8-3",],
	// ] : null;

	// const  RNA_sample_name_list = merge_RNA_TPM_repeat_methyl ? [
	// 	"QM6a",
	// 	"CBS_1-1",
	// 	"WT-D1",
	// 	"WT-D2",
	// 	"WT-D3",
	// 	"WT-D4",
	// 	"WT-D5",
	// 	"WT-D6",
	// 	"WT-D7",
	// 	"WT-D8",
	// 	"D1","D2","D3","D4","D5",
	// ] : null;
	
	// /**
	//  * @param {string} file_path
	//  * @returns {Map<string, { [sampleName:string]:{ stdev:number; average:number; } }>}
	//  */
	// function load_TPM_table(file_path) {
	// 	const text = fs.readFileSync(file_path).toString();
		
	// 	// remove rid*
	// 	const rows = tsv_parse(text);//.map(row => row.slice(0, -30));
		
	// 	const table = table_to_object_list(rows, 0, { header_map: function (rawHead) {
	// 		if (rawHead == "Geneid") {
	// 			return "geneID";
	// 		}
	// 		else {
	// 			return rawHead.slice(0, -("_Aligned.out.sam".length));
	// 		}
	// 	}});

	// 	const map = new Map();

	// 	table.forEach(row => {
	// 		/** @type {{ [sampleName:string]:{ stdev:number; average:number; } }} */
	// 		const results = {};
	// 		_RNA_sample_name_list.map((key_list, idx) => {
	// 			const s_name = RNA_sample_name_list[idx];
	// 			const val_list = key_list.map(key => Number(row[key]));
	// 			results[s_name] = average_stdev(val_list);
	// 		});
	// 		map.set(row.geneID, results);
	// 	});

	// 	return map;
	// }

	// const RNA_TPM = merge_RNA_TPM_repeat_methyl ? {
	// 	[dataset.parental_list[0]]: load_TPM_table(`${dataset.output_path}/TPM/20200909_QM6a_mix_CC_TPM.txt`),
	// 	[dataset.parental_list[1]]: load_TPM_table(`${dataset.output_path}/TPM/20200909_CBS1-1_TPM.txt`),
	// } : null;
	
	// /**
	//  * @param {RegExpMatchArray} m
	//  */
	// function sss_mapto_segment(m) {
	// 	return {
	// 		start: m.index,
	// 		end: m.index + m[0].length,
	// 	};
	// }
	
	// // const sss_50 = get_strain_specific_region_by_length(seq_list, 50);//{ qqq_sss_t_len, ccc_sss_t_len, qqq_sss, ccc_sss }
	// // const sss_100 = get_strain_specific_region_by_length(seq_list, 100);
	// const sss_60 = get_strain_specific_region_by_length(seq_list, 60);
	// const strain_specific = {
	// 	// preset: ["50", "100"],
	// 	profiles: ["60"],
	// 	[dataset.parental_list[0]]: {
	// 		// "50": sss_50.qqq_sss,//.map(sss_mapto_segment),
	// 		// "100": sss_100.qqq_sss,//.map(sss_mapto_segment),
	// 		"60": sss_60.qqq_sss,
	// 	},
	// 	[dataset.parental_list[1]]: {
	// 		// "50": sss_50.ccc_sss,//.map(sss_mapto_segment),
	// 		// "100": sss_100.ccc_sss,//.map(sss_mapto_segment),
	// 		"60": sss_60.ccc_sss,
	// 	},
	// };

	let analysis_options = {
		get mode() { return dataset.mode; },
		get nChr() { return nChr },
		
		get point_filter() { return argv_point_filter; },
		get telomere() { return dataset.telomere[nChr]; },

		get rDNA_info() { return dataset.rDNA_info; },
		// get rDNA_info() {
		// 	return {
		// 		chr: dataset.rDNA.nChr,
		// 		region_start: dataset.rDNA.region[0],
		// 		region_end: dataset.rDNA.region[1],
		// 	};
		// },

		// get show_rDNA_snp() { return false; },
		get co_list() { return co_list; },//in this stage
		get nco_list() { return nco_list; },//in this stage
		get fill_prev_color() { return true; },

		/** @type {ChromosomeData[]} - gemome[*].chrInfo */
		get chrInfo_list() { return genome_info_list.map(gInfo => gInfo.chr_list[analysis_options.nChr - 1]); },
		
		get show_rDNA_snp() { return false; },
		get show_rDNA_non_InDel() { return true; },

		// // 20200818 // get find_rip() { return true; },

		// rip_gff: false && gff_files.length > 0,
		// rip_repeat_methyl_gff: true && gff_files.length > 0,

		// repeat_segment: merge_RNA_TPM_repeat_methyl ? repeat_segment : null,
		// // strain_specific: strain_specific,
		
		// methratio: merge_RNA_TPM_repeat_methyl ? methratio : null,
		// methyl_left: merge_RNA_TPM_repeat_methyl ? 10 : null,
		// methyl_right: merge_RNA_TPM_repeat_methyl ? 10 : null,

		// RNA_TPM: merge_RNA_TPM_repeat_methyl ? RNA_TPM : null,

		// // gff_sIdx: gff_sIdx,
		// gff: merge_RNA_TPM_repeat_methyl ? gff : null,
		// gene_rip_color: merge_RNA_TPM_repeat_methyl ? {
		// 	"gene": "lime",
		// 	"mRNA": "orange",
		// 	"CDS": "purble",
		// 	"exon": "red",
		// } : null,

		ref_snp_list: true,
	};

	initAnalyser(analysis_options);
	let data = loadCmpData(seq_list, analysis_options);

	// analyser_gff_methyl(gff, data, analysis_options, seq_list);

	if (nco_list && argv_nco_list) {
		// fs.writeFileSync("plot_nco_ch" + (analysis_options.nChr) + ".txt", "");// clear

		let list = nco_list.splice(0);
		/**
		 * @param {{ rip_snv22_markers: any; maybe_rip_snv40_markers: any; type: string; }} nco
		 */
		list.forEach(nco => {
			delete nco.rip_snv22_markers;//20200629
			delete nco.maybe_rip_snv40_markers;;//20200629

			if (nco.type != "2NCO") {
				return true;
			}
			const _2nco = nco;
			// let output_plot = false;

			// if (nco.rip_snv22_markers) {
			// 	nco.why_remove = nco.rip_snv22_markers.length + " rip* snv22";
			// 	output_plot = true;
			// }

			// if (!nco.is_rip && nco.maybe_rip_snv40_markers) {
			// 	nco.why_remove = nco.maybe_rip_snv40_markers.length + " RIP^ marker";
			// 	// throw new Error(nco.why_remove);
			// 	output_plot = true;
			// }

			// if (nco.n_rip_snv40) {
			// 	nco.why_remove = nco.n_rip_snv40 + " rip* snv40";
			// 	output_plot = true;
			// }

			// delete _2nco.rip_snv22_markers;//20200629
			// delete _2nco.maybe_rip_snv40_markers;//20200629
			if (_2nco.is_rip) {
				_2nco.why_remove = "all marker are RIP";
			}

			// if (output_plot) {
			// 	plot_nco(nco);
			// }

			// if (_2nco.is_rip) {
			// 	return false;
			// }
			// else {
			// 	return true;
			// }
		});
		nco_list.push(...list);
		fs.writeFile(argv_nco_list, JSON.stringify(nco_list), function (err) {
			if (err) {
				console.error(err);
			}
		});
	}
	/**
	 * @param {{ snp_start_out: any; snp_start_in: any; snp_end_in: any; snp_end_out: any; why_remove: any; }} nco
	 */
	function plot_nco(nco) {
		const snp_start_out = data.ref1_pos_uint32array[Number(nco.snp_start_out)];
		const snp_start_in = data.ref1_pos_uint32array[Number(nco.snp_start_in)];
		const snp_end_in = data.ref1_pos_uint32array[Number(nco.snp_end_in)];
		const snp_end_out = data.ref1_pos_uint32array[Number(nco.snp_end_out)];
		
		const m1 = snp_start_in - snp_start_out - 1;
		const m2 = snp_end_in - snp_start_out - (m1 + 1) - 1;
		
		const local_seq = seq_list.map(sss => sss.slice(snp_start_out, snp_end_out)).join("\n");
		const marker = (m1 >= 0 ? " ".repeat(m1) : "") + "^" + (m2 >= 0 ? " ".repeat(m2) : "") + "^";
		const out_text = `${nco.why_remove} [${snp_start_out}, ${snp_start_in}, ${snp_end_in}, ${snp_end_out}]\n${local_seq}\n${marker}\n\n`;

		fs.writeFileSync("plot_nco_ch" + analysis_options.nChr + ".txt", out_text, { flag: "a" });
	}
	
	const output_path = `${dataset.tmp_path}/table`;
	if (!fs.existsSync(output_path)) {
		fs.mkdirSync(output_path, { recursive: true, });
	}

	const ref_snp = data.ref_snp_list;
	// const ref_snp = [...seq_list[0]].map((q, i) => [i + 1, q, seq_list[1][i]]).filter(([pos, q, c]) => !point_filter(Number(pos)) && q != c);
	const ref_snv = ref_snp.filter(([pos, q, c]) => q != c && q != "-" && c != "-");
	const ref_snp_indel = ref_snp.filter(([pos, q, c]) => q != c && ((q == "-" && c != "-") || (q != "-" && c == "-")));//q is del or c is del
	// const ref_snp = [...seq_list[0]].map((q, i) => [data.pos_ref1_uint32array[i], q, seq_list[1][i], data.pos_ref2_uint32array[i], i]).filter(([ref1_pos, q, c, ref2_pos, i]) => !point_filteri) && q != c);
	// const ref_snv = ref_snp.filter(([ref1_pos, q, c, ref2_pos, i]) => q != c && q != "-" && c != "-");
	// const ref_snp_indel = ref_snp.filter(([ref1_pos, q, c, ref2_pos, i]) => q != c && ((q == "-" && c != "-") || (q != "-" && c == "-")));//q is del or c is del

	// if (0) {
	// 	ref_snp = [...seq_list[0]].map((q, i) => q != seq_list[1][i] ? [i, q, seq_list[1][i]] : null).filter(a => a);
	// 	ref_snv = ref_snp.filter(([i, q, c]) => q != "-" && c != "-");
	// 	ref_snp_indel = ref_snp.filter(([i, q, c]) => (q == "-" && c != "-") || (q != "-" && c == "-"));//q is del or c is del
	// }

	fs.writeFileSync(`${output_path}/${argv_output_prefix}_snp.json`, JSON.stringify(ref_snp));
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_snv.json`, JSON.stringify(ref_snv));
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_snp_indel.json`, JSON.stringify(ref_snp_indel));

	// const mut_ref1 = data.spore_cmp.rip.values.filter(rip => rip.mut_ref == 1);
	// const mut_ref2 = data.spore_cmp.rip.values.filter(rip => rip.mut_ref == 2);
	const mut_ref1 = data.spore_cmp.rip_Q.values;
	const mut_ref2 = data.spore_cmp.rip_C.values;
	const rip_QC = data.spore_cmp.rip_QC.values;
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_ref1_rip.json`, JSON.stringify(mut_ref1));
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_ref2_rip.json`, JSON.stringify(mut_ref2));
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_QC_rip.json`, JSON.stringify(rip_QC));

	// const rip_2_ref1 = data.spore_cmp.rip_2.values.filter(rip => rip.mut_ref == 1);
	// const rip_2_ref2 = data.spore_cmp.rip_2.values.filter(rip => rip.mut_ref == 2);
	const rip_2_ref1 = data.spore_cmp.rip_2_Q.values;
	const rip_2_ref2 = data.spore_cmp.rip_2_C.values;
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_ref1_2_rip.json`, JSON.stringify(rip_2_ref1));
	fs.writeFileSync(`${output_path}/${argv_output_prefix}_ref2_2_rip.json`, JSON.stringify(rip_2_ref2));

	[].concat(mut_ref1, mut_ref2, rip_QC, rip_2_ref1, rip_2_ref2).forEach(rip => {
		if (rip.is_2nco) {
			debugger;
		}
	});

	{// 20200709
		// && seq_list.slice(2).filter(ss => ss[i] != "-").length == 2

		const output_map = {};

		const { qqq_sss_t_len, ccc_sss_t_len, qqq_sss, ccc_sss } = get_strain_specific_region_by_length(seq_list);

		output_map.Chromosome = nChr;
			
		output_map.QC_SNV = ref_snv.length;
		output_map.QC_InDel = ref_snp_indel.length;
		// output_map.QC_SNP = ref_snp.length;

		ref_snp_indel.join()

		output_map["QM6a specific"] = qqq_sss_t_len;
		output_map["CBS1-1 specific"] = ccc_sss_t_len;

		const s_tbl_qc_bp_len = Object.keys(output_map).map(k => `${output_map[k]}`).join("\t");
		fs.writeFile(`${output_path}/${argv_output_prefix}_qc_sss_bp.txt`, s_tbl_qc_bp_len, (err) => err ? console.error(err) : void 0);

		const s_tbl_qc_sss = [
			//                                     pos,   ch,                   strain start pos,      length
			qqq_sss.map(m => [dataset.parental_list[0], nChr, data.pos_ref1_uint32array[m.index], m[0].length].join("\t")).join("\n"),
			ccc_sss.map(m => [dataset.parental_list[1], nChr, data.pos_ref2_uint32array[m.index], m[0].length].join("\t")).join("\n")
		].join("\n");

		fs.writeFile(`${output_path}/${argv_output_prefix}_qc_sss.txt`, s_tbl_qc_sss, (err) => err ? console.error(err) : void 0);
	}

	//output table 1
	{
		// const output_head = [
		// 	"Chromosome",
		// 	"simple CO", "CO(NCO)",
		// 	"Q/C SNV", "Q/C SNP", "Q/C InDel",
		// 	"RIP Q", "RIP C",
		// 	"illegitimate mutation",
		// 	"SNV 2:2", "NCO 3:1", "NCO 4:0",
		// 	"1n:3", "2n:2", "3n:1", "4n:0",
			
		// 	"Strain-specific sequences 2:2",
		// 	"Strain-specific sequences 3:1",
		// 	"Strain-specific sequences 4:0",

		// 	"illegitimate mutation (not 3:1)",
		// 	"illegitimate mutation (not 4:0)",
		// ];

		// /**
		//  * @type {FinalTableRow}
		//  */
		// const output_map = Object.assign(new FinalTableRow(), {
		// 	"Chromosome": nChr,
		// 	"simple CO": analysis_options.co_list.filter(co => co.type == "CO").length,
		// 	"CO(NCO)": analysis_options.co_list.filter(co => co.type == "CO(NCO)").length,
			
		// 	"Q/C SNV": ref_snv.length,
		// 	"Q/C SNP": ref_snp.length,
		// 	"Q/C InDel": ref_snp_indel.length,
			
		// 	"RIP Q": mut_ref1.length,
		// 	"RIP C": mut_ref2.length,

		// 	"illegitimate mutation": data.spore_cmp.illegitimate_mutation_list.values.length,
			
		// 	"SNV 2:2": data.spore_cmp.s22.values.length,
		// 	"NCO 3:1": data.spore_cmp.s31.values.length,
		// 	"NCO 4:0": data.spore_cmp.s40.values.length,
			
		// 	"1n:3": data.spore_cmp.s1n3.values.length,
		// 	"2n:2": data.spore_cmp.s2n2.values.length,
		// 	"3n:1": data.spore_cmp.s3n1.values.length,
		// 	"4n:0": data.spore_cmp.s4n0.values.length,
		
		// 	"Strain-specific sequences 2:2": data.spore_cmp.sss_22.values.length,
		// 	"Strain-specific sequences 3:1": data.spore_cmp.sss_31.values.length,
		// 	"Strain-specific sequences 4:0": data.spore_cmp.sss_40.values.length,

		// 	"illegitimate mutation (not 3:1)": data.spore_cmp.illegitimate_mutation_31.values.length,
		// 	"illegitimate mutation (not 4:0)": data.spore_cmp.illegitimate_mutation_40.values.length,
		// });

		/** @type {FinalTableRow<number|string>} */
		const output_map = new FinalTableRow();
		
		output_map.Chromosome = nChr;
		output_map.ref1_len = genome_info_list[0].chr_list[nChr - 1].length;
		output_map.ref2_len = genome_info_list[1].chr_list[nChr - 1].length;
		
		const filtered_co = analysis_options.point_filter ? analysis_options.co_list.filter(co => !co.why_remove) : analysis_options.co_list;
		output_map.simple_CO = filtered_co.filter(co => co.type == "CO").length;
		output_map.CO_NCO =  filtered_co.filter(co => co.type == "CO(NCO)").length;

		const filtered_nco = analysis_options.point_filter ? analysis_options.nco_list.filter(co => !co.why_remove) : analysis_options.nco_list;
		output_map.NCO =  filtered_nco.length;
		output_map.NCO_2p =  filtered_nco.filter(a => a.GCasso_marker >= 2).length;
		
		output_map.NCO1 =  filtered_nco.filter(a => a.nco_type == "NCO1").length;
		output_map.NCO2 =  filtered_nco.filter(a => a.nco_type == "NCO2").length;
		output_map.NCO1_2p =  filtered_nco.filter(a => a.nco_type == "NCO1" && a.GCasso_marker >= 2).length;
		output_map.NCO2_2p =  filtered_nco.filter(a => a.nco_type == "NCO2" && a.GCasso_marker >= 2).length;

		const num_of_CO = Number(output_map.simple_CO) + Number(output_map.CO_NCO);
		const num_of_CO_and_NCO = num_of_CO + Number(output_map.NCO_2p);
		output_map.CO_div_nCO = (num_of_CO / num_of_CO_and_NCO).toFixed(2);

		output_map.QC_SNV = ref_snv.length;
		output_map.QC_SNP = ref_snp.length;
		output_map.QC_InDel = ref_snp_indel.length;
		// {
		// 	function aaa() {
		// 		return [...seq_list[0]].map((q, i) => [
		// 			data.pos_ref1_uint32array[i],
		// 			q, seq_list[1][i],
		// 			data.pos_ref2_uint32array[i],
		// 			i
		// 		]);
		// 	}
		//
		// 	function bbb() {
		// 		return aaa().map(([ref1_pos, q, c, ref2_pos, i]) => [
		// 			ref1_pos,
		// 			q, c,
		// 			ref2_pos,
		// 			i,
		// 			!point_filteri) && q != c,
		// 			!point_filteri) && q != c && ((q == "-" && c != "-") || (q != "-" && c == "-"))//q is del or c is del
		// 		]);
		// 	}
		//
		// 	function ccc(bbb) {
		// 		return bbb.map(([ref1_pos, q, c, ref2_pos, i, isSNP, isIndel]) => isIndel ? 1 : 0).join("");
		// 	}
		//
		// 	/** @returns {RegExpMatchArray[]} */
		// 	function mmm(bbb) {
		// 		return [...ccc(bbb).matchAll(/1+/g)];
		// 	}
		//
		// 	const rrr = (function () {
		// 		return mmm(bbb()).map(a => a[0].length).sort((a, b) => a - b);
		// 	})();
		//
		// 	let indel_map = {};
		// 	rrr.forEach(n_indel => indel_map[n_indel] = (indel_map[n_indel] | 0) + 1);
		// 	// fs.writeFileSync(`${output_path}/ch${argv_output_chr}_snp_map.json`, JSON.stringify(rrr));
		//
		// 	const out_text = Object.keys(indel_map).map(n_indel => [n_indel, indel_map[n_indel]]).filter(([n_indel, n]) => n > 0).map(row => row.join("\t")).join("\n");
		// 	fs.writeFileSync(`${output_path}/ch${argv_output_chr}_snp_map.txt`, `length(bp)\tn\n${out_text}`);
		// }

		output_map.SNV_22 = data.spore_cmp.s22.values.length;
		output_map.NCO_31 = data.spore_cmp.s31.values.length;
		output_map.NCO_40 = data.spore_cmp.s40.values.length;
		
		// if (dataset.options.includes("RIP")) {
			output_map.RIP_Q = mut_ref1.length;
			output_map.RIP_C = mut_ref2.length;
			output_map.RIP_QC = rip_QC.length;
		// }
		
		output_map.del_1n3 = data.spore_cmp.s1n3.values.length;
		output_map.del_2n2 = data.spore_cmp.s2n2.values.length;
		output_map.del_3n1 = data.spore_cmp.s3n1.values.length;
		output_map.del_4n0 = data.spore_cmp.s4n0.values.length;
		
		// if (dataset.options.includes("RIP")) {
			output_map.RIP_2_Q = rip_2_ref1.length;
			output_map.RIP_2_C = rip_2_ref2.length;
		// }
		
		output_map.sss22 = data.spore_cmp.sss_22.values.length;
		output_map.sss13 = data.spore_cmp.sss_13.values.length;
		output_map.sss04 = data.spore_cmp.sss_04.values.length;
		output_map.sss31 = data.spore_cmp.sss_31.values.length;
		output_map.sss40 = data.spore_cmp.sss_40.values.length;

		debugger
		
		output_map.IM1 = data.spore_cmp.illegitimate_mutation_list.values.length;
		// output_map.IM2 = data.spore_cmp.illegitimate_mutation_31.values.length;
		// output_map.IM3 = data.spore_cmp.illegitimate_mutation_40.values.length;
		// output_map.IM4 = data.spore_cmp.illegitimate_mutation_deletion.values.length;
		output_map.IM2 = data.spore_cmp.illegitimate_mutation_indel.values.length;
		output_map.IM3 = data.spore_cmp.illegitimate_mutation_deletion.values.length;

		// const chr = nChr;

		// const simpleCO_list = analysis_options.co_list.filter(co => co.type == "CO");
		// const CO_NCO_list = analysis_options.co_list.filter(co => co.type == "CO(NCO)");

		// const snv = ref_snv;
		// const snp = ref_snp;
		// const snp_indel = ref_snp_indel;

		// const illeg = data.spore_cmp.illegitimate_mutation_list;

		// const s22 = data.spore_cmp.s22;
		// const s31 = data.spore_cmp.s31;
		// const s40 = data.spore_cmp.s40;

		// const s1n3 = data.spore_cmp.s1n3;
		// const s2n2 = data.spore_cmp.s2n2;
		// const s3n1 = data.spore_cmp.s3n1;
		// const s4n0 = data.spore_cmp.s4n0;

		// Strain-specific sequences
		// const sss_22 = data.spore_cmp.sss_22;
		// const sss_31 = data.spore_cmp.sss_31;
		// const sss_40 = data.spore_cmp.sss_40;
		// const illegitimate_mutation_31 = data.spore_cmp.illegitimate_mutation_31;
		// const illegitimate_mutation_40 = data.spore_cmp.illegitimate_mutation_40;

		let output_text = "";

		console.table(output_map);

		// output_text += output_head.join("\t") + "\n";
		// output_text += [
		// 	chr,
		// 	simpleCO_list.length, CO_NCO_list.length,
		// 	snv.length, snp.length, snp_indel.length,
		// 	mut_ref1.length, mut_ref2.length,
		// 	illeg.values.length,
		// 	s22.values.length, s31.values.length, s40.values.length,
		// 	s1n3.values.length, s2n2.values.length, s3n1.values.length, s4n0.values.length,

		// 	// Strain-specific sequences
		// 	sss_22.values.length, sss_31.values.length, sss_40.values.length,
		// 	illegitimate_mutation_31.values.length, illegitimate_mutation_40.values.length,
		// ].join("\t") + "\n";

		// output_text += Object.keys(output_map).join("\t") + "\n";
		// output_text += Object.keys(output_map).map(k => output_map[k]).join("\t") + "\n";

		output_text += output_map.colNameList(dataset.parental_list[0], dataset.parental_list[1]).join("\t") + "\n";
		output_text += Object.values(output_map.outputTableRow(dataset.parental_list[0], dataset.parental_list[1])).join("\t") + "\n";

		fs.writeFileSync(`${output_path}/${argv_output_prefix}_summary.txt`, output_text);
	}

	if (argv_snp_count_window != null && argv_snp_count_window > 0) {
		[
			[ref_snp, "snp"],
			[ref_snv, "snv"],
			[ref_snp_indel, "indel"]
		].forEach(function (sn_info) {
			[
				[0, "Q"],
				[3, "C"]
			].forEach(function (pos_info) {
				let output_dir = `${output_path}/${sn_info[1]}_per_${argv_snp_count_window}bp/${pos_info[1]}`;
				fs.mkdirSync(output_dir, { recursive: true });

				let count_list = [];
	
				/**
				 * @param {{ [x: string]: any; }} snp
				 */
				sn_info[0].forEach(snp => {
					// const [r1, r2] = snp;
					const pos = snp[pos_info[0]];
					let wId = Math.trunc(Number(pos) / argv_snp_count_window);
					count_list[wId] = (count_list[wId] | 0) + 1;
				});
	
				const output_path = `${output_dir}/${argv_output_prefix}_${pos_info[1]}_${sn_info[1]}_per_${argv_snp_count_window}bp_ch${nChr}.txt`;
				const ws = fs.createWriteStream(output_path);
				ws.write(["chr", "start", "end", "count"].join("\t") + "\n");
				count_list.forEach((count, wId) => {
					const start = wId * argv_snp_count_window + 1;
					const end = start - 1 + argv_snp_count_window;
					let row = [nChr, start, end, count].join("\t") + "\n";
					ws.write(row);
				});
				ws.end();
			});
		});
	}
	
	fs.writeFile(`${output_path}/tmp_ch${analysis_options.nChr}.json`, JSON.stringify({
		marker_map: data.allMarker.map,
		pos_ref_map: [
			Array.from(data.pos_ref1_uint32array),
				Array.from(data.pos_ref2_uint32array),
		],
		ref_pos_map: [
			Array.from(data.ref1_pos_uint32array),
			Array.from(data.ref2_pos_uint32array),
		],
	}), function (err) {
		if (err) console.error(err);
	});
}

function analyser_gff_methyl(gff, data, analysis_options, seq_list) {
	if (gff) {
		dataset.parental_list.forEach((gff_ref, gff_sIdx) => {
			const pos_ref_map = [
				data.pos_ref1_uint32array,
				data.pos_ref2_uint32array,
			][gff_sIdx];

			const ref_pos_map = [
				data.ref1_pos_uint32array,
				data.ref2_pos_uint32array,
			][gff_sIdx];

			console.log("gff_sIdx", gff_sIdx, gff_ref);

			const out_path = `${dataset.output_path}/rip_gff`;
			if (!fs.existsSync(out_path)) {
				fs.mkdirSync(out_path);
			}
			const out_fa_dir_path = `${dataset.output_path}/rip_gff/fa`;
			if (!fs.existsSync(out_fa_dir_path)) {
				fs.mkdirSync(out_fa_dir_path);
			}

			console.log({
				"typeof gff": typeof gff,
			});

			if (analysis_options.rip_gff) {
				// const gene_type = "gene";
				// const rip_gff = find_rip_gff(analysis_options.nChr, seq_list, gff, data, gene_type);
				// const gene_type_map = {
				// 	"gene": "gene",
				// 	"CDS": "CDS",
				// };
				const rip_gff = data.rip_gff;

				console.log({
					"rip_gff.length": rip_gff.length,
				});

				/** @type {Set<GFF_ROW>} */
				const all_RIP_gene = new Set();

				/** @type {[ Set<GFF_ROW>, Set<GFF_ROW>, Set<GFF_ROW> ]} */
				const all_RIP_ref = [new Set(), new Set(), new Set()];

				const ref_rip_header = [
					"chr",
					"TSETA pos",
					`${gff_ref} pos`,
					"gene ID",
					"gene/mRNA/CDS/exon",
					"gene start",
					"gene end",
					"rip type",
					"rip ref",
				].join("\t");
				/** @type {string[][]} */
				const ref_rip_list = [
					[],
					[],
					[], // CBS
				];

				rip_gff.forEach(({ pos, marker, gene, rip_type, gc_value, AT_island }) => {
					// const pos = marker.pos;
					const ref1_pos = pos_ref_map[marker.pos];
					const chr_name = gene.seqid;

					/** @type {string} */
					const gene_id = gene.attributes.ID;
					/** @type {string} */
					const gene_type = gene.type;

					/** @type {number} */
					const gene_start = gene.start;
					/** @type {number} */
					const gene_end = gene.end;

					/** @type {"-"|"QM6a"|"CBS1-1"} */
					const mut_ref_name = (["-", ...dataset.parental_list][marker.mut_ref]);

					if (!all_RIP_gene.has(gene)) {
						all_RIP_gene.add(gene);
					}
					if (!all_RIP_ref[marker.mut_ref].has(gene)) {
						all_RIP_ref[marker.mut_ref].add(gene);
					}

					const line = [
						chr_name,
						pos,
						pos_ref_map[pos],

						gene_id,
						gene_type,
						gene_start,
						gene_end,

						gc_value,
						AT_island,

						rip_type,
						mut_ref_name,

						seq_list.map(ss => ss[marker.pos]),
						// endl
					].join("\t");

					ref_rip_list[marker.mut_ref].push(line);

					try {
						const gene_name = `ch${analysis_options.nChr}_${gene_id}`;
						const out_fa_file = `${out_fa_dir_path}/${gff_ref}_${mut_ref_name}_${gene_name}.fa`;
						const fa = {
							[gene_name]: seq_list[gff_sIdx].slice(gene.$start, gene.$end + 1).replace(/-/g, ""),
						};
						saveFasta(out_fa_file, fa);
					}
					catch (ex) {
					}
				});

				all_RIP_ref.forEach((gene_set, mut_ref) => {
					if (mut_ref == 0) {
						return;
					}
					/** @type {"-"|"QM6a"|"CBS1-1"} */
					const mut_ref_name = (["-", ...dataset.parental_list][mut_ref]);

					const out_file = `${out_path}/${gff_ref}_gene_${mut_ref_name}_rip_ch${analysis_options.nChr}.txt`;

					fs.writeFileSync(out_file, [
						`${mut_ref_name} RIP`, gene_set.size, [...gene_set].map(gene => gene.attributes.ID).join("\t"),
					].join("\t") + "\n", { flag: "a" });

					fs.writeFileSync(out_file, "\n\n", { flag: "a" });

					fs.writeFileSync(out_file, ref_rip_header + "\n" + ref_rip_list[mut_ref].join("\n"));
				});

				// fs.writeFileSync(out_file, [
				// 	"gene", all_RIP_gene.size, [...all_RIP_gene].map(gene => gene.attributes.ID).join("\t"),
				// ].join("\t") + "\n", { flag: "a" });
				// console.log({
				// 	"all_RIP_gene.size": all_RIP_gene.size,
				// });
				global.all_RIP_ref = all_RIP_ref;
				global.all_RIP_gene = all_RIP_gene;
			}

			if (analysis_options.rip_repeat_methyl_gff) {
				const IM_out_file = `${out_path}/${gff_ref}_gene_RIP_IM_ch${analysis_options.nChr}.txt`;

				class StructMapping {
					constructor() {
						this.ch = 0;

						/** @type {number} TSETA pos */
						this.pos = undefined;
						/** @type {number} QM6a pos */
						this.pos_ref_1 = undefined;
						/** @type {number} CBS1-1 pos */
						this.pos_ref_2 = undefined;

						/** @type {string} mutated referece */
						this.mut_ref = undefined;

						/** @type {"SNV"|"InDel"} */
						this.snv_indel = undefined;

						// A -> N
						this.mutation_type = undefined;

						/** @type {"intragenic"|"extragenic"} */
						this.inout_gene = undefined;
						/** @type {string} */
						this.gene_ID = undefined;
						/** @type {"gene"|"mRNA"|"exon"|"CDS"} */
						this.struct_type = undefined;

						/** @type {number} */
						this.gc_value = undefined;
						/** @type {"+"|"-"} */
						this.AT_island = undefined;

						// /** @type {"+"|"-"} in repeat seq ? "+" : "-" */
						/** @type {number} */
						this.repeat = undefined;

						// /** @type {"+"|"-"} in repeat seq ? "+" : "-" */
						/** @type {number} */
						this.repeat_60 = undefined;

						// // /** @type {"+"|"-"} in repeat seq ? "+" : "-" */
						// /** @type {number} */
						// this.repeat_70 = undefined;
						// // /** @type {"+"|"-"} in repeat seq ? "+" : "-" */
						// /** @type {number} */
						// this.repeat_65 = undefined;
						// /** @type {"+"|"-"} in strain specific seq ? "+" : "-" */
						// this.sss_50bp = undefined;
						// /** @type {"+"|"-"} in strain specific seq ? "+" : "-" */
						// this.sss_100bp = undefined;
						/** @type {"+"|"-"} in strain specific seq ? "+" : "-" */
						this.sss_60bp = undefined;

						this.methyl_ratio_l10 = undefined;
						this.methyl_ratio_l9 = undefined;
						this.methyl_ratio_l8 = undefined;
						this.methyl_ratio_l7 = undefined;
						this.methyl_ratio_l6 = undefined;
						this.methyl_ratio_l5 = undefined;
						this.methyl_ratio_l4 = undefined;
						this.methyl_ratio_l3 = undefined;
						this.methyl_ratio_l2 = undefined;
						this.methyl_ratio_l1 = undefined;
						this.methyl_ratio_r0 = undefined;
						this.methyl_ratio_r1 = undefined;
						this.methyl_ratio_r2 = undefined;
						this.methyl_ratio_r3 = undefined;
						this.methyl_ratio_r4 = undefined;
						this.methyl_ratio_r5 = undefined;
						this.methyl_ratio_r6 = undefined;
						this.methyl_ratio_r7 = undefined;
						this.methyl_ratio_r8 = undefined;
						this.methyl_ratio_r9 = undefined;
						this.methyl_ratio_r10 = undefined;

						// this.methyl_ratio_4l10 = undefined;
						// this.methyl_ratio_4l9 = undefined;
						// this.methyl_ratio_4l8 = undefined;
						// this.methyl_ratio_4l7 = undefined;
						// this.methyl_ratio_4l6 = undefined;
						// this.methyl_ratio_4l5 = undefined;
						// this.methyl_ratio_4l4 = undefined;
						// this.methyl_ratio_4l3 = undefined;
						// this.methyl_ratio_4l2 = undefined;
						// this.methyl_ratio_4l1 = undefined;
						// this.methyl_ratio_4r0 = undefined;
						// this.methyl_ratio_4r1 = undefined;
						// this.methyl_ratio_4r2 = undefined;
						// this.methyl_ratio_4r3 = undefined;
						// this.methyl_ratio_4r4 = undefined;
						// this.methyl_ratio_4r5 = undefined;
						// this.methyl_ratio_4r6 = undefined;
						// this.methyl_ratio_4r7 = undefined;
						// this.methyl_ratio_4r8 = undefined;
						// this.methyl_ratio_4r9 = undefined;
						// this.methyl_ratio_4r10 = undefined;
						// this.methyl_ratio_8l10 = undefined;
						// this.methyl_ratio_8l9 = undefined;
						// this.methyl_ratio_8l8 = undefined;
						// this.methyl_ratio_8l7 = undefined;
						// this.methyl_ratio_8l6 = undefined;
						// this.methyl_ratio_8l5 = undefined;
						// this.methyl_ratio_8l4 = undefined;
						// this.methyl_ratio_8l3 = undefined;
						// this.methyl_ratio_8l2 = undefined;
						// this.methyl_ratio_8l1 = undefined;
						// this.methyl_ratio_8r0 = undefined;
						// this.methyl_ratio_8r1 = undefined;
						// this.methyl_ratio_8r2 = undefined;
						// this.methyl_ratio_8r3 = undefined;
						// this.methyl_ratio_8r4 = undefined;
						// this.methyl_ratio_8r5 = undefined;
						// this.methyl_ratio_8r6 = undefined;
						// this.methyl_ratio_8r7 = undefined;
						// this.methyl_ratio_8r8 = undefined;
						// this.methyl_ratio_8r9 = undefined;
						// this.methyl_ratio_8r10 = undefined;
						this["QM6a"] = undefined;
						this["CBS_1-1"] = undefined;
						this["WT-D1"] = undefined;
						this["WT-D2"] = undefined;
						this["WT-D3"] = undefined;
						this["WT-D4"] = undefined;
						this["WT-D5"] = undefined;
						this["WT-D6"] = undefined;
						this["WT-D7"] = undefined;
						this["WT-D8"] = undefined;
						this["D1"] = undefined;
						this["D2"] = undefined;
						this["D3"] = undefined;
						this["D4"] = undefined;
						this["D5"] = undefined;

						/** @type {string} seq.length = 1 */
						this.s = undefined;

						/** @type {string} peek_seq_range = 10; seq.slice(pos - peek_seq_range, pos + peek_seq_range + 1) */
						this.peek = undefined;
					}
				}

				/**
				 * @template T
				 * @template V
				 * @typedef TypeMapping<T,V>
				 * @type {{ [P in keyof T]: V; }}
				 */
				// const o_methyl_ratio = {};
				// for (let i = analysis_options.methyl_left; i > 0; --i) {
				// 	o_methyl_ratio[`methyl_ratio_l${i}`] = `me % -${i}`;
				// }
				// for (let i = 0; i <= analysis_options.methyl_right; ++i) {
				// 	o_methyl_ratio[`methyl_ratio_r${i}`] = `me % +${i}`;
				// }
				/**
				 * @type {TypeMapping<StructMapping, string>}
				 */
				const output_header = Object.assign(new StructMapping(), {
					ch: "ch",
					pos: "TSETA pos",
					pos_ref_1: `${dataset.parental_list[0]} pos`,
					pos_ref_2: `${dataset.parental_list[1]} pos`,
					mut_ref: "mutated referece",
					snv_indel: "SNV/InDel",
					mutation_type: "mutation",
					inout_gene: "intragenic/extragenic",
					gene_ID: "gene ID",
					struct_type: "CDS/exon/intron",
					gc_value: "GC content %",
					AT_island: "AT island",
					// repeat: "repeat (identity >= 70)",
					// repeat_70: "repeat (65 <= identity < 70)",
					// repeat_65: "repeat (identity < 65)",
					repeat: "repeat (identity >= 65)",
					repeat_60: "repeat (identity >= 60)",
					// sss_50bp: "strain specific 50bp",
					// sss_100bp: "strain specific 100bp",
					sss_60bp: "strain specific 60bp",

					methyl_ratio_l10: "me % -10",
					methyl_ratio_l9: "me % -9",
					methyl_ratio_l8: "me % -8",
					methyl_ratio_l7: "me % -7",
					methyl_ratio_l6: "me % -6",
					methyl_ratio_l5: "me % -5",
					methyl_ratio_l4: "me % -4",
					methyl_ratio_l3: "me % -3",
					methyl_ratio_l2: "me % -2",
					methyl_ratio_l1: "me % -1",
					methyl_ratio_r0: "me % +0",
					methyl_ratio_r1: "me % +1",
					methyl_ratio_r2: "me % +2",
					methyl_ratio_r3: "me % +3",
					methyl_ratio_r4: "me % +4",
					methyl_ratio_r5: "me % +5",
					methyl_ratio_r6: "me % +6",
					methyl_ratio_r7: "me % +7",
					methyl_ratio_r8: "me % +8",
					methyl_ratio_r9: "me % +9",
					methyl_ratio_r10: "me % +10",

					// methyl_ratio_4l10: "D4 me % -10",
					// methyl_ratio_4l9: "D4 me % -9",
					// methyl_ratio_4l8: "D4 me % -8",
					// methyl_ratio_4l7: "D4 me % -7",
					// methyl_ratio_4l6: "D4 me % -6",
					// methyl_ratio_4l5: "D4 me % -5",
					// methyl_ratio_4l4: "D4 me % -4",
					// methyl_ratio_4l3: "D4 me % -3",
					// methyl_ratio_4l2: "D4 me % -2",
					// methyl_ratio_4l1: "D4 me % -1",
					// methyl_ratio_4r0: "D4 me % +0",
					// methyl_ratio_4r1: "D4 me % +1",
					// methyl_ratio_4r2: "D4 me % +2",
					// methyl_ratio_4r3: "D4 me % +3",
					// methyl_ratio_4r4: "D4 me % +4",
					// methyl_ratio_4r5: "D4 me % +5",
					// methyl_ratio_4r6: "D4 me % +6",
					// methyl_ratio_4r7: "D4 me % +7",
					// methyl_ratio_4r8: "D4 me % +8",
					// methyl_ratio_4r9: "D4 me % +9",
					// methyl_ratio_4r10: "D4 me % +10",
					// methyl_ratio_8l10: "D8 me % -10",
					// methyl_ratio_8l9: "D8 me % -9",
					// methyl_ratio_8l8: "D8 me % -8",
					// methyl_ratio_8l7: "D8 me % -7",
					// methyl_ratio_8l6: "D8 me % -6",
					// methyl_ratio_8l5: "D8 me % -5",
					// methyl_ratio_8l4: "D8 me % -4",
					// methyl_ratio_8l3: "D8 me % -3",
					// methyl_ratio_8l2: "D8 me % -2",
					// methyl_ratio_8l1: "D8 me % -1",
					// methyl_ratio_8r0: "D8 me % +0",
					// methyl_ratio_8r1: "D8 me % +1",
					// methyl_ratio_8r2: "D8 me % +2",
					// methyl_ratio_8r3: "D8 me % +3",
					// methyl_ratio_8r4: "D8 me % +4",
					// methyl_ratio_8r5: "D8 me % +5",
					// methyl_ratio_8r6: "D8 me % +6",
					// methyl_ratio_8r7: "D8 me % +7",
					// methyl_ratio_8r8: "D8 me % +8",
					// methyl_ratio_8r9: "D8 me % +9",
					// methyl_ratio_8r10: "D8 me % +10",
					// ...o_methyl_ratio,
					["QM6a"]: "QM6a",
					["CBS_1-1"]: "CBS 1-1",
					["WT-D1"]: "WT-D1",
					["WT-D2"]: "WT-D2",
					["WT-D3"]: "WT-D3",
					["WT-D4"]: "WT-D4",
					["WT-D5"]: "WT-D5",
					["WT-D6"]: "WT-D6",
					["WT-D7"]: "WT-D7",
					["WT-D8"]: "WT-D8",
					["D1"]: "CC D1",
					["D2"]: "CC D2",
					["D3"]: "CC D3",
					["D4"]: "CC D4",
					["D5"]: "CC D5",

					s: "s",
					peek: "pos-10 ~ pos+10",
				});
				const output_keys = Object.keys(output_header);

				const deletetion_sign = "Del";

				const IM_out_text = data.IM_gff[gff_ref].map(_row => {
					/** @type {StructMapping} */
					const row = _row;

					if (analysis_options.nChr == 6 && row.mut_ref == 2 && row.pos <= 9116) {
						return null;
					}

					/** @type {"-"|"QM6a"|"CBS1-1"} */
					const mut_ref_name = (["-", ...dataset.parental_list][row.mut_ref]);

					row.ch = analysis_options.nChr;
					row.pos_ref_1 = data.pos_ref1_uint32array[row.pos] + 1;
					row.pos_ref_2 = data.pos_ref2_uint32array[row.pos] + 1;
					row.mut_ref = mut_ref_name;

					// row.mutation_type = `"` + row.mutation_type.replace(/-🢂/g, "-->").replace(/🢂-/g, "->-").replace(/🢂/g, "->") + `"`;
					// row.mutation_type = `"` + row.mutation_type.replace(/-🢂/g, "N->").replace(/🢂-/g, "->N").replace(/🢂/g, "->") + `"`;
					row.mutation_type = `"` + row.mutation_type.replace(/-🢂/g, `${deletetion_sign}->`).replace(/🢂-/g, `->${deletetion_sign}`).replace(/🢂/g, "->") + `"`;
					row.peek = `"${row.peek}"`;

					const line = output_keys.map(k => row[k]).join("\t");

					return line;
				}).filter(a => a).join("\n");

				// `${gff_ref} pos`,// ref pos
				fs.writeFileSync(IM_out_file, Object.values(output_header).join("\t") + "\n" + IM_out_text);
			}

		});
	}

	// if (data.rip_seq_tuple) {
	// 	const out_path = `${dataset.output_path}/rip`;
	// 	if (!fs.existsSync(out_path)) {
	// 		fs.mkdirSync(out_path);
	// 	}
	// 	const path_prefix = `${out_path}/rip_ch${analysis_options.nChr}`;

	// 	// console.log({
	// 	// 	"data.rip_seq_tuple": data.rip_seq_tuple.length,
	// 	// 	"data.rip_seq_tuple_map": data.rip_seq_tuple.map(a => a.length).join("-"),
	// 	// });

	// 	const fin_rip = data.rip_seq_tuple.map(grp => {
	// 		return grp.map(tuple => {
	// 			const [info, seq, pos, table_row] = tuple;

	// 			// console.log("ss", analysis_options.nChr, pos);

	// 			const input_path = `${path_prefix}_${pos}.fa`;
	// 			const re_path = `${path_prefix}_${pos}_re.txt`;

	// 			fs.writeFileSync(input_path, seq);

	// 			const num_thread = 1;
	// 			const _algorithm = "";
	// 			const maxiterate = 10000;

	// 			const arg_algorithm = _algorithm ? `--${_algorithm}` : "";
	// 			const arg_thread = num_thread >= 1 ? `--thread ${20}` : "";
	// 			const arg_maxiterate = maxiterate >= 0 ? `--maxiterate ${maxiterate}` : "";

	// 			const output_file = `mafft_${input_path}`;

	// 			const mafft_cmd = `mafft --quiet ${arg_thread} ${arg_algorithm} ${arg_maxiterate} ${input_path} > ${output_file}`;
				
	// 			child_process.execSync(mafft_cmd);

	// 			const re_cmd = `node tetrad_rip.js -dataset ${argv_dataset_path} -chr ${analysis_options.nChr} -seq ${output_file} -output ${re_path}`;
	// 			child_process.execSync(re_cmd);

	// 			return info;
	// 		}).join("\n" + "-".repeat(60) + "\n");
	// 	}).join("\n" + "=".repeat(60) + "\n");
		
	// 	fs.writeFileSync(`${path_prefix}_info.txt`, fin_rip);
	// }
}

/**
 * @param {string[]} seq_list
 * @param {length} [min_specific_length] default = 10
 */
function get_strain_specific_region_by_length(seq_list, min_specific_length = 10) {
	const regexp = new RegExp(`1{${min_specific_length},}`, "g");

	const qqq = [...seq_list[0]].map((a, i) => a != "-" && seq_list[1][i] == "-" ? "1" : "0").join("");
	/** @type {RegExpMatchArray[]} */
	// const qqq_sss = [...qqq.matchAll(/1{10,}/g)];
	const qqq_sss = [...qqq.matchAll(regexp)];
	const qqq_sss_t_len = qqq_sss.reduce((t, v) => t + v[0].length, 0);

	const ccc = [...seq_list[0]].map((a, i) => a == "-" && seq_list[1][i] != "-" ? "1" : "0").join("");
	/** @type {RegExpMatchArray[]} */
	// const ccc_sss = [...ccc.matchAll(/1{10,}/g)];
	const ccc_sss = [...ccc.matchAll(regexp)];
	const ccc_sss_t_len = ccc_sss.reduce((t, v) => t + v[0].length, 0);

	return { qqq_sss_t_len, ccc_sss_t_len, qqq_sss, ccc_sss };
}

/*
node .\src\tetrad_summary.js -dataset QCt.json -min-co 5000 --verbose --snp-count-window 500
node .\src\tetrad_summary.js -dataset QCt.json -min-co 5000 --verbose --snp-count-window 3000
node .\src\tetrad_summary.js -dataset QCt.json -min-co 5000 --verbose --snp-count-window 5000
*/

/**
 * @param {number} sIdx
 * @param {string} file_path
 * @param {number} nChr
 * param {Uint32Array|number[]} pos_map
 */
function load_gff(sIdx, file_path, nChr) {
	// const chr_name_list = dataset.genomeNameList.map((n, si) => genome_info_list[si].chr_list[nChr - 1]).map((a, si) => a.chr.replace(dataset.genomeNameList[si], "").slice(1).replace(/_unitig_.*consensus/, ""));
	const chr_name_list = dataset.genomeNameList.map((n, si) => genome_info_list[si].chr_list[nChr - 1]).map((a, si) => a.symbol);

	// const pos_map_map = {
	// 	"Q": {
	// 		pos_map: ref1_pos_uint32array,
	// 		file: "QM6a.gff3",
	// 	},
	// 	"C": {
	// 		pos_map: ref2_pos_uint32array,
	// 		file: "Trichoderma_reesei_CBS.gff3",
	// 	},
	// };

	const text = fs.readFileSync(file_path).toString();

	const gff_data = parseGFF(text);
	console.log({
		gff_keys: Object.keys(gff_data),
	});

	// Object.values(gff_data).forEach(gene_list => {
	// 	// console.log({
	// 	// 	"gene_list.length": gene_list.length,
	// 	// });
	// 	gene_list.forEach(gene => {
	// 		gene.$start = pos_map[gene.start - 1];
	// 		gene.$length = pos_map[gene.end - 1] - gene.$start + 1;
	// 	});
	// });

	// row.seqid
	// row.type == "gene";
	// row.type == "CDS";

	// console.log(chr_name_list[0]);

	return gff_data[chr_name_list[sIdx]];
}

class RepeatSegment {
	constructor() {
		/** TSETA pos */
		this.start = 0;

		/** TSETA pos */
		this.end = 0;
		
		this.q_start = 0;
		this.q_end = 0;

		this.identity = 0;
		
		this.alignment = 0;
		
		this.q_len = 0;

		this.s_len = 0;
	}
}

/**
 * @param {string} url
 * @param {number} chr_length
 */
function _loadRepeatSegment(url, chr_length) {
	/** @type {string} */
	const text = fs.readFileSync(url).toString();
	const rows = text.split("\n").map(line => {
		const row = line.trim().split("\t");
		if (row.length <= 11) {
			return;
		}
		// return [
		// 	// row[0], // q
		// 	// row[1], // s
		// 	// row[2], // identity
		// 	// row[3], // alignment length
		// 	// row[4], // mismatch
		// 	// row[5], // gap
		// 	row[6], // q.start
		// 	row[7], // q.end
		// 	// row[8], // s.start
		// 	// row[9], // s.end
		// 	// row[10], // evalue
		// 	// row[11], // score
		// ];

		const q_start = Number(row[6].trim());
		const q_end = Number(row[7].trim());
		const s_start = Number(row[8].trim());
		const s_end = Number(row[9].trim());

		const q_len = Math.abs(q_end - q_start) + 1;
		const s_len = Math.abs(s_end - s_start) + 1;

		// const start = ref_pos_map[q_start - 1];
		// const end = ref_pos_map[q_end - 1];
		const identity = Number(row[2].trim());
		const alignment = Number(row[3].trim());// alignment (end - start + 1) >= 100

		if (alignment < chr_length) {
			// if (identity >= 65) {
			// if (alignment >= 100 && q_len >= 100 && s_len >= 100) {
			if (alignment >= 100) {// || q_len >= 100 || s_len >= 100
				return Object.assign(new RepeatSegment(), {
					// start: start,
					// end: end,
					q_start: q_start,
					q_end: q_end,
					identity: identity,
					alignment: alignment,
					q_len: q_len,
					s_len: s_len,
				});
			}
		// }
		}
	}).filter(a => a);
	return rows;
}

/**
 * @param {string} ref
 * @param {number} nChr
 */
function LoadRepeatSegment(ref, nChr) {
	// const url = `repeat_segment/${ref}_ch${nChr}_repeat_${in_or_out}.txt`;
	// const url = `repeat_segment/${ref}_ch${nChr}_repeat_i65_${["inout", "in", "out"][in_or_out]}.txt`;
	const url = `${dataset.output_path}/repeat_segment/${ref}_ch${nChr}_repeat_blastn.txt`;
	// const chr_length = analysis_options.chrInfo_list[dataset.parental_list.indexOf(ref)].length;
	const chr_length = genome_info_list[dataset.parental_list.indexOf(ref)].chr_list[nChr - 1].length;
	// console.log(genome_info_list[dataset.parental_list.indexOf(ref)].chr_list[nChr - 1].symbol);
	return _loadRepeatSegment(url, chr_length);
}

/**
 * @param {number} gff_idx
 * @param {number} nChr
 * @param {string} file_path
 */
function LoadMethylRatio(gff_idx, nChr, file_path) {
	const buffer = fs.readFileSync(file_path).buffer;

	/** @type {number[]} */
	const chrLengthList = genome_info_list[gff_idx].chr_list.map(a => a.length);

	const chrOffsetInByte = [];
	chrOffsetInByte[0] = 0;
	chrLengthList.reduce((pos, len, i) => {
		let offset = pos + (len * 4);
		chrOffsetInByte[i + 1] = offset;
		return offset;
	}, 0);

	const chrIdx = nChr - 1;
	const offsetStart = chrOffsetInByte[chrIdx];
	
	// console.log({
	// 	file_path,
	// 	gff_idx,
	// 	chrIdx,
	// 	symbol: genome_info_list[gff_idx].chr_list[chrIdx].symbol,
	// 	"buffer.byteLength": buffer.byteLength,
	// 	offsetStart,
	// 	"chrLengthList[chrIdx]": chrLengthList[chrIdx],
	// });

	/** array view */
	const f32a = new Float32Array(buffer, offsetStart, chrLengthList[chrIdx]);

	return f32a;
}


main();
