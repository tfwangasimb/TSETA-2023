// @ts-check

// 2021/08/31 11:17
// 2021/10/26 11:44

// conda install -c bioconda ucsc-bigwigtobedgraph
// conda install -c bioconda ucsc-bigwigtowig

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;

if (DEBUG) {
	process.on("unhandledRejection", (reason, p) => {
		console.log("Unhandled Rejection at: Promise", p, "reason:", reason);
		writeLog(p.toString());
		writeLog(reason.toString());
		writeLog("\n");
		process.exit(1);
	});
}

const Path = require("path");
if (process.env.PATH && process.env.PATH.indexOf("bedops") < 0) {
	throw new Error("bedops is required")
}

const os = require("os");
const fs = require("fs");
const child_process = require("child_process");
const cluster = require("cluster");
// const { Worker, workerData, isMainThread, parentPort, threadId } = require("worker_threads");

const samtools = require("../../samtools.js");
const process_utils = require("../../process_utils.js");
// const trimmomatic = require("./run_trimmomatic.js");

const TAG_DATE = (function () {
	const date = new Date();
	return [
		date.getFullYear(),
		date.getMonth() + 1,
		date.getDate(),
	].join("");
})();

// // ssDNA remove unpaired reads
// // remove unpaired reads
// const RUN_TRIMMOMATIC = false;


// const trimmomaticOptions = new trimmomatic.TrimmomaticOptions();
// trimmomatic.TrimmomaticOptions.setBinPath("/opt/app/Trimmomatic-0.32/trimmomatic-0.32.jar");

// const OUTPUT_TRIM_PATH = "./trim";

// const output_results_dir = "./results";

const MAX_THREAD = 16;

// /**
//  * @type {{ refId: string; genome_path: string; gff_path: string; samples: { output_tag: string; reads: string[]; }[]; }[]}
//  */
// let dataset_list;

/**
 * @type {{ param: { binSize: number; min_lg: number; max_lg: number; min_peak_len: number; pValue: number; }; cmp: { output_tag: string; sample_1: string; sample_2: string; }[]; sample: { [sample: string]: string[]; }; }}
 */
const cmp_bam = {
	"param": {
		"binSize": 10,
		"min_lg": 2,
		"max_lg": 4,
		"min_peak_len": 11,
		"pValue": 0.05,
	},
	"cmp": [
		// {
		// 	"output_tag": "CBS1-1_rad51_vs_sae2",
		// 	"sample_1": "CBS1-1_SS_Rad51",
		// 	"sample_2": "CBS1-1_SS_Sae2",
		// },
		// {
		// 	"output_tag": "CBS1-1_spo11rad51_vs_sae2",
		// 	"sample_1": "CBS1-1_ss_spo11_rad51",
		// 	"sample_2": "CBS1-1_SS_Sae2"
		// },
		// {
		// 	"output_tag": "CBS1-1_spo11rad51_vs_spo11sae2",
		// 	"sample_1": "CBS1-1_ss_spo11_rad51",
		// 	"sample_2": "CBS1-1_ss_spo11_sae2"
		// },
		// {
		// 	"output_tag": "CBS1-1_spo11sae2_vs_sae2",
		// 	"sample_1": "CBS1-1_ss_spo11_sae2",
		// 	"sample_2": "CBS1-1_SS_Sae2"
		// },
		// {
		// 	"output_tag": "rad51_vs_sae2",
		// 	"sample_1": "QM6a_SS_Rad51",
		// 	"sample_2": "QM6a_SS_Sae2"
		// },
		// {
		// 	"output_tag": "spo11rad51_vs_sae2",
		// 	"sample_1": "QM6a_ss_spo11_rad51",
		// 	"sample_2": "QM6a_SS_Sae2"
		// },
		// {
		// 	"output_tag": "spo11rad51_vs_spo11sae2",
		// 	"sample_1": "QM6a_ss_spo11_rad51",
		// 	"sample_2": "QM6a_ss_spo11_sae2"
		// },
		// {
		// 	"output_tag": "spo11sae2_vs_sae2",
		// 	"sample_1": "QM6a_ss_spo11_sae2",
		// 	"sample_2": "QM6a_SS_Sae2"
		// },
	],
	"sample": {
		"QM6a_rad51": [
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_1.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_2.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_3.sort.bam",
		],
		// "QM6a_rad51_1": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_1.sort.bam",
		// ],
		// "QM6a_rad51_2": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_2.sort.bam",
		// ],
		// "QM6a_rad51_3": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_3.sort.bam",
		// ],
		"QM6a_spo11rad51": [
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_1.sort.bam",
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_2.sort.bam",
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_3.sort.bam",
		],
		// "QM6a_spo11rad51_1": [
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_1.sort.bam",
		// ],
		// "QM6a_spo11rad51_2": [
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_2.sort.bam",
		// ],
		// "QM6a_spo11rad51_3": [
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_3.sort.bam",
		// ],
		"QM6a_sae2": [
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_1.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_2.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_3.sort.bam",
		],
		"QM6a_spo11sae2": [
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_1.sort.bam",
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_2.sort.bam",
			"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_3.sort.bam",
		],
		
		"CBS1-1_rad51": [
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_1.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_2.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_3.sort.bam",
		],
		// "CBS1-1_rad51_1": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_1.sort.bam",
		// ],
		// "CBS1-1_rad51_2": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_2.sort.bam",
		// ],
		// "CBS1-1_rad51_3": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_3.sort.bam",
		// ],
		"CBS1-1_spo11rad51": [
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_1.sort.bam",
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_2.sort.bam",
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_3.sort.bam",
		],
		// "CBS1-1_spo11rad51_1": [
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_1.sort.bam",
		// ],
		// "CBS1-1_spo11rad51_2": [
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_2.sort.bam",
		// ],
		// "CBS1-1_spo11rad51_3": [
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_3.sort.bam",
		// ],
		"CBS1-1_sae2": [
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_1.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_2.sort.bam",
			"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_3.sort.bam",
		],
		"CBS1-1_spo11sae2": [
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_1.sort.bam",
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_2.sort.bam",
			"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_3.sort.bam"
		],

		// "QM6a_SS_Rad51": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_1.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_2.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_3.sort.bam",
		// ],
		// // "../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Rad51_merge.bam",
		// "QM6a_SS_Sae2": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_1.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_2.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_3.sort.bam",
		// ],
		// // "../20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_SS_Sae2_merge.bam",
		// "QM6a_ss_spo11_rad51": [
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_1.sort.bam",
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_2.sort.bam",
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_rad51_3.sort.bam",
		// ],
		// "QM6a_ss_spo11_sae2": [
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_1.sort.bam",
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_2.sort.bam",
		// 	"../20210503_ssDNA_mapto_QM6a/QM6a_ss_spo11_sae2_3.sort.bam",
		// ],

		// "CBS1-1_SS_Rad51": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_1.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_2.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Rad51_3.sort.bam"
		// ],
		// "CBS1-1_ss_spo11_rad51": [
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_1.sort.bam",
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_2.sort.bam",
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_rad51_3.sort.bam"
		// ],
		// "CBS1-1_SS_Sae2": [
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_1.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_2.sort.bam",
		// 	"../20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_SS_Sae2_3.sort.bam"
		// ],
		// "CBS1-1_ss_spo11_sae2": [
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_1.sort.bam",
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_2.sort.bam",
		// 	"../20211026_ssDNA_mapto_CBS1-1/CBS1-1_ss_spo11_sae2_3.sort.bam"
		// ],

		// "../20210623_ssDNA_mapto/rad51_f_rad51_1.sort.bam",
		// "../20210623_ssDNA_mapto/rad51_f_rad51_2.sort.bam",
		// "../20210623_ssDNA_mapto/rad51_f_rad51_3.sort.bam",
		// "../20210623_ssDNA_mapto/rad51_m_rad51_1.sort.bam",
		// "../20210623_ssDNA_mapto/rad51_m_rad51_2.sort.bam",
		// "../20210623_ssDNA_mapto/rad51_m_rad51_3.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_f_sae2_1.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_f_sae2_2.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_f_sae2_3.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_m_sae2_1.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_m_sae2_2.sort.bam",
		// "../20210623_ssDNA_mapto/sae2_m_sae2_3.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_1.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_2.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_3.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_1.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_2.sort.bam",
		// "../20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_3.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_1.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_2.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_3.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_1.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_2.sort.bam",
		// "../20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_3.sort.bam",
	},
};

function each_trip() {
	const rep_lo = 1;// 1
	const rep_hi = 3;// 3
	[
		"QM6a",
		"CBS1-1",
	].forEach(ref => {
		[
			"rad51",
			"spo11rad51",
		].forEach(target => {
			for (let i = rep_lo; i <= rep_hi; ++i) {
				const output_tag = `${ref}_${target}_${i}_vs_sae2`;

				// if ([
				// 	"CBS1-1_rad51_3_vs_sae2",// CBS1-1_rad51_3_vs_sae2/CBS1-1_sae2_all.sort.bam.scaled.sum_mean.tsv
				// 	// "QM6a_spo11rad51_1_vs_sae2",// QM6a_spo11rad51_1_vs_sae2/QM6a_sae2_all.sort.bam.scaled.sum_mean.tsv
				// 	// "QM6a_spo11rad51_3_vs_sae2",// QM6a_spo11rad51_3_vs_sae2/QM6a_sae2_all.sort.bam.scaled.sum_mean.tsv
				// ].includes(output_tag)) {
					cmp_bam.cmp.push({
						"output_tag": output_tag,
						"sample_1": `${ref}_${target}_${i}`,
						"sample_2": `${ref}_sae2`,
					});
				// }//output_tag
			}
		});
	});
}
/**
 * 
 * @param {("rad51"|"spo11rad51")[]} target
 * @param {"sae2"|"spo11sae2"} control
 */
function each_cmp(target, control) {
	[
		"QM6a",
		"CBS1-1",
	].forEach(ref => {
		target.forEach(target => {
			const output_tag = `${ref}_${target}_vs_${control}`;

			cmp_bam.cmp.push({
				"output_tag": output_tag,
				"sample_1": `${ref}_${target}`,
				"sample_2": `${ref}_${control}`,
			});
		});
	});
}
// console.log(cmp_bam.cmp);


// [
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_rad51_1_1_bin10_scaled/QM6a_SS_Rad51_1.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_sae2_1_bin10_scaled/QM6a_SS_Sae2_1.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_sae2_2_bin10_scaled/QM6a_SS_Sae2_2.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_sae2_3_bin10_scaled/QM6a_SS_Sae2_3.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_rad51_2_1_bin10_scaled/QM6a_SS_Rad51_2.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/QM6a_rad51_3_1_bin10_scaled/QM6a_SS_Rad51_3.sort.bam.scaled.bed",
// 	"/ssDNA/20210503_ssDNA_mapto_QM6a/QM6a_spo11rad51_1_1_bin10_scaled/QM6a_ss_spo11_rad51_1.sort.bam.scaled.bed",
// 	"/ssDNA/20210503_ssDNA_mapto_QM6a/QM6a_spo11rad51_2_1_bin10_scaled/QM6a_ss_spo11_rad51_2.sort.bam.scaled.bed",
// 	"/ssDNA/20210503_ssDNA_mapto_QM6a/QM6a_spo11rad51_3_1_bin10_scaled/QM6a_ss_spo11_rad51_3.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_rad51_1_1_bin10_scaled/CBS1-1_SS_Rad51_1.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_sae2_1_bin10_scaled/CBS1-1_SS_Sae2_1.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_sae2_2_bin10_scaled/CBS1-1_SS_Sae2_2.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_sae2_3_bin10_scaled/CBS1-1_SS_Sae2_3.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_rad51_2_1_bin10_scaled/CBS1-1_SS_Rad51_2.sort.bam.scaled.bed",
// 	"/ssDNA/20210324_ssDNA_mapto_QM6a_CBS1-1/CBS1-1_rad51_3_1_bin10_scaled/CBS1-1_SS_Rad51_3.sort.bam.scaled.bed",
// 	"/ssDNA/20211026_ssDNA_mapto_CBS1-1/CBS1-1_spo11rad51_1_1_bin10_scaled/CBS1-1_ss_spo11_rad51_1.sort.bam.scaled.bed",
// 	"/ssDNA/20211026_ssDNA_mapto_CBS1-1/CBS1-1_spo11rad51_2_1_bin10_scaled/CBS1-1_ss_spo11_rad51_2.sort.bam.scaled.bed",
// 	"/ssDNA/20211026_ssDNA_mapto_CBS1-1/CBS1-1_spo11rad51_3_1_bin10_scaled/CBS1-1_ss_spo11_rad51_3.sort.bam.scaled.bed",
// ].forEach(a => {
// 	// document.write(`<a href="${a}">${a}</a><br>`)
// })

/*

awk '{print $5}' QM6a_SS_Rad51_1.sort.bam.scaled.bed > QM6a_rad51_1.csv
awk '{print $5}' QM6a_SS_Rad51_2.sort.bam.scaled.bed > QM6a_rad51_2.csv
awk '{print $5}' QM6a_SS_Rad51_3.sort.bam.scaled.bed > QM6a_rad51_3.csv
awk '{print $5}' QM6a_SS_Sae2_1.sort.bam.scaled.bed > QM6a_sae2_1.csv
awk '{print $5}' QM6a_SS_Sae2_2.sort.bam.scaled.bed > QM6a_sae2_2.csv
awk '{print $5}' QM6a_SS_Sae2_3.sort.bam.scaled.bed > QM6a_sae2_3.csv
awk '{print $5}' QM6a_ss_spo11_rad51_1.sort.bam.scaled.bed > QM6a_spo11rad51_1.csv
awk '{print $5}' QM6a_ss_spo11_rad51_2.sort.bam.scaled.bed > QM6a_spo11rad51_2.csv
awk '{print $5}' QM6a_ss_spo11_rad51_3.sort.bam.scaled.bed > QM6a_spo11rad51_3.csv
awk '{print $5}' CBS1-1_SS_Rad51_1.sort.bam.scaled.bed > CBS1-1_rad51_1.csv
awk '{print $5}' CBS1-1_SS_Rad51_2.sort.bam.scaled.bed > CBS1-1_rad51_2.csv
awk '{print $5}' CBS1-1_SS_Rad51_3.sort.bam.scaled.bed > CBS1-1_rad51_3.csv
awk '{print $5}' CBS1-1_SS_Sae2_1.sort.bam.scaled.bed > CBS1-1_spo11rad51_1.csv
awk '{print $5}' CBS1-1_SS_Sae2_2.sort.bam.scaled.bed > CBS1-1_spo11rad51_2.csv
awk '{print $5}' CBS1-1_SS_Sae2_3.sort.bam.scaled.bed > CBS1-1_spo11rad51_3.csv
awk '{print $5}' CBS1-1_ss_spo11_rad51_1.sort.bam.scaled.bed > CBS1-1_sae2_1.csv
awk '{print $5}' CBS1-1_ss_spo11_rad51_2.sort.bam.scaled.bed > CBS1-1_sae2_2.csv
awk '{print $5}' CBS1-1_ss_spo11_rad51_3.sort.bam.scaled.bed > CBS1-1_sae2_3.csv

*/

function bin_scaledCoverage_dist() {
	[
		"QM6a_rad51",
		"QM6a_spo11rad51",
		"QM6a_sae2",
		"CBS1-1_rad51",
		"CBS1-1_spo11rad51",
		"CBS1-1_sae2",
	].forEach(fname => {
		const mm = [
		];

		const high_marker = 50;

		for (let i = 0; i <= high_marker; ++i) {
			const lower = i / 10;
			// const upper = i / 10 + 1;
			mm[i] = [0, 0, 0];
		}

		for (let i = 0; i < 3; ++i) {
			console.log(fname, i);

			const vec = fs.readFileSync(`bed/${fname}_${(i + 1)}.csv`).toString().split("\n").map(a => Number(a));
			vec.forEach(a => {
				const aa = Math.trunc(a * 10);

				if (aa < high_marker) {
					mm[aa][i] += 1;
				}
				else {
					mm[high_marker][i] += 1;
				}
			});
		}

		fs.writeFileSync(`bed/${fname}.txt`, mm.map((v, i) => [i, ...v].join("\t")).join("\n"));
	});
}

/**
 * @param {string[]} cmp_list
 */
function diff_bin_avgCoverage_dist(cmp_list) {
	if (!fs.existsSync("./diff_bin_avgCoverage_dist")) {
		fs.mkdirSync("./diff_bin_avgCoverage_dist");
	}
	
	cmp_list.forEach(fname => {
		const mm = [
		];

		const high_marker = 50;

		for (let i = -high_marker; i <= high_marker; ++i) {
			const lower = i / 10;
			// const upper = i / 10 + 1;
			mm[high_marker + i] = 0;
		}

		console.log(fname);

		const vec = fs.readFileSync(`${fname}/diff_bins_scaled_${fname}.tsv`).toString().split("\n").map(line => {
			const [
				strChr, strStart, strEnd, strValue,
			] = line.split("\t");

			return Number(strValue);
		});
		vec.forEach(a => {
			const aa = Math.floor(high_marker + a * 10);

			if (0 < aa && aa < (high_marker * 2)) {
				mm[aa] += 1;
			}
			else {
				mm[a >= 0 ? (high_marker + high_marker) : (-high_marker + high_marker)] += 1;
			}
		});

		fs.writeFileSync(`diff_bin_avgCoverage_dist/${fname}.txt`, mm.map((v, i) => [i - high_marker, v].join("\t")).join("\n"));
	});
}

if (cluster.isMaster) {
	(async function () {
		each_cmp([
			"rad51",
			"spo11rad51"
		], "sae2");
		each_cmp([
			"spo11rad51"
		], "spo11sae2");
		
		await main();

		// diff_bin_avgCoverage_dist(Object.values(cmp_bam.cmp).map(a => a.output_tag));
	})();
}
else {
	process.on("message", async function(msg) {
		const args = msg;
		const {
			output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue,
		} = args;

		if (!fs.existsSync(output_tag)) {
			fs.mkdirSync(output_tag);
		}
		process.chdir(output_tag);

		try {
			await _cmp_task(output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue);
		}
		catch (ex) {
			console.error({
				error: ex,
				output_tag,
			});
		}

		process.send(null);
	});
}

async function main() {
	const {
		binSize,// 10;
		min_lg,// = 2;
		max_lg,// = 4;
		min_peak_len,// = 11;// >10
		pValue,// = 0.05;
	} = cmp_bam.param;

	await Promise.all(cmp_bam.cmp.map(async (cmp, cmp_idx) => {
		const { sample_1, sample_2, output_tag } = cmp;
		await run_cmp_task(output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue);
	}));
}

/**
 * @param {number} binSize
 */
async function run_all_BAMScale(binSize) {
	const sample_name_list = enum_all_bam();
	await Promise.all(sample_name_list.map(async (sample) => {
		await scaleSample(sample, binSize, { run: true, });
	}));
}

/**
 * @returns {string[]}
 */
function enum_all_bam() {
	const sample_set = new Set();
	cmp_bam.cmp.forEach((cmp, cmp_idx) => {
		const { sample_1, sample_2, output_tag } = cmp;
		sample_set.add(sample_1);
		sample_set.add(sample_2);
	});
	return [...sample_set.values()];
}

/**
 * @param {string} output_tag
 * @param {string} sample_1
 * @param {string} sample_2
 * @param {number} binSize
 * @param {number} min_lg
 * @param {number} max_lg
 * @param {number} min_peak_len
 * @param {number} pValue
 */
function run_cmp_task(output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue) {
	const task = new Promise(function (resolve) {
		const worker = cluster.fork();

		worker.send({
			output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue,
		});

		worker.on("message", async function() {
			worker.disconnect();// disconnect IPC
			resolve();
		});
	});

	return task;
}

/**
 * @param {string} output_tag
 * @param {string} sample_1
 * @param {string} sample_2
 * @param {number} binSize
 * @param {number} min_lg
 * @param {number} max_lg
 * @param {number} min_peak_len
 * @param {number} pValue
 */
async function _cmp_task(output_tag, sample_1, sample_2, binSize, min_lg, max_lg, min_peak_len, pValue) {
	const {
		chr_nameList,
		chr_lengthList,
	} = samtools.getChrInfo(cmp_bam.sample[sample_1][0]);

	// const fin = (function () {
	// 	for (let i = min_lg; i <= max_lg; ++i) {
	// 		const fin_out = `Diff_peaks_lg${i}_${output_tag}_KStest_p${pValue}.txt`;
	// 		if (!fs.existsSync(fin_out)) {
	// 			return false;
	// 		}
	// 		else {
	// 			console.log("found final result:", fin_out);
	// 			const stat = fs.statSync(fin_out);
	// 			console.log({
	// 				size: stat.size,
	// 				mtime: stat.mtime,
	// 				birthtime: stat.birthtime,
	// 			});
	// 		}
	// 	}
	// 	return true;
	// })();
	// if (fin) {
	// 	return;
	// }

	const [
		sum_mean_1, sum_mean_2
	] = await Promise.all([
		runSampleSumMean(sample_1, binSize, { run: false, }),
		runSampleSumMean(sample_2, binSize, { run: false, }),
	]);

	if (sum_mean_1 == null || sum_mean_2 == null) {
		writeLog(JSON.stringify({
			output_tag,
			sum_mean_1,
			sum_mean_2,
			error: new Error("if (sum_mean_1 == null || sum_mean_2 == null) {").stack,
		}, null, "\t") + "\n");
		return;
	}
	
	// await Promise.all([
	// 	await rep_corr_scaled(sum_mean_1, chr_nameList, sample_1),
	// 	await rep_corr_scaled(sum_mean_2, chr_nameList, sample_2),
	// ]);

	await Diff_bin_calculate(sum_mean_1, sum_mean_2, chr_nameList, min_lg, max_lg, output_tag);
	// output: QM6a_rad51_vs_sae2/diff_bins_scaled_rad51_vs_sae2.tsv

	// cutoff = 2, 3, 4
	const cmd_Diff_peak_cutoff_lg = [
		"python ../src/Diff_peak_cutoff_lg.py",
		output_tag,
		Array.from(Array(max_lg + 1)).map((_, i) => i).slice(min_lg, max_lg + 1).join(",")
	].join(" ");
	const cp_cf_list = await execAsync("Diff_peak_cutoff_lg", cmd_Diff_peak_cutoff_lg, true);

	// A - B > 2
	const tasks = Array.from(new Array(max_lg - min_lg + 1)).map((_, i) => min_lg + i).map(async lg => {
		const cmd_Diff_peaks_frag_pvalue_lg = [
			`Rscript ../src/Diff_peaks_frag_pvalue_lg.R`,
			output_tag,
			sum_mean_1,
			sum_mean_2,
			lg,
			lg,
			min_peak_len,
			pValue,
		].join(" ");
		await execAsync("Diff_peaks_frag_pvalue_lg", cmd_Diff_peaks_frag_pvalue_lg, true);
	});
	await Promise.all(tasks);
}

// `Rscript ../src/rep_corr_scaled.R`

/**
 * 
 * @param {string} sum_mean `${sample}_all.sort.bam.scaled.sum_mean.tsv`
 * @param {string[]} chr_list
 * @param {string} output_name AAA
 */
async function rep_corr_scaled(sum_mean, chr_list, output_name) {
	const cmd =[
		`Rscript ../src/rep_corr_scaled.R`,
		sum_mean,
		chr_list.join(","),
		output_name
	].join(" ");
	await execAsync("Diff_bin_calculate", cmd, true);
}

/**
 * 
 * @param {string} sum_mean_1 `${sample_1}_all.sort.bam.scaled.sum_mean.tsv`
 * @param {string} sum_mean_2 `${sample_2}_all.sort.bam.scaled.sum_mean.tsv`
 * @param {string[]} chr_list
 * @param {number} min_lg
 * @param {number} max_lg
 * @param {string} output_name AAA_vs_BBB
 */
async function Diff_bin_calculate(sum_mean_1, sum_mean_2, chr_list, min_lg, max_lg, output_name) {
	try {
		const cmd =[
			`Rscript ../src/Diff_bin_calculate.R`,
			sum_mean_1, sum_mean_2,
			chr_list.join(","),
			min_lg, max_lg,
			output_name
		].join(" ");
		
		await execAsync("Diff_bin_calculate", cmd, true);

		const cutoff_list = Array.from(Array(max_lg + 1)).map((_, i) => i).slice(min_lg, max_lg + 1);
		
		return {
			all: `diff_bins_scaled_${output_name}.tsv`,
			cutoff_list: cutoff_list.map(cf => `diff_bins_scaled_lg${cf}_${output_name}.bedgraph`),
		};
	}
	catch (ex) {
		writeLog(ex.stack + "\n");
	}
}

/**
 * @param {string} multi_sample
 * @param {number} binSize
 * @param {{ run: boolean }} _ if run == false then skip BAMscale
 * @returns {Promise<string[]>} .bed
 */
async function scaleSample(multi_sample, binSize, { run = false }) {
	const task_sample_1_list = cmp_bam.sample[multi_sample].map(async (bam, bam_idx, bam_list) => {
		const sampleName = (function () {
			if (bam_list.length > 1) {
				return `${multi_sample}_${Number(bam_idx) + 1}`;
			}
			else {
				return `${multi_sample}`;
			}
		})();
		
		const scale_bw = await run_BAMscale(sampleName, bam, "scaled", binSize, { run: run, });
		return await run_bigWigToBed(sampleName, scale_bw, { run: run, });
	});
	const sample_1_list = await Promise.all(task_sample_1_list);
	
	if (sample_1_list.length == 0) {
		console.error("sample_1_list.length = 0,", multi_sample);
	}

	return sample_1_list;
}

/**
 * @param {string} multi_sample
 * @param {number} binSize
 * @param {{ run: boolean }}
 * @returns {Promise<string>} `${multi_sample}_all.sort.bam.scaled.sum_mean.tsv`
 */
async function runSampleSumMean(multi_sample, binSize, { run = false }) {
	try {
		const sum_mean = `${multi_sample}_all.sort.bam.scaled.sum_mean.tsv`;

		const sample_1_list = await scaleSample(multi_sample, binSize, { run: run, });
		
		const sample_1_all = await combineRepeat(multi_sample, sample_1_list);
		if (sample_1_all == null) {
			console.error({
				err: "runSampleSumMean -> combineRepeat",
				multi_sample,
				sample_1_list,
			})
			throw new Error("combineRepeat no return");
		}

		const cmd_sum_mean = `awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$10"\t"$15"\t"$5+$10+$15"\t"($5+$10+$15)/3}' ${sample_1_all} > ${sum_mean}`;
		await execAsync("sum_mean", cmd_sum_mean, false);

		console.log([sum_mean, fs.statSync(sum_mean).size].join("\t"));

		return sum_mean;
	}
	catch (ex) {
		writeLog(ex.stack + "\n");
	}
}

/**
 * @param {string} sample
 * @param {string[]} input_list
 */
async function combineRepeat(sample, input_list) {
	try {
		const out_file = `${sample}_all.sort.bam.scaled.bed`;
		const list = (function () {
			if (input_list.length == 1) {
				return [
					input_list[0],
					input_list[0],
					input_list[0],
				];
			}
			else if (input_list.length == 3) {
				return input_list;
			}
			else {
				const err = new Error("require 1 or 3 sample(s)");
				console.error(err);
				writeLog(err.stack + "\n");
				throw err;
			}
		})();

		// console.log(input_list);
		// console.log(list);

		const cmd = `paste ${list.join(" ")} > ${out_file}`;
		await execAsync("paste", cmd, false);
		return out_file;
	}
	catch (ex) {
		console.error(ex);
		writeLog(ex.stack + "\n");
	}
}

/**
 * @param {string} sampleName
 * @param {string} scale_bw
 * @param {{ run: boolean }}
 * @returns {Promise<string>} out_file
 */
async function run_bigWigToBed(sampleName, scale_bw, { run = false }) {
	try {
		// const out_file = `${sampleName}.sort.bam.scaled.bed`;
		const out_file = scale_bw.replace(/\.bw$/, ".bed");

		if (run) {
			const bedGraph = scale_bw.replace(/\.bw$/, ".bedGraph");
			const wig = scale_bw.replace(/\.bw$/, ".wig");

			await execAsync("bigWigToBed", `bigWigToBedGraph ${scale_bw} ${bedGraph}`, false);
			
			await execAsync("bigWigToBed", `bigWigToWig ${scale_bw} ${wig}`, true);
			await execAsync("bigWigToBed", `wig2bed < ${wig} > ${out_file}`, true);
		}
		
		return out_file;
	}
	catch (ex) {
		console.log("error", sampleName);
		console.log("error cwd", process.cwd());
		console.log(ex);
		writeLog(ex.stack + "\n");
	}
}

/**
 * 
 * @param {string} sample_name
 * @param {string} bam_path
 * @param {string} operation default: scaled
 * @param {number} binSize default: 10
 * @param {{ run: boolean }}
 * @returns {Promise<string>} output file path ${input_file_dir}/${BAMscale args...}/${input_file_name}.scaled.bw
 */
async function run_BAMscale(sample_name, bam_path, operation, binSize, { run = false }) {
	try {
		const input_file_dir = Path.dirname(bam_path);
		const input_file_name = Path.basename(bam_path);
		const outdir = `${input_file_dir}/${sample_name}_bin${binSize}_${operation}`;

		if (run) {
			const bai = `${bam_path}.bai`;
			if (!fs.existsSync(bai)) {
				await execAsync("samtools.index", `samtools index ${bam_path}`, true);
			}

			const cmd = `BAMscale scale -t ${MAX_THREAD} --operation ${operation} --binsize ${binSize} --outdir ${outdir} --bam ${bam_path}`;
			await execAsync("BAMscale", cmd, true);
		}

		return `${outdir}/${input_file_name}.scaled.bw`;
	}
	catch (ex) {
		writeLog(ex.stack + "\n");
	}
}

/**
 * @param {string} content
 */
function writeLog(content) {
	fs.writeFileSync(`run.${TAG_DATE}.${process.pid}.log`, content, { flag: "a" });
}

/**
 * @param {string} taskTag
 * @param {string} cmd
 * @param {boolean} writeLogFile
 * @returns {Promise<{ err?: Error; stdout?: string; stderr?: string; }>}
 */
async function execAsync(taskTag, cmd, writeLogFile) {
	const promise = new Promise((resolve, reject) => {
		writeLog(`${TAG_DATE} [${process.pid}]$ ${cmd}\n`);

		console.log(cmd);

		const proc = child_process.exec(cmd, function (err, stdout, stderr) {
			if (writeLogFile) {
				if (stdout.length) {
					_writeLog("stdout", stdout);
				}
				if (stderr.length) {
					_writeLog("stderr", stderr);
				}
			}
			if (err) {
				const ro = {
					err,
					stdout,
					stderr,
					pid: process.pid,
					wd: process.cwd(),
				};
				writeLog(`${TAG_DATE} [${process.pid}]$ abort.\n${JSON.stringify(ro, null, "\t")}\n`);
				reject(ro);
			}
			else {
				writeLog(`${TAG_DATE} [${process.pid}]$ done.\n`);
				resolve({
					stdout,
					stderr,
				});
			}
		});
		os.setPriority(proc.pid, os.constants.priority.PRIORITY_BELOW_NORMAL);
	});

	return await promise;

	/**
	 * @param {string} className
	 * @param {string} content
	 */
	function _writeLog(className, content) {
		fs.writeFileSync(`${taskTag}.${className}.txt`, `${TAG_DATE} [${process.pid}]$ ${cmd}\n`, { flag: "a" });
		fs.writeFileSync(`${taskTag}.${className}.txt`, content, { flag: "a" });
		fs.writeFileSync(`${taskTag}.${className}.txt`, "\n\n", { flag: "a" });
	}
}

