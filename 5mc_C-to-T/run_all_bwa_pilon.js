// @ts-check

// 2021/08/23 change:

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;

const fs = require("fs");
const Path = require("path");

const process_utils = require("./process_utils.js");
const trimmomatic = require("./run_trimmomatic.js");
const mapping = require("./run_mapping.js");
const polish = require("./run_pilon.js");

const { GenomeInfo } = require("./genome_info.js");


const trimmomaticOptions = new trimmomatic.TrimmomaticOptions();
//
trimmomatic.TrimmomaticOptions.setBinPath("/opt/app/Trimmomatic-0.32/trimmomatic-0.32.jar");

const mappingOptions = new mapping.MappingOptions();
mapping.MappingOptions.setBWABinPath("bwa");
mapping.MappingOptions.setSAMtoolsBinPath("samtools");

const pilonOptions = new polish.PilonOptions();
//
polish.PilonOptions.setBinPath("/opt/app/pilon-1.23.jar");

console.error({
	bwa_bin_path: mapping.MappingOptions.bwa_bin_path,
	samtools_bin_path: mapping.MappingOptions.samtools_bin_path,
});

const output_dir = process.argv[3] || ".";
// const OUTPUT_TRIM_PATH = `${output_dir}/trim`;
const OUTPUT_POLISH_PATH = `${output_dir}/polish`;
const OUTPUT_MAPPING_PATH = `${output_dir}/mapping`;

class NGS_MapTo_Ref_dataset_GenomeGroup {
	constructor() {
		this.ref_parental_genome = "";
		this.raw_template_genome_path = "";
		
		this.genome_id = "";
		this.template_genome_path = "";
		/** @type {NGS_MapTo_Ref_dataset_Sample[]} */
		this.samples = [];
	}
}
class NGS_MapTo_Ref_dataset_Sample {
	constructor() {
		this.output_tag = "";

		/** @type {[R1: string[], R2: string[]]} */
		this.reads = [
			[],
			[],
		];
	}
}

/**
 * @type {NGS_MapTo_Ref_dataset_GenomeGroup[]}
 */
const dataset = JSON.parse(fs.readFileSync(process.argv[2]).toString());

async function main() {
	if (!fs.existsSync("./logfiles")) {
		fs.mkdirSync("./logfiles");
	}
	if (!fs.existsSync(OUTPUT_POLISH_PATH)) {
		fs.mkdirSync(OUTPUT_POLISH_PATH);
	}
	if (!fs.existsSync(OUTPUT_MAPPING_PATH)) {
		fs.mkdirSync(OUTPUT_MAPPING_PATH);
	}

	const tasks = dataset.map(async function (gemome_group) {
		if (!gemome_group.template_genome_path) {
			gemome_group.template_genome_path = Path.resolve("./", `${gemome_group.genome_id}.fa`);

			const cmd_det_chr = [
				"node",
				`../src/det_chr.js`,
				`-r ${gemome_group.ref_parental_genome}`,
				`-i ${gemome_group.raw_template_genome_path}`,
				`-o ${gemome_group.template_genome_path}`,
				`--chr_name_template ` + "'ch_${(chrIdx + 1)}'"
			].join(" ");

			await process_utils.exec(cmd_det_chr, "./logfiles", "det_chr", "det_chr");
		}

		if (!fs.existsSync(gemome_group.template_genome_path)) {
			throw new Error("Not found template genome: " + gemome_group.template_genome_path);
		}
		
		await make_bwa_index(gemome_group.template_genome_path);

		const sub_task = gemome_group.samples.map(async function (sample) {
			if (sample.reads[0].length != 1 || sample.reads[1].length != 1) {
				throw new TypeError("no support multiple R1, R2 reads");
			}

			// const trim_results = await trimmomatic.run_trimmomatic(progeny_name, raw_r1, raw_r2, OUTPUT_TRIM_PATH, trimmomaticOptions);
			// const [trim_r1, trim_r2] = [trim_results.output_R1, trim_results.output_R2];
			const [trim_r1, trim_r2] = [sample.reads[0][0], sample.reads[1][0]];

			const ref_info = new GenomeInfo(gemome_group.genome_id, gemome_group.template_genome_path, null);

			try {
				await mapping.run_mapping(ref_info, sample.output_tag, trim_r1, trim_r2, OUTPUT_MAPPING_PATH, mappingOptions);
			}
			catch (ex) {
				console.error("error run_mapping", sample.output_tag, ex);
			}

			try {
				await polish.run_pilon(ref_info, sample.output_tag, OUTPUT_MAPPING_PATH, OUTPUT_POLISH_PATH, pilonOptions);
			}
			catch (ex) {
				console.error("error run_pilon", sample.output_tag, ex);
			}
		});
		await Promise.all(sub_task);
	});
	await Promise.all(tasks);
}


/**
 * @param {string} fasta_path 
 */
function check_bwa_index(fasta_path) {
	const amb = fs.existsSync(`${fasta_path}.amb`);
	const ann = fs.existsSync(`${fasta_path}.ann`);
	const bwt = fs.existsSync(`${fasta_path}.bwt`);
	const pac = fs.existsSync(`${fasta_path}.pac`);
	const sa = fs.existsSync(`${fasta_path}.sa`);
	// console.log({
	// 	amb, ann, bwt, pac, sa,
	// });
	return amb && ann && bwt && pac && sa;
}

async function make_bwa_index(filepath) {
	if (!check_bwa_index(filepath)) {
		const bam_idx_cmd = `${mapping.MappingOptions.bwa_bin_path} index ${filepath}`;
		await process_utils.exec(bam_idx_cmd, "./logfiles", "bwa_index", "bwa_index");
	}
}

main();
