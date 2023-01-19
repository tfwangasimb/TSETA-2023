// @ts-check

const fs = require("fs");
const Path = require("path");

const process_utils = require("./process_utils.js");
const { GenomeInfo } = require("./genome_info.js");

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;


class MappingOptions {
	constructor() {
		this.max_memory = "16G";
	}

	/**
	 * @param {string} bin_path
	 */
	static setBWABinPath(bin_path) {
		MappingOptions.bwa_bin_path = bin_path;
	}

	/**
	 * @param {string} bin_path
	 */
	static setSAMtoolsBinPath(bin_path) {
		MappingOptions.samtools_bin_path = bin_path;
	}
}
MappingOptions.bwa_bin_path = "";
MappingOptions.samtools_bin_path = "";

/**
 * @param {GenomeInfo} refInfo
 * @param {string} sample_id
 * @param {string} input_R1
 * @param {string} input_R2
 * @param {string} output_dir - dir path
 * @param {MappingOptions} options
 */
async function run_mapping(refInfo, sample_id, input_R1, input_R2, output_dir, options) {
	if (!(options instanceof MappingOptions)) {
		throw new TypeError("run_mapping(..., options:BWAOptions)");
	}
	
	const output_bam = `${output_dir}/${sample_id}.sort.bam`;

	if (!DEBUG && fs.existsSync(output_bam)) {
		console.log("found:", {
			output_bam,
		});
		return null;
	}

	if (!refInfo.file || refInfo.file == "") {
		console.error(refInfo);
		throw new TypeError("refInfo.file");
	}

	await process_utils.exec(`bwa index ${refInfo.file}`, output_dir, sample_id, "bwa_mem");

	const params = [
		`mem`,
		"-t", `${8}`,
		refInfo.file,
		`${input_R1}`,
		`${input_R2}`,
	];
	const cmd = `${MappingOptions.bwa_bin_path} ${params.join(" ")} | ${MappingOptions.samtools_bin_path} sort -@ 8 -l 9 -o ${output_bam}`;
	console.log(cmd);
	await process_utils.exec(cmd, output_dir, sample_id, "bwa_mem");

	const bam_idx_cmd = `${MappingOptions.samtools_bin_path} index ${output_bam}`;
	await process_utils.exec(bam_idx_cmd, output_dir, sample_id, "samtools_index");

	return output_bam;
}


module.exports.MappingOptions = MappingOptions;
module.exports.run_mapping = run_mapping;

