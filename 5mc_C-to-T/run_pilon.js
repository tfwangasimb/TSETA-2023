// @ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const process_utils = require("./process_utils.js");
const { GenomeInfo } = require("./genome_info.js");

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;


class PilonOptions {
	constructor() {
		const max_memory_GB = 16;

		/**
		 * @type {number} - memory in GB
		 */
		this.max_memory = max_memory_GB;
	}

	/**
	 * @param {number} number_of_task
	 */
	setPilonMaxMemory(number_of_task) {
		try {
			const max_memory_GB = Math.max(1, Math.trunc((parseInt(child_process.execSync(`free -g | awk '{print $7}' | grep -E "[0-9]+"`).toString().trim()) - 10) / number_of_task));
			this.max_memory = max_memory_GB;
			console.log("pilon max memory:", this.max_memory, "GB");
		}
		catch (ex) {
		}
	}

	/**
	 * @param {string} bin_path
	 */
	static setBinPath(bin_path) {
		PilonOptions.bin_path = bin_path;
	}
}
PilonOptions.bin_path = "";

/**
 * @param {GenomeInfo} refInfo
 * @param {string} sample_id - input prefix; *.sort.bam
 * @param {string} input_dir - dir path
 * @param {string} output_dir - dir path
 * @param {PilonOptions} options
 */
async function run_pilon(refInfo, sample_id, input_dir, output_dir, options) {
	if (!(options instanceof PilonOptions)) {
		throw new TypeError("run_pilon(..., options:PilonOptions)");
	}
	const input_bam = `${input_dir}/${sample_id}.sort.bam`;

	const output_prifix = `${sample_id}`;
	const _output_fasta = `${output_dir}/${output_prifix}.fasta`;
	const _output_vcf = `${output_dir}/${output_prifix}.vcf`;
	if (!DEBUG && fs.existsSync(_output_fasta) && fs.existsSync(_output_vcf)) {
		console.log("found:", {
			_output_fasta,
			_output_vcf,
		});
		return {
			fasta: _output_fasta,
			vcf: _output_vcf,
		};
	}

	const params = [
		`-Xmx${options.max_memory}G`,
		"-jar", PilonOptions.bin_path,
		"--genome", refInfo.file,
		"--frags", input_bam,
		"--output", `${output_prifix}`,
		"--outdir", `${output_dir}`,
		"--vcf",
		"--fix", "snps,indels",
		// "--threads", `${8}`, // don not use this argument, will cause error
	];
	await process_utils.spawn("java", params, output_dir, sample_id, "pilon");

	return {
		fasta: `${output_dir}/${output_prifix}.fasta`,
		vcf: `${output_dir}/${output_prifix}.vcf`,
	}
}


module.exports.PilonOptions = PilonOptions;
module.exports.run_pilon = run_pilon;

