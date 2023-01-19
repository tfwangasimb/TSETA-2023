// @ts-check

const fs = require("fs");
const Path = require("path");

const process_utils = require("./process_utils.js");

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;


class TrimmomaticOptions {
	constructor() {
		this.ILLUMINACLIP = "TruSeq3-PE.fa:2:30:10";
		this.LEADING = 3;
		this.TRAILING = 3;
		this.SLIDINGWINDOW = [4, 20];
		this.MINLEN = 100;
	}
	getParams() {
		return [
			`ILLUMINACLIP:${this.ILLUMINACLIP}`,
			`LEADING:${this.LEADING}`,
			`TRAILING:${this.TRAILING}`,
			`SLIDINGWINDOW:${this.SLIDINGWINDOW.join(":")}`,
			`MINLEN:${this.MINLEN}`,
		];
	}
	toString() {
		return this.getParams().join(" ");
	}

	/**
	 * @param {string} bin_path
	 */
	static setBinPath(bin_path) {
		TrimmomaticOptions.bin_path = bin_path;
	}
}
TrimmomaticOptions.bin_path = "";

/**
 * @param {string} id
 * @param {string} input_R1 - input_forward.fq.gz
 * @param {string} input_R2 - input_reverse.fq.gz
 * @param {string} output_dir - dir path
 * @param {TrimmomaticOptions} options
 */
async function run_trimmomatic(id, input_R1, input_R2, output_dir, options, skip = false) {
	if (!(options instanceof TrimmomaticOptions)) {
		throw new TypeError("run_trimmomatic(..., options:TrimmomaticOptions)");
	}
	const output_R1 = Path.join(output_dir, `${id}_R1_paired.fastq.gz`);
	const output_R2 = Path.join(output_dir, `${id}_R2_paired.fastq.gz`);
	const output_unpaired_R1 = Path.join(output_dir, `${id}_R1_unpaired.fastq.gz`);
	const output_unpaired_R2 = Path.join(output_dir, `${id}_R2_unpaired.fastq.gz`);

	const results = {
		output_R1, output_R2,
		output_unpaired_R1, output_unpaired_R2,
	};

	if (skip) {
		return results;
	}

	if (!DEBUG && fs.existsSync(output_R1) && fs.existsSync(output_R2)) {
		console.log("found:", {
			output_R1,
			output_R2,
		});
		results.found_exist = true;
		return results;
	}
	const params = [
		"-jar", TrimmomaticOptions.bin_path,
		"PE",
		input_R1, input_R2,
		output_R1, output_unpaired_R1,
		output_R2, output_unpaired_R2,
		...options.getParams(),
	];
	if (DEBUG) {
		console.log("DEBUG>", ["java", ...params].join(" "));
	}
	else {
		await process_utils.spawn("java", params, output_dir, id, "trimming");
	}

	return results;
}


module.exports.TrimmomaticOptions = TrimmomaticOptions;
module.exports.run_trimmomatic = run_trimmomatic;

