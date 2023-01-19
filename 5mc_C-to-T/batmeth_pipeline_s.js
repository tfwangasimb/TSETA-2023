// @ts-check

// 2021/08/31 11:17
// 2022/09/23 20220923: upgrade BatMeth2 to @af3df5f

const os = require("os");
const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");
const worker_threads = require("worker_threads");
// const { Worker, workerData, isMainThread, parentPort, threadId } = require("worker_threads");

const process_utils = require("./process_utils.js");
const trimmomatic = require("./run_trimmomatic.js");
const { ArgvParser } = require("../src/tools/ArgvParser.js");

// ssDNA remove unpaired reads
// remove unpaired reads
const RUN_TRIMMOMATIC = false;

const trimmomaticOptions = new trimmomatic.TrimmomaticOptions();
trimmomatic.TrimmomaticOptions.setBinPath("/opt/app/Trimmomatic-0.32/trimmomatic-0.32.jar");

const program_bin_root = "./BatMeth2-af3df5f/bin/";
const MAX_THREAD = 8;
const g_ReDo = false;

const OUTPUT_TRIM_PATH = "./trim";

class NGS_SampleData {
	output_dir = "";

	/**
	 * @example [refId, mutName, day, seqType].join("_")
	 */
	output_tag = "";

	skip = false;

	refId = "";
	
	genome_path = "";

	/** @type {{ R1: string; R2: string; }[]} */
	reads = [];
}

const cmd_log = `log_cmd_${process.pid}_${makeDateString()}.txt`;
console.log(`[${new Date().toLocaleString()}] write log file: ${cmd_log}`);
function writeln_log(text) {
	fs.writeFileSync(cmd_log, `${text}\n`, { flag: "a" });
}

main();

function makeDateString() {
	const t = new Date();

	const d = [
		t.getFullYear().toFixed(0).slice(2),
		(t.getMonth() + 1).toFixed(0).padStart(2, "0"),
		t.getDate().toFixed(0).padStart(2, "0"),
	].join("");

	const s = [
		t.getHours().toFixed(0).padStart(2, "0"),
		t.getMinutes().toFixed(0).padStart(2, "0"),
		t.getSeconds().toFixed(0).padStart(2, "0"),
	].join("");

	return d + "_" + s;
}

async function main() {
	process.on("uncaughtException", err => {
		writeln_log(err.stack);
	});

	// /** @type {NGS_SampleData} */
	// const sample_data = JSON.parse(fs.readFileSync(process.argv[2]).toString());

	const argv = new ArgvParser(process.argv);

	const sample_data = new NGS_SampleData();
	sample_data.output_dir = argv.get(/--output_dir=(.*)/, (arg, args) => args[0]) + "";
	sample_data.output_tag = argv.get(/--output_tag=(.*)/, (arg, args) => args[0]) + "";//[refId, mutName, day, seqType].join("_")
	// sample_data.skip = !!argv.get(/--skip=(.*)/, (arg, args) => args[0]);
	sample_data.refId = argv.get(/--refId=(.*)/, (arg, args) => args[0]) + "";
	sample_data.genome_path = argv.get(/--genome_path=(.*)/, (arg, args) => args[0]) + "";
	sample_data.reads = [{ R1: "", R2: "", }];
	sample_data.reads[0].R1 = argv.get(/--R1=(.*)/, (arg, args) => args[0]) + "";
	sample_data.reads[0].R2 = argv.get(/--R2=(.*)/, (arg, args) => args[0]) + "";

	if (![
		sample_data.output_dir,
		sample_data.output_tag,
		sample_data.refId,
		sample_data.genome_path,
		sample_data.reads[0],
		sample_data.reads[0].R1,
		sample_data.reads[0].R2,
	].every(a => !!a)){
		throw new Error("argv");
	}
	
	// let dirty = false;

	if (RUN_TRIMMOMATIC) {
		if (!fs.existsSync(OUTPUT_TRIM_PATH)) {
			fs.mkdirSync(OUTPUT_TRIM_PATH);
		}
	}

	// if (!fs.existsSync(sample_data.output_dir)) {
	// 	fs.mkdirSync(sample_data.output_dir);
	// }
	sample_data.output_dir = fs.realpathSync(sample_data.output_dir);

	const genome_path = Path.basename(sample_data.genome_path);
	fs.copyFileSync(sample_data.genome_path, genome_path);

	sample_data.genome_path = genome_path;

	await BatMeth21_index(genome_path);

	if (!fs.existsSync(`${genome_path}.fai`)) {
		const { cmd } = await process_utils.exec(`samtools faidx ${genome_path}`, ".", "samtools", "faidx");
		writeln_log(cmd);
	}
	else {
		writeln_log(`skip: ${genome_path}.fai`);
	}

	if (!fs.existsSync(`${genome_path}.len`)) {
		const cmd = `cat ${genome_path}.fai  | awk '{print $1 "\t" $2}' > ${genome_path}.len`;
		child_process.execSync(cmd);
		writeln_log(cmd);
	}
	else {
		writeln_log(`skip: ${genome_path}.len`);
	}

	if (sample_data.skip === true) {
		return console.log("skip:", sample_data.skip, sample_data.output_tag);
	}

	// const sam = `${sample_data..output_tag}.sam`;
	const bam = `${sample_data.output_tag}.sort.bam`;
	const batmeth_file = `${sample_data.output_tag}.methratio.txt`;

	if (sample_data.reads.length != 1) {
		throw new Error(`too more reads file, sample_data.reads.length = ${sample_data.reads.length}`);// deprecated
	}

	if (isFileValid(bam) == false) {
		/** @see {@link https://github.com/GuoliangLi-HZAU/BatMeth2/issues/32#issuecomment-1123409625} */
		const cmd_align = `${program_bin_root}BatMeth2 align -g ${genome_path} -1 ${sample_data.reads[0].R1} -2 ${sample_data.reads[0].R2} -o ${sample_data.output_tag} -p ${MAX_THREAD}`;
		writeln_log(cmd_align);

		await process_utils.exec(cmd_align, ".", "BatMeth2_pipel_align", sample_data.output_tag);
	}
	
	if (isFileValid(bam)) {
		// await BatMeth2_calmeth();
		await BatMeth21_calmeth();
	}
	else {
		throw new Error(`invalid bam: bam`);
	}

	// async function BatMeth2_calmeth() {
	// 	await BatMeth2_index(genome_path);
	// 	if (isFileValid(batmeth_file) == false) {
	// 		const calmeth = `/opt/app/_BatMeth2/bin/calmeth`;
	// 		// const calmeth = `${program_bin_root}calmeth`;
	// 		const cmd_calmeth = `${calmeth} -g ${genome_path} -n 0.4 -b ${bam} -m ./${sample_data.output_tag}`;
	// 		/**
	// 		 * -n: Number of mismatches, default 0.06 percentage of read length. [0-1]
	// 		 */
	// 		const { cmd } = await process_utils.exec(cmd_calmeth, ".", "BatMeth2_calmeth", sample_data.output_tag);// *.log.txt
	// 		writeln_log(cmd);
	// 	}
	// 	else {
	// 		isFileValid(batmeth_file);
	// 	}
	// }
	async function BatMeth21_calmeth() {
		await BatMeth21_index(genome_path);
		
		if (isFileValid(batmeth_file) == false) {
			const calmeth = `${program_bin_root}calmeth`;
			const cmd_calmeth = `${calmeth} -g ${genome_path} -n 0.4 -b ${bam} -m ./${sample_data.output_tag}`;

			/**
			 * -n: Number of mismatches, default 0.06 percentage of read length. [0-1]
			 */
			const { cmd } = await process_utils.exec(cmd_calmeth, ".", "BatMeth2_calmeth", sample_data.output_tag);// *.log.txt
			writeln_log(cmd);
		}
		else {
			isFileValid(batmeth_file);
		}
	}

	if (isFileValid(bam)) {// is results valid ?
		if (isFileValid(batmeth_file)) {// is results valid ?
			await Promise.all([
				(async () => {
					if (isFileValid(`${bam}.bai`) == false) {
						await execAsync("samtools", [
							"index", bam
						]);
					}
					else {
						console.log("found", `${bam}.bai`);
					}
				})(),
				call_coverage_to_uint32array(bam, sample_data, { re_do: g_ReDo, }),
				call_convert_sample_to_float32array(batmeth_file, sample_data, { re_do: g_ReDo, }),
			]);

			// bam
			const meth_valid = validate_calmeth(`${sample_data.output_tag}.methlog.txt`);
			if (meth_valid) {
				console.log("done methratio.float32 and depth.uint32:", batmeth_file);
			}
			else {
				throw new Error("calmeth: maybe some error: " + batmeth_file);
			}
		}
		else {
			throw new Error("err not found BatMeth2 result: " + sample_data.output_tag + "; " + batmeth_file);
		}
	}
	else {
		throw new Error(`not found bam: ${bam}`);
	}

	console.log("exit");
}

/**
 * @param {string} output_tag
 * @param {string} raw_r1
 * @param {string} raw_r2
 * @returns {Promise<[string, string]>}
 */
async function run_trim(output_tag, raw_r1, raw_r2) {
	if (RUN_TRIMMOMATIC) {
		// let trim_results;
		// try {
			const trim_results = await trimmomatic.run_trimmomatic(output_tag, raw_r1, raw_r2, OUTPUT_TRIM_PATH, trimmomaticOptions);
			return [trim_results.output_R1, trim_results.output_R2];
		// } catch (ex) {
		// 	console.error(ex);
		// 	writeln_log(`err trimmomatic: ${output_tag} <- genome[${data_idx}].sample[${sample_idx}] ${ex.stack}`);
		// 	// throw ["err", "trimmomatic", data_idx, sample_idx, output_tag].join(" ");
		// }
	}
	else {
		return [raw_r1, raw_r2];
	}
}

async function BatMeth2_index(genome_path) {
	const fa_ann_location = `${genome_path}.ann.location`;
	// const fa_bin = `${genome_path}.bin`;
	// if (isFileValid(`${genome_path}.sort`) == false || isFileValid(fa_bin) == false) {
	if (isFileValid(fa_ann_location) == false) {
		try {
			await process_utils.exec(`/opt/app/_BatMeth2/bin/BatMeth2 build_index ${genome_path}`, ".", "BatMeth2", "build_index");
		}
		catch (ex) {
			console.error(ex);
			throw ["err", "BatMeth2", "build_index"].join(" ");
		}
	}
}

/**
 * @param {string} genome_path
 */
async function BatMeth21_index(genome_path) {
	if (hasGenomeSeqIndex(genome_path)) {
	}
	else {
		try {
			await process_utils.exec(`${program_bin_root}BatMeth2 index -g ${genome_path}`, ".", "BatMeth2", "build_index");
		} catch (ex) {
			console.error(ex);
			throw ["err", "BatMeth2", "build_index"].join(" ");
		}
	}
}

/**
 * @param {string} fa
 */
function hasGenomeSeqIndex(fa) {
	return [
		`${fa}.batmeth2.fa`,
		`${fa}.batmeth2.fa.amb`,
		`${fa}.batmeth2.fa.ann`,
		`${fa}.batmeth2.fa.bwt`,
		`${fa}.batmeth2.fa.pac`,
		`${fa}.batmeth2.fa.sa`,
		`${fa}.bin`,
		`${fa}.fai`,
		// `${fa}.len`,
	].every(file_path => {
		const flag = isFileValid(file_path);
		if (!flag) {
			console.log("hasGenomeSeqIndex >> 404:", file_path);
		}
		return flag;
	});
}

/** @param {string} file_path */
function isFileValid(file_path) {
	if (fs.existsSync(file_path)) {
		const stat = fs.statSync(file_path);
		if (stat.isDirectory() == false && stat.size > 0) {
			return true;
		}
		else {
			console.log({
				path: file_path,
				full_path: Path.resolve(file_path),

				mode: stat.mode,
				size: stat.size,
				atime: stat.atime,
				ctime: stat.ctime,
				birthtime: stat.birthtime,

				// isFile: stat.isFile(),
				isDirectory: stat.isDirectory(),
				// isSymbolicLink: stat.isSymbolicLink(),
				// isBlockDevice: stat.isBlockDevice(),
				// isCharacterDevice: stat.isCharacterDevice(),
				// isFIFO: stat.isFIFO(),
				// isSocket: stat.isSocket(),
				stack: new Error().stack,
			});
		}
	}
	return false;
}

/**
 * @param {string} input_file
 * @param {NGS_SampleData} sample_data
 */
async function call_coverage_to_uint32array(input_file, sample_data, { re_do = false, }) {
	const output_file = Path.join(sample_data.output_dir, `${sample_data.output_tag}.depth.uint32`);
	if (re_do == false && fs.existsSync(output_file)) {
		console.log("found", output_file);
		return;
	}
	else {
		if (re_do) {
			console.log("redo:", output_file);
		}
		return await execAsync("node", [
			"../tools/uint32array_coverage.js",
			input_file,
			output_file,
		], [
			`${sample_data.output_tag}.depth.stdout.txt`,
			`${sample_data.output_tag}.depth.stderr.txt`,
		]);
	}
}

/**
 * @param {string} input_file
 * @param {NGS_SampleData} sample_data
 */
async function call_convert_sample_to_float32array(input_file, sample_data, { re_do = false, }) {
	const output_file = Path.join(sample_data.output_dir, `${sample_data.output_tag}.methratio.float32`);
	if (re_do == false && fs.existsSync(output_file)) {
		console.log("found", output_file);
		return;
	}
	else {
		if (re_do) {
			console.log("redo:", output_file);
		}
		return await execAsync("node", [
			`${__dirname}/methratio_to_float32array_v2.js`,
			`--input=${input_file}`,
			`--genome=${sample_data.genome_path}`,
			`--sample_name=${sample_data.output_tag}`,
			`--output-dir=${sample_data.output_dir}`,
		], [
			`${sample_data.output_tag}.methratio.log.txt`,
			`${sample_data.output_tag}.methratio.log.txt`,
		]);
	}
}

/**
 * @param {any} value
 * @returns {boolean} !Number.isNaN(value) && Number.isFinite(value)
 */
 function validateNumber(value) {
	return !Number.isNaN(value) && Number.isFinite(value);
}

/**
 * @param {any} value
 * @returns {boolean} !Number.isNaN(value) && Number.isFinite(value) && Number.isSafeInteger(value)
 */
function validateInteger(value) {
	return !Number.isNaN(value) && Number.isFinite(value) && Number.isSafeInteger(value);
}

/**
 * @see {@link ./validate_calmeth.js}
 * @param {string} log_file
 */
function validate_calmeth(log_file) {
	const str_log_text = fs.readFileSync(log_file).toString();
	return str_log_text.split("\n").map(line => {
		const [k, v] = line.split("\t").map(a => a.trim());
		switch (k) {
			case "Raw count of Met_C in CG:":
			case "Raw count of Non_Met_C in CG:":
			case "Raw count of Met_C in CHG:":
			case "Raw count of Non_Met_C in CHG:":
			case "Raw count of Met_C in CHH:":
			case "Raw count of Non_Met_C in CHH:":
				let numVal = parseInt(v, 10);
				// console.log(k, numVal);
				return validateInteger(numVal);
			case "[CpG]":
			case "[mC]":
				let numArr = [...v.matchAll(/\d+/g)].map(a => parseInt(a[0], 10));
				// let [M, Mh, H, hU, U] = numArr;
				// console.log(k, numArr);
				return numArr.every(a => validateInteger(a));
			case "mC/(C+T)":
			case "mCG/(CG+TG)":
			case "mCHG/(CHG+THG)":
			case "mCHH/(CHH+THH)":
				let va = v.match(/{(.*) \/ (.*)} = (.*)%/);
				// console.log([k, v, va]);
				return va && validateInteger(parseInt(va[1], 10)) && validateInteger(parseInt(va[2], 10)) && validateNumber(parseFloat(va[3]));
			// default: 
			//   console.log(k);
			//   break;
		}
		return k ? k : true;// skip this line
	});
}

/**
 * @param {string} cmd
 * @param  {string[]} args
 * @param  {[string|null, string|null]} logfile
 */
async function execAsync(cmd, args, logfile = [null, null]) {
	const task = new Promise((resolve, reject) => {
		console.log(cmd, ...args);

		const proc = child_process.spawn(cmd, args, {
				stdio: [
					"ignore",
					logfile && logfile[0] ? "pipe" : "ignore",
					logfile && logfile[1] ? "pipe" : "ignore",
				],
			}
		);

		if (logfile && logfile[0]) {
			proc.stdout?.pipe?.(fs.createWriteStream(logfile[0]));
		}
		if (logfile && logfile[1]) {
			proc.stderr?.pipe?.(fs.createWriteStream(logfile[1]));
		}

		// os.setPriority(proc.pid, os.constants.priority.PRIORITY_BELOW_NORMAL);

		proc.on("exit", (code, signal) => {
			if (code || signal) {
				reject({
					"err": cmd, code, signal,
				});
			}
			else {
				resolve({
					cmd, code, signal,
				});
			}
		});
	});

	return await task;
}
