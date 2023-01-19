// @ts-check

const os = require("os");
const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { ArgvParser } = require("../src/tools/ArgvParser.js");

const cmd_log = `log_cmd_${process.pid}_${makeDateString()}.txt`;
console.log(`[${new Date().toLocaleString()}] write log file: '${Path.resolve(cmd_log)}'`);
function writeln_log(text) {
	fs.writeFileSync(cmd_log, `${text}\n`, { flag: "a" });
}

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

/**
 * @example
 * node pipeline_batch.js ../BS-seq/dim2_rid1_dimrid/Batmeth2_20220923_ddr/Batmeth2_20220923_ddr.json
 * node pipeline_batch.js ../BS-seq/dim2_rid1_dimrid/Batmeth2_20221005_dsiplay-2/data_list.json
 */
async function main() {
	class SampleData {
		/**
		 * @param {string} output_tag
		 * @param {string} refId
		 * @param {string} genome_path
		 * @param {string} R1
		 * @param {string} R2
		 */
		constructor(output_tag, refId, genome_path, R1, R2) {
			this.output_dir = Path.join(process.cwd(), "results", "/");
			this.output_tag = output_tag;
			this.refId = refId;
			this.genome_path = genome_path;
			this.R1 = R1;
			this.R2 = R2;
			if ([this.output_tag, this.refId, this.genome_path, this.R1, this.R2].some(s => typeof s !== "string")) {
				throw new TypeError(`Bad args type: constructor(${output_tag}, ${refId}, ${genome_path}, [...]/${Path.basename(R1)}, [...]/${Path.basename(R2)})`);
			}
		}

		pipline_args() {
			return [
				`--output_dir=${this.output_dir}`,
				`--output_tag=${this.output_tag}`,
				`--refId=${this.refId}`,
				`--genome_path=${this.genome_path}`,
				`--R1=${this.R1}`,
				`--R2=${this.R2}`,
			];
		}
	}

	/** @type {SampleData[]} */
	const data_list = [];
	
	const old_data_list = JSON.parse(fs.readFileSync(process.argv[2]).toString());

	old_data_list.forEach(arr => {
		if (Array.isArray(arr.samples)) {
			const refGroup = arr;
			const samples = refGroup.samples.filter(a => a.skip === true);
			samples.forEach(smp => {
				data_list.push(new SampleData(
					smp.output_tag,
					refGroup.refId,
					refGroup.genome_path,
					smp.reads[0][0],
					smp.reads[1][0],
				));
			});
		}
		else {
			/** @type {SampleData} */
			const od = arr;
			data_list.push(new SampleData(od.output_tag, od.refId, od.genome_path, od.R1, od.R2));
		}
	});

	// return console.log(JSON.stringify(data_list, null, "\t"));
	writeln_log(JSON.stringify(data_list, null, "\t"));
	
	const promise_list = data_list.map(ds => {
		const cwd = Path.join(process.cwd(), ds.output_tag);
		if (!fs.existsSync(cwd)) {
			fs.mkdirSync(cwd);
		}

		const proc = child_process.spawn("node", [
			"./batmeth_pipeline_s.js",
			...ds.pipline_args(),
		], {
			cwd: cwd,
		});

		writeln_log(proc.spawnargs.join(" "));

		proc.stdout.on("data", chunk => writeln_log([ds.output_tag, "stdout", chunk].join("\t") + "\n"));
		proc.stderr.on("data", chunk => writeln_log([ds.output_tag, "stderr", chunk].join("\t") + "\n"));

		return new Promise((resolve, reject) => {
			proc.on("exit", (code, signal) => {
				resolve({
					output_tag: ds.output_tag,
					proc,
					code,
					signal,
					fail: (code != null && code != 0) || signal != null,
				});
			})
		});
	});

	const proc_list = await Promise.all(promise_list);

	writeln_log("fin: " + JSON.stringify(proc_list.filter(a => !a.fail).map(a => a.output_tag), null, "\t"));
	writeln_log("fail: " + JSON.stringify(proc_list.filter(a => a.fail).map(a => a.output_tag), null, "\t"));
}

main();
