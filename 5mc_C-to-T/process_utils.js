
// @ts-check

const os = require("os");
const fs = require("fs");
const child_process = require("child_process");

const DEBUG = process.argv.includes("--DEBUG") || process.env.DEBUG;


/**
 * @param {string} cmd
 * @param {string} output_dir
 * @param {string} output_prifix
 * @param {string} tag
 * @returns {Promise<{ code?: string; signal?: string; cmd: string; }>}
 */
async function exec(cmd, output_dir, output_prifix, tag) {
	const pid = process.pid;
	console.log(pid, "$", cmd);

	if (!DEBUG) {
		const promise = new Promise(function (resolve, reject) {
			const date = new Date();
			const key = date.toLocaleString();

			fs.writeFileSync(`${output_dir}/log_${pid}_${output_prifix}_${tag}.stderr.txt`, `${key}\n${pid} $ ${cmd}\n`, { flag: "a" });
			fs.writeFileSync(`${output_dir}/log_${pid}_${output_prifix}_${tag}.stdout.txt`, `${key}\n${pid} $ ${cmd}\n`, { flag: "a" });
			
			const proc = child_process.exec(cmd);

			os.setPriority(proc.pid, os.constants.priority.PRIORITY_BELOW_NORMAL);

			proc.on("error", function (err) {
				fs.writeFileSync(`${output_dir}/log_${pid}_${tag}_status.txt`, `${output_prifix}\t${err.name}\t"${err.stack}"\t${key}\t${cmd}\n`, { flag: "a" });
			});

			proc.stderr.on("data", function (chunk) {
				fs.writeFileSync(`${output_dir}/log_${pid}_${output_prifix}_${tag}.stderr.txt`, chunk, { flag: "a" });
			});
			proc.stdout.on("data", function (chunk) {
				fs.writeFileSync(`${output_dir}/log_${pid}_${output_prifix}_${tag}.stdout.txt`, chunk, { flag: "a" });
			});
			proc.once("exit", function (code, signal) {
				fs.writeFileSync(`${output_dir}/log_${pid}_${tag}_status.txt`, `${output_prifix}\t${code}\t${signal}\t${key}\t${cmd}\n`, { flag: "a" });
				resolve({
					code, signal,
					cmd,
				});
			});
		});

		return await promise;
	}
	
	return {
		cmd,
	};
}

/**
 * @param {string} cmd
 * @param {string[]} params
 * @param {string} output_dir
 * @param {string} output_prifix
 * @param {string} tag
 */
async function spawn(cmd, params, output_dir, output_prifix, tag) {
	const full_cmd = [cmd, ...params].join(" ");
	//
	console.log(full_cmd);

	if (!DEBUG) {
		const promise = new Promise(function (resolve, reject) {
			const date = new Date();
			const time_1 = date.toUTCString();

			fs.writeFileSync(`${output_dir}/log_${output_prifix}_${tag}.stderr.txt`, time_1 + "\n", { flag: "a" });
			fs.writeFileSync(`${output_dir}/log_${output_prifix}_${tag}.stdout.txt`, time_1 + "\n", { flag: "a" });
			
			const proc = child_process.spawn(cmd, params);

			os.setPriority(proc.pid, os.constants.priority.PRIORITY_BELOW_NORMAL);

			proc.on("error", function (err) {
				fs.writeFileSync(`${output_dir}/log_${tag}_status.txt`, `${output_prifix}\t${err.name}\t"${err.stack}"\t${time_1}\t${full_cmd}\n`, { flag: "a" });
			});

			proc.stderr.on("data", function (chunk) {
				fs.writeFileSync(`${output_dir}/log_${output_prifix}_${tag}.stderr.txt`, chunk, { flag: "a" });
			});
			proc.stdout.on("data", function (chunk) {
				fs.writeFileSync(`${output_dir}/log_${output_prifix}_${tag}.stdout.txt`, chunk, { flag: "a" });
			});
			proc.once("exit", function (code, signal) {
				const date_2 = new Date();
				const time_2 = date_2.toUTCString();
				fs.writeFileSync(`${output_dir}/log_${tag}_status.txt`, `${output_prifix}\t${code}\t${signal}\t${time_1}\t${time_2}\t${full_cmd}\n`, { flag: "a" });
				resolve({
					code, signal,
					params,
				});
			});
		});

		await promise;
	}
}


module.exports.exec = exec;
module.exports.spawn = spawn;

