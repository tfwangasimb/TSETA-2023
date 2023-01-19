// @ts-check

const child_process = require("child_process");
const fs = require("fs");
const Path = require("path");

const { loadSetting } = require("./setting.js");
const {
	BlastnCoord,
	parseBlastnResults,
	isCollide, removeOverlap, groupByOverlap,
} = require("./tools/BlastnCoord.js");

const {
	parse_blastn_results,
} = require("./tools/blastn_stream_parser.js");

const VERBOSE = process.argv.indexOf("--verbose") >= 0;

const setting = loadSetting();

const WRITE_VERBOSE_LOG = process.env.WRITE_VERBOSE_LOG;


/**
 * @param {string} cmd
 * @param {boolean} output_stdout
 * @param {boolean} output_stderr
 */
function execAsync(cmd, output_stdout, output_stderr) {
	// return new Promise(function (resolve, reject) {
	// 	let proc = child_process.spawn(cmd);
		
	// 	proc.on("error", function(err) {
	// 		reject(err);
	// 	});

	// 	let stdout_buf_list = [];
	// 	proc.stdout.on("data", function(buffer) {
	// 		if (output_stdout) {
	// 			stdout_buf_list.push(buffer);
	// 		}
	// 	});
		
	// 	let stderr_buf_list = [];
	// 	proc.stderr.on("data", function(buffer) {
	// 		if (output_stderr) {
	// 			stderr_buf_list.push(buffer);
	// 		}
	// 	});
		
	// 	proc.on("exit", function(code, signal) {
	// 		let stdout, stderr;
	// 		if (output_stdout) {
	// 			stdout = Buffer.concat(stdout_buf_list);
	// 		}
	// 		if (output_stderr) {
	// 			stderr = Buffer.concat(stderr_buf_list);
	// 		}
	// 		resolve({
	// 			code, signal, stdout, stderr
	// 		});
	// 	});
	// });

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, {
			maxBuffer: 1024 * 1024 * 128,
		}, function (err, stdout, stderr) {
			if (output_stderr && stderr) {
				console.error(stderr.toString());
			}
		//
			if (err) {
				reject(err);
			}
			else if (output_stdout) {
				resolve(stdout.toString());
			}
			else {
				resolve();
			}
		});
	});
}

/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} [qstart]
 * @param {number} [qend]
 * @param {number} [sstart]
 * @param {number} [send]
 * @param {string} [_task_name]
 * @returns {Promise<string>}
 */
function exec_blastn(query_file, subject_file, qstart, qend, sstart, send, _task_name) {
	let query_loc = Number.isSafeInteger(qstart) && Number.isSafeInteger(qend) ? `-query_loc ${qstart}-${qend}` : "";
	let subject_loc = Number.isSafeInteger(sstart) && Number.isSafeInteger(send) ? `-subject_loc ${sstart}-${send}` : "";
	let cmd = `${setting.blastn_bin} -query ${query_file} ${query_loc} -subject ${subject_file} ${subject_loc} -outfmt 6`;

	_task_name = _task_name || "default";

	if (VERBOSE) {
		console.log(cmd);
	}

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, {
			maxBuffer: 1024 * 1024 * 128,
		}, function (err, stdout, stderr) {
			if (stderr) {
				console.error(cmd);
				console.error(stderr.toString());
			}

			// @ts-ignore
			const tmp_path = global.dataset.tmp_path;

			if (err) {
				console.error(err, { query_file, subject_file, qstart, qend, sstart, send, _task_name });
				
				fs.writeFileSync(`${tmp_path}/ma_util_blastn/ma_util_blastn_${_task_name}.txt`, JSON.stringify(err) + "\n", { flag: "a" });
				reject(err);
			}
			else {
				// resolve({
				// 	stdout: stdout.toString(),
				// 	//stderr: stderr.toString(),
				// });
				let text = stdout.toString();

				if (WRITE_VERBOSE_LOG) {
					fs.writeFileSync(`${tmp_path}/ma_util_blastn/ma_util_blastn_${_task_name}.txt`, cmd + "\n" + text + "\n", { flag: "a" });
				}
				
				resolve(text);
			}
		});
	});
}
/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} [qstart]
 * @param {number} [qend]
 * @param {number} [sstart]
 * @param {number} [send]
 * @param {string} [args]
 * @returns {Promise<string>}
 */
function exec_blastn_Ex(query_file, subject_file, qstart, qend, sstart, send, args) {
	let query_loc = Number.isSafeInteger(qstart) && Number.isSafeInteger(qend) ? `-query_loc ${qstart}-${qend}` : "";
	let subject_loc = Number.isSafeInteger(sstart) && Number.isSafeInteger(send) ? `-subject_loc ${sstart}-${send}` : "";
	let cmd = `${setting.blastn_bin} -query ${query_file} ${query_loc} -subject ${subject_file} ${subject_loc} ${args} -outfmt 6`;

	if (VERBOSE) {
		console.log(cmd);
	}

	return new Promise(function (resolve, reject) {
		child_process.exec(cmd, {
			maxBuffer: 1024 * 1024 * 128,
		}, function (err, stdout, stderr) {
			if (stderr) {
				console.error(cmd);
				console.error(stderr.toString());
			}

			if (err) {
				console.error(err, { query_file, subject_file, qstart, qend, sstart, send });
				reject(err);
			}
			else {
				let text = stdout.toString();

				resolve(text);
			}
		});
	});
}

/**
 * @param {string} query_file
 * @param {string} subject_file
 * @param {number} qstart
 * @param {number} qend
 * @param {number} sstart
 * @param {number} send
 * @param {function(Partial<BlastnCoord>,any[]):boolean} filter
 * @param {any[]} filter_params
 * @param {string} _task_name
 * @returns {Promise<BlastnCoord[]>}
 */
async function blastn_coord(query_file, subject_file, qstart, qend, sstart, send, filter, filter_params, _task_name) {
	if (qstart >= qend || sstart >= send) {
		debugger;
		return null;
	}
	let text = await exec_blastn(query_file, subject_file, qstart, qend, sstart, send, _task_name);

	let _table = parseBlastnResults(text);

	// _table.forEach(row => {
	// 	if (row.send == 2641826) {//2004187
	// 		debugger;
	// 	}
	// });
	
	//_table = _table.sort((a, b) => b.score - a.score);
	_table = _table.sort((a, b) => a.send - b.send);

	// let table = _table.filter(row => {
	// 	return (
	// 		row.qstart < row.qend &&
	// 		row.sstart < row.send &&
	// 		//row.sstart >= next_start
	// 		row.send >= next_start
	// 	);
	// });

	let table = removeOverlap(_table).filter(row => filter(row, filter_params));
	//let table = _table.filter(row => filter(row, filter_params));

	return table;

	// //console.log(table);
	// if (table.length) {
	// 	return table[0];
	// }
	// //console.log(next_start, _table, query_file, subject_file, qstart, qend, sstart, send, next_start);
	// return null;
}


module.exports.execAsync = execAsync;
module.exports.exec_blastn = exec_blastn;
module.exports.exec_blastn_Ex = exec_blastn_Ex;
module.exports.blastn_coord = blastn_coord;

module.exports.BlastnCoord = BlastnCoord;

module.exports.parseBlastnResults = parseBlastnResults;
module.exports.isCollide = isCollide;
module.exports.removeOverlap = removeOverlap;
module.exports.groupByOverlap = groupByOverlap;
module.exports.parse_blastn_results = parse_blastn_results;

