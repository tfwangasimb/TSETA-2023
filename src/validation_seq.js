
// @ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { argv_parse } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { Dataset } = require("./dataset.js");

const argv = argv_parse(process.argv);

const argv_input_path = String(argv["-i"] || "");
const argv_dataset_patha = String(argv["-dataset"] || "");

const argv_find_border = String(argv["--find-border"] || "");

const dataset = Dataset.loadFromFile(argv_dataset_patha);
const genome_info_list = dataset.loadGenomeInfoList();


class ErrInfo {
	constructor() {
		this.start = -1;
		this.end = -1;

		this.min_raw = null;
		this.min_new = null;

		this.full_raw = null;
		this.full_new = null;
		this.old_fa = null;
	}
}

let VERBOSE = false;
if (fs.realpathSync(process.argv[1]) == __filename) {
	VERBOSE = true;
	main();
}
else {
	VERBOSE = process.argv.indexOf("--verbose") >= 0;
}

function main() {
	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		validation_chr(nChr, dataset.tmp_path, false);
	}
}

/**
 * @param {number} nChr
 * @param {string} directory_path
 * @param {boolean} skip_if_not_exist
 */
function validation_chr(nChr, directory_path, skip_if_not_exist) {
	const input_path = `${directory_path}/mafft_ch${nChr}.fa`;
	let chr_has_error = false;
	if (fs.existsSync(input_path)) {
		const input_fasta = readFasta(Path.join(argv_input_path, input_path));

		if (VERBOSE) {
			console.log("input_path", input_path);
		}
			
		const seq_id_list = Object.keys(input_fasta);

		let seq_has_error = validation_seq(seq_id_list, input_fasta, nChr);
		chr_has_error = chr_has_error || seq_has_error;

		if (argv_find_border) {
			find_border(seq_id_list, input_fasta, nChr);
		}
	}
	else if (!skip_if_not_exist) {
		console.log("skip:", input_path);
	}
	
	return chr_has_error;
}

function validation_seq(seq_id_list, fa, nChr) {
	let seq_has_error = false;
	let raw_seq = genome_info_list.map((genome_info, i) => {
		try {
			return genome_info.chr_map[seq_id_list[i]].loadSeq();
		}
		catch (ex) {
			try {
				return genome_info_list[i].chr_list[nChr - 1].loadSeq();
			}
			catch (err) {
				console.error(err);
			}
		}
	});
	
	seq_id_list.forEach((fa_seq_name, idx) => {
		const ma_to_raw_seq = fa[fa_seq_name].replace(/-/g, "");
		if (ma_to_raw_seq.length != raw_seq[idx].length) {
			console.log("ch", nChr, fa_seq_name, "in  len:", ma_to_raw_seq.length);
			console.log("ch", nChr, fa_seq_name, "raw len:", raw_seq[idx].length);
		}
		let s_name = seq_id_list[idx];
		let has_error = false;
		let to_raw_pos = 0;
		for (let i = 0; i < fa[fa_seq_name].length; ++i) {
			if (fa[fa_seq_name][i] != "-") {
				if (fa[fa_seq_name][i] == raw_seq[idx][to_raw_pos]) {
				}
				else {
					const ext = 10;
					console.log(s_name, "error in:", (i + 1), "raw:", to_raw_pos);
					console.log("in ", fa[fa_seq_name].slice(i - ext, i + ext * 2));
					console.log("raw", raw_seq[idx].slice(to_raw_pos - ext, to_raw_pos + ext * 2));
					let mark = [...fa[fa_seq_name].slice(i - ext, i)].fill(".").join("") + "^";
					console.log("mrk", mark);

					has_error = true;
					break;
				}
				++to_raw_pos;
			}
		}
		if (!has_error) {
			seq_has_error = seq_has_error || has_error;

			if (VERBOSE) {
				console.log(s_name, "ok");
			}
		}
	});
	
	return seq_has_error;
}

function find_border(seq_id_list, fa, nChr) {
	console.log("validation seq");

	let seq_list = [
		"", "",
		"", "", "", "",
	];
	seq_id_list.forEach((id, i) => seq_list[i] = fa[id]);

	if (!fs.existsSync("err-ma")) {
		fs.mkdirSync("err-ma");
	}

	/** @type {ErrInfo[]} */
	let err_info = [];

	let count = 0;

	let ref1_pos = 1;
	for (let i = 0; i < seq_list[0].length - 1; ++i) {
		if (!dataset.isIn_rDNA(nChr, i)) {
			if (seq_list[0][i] != "-") {
				check(ref1_pos);
				++ref1_pos;
			}
			else if (ref1_pos - 1) {//indel
				check(ref1_pos - 1);
			}
		}
		function check(ref1_pos) {
			let valid = seq_list.length;
			let valid_2 = seq_list.length;
			if (dataset.isInCentromere(nChr, ref1_pos)) {
				//skip centromere
			}
			else {
				for (let j = 0; j < seq_list.length; ++j) {
					let current = seq_list[j][i];
					let next = seq_list[j][i + 1];
					if (current == "-" && next != "-") {
						--valid;
					}
					else if (current != "-" && next == "-") {
						--valid;
					}
					if (current == "-" && next == "-") {
						--valid_2;
					}
				}
				if (valid == 0 || (valid < seq_list.length && (seq_list.length - valid) == valid_2)) {
					console.log("error at:", i);
					let max_lower = 0;
					let max_upper = 0;
					for (let j = 0; j < seq_list.length; ++j) {
						for (let k = i - 1; k >= 0; --k) {
							let current = seq_list[j][k];
							if (current != "-") {
								max_lower = Math.max(max_lower, i - k);
								break;
							}
						}
						for (let k = i + 1; k < seq_list[j].length; ++k) {
							let current = seq_list[j][k];
							if (current != "-") {
								max_upper = Math.max(max_upper, k - i);
								break;
							}
						}
					}
					const show_more = 0;
					const end_200207 = 1;
					let min = Math.min(max_lower, max_upper);
					let max = Math.max(max_lower, max_upper);
					let err_start = i - max_lower - show_more + 1;
					let err_end = i + max_upper + show_more + end_200207;
					let err_min_start = i - min - show_more + 1;
					let err_min_end = i + min + show_more + end_200207;
					let err_min_fa = {};
					let err_fa = {};
					let old_fa = {};
					for (let j = 0; j < seq_list.length; ++j) {
						if (Math.abs(max_upper - max_lower) < 100) {
							let err_seq = seq_list[j].slice(err_start, err_end);
							console.log("seq" + j, max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, err_seq);
							old_fa[seq_id_list[j]] = err_seq;
						}
						else {
							let err_seq;
							if (min < 100) {
								err_seq = seq_list[j].slice(i - min - show_more + 1, i + min + show_more);
								console.log("seq", max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, err_seq, "...more", Math.abs(max_upper - max_lower));
								old_fa[seq_id_list[j]] = err_seq;
							}
							else {
								console.log("seq", max_lower, max_lower < max_upper ? "<" : ">", max_upper, min, "...more", Math.abs(max_upper - max_lower));
							}
						}
						err_min_fa[seq_id_list[j]] = seq_list[j].slice(err_min_start, err_min_end).replace(/-/g, "");
						err_fa[seq_id_list[j]] = seq_list[j].slice(err_start, err_end).replace(/-/g, "");
					}
					
					let ei = new ErrInfo();
					ei.start = err_min_start;
					ei.end = err_min_end;
					
					// let all_match = seq_id_list.slice(1).every(sname => err_min_fa[seq_id_list[0]] == err_min_fa[sname]);
					// if (all_match) {
					// }
					// else {
						const min_ma_name = `err-ch${nChr}-min-${err_min_start}-${err_min_end}.fa`;
						const min_ma_path = `err-ma/${min_ma_name}`;
						const min_output_path = `err-ma/mafft-${min_ma_name}`;
						
						const full_ma_name = `err-ch${nChr}-full-${err_start}-${err_end}.fa`;
						const full_ma_path = `err-ma/${full_ma_name}`;
						const full_output_path = `err-ma/mafft-${full_ma_name}`;

						ei.min_raw = err_min_fa;
						ei.full_raw = err_fa;
						ei.old_fa = old_fa;

						// saveFasta(min_ma_path, err_min_fa);
						// saveFasta(full_ma_path, err_fa);

						// let min_mafft_cmd = `mafft --quiet --thread ${num_thread} ${algorithm} --maxiterate ${maxiterate} ${min_ma_path} > ${min_output_path}`;
						// let full_mafft_cmd = `mafft --quiet --thread ${num_thread} ${algorithm} --maxiterate ${maxiterate} ${full_ma_path} > ${full_output_path}`;
						// try {
						// 	if (min < 100) {
						// 		child_process.execSync(min_mafft_cmd);
						// 		let ma_min_new = readFasta(min_output_path);
						// 		ei.min_new = ma_min_new;
						// 	}
						// 	else {
						// 		ei.min_new = "...more";
						// 	}
						// 	if (max < 100) {
						// 		child_process.execSync(full_mafft_cmd);
						// 		let ma_full_new = readFasta(full_output_path);
						// 		ei.full_new = ma_full_new;
						// 	}
						// 	else {
						// 		ei.full_new = "...more";
						// 	}
						// }
						// catch (ex) {
						// 	console.error(ex);
						// }
					// }

					err_info.push(ei);

					++count;
				}
			}
		}
	}

	//fs.writeFileSync(`err-ch${nChr}.json`, JSON.stringify(err_info, null, "\t"));

	console.log("count", count);
}


module.exports.validation_chr = validation_chr;

