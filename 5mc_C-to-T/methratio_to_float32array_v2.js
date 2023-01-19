// @ts-check

const fs = require("fs");
const Path = require("path");
const { readFasta, reverseComplement } = require("../src/fasta_util.js");

/**
 * @param {RegExp} regexp
 */
function get_arg(regexp) {
	const m = process.argv.map(a => a.match(regexp)).filter(a => a)[0];
	if (m) {
		return m[1]
	}
}

const input_file = get_arg(/--input=(.*)/);
const genome_path = get_arg(/--genome=(.*)/);
const sample_name = get_arg(/--sample_name=(.*)/);
const output_results_dir = get_arg(/--output-dir=(.*)/);

if (input_file && genome_path && sample_name && output_results_dir) {
	_convert_sample_to_float32array(input_file, genome_path, sample_name, output_results_dir);
}
else {
	console.error("process.argv", process.argv);
	throw new Error("process.argv");
}

/**
 * @param {string} input_file
 * @param {string} genome_path
 * @param {string} sample_name
 * @param {string} output_results_dir
 */
async function _convert_sample_to_float32array(input_file, genome_path, sample_name, output_results_dir) {
	if (!fs.existsSync(output_results_dir)) {
		fs.mkdirSync(output_results_dir);
	}
	const output_file = Path.join(output_results_dir, `${sample_name}.methratio.float32`);
	// if (fs.existsSync(output_file) && fs.statSync(output_file).size > 0) {
	// 	return;
	// }
	console.log(`fs.existsSync("${output_file}") ==`, fs.existsSync(output_file));

	/** @type {Float32Array[]} */
	const output_chr_map = [];

	/** @type {number[]} */
	const mc_cnt_map = [];

	/** @type {number[]} */
	const non_mc_cnt_map = [];
	
	// const len_info = fs.readFileSync(`${genome_path}.len`).toString().trim().split("\n").map(a => {
	// 	const [name, len] = a.trim().split("\t");
	// 	return {
	// 		name,
	// 		len: Number(len),
	// 	};
	// });

	const fa = readFasta(genome_path);
	const chrInfoMap = Object.fromEntries(Object.entries(fa).map(([name, seq], chrIdx) => {
		// const seq = fa[name];
		const len = seq.length;
		console.log("ch", chrIdx + 1, name, len);
		return [
			name,
			{
				name: name,
				chrIdx: chrIdx,
				len: Number(len),
				seq: seq,
			},
		];
	}));

	Object.keys(fa).forEach((strChr, chrIdx) => {
		const arr = new Float32Array(fa[strChr].length);
		arr.fill(NaN);
		output_chr_map[chrIdx] = arr;

		mc_cnt_map[chrIdx] = 0;
		non_mc_cnt_map[chrIdx] = 0;
	});
	
	await readTableFromStream(input_file, function rowParser(cols, index) {
		// skip comment
		if (cols[0].startsWith("#")) {
			console.log(input_file, cols[0]);
			return;
		}

		const {
			[1 - 1]: strChr,
			[2 - 1]: strPos,
			[3 - 1]: strand,
			// [4 - 1]: context,// "CG" | "CHG" | "CHH"
			[5 - 1]: str_C_count,
			[6 - 1]: str_CT_count,
			// [7 - 1]: strMethRatio
			// [8 - 1]: str_eff_CT_count, // fwd or rev
			// [9 - 1]: str_rev_G_count,
			// [10 - 1]: str_rev_GA_count,
			// [11 - 1]: str_MethContext,// level: "M" | "Mh" | "H" | "hU" | "U"
			[12 - 1]: str_5context, // seq.slice()
		} = cols;

		// 	"12345".slice((3 - 1) - 2, 3 + 0); // 123
		// 	"12345".slice((3 - 1) - 0, 3 + 2); //   345
		// 	"12345".slice((3 - 1) - 2, 3 + 2); // 12345

		const nPos = Number(strPos.trim()) - 1;
		// const nMethRatio = Number(strMethRatio);
		const nC_count = Number(str_C_count.trim());
		const nCT_count = Number(str_CT_count.trim());
		const nMethRatio = nC_count / nCT_count;

		const {
			chrIdx,
			seq: chrSeq,
		} = chrInfoMap[strChr];

		if (strChr == "ch6_QM6a_NP43_pilon" && nPos >= 1366205 && nPos <= 1366312) {
			console.log(strChr, strPos, strand, str_5context, [
				chrSeq[nPos - 2],
				chrSeq[nPos - 1],
				chrSeq[nPos + 0],
				chrSeq[nPos + 1],
				chrSeq[nPos + 2],
			].join(""));
		}

		const mathArr = output_chr_map[chrIdx];
		if (!mathArr) {
			throw new Error(`not found chr '${strChr}' in ${JSON.stringify(output_chr_map)}`);
		}

		if (Number.isFinite(nPos) && !Number.isNaN(nPos) && Number.isSafeInteger(nPos) &&
			Number.isFinite(nC_count) && !Number.isNaN(nC_count) &&
			Number.isFinite(nCT_count) && !Number.isNaN(nCT_count)
		) {
			if (strand == "+" && (
					(chrSeq[nPos - 2] != str_5context[0] && str_5context[0] != "N") ||
					(chrSeq[nPos - 1] != str_5context[1] && str_5context[1] != "N") ||
					(chrSeq[nPos + 0] != str_5context[2] && str_5context[2] != "N") ||
					(chrSeq[nPos + 1] != str_5context[3] && str_5context[3] != "N") ||
					(chrSeq[nPos + 2] != str_5context[4] && str_5context[4] != "N")
			)) {
				console.error({
					err: "diff seq", input_file, index, strChr, strPos, strand,
					// strMethRatio,
					str_C_count,
					str_CT_count,
					nPos: nPos,
					[`mathArr[${nPos}]`]: mathArr[nPos],
					ref: chrSeq.slice(nPos - 2, nPos + 3),
					m5: str_5context,
				});
				throw new Error("diff seq");
			}
			else if (strand == "-") {
				const [r1, r2, r3, r4, r5] = chrSeq.slice(nPos - 2, nPos + 3);
				const [b1, b2, b3, b4, b5] = reverseComplement(str_5context);
				if (
					(r1 != b1 && b1 != "N") ||
					(r2 != b2 && b2 != "N") ||
					(r3 != b3 && b3 != "N") ||
					(r4 != b4 && b4 != "N") ||
					(r5 != b5 && b5 != "N")
				) {
					console.error({
						err: "diff seq", input_file, index, strChr, strPos, strand,
						// strMethRatio,
						str_C_count,
						str_CT_count,
						nPos: nPos,
						[`mathArr[${nPos}]`]: mathArr[nPos],
						ref: chrSeq.slice(nPos - 2, nPos + 3),
						m5: reverseComplement(str_5context),
					});
					throw new Error("diff seq");
				}
			}

			if (!Number.isNaN(mathArr[nPos])) {
				console.error({
					err: "dup pos", input_file, index, strChr, strPos, strand,
					// strMethRatio,
					str_C_count,
					str_CT_count,
					nPos: nPos,
					[`mathArr[${nPos}]`]: mathArr[nPos],
				});
				if (nMethRatio <= mathArr[nPos]) {
					return;
				}
			}

			if (nMethRatio == 0) {
				non_mc_cnt_map[chrIdx] += 1;
			}
			else {
				mc_cnt_map[chrIdx] += 1;
			}

			if (strand == "+") {
				mathArr[nPos] = nMethRatio; // 1 * 0 => 0
				// ui4a[0] >>>  31 => 0
			}
			else if (strand == "-") {
				mathArr[nPos] = -nMethRatio; // -1 * 0 => -0
				// ui4a[0] >>>  31 => 1
			}
			else {
				throw new Error();
			}
		}
		else {
			console.error({
				input_file, index, strChr, strPos, strand,
				// strMethRatio,
				str_C_count,
				str_CT_count,
			});
			throw "error table format: cols[" + index + "]";
		}
	});
	if (globalThis.gc) {
		globalThis.gc();
		//--expose-gc
	}

	// for (let i = -5; i < 10; ++i) {
	// 	console.log("ChVI_QM6a", `[341675 + ${i}]`, output_chr_map["ChVI_QM6a"][341675 + i]);
	// }
		
	console.log("loaded:", input_file);

	// console.log("begin write file", output_file, Object.keys(output_chr_map), Object.values(output_chr_map).map(a => a.length));

	const ws = fs.createWriteStream(output_file);
	ws.once("close", function () {
		console.log("closed:", output_file);
	});

	const resultInfo_list = Object.values(chrInfoMap).map(({ name: strChr, len }, chrIdx) => {
		ws.write(Buffer.from(output_chr_map[chrIdx].buffer));
		
		console.log("write:", output_file, strChr);
		
		const chr_len = output_chr_map[chrIdx].length;

		const cnt_data = output_chr_map[chrIdx].reduce((acc, v) => acc + (Number.isNaN(v) ? 0 : 1), 0);
		console.log("cnt_data", strChr, cnt_data, (mc_cnt_map[chrIdx] + non_mc_cnt_map[chrIdx]));
		
		if (chr_len == len) {
			non_mc_cnt_map[chrIdx];
			mc_cnt_map[chrIdx];

			const info = {
				input_file,
				name: strChr,

				chr_idx: chrIdx,
				chr_len,

				mc_cnt: mc_cnt_map[chrIdx],
				non_mc_cnt: non_mc_cnt_map[chrIdx],
				
				mc_rate: (mc_cnt_map[chrIdx] / chr_len).toFixed(2),
				non_mc_rate: (non_mc_cnt_map[chrIdx] / chr_len).toFixed(2),
				
				data_rate: ((mc_cnt_map[chrIdx] + non_mc_cnt_map[chrIdx]) / chr_len).toFixed(2),

				mc_cnt_map,
				non_mc_cnt_map,
			};
			console.table(info);

			// fs.writeFileSync(`${sample_name}.info.json`, JSON.stringify(info, null, "\t"));
			return info;
		}
		else {
			throw new Error(`chr.len: ${chr_len} != ${len}`);
		}
	});
	ws.end();

	fs.writeFileSync(`${sample_name}.info.json`, JSON.stringify(resultInfo_list, null, "\t"));
	console.log(`save info: ${sample_name}.info.json`);
}

/**
 * @template T
 * @param {string} filePath
 * @param {(cols:string[],index:number) => T} rowParser
 */
 async function readTableFromStream(filePath, rowParser) {
	const stream = fs.createReadStream(filePath, {
		highWaterMark: 128 * 1024 * 1024,// 128Mib
	});

	const promise = new Promise(function (resolve, reject) {
		let last_str = "";

		stream.on("error", function (err) {
			reject(err);
		});

		stream.on("data", async function (chunk) {
			const buf = chunk.toString();
			
			const str_1 = buf.slice(0, buf.lastIndexOf("\n"));
	
			const str_part = last_str + str_1;
			
			last_str = buf.slice(buf.lastIndexOf("\n") + 1);

			let _lines = str_part.split("\n");
			
			// _lines.forEach((line, line_idx) => {
			for (let line_idx = 0; line_idx < _lines.length; ++line_idx) {
				const line = _lines[line_idx];
				const cols = line.split("\t").map(a => a.trim());
				// try {
					rowParser(cols, line_idx);
				// }
				// catch (ex) {
				// 	console.error(ex);
				// }
			}
			// });
		});

		stream.on("close", function() {
			resolve();
		});
	});

	await promise;
}


