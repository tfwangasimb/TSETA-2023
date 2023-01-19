// @ts-check

// 2022/10/21: allow input align.sam + genome.fasta

const fs = require("node:fs");
const child_process = require("node:child_process");
const stream = require("node:stream");
const assert = require("node:assert");
const node_path = require("node:path");

// samtools depth QM6a-1_Aligned.sortedByCoord.out.bam > QM6a-1_coverage &
// node uint32array_coverage.js rid-D81_coverage C_d8r_coverage.uint32 &

const input_bam_file = process.argv[2];
const output_file = process.argv[3];


async function main() {
	if (input_bam_file == null || output_file == null) {
		console.log("cmd:", "node", __filename, "<input_bam_file>", "<output_file>");
		return;
	}
	
	// const chr_map = Object.fromEntries(chr_nameList.map((name, i) => [name, new Uint32Array(chr_lengthList[i])]));
	/** @type {string[]} */
	const chr_nameList = [];
	/** @type {number[]} */
	const chr_lengthList = [];
	/** @type {{ [chrName:string]:Uint32Array }} */
	const chr_map = {
	};

	if (node_path.extname(input_bam_file) == ".bam") {
		if (!fs.existsSync(`${input_bam_file}.bai`)) {
			try {
				child_process.execSync(`samtools index -@ 32 ${input_bam_file}`);
			}
			catch (ex) {
				console.error(ex);
				throw ex;
			}
		}

		try {
			const text = child_process.execSync(`samtools idxstats ${input_bam_file}`).toString();
			text.split("\n").map(line => line.split("\t")).map((cols, idx) => {
				const name = cols[0];
				const len = Number(cols[1]);
				
				if (name && len &&
					name != "*" &&
					!Number.isNaN(len) && Number.isSafeInteger(len)
				) {
					chr_map[name] = new Uint32Array(len);

					chr_nameList.push(name);
					chr_lengthList.push(len);

					console.log([Number(idx) + 1, name, len].join("\t"));
				}
			});
		}
		catch (ex) {
			console.error(ex);
			throw ex;
		}
	}
	else {
		const { readFasta } = require("../fasta_util.js");
		const { ArgvParser } = require("./ArgvParser.js");

		const argv_parser = new ArgvParser(process.argv);
		const genome_path = argv_parser.get(/--genome=(.+)/, (arg, args) => args[0]);
		if (genome_path && fs.existsSync(genome_path)) {
			const fa = readFasta(genome_path);
			Object.entries(fa).forEach(([name, seq], chrIdx) => {
				// const seq = fa[name];
				const len = seq.length;
				const nChr = chrIdx + 1;

				chr_nameList.push(name);
				chr_lengthList.push(len);
				chr_map[name] = new Uint32Array(len);
			});
		}
		else {
			throw new Error(`404: genome_path=${genome_path}`);
		}
	}

	if (Object.keys(chr_map).length <= 0) {
		throw new Error("chr_map");
	}

	const samtools = "samtools";
	const args = [
		"depth",
		input_bam_file,
	];
	const cmd = [samtools, ...args].join(" ");
	console.log(cmd);

	const proc_samtools_depth = child_process.spawn(samtools, args);

	// const e_UNMAP = 0x4;
	// const e_SECONDARY = 0x100;
	// const e_QCFAIL = 0x200;
	// const e_DUP = 0x400;
	// const gG_flags = `0x${(e_UNMAP | e_QCFAIL | e_DUP).toString(16)}`;

	// const g_flags = 0;// gG_flags;
	// const G_flags = 0;// gG_flags;

	// const proc_samtools_depth = child_process.exec(`samtools depth ${input_bam_file}`, {
	// 	maxBuffer: 1024 * 1024 * 1024, // 1GB
	// });

	let number_of_non_zero_value = 0;

	console.log("<readTableFromStream>");
	await readTableFromStream(proc_samtools_depth.stdout, function rowParser(cols) {
		// const obj = {};
		// header.reduce((prop, i) => obj[prop] = cols[i]);
		// rows.push(obj);

		const [strChr, strPos, value] = cols;
		try {
			const nVal = Number(value);
			chr_map[strChr][Number(strPos)] = nVal;

			if (nVal) {
				++number_of_non_zero_value;
			}
		}
		catch (ex) {
			console.error({ input_bam_file, output_file, strChr, strPos, value });
			throw ex;
		}
	});
	console.log("</readTableFromStream>");

	console.log({ number_of_non_zero_value });

	const ws = fs.createWriteStream(output_file);
	chr_nameList.forEach(name => {
		console.log("<" + name + "=" + chr_map[name].length + " />");
		ws.write(Buffer.from(chr_map[name].buffer));
	});

	if (typeof ws != "undefined") {
		ws.end(() => {
			const buffer = fs.readFileSync(output_file);
			const ui32a = new Uint32Array(buffer.buffer);
			chr_lengthList.reduce((pos, len, i) => {
				let j = 1;
				for (; j <= len; ++j) {
					if (ui32a[pos + j] > 0) {
						console.log([
							chr_nameList[i],
							pos,			// chr.offset
							j,				// chr.cursor
							ui32a[pos + j]	// chr[j].val
						].join("\t"));
						break;
					}
				}
				
				const aaa = new Uint32Array(buffer.buffer, pos * 4, len);
				assert.strictEqual(aaa.length, len, `aaa.length == len (${aaa.length} == ${len})`);
				assert.strictEqual(aaa[j], ui32a[pos + j], `aaa[${j}] == ui32a[${pos + j}] ${ui32a[pos + j]} (${aaa[j]} == ${ui32a[pos + j]})`);

				for (let k = 0; k < 10; ++k) {
					j = Math.trunc((0.25 + Math.random() * 0.5) * len);
					assert.strictEqual(aaa[j], ui32a[pos + j], `aaa[${j}] == ui32a[${pos + j}] ${ui32a[pos + j]} (${aaa[j]} == ${ui32a[pos + j]})`);
				}

				// console.log(0, ui32a[0]);
				// console.log(len - 1, ui32a[len - 1]);
				return pos + len;
			}, 0);
		});
	}
}

/**
 * @template T
 * @param {stream.Readable} read_stream
 * @param {(cols:string[]) => T} rowParser
 * @returns {Promise<void>}
 */
async function readTableFromStream(read_stream, rowParser) {
	// const stream = fs.createReadStream(filePath);

	const promise = new Promise(function (resolve, reject) {
		let last_str = "";

		read_stream.on("error", function (err) {
			reject(err);
		});

		read_stream.on("data", async function (chunk) {
			// console.log(chunk.length);

			const buf = chunk.toString();
			
			const str_1 = buf.slice(0, buf.lastIndexOf("\n"));
	
			const str_part = last_str + str_1;
			
			last_str = buf.slice(buf.lastIndexOf("\n") + 1);

			let _lines = str_part.split("\n");
			
			_lines.forEach((line, line_idx) => {
				const cols = line.split("\t").map(a => a.trim());
				rowParser(cols)
			});
		});

		read_stream.on("close", function() {
			resolve();
		});
	});

	await promise;
}

/**
 * @template T
 * @param {string} filePath
 * @param {(cols:string[]) => T} rowParser
 * @returns {Promise<T[]>}
 */
async function parseTableFromStream(filePath, rowParser) {
	const stream = fs.createReadStream(filePath);

	/** @type {Promise<T[]>} */
	const promise = new Promise(function (resolve, reject) {
		/** @type {T[]} */
		const rows = [];
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
			
			_lines.forEach((line, line_idx) => {
				const cols = line.split("\t").map(a => a.trim());
				rows.push(rowParser(cols));
			});
		});

		stream.on("close", function() {
			resolve(rows);
		});
	});

	const table = await promise;
	
	return table;
}

main();

class Row {
	constructor() {
		this.chr = "";
		this.pos = 0;
		this.value = 0;
	}
}
const RowDesc = {
	chr: String,
	pos: Number,
	value: Number,
};
const header = Object.keys(RowDesc);

// const chr_nameList = [
// 	"ChI_QM6a",
// 	"ChII_QM6a",
// 	"ChIII_QM6a",
// 	"ChIV_QM6a",
// 	"ChV_QM6a",
// 	"ChVI_QM6a",
// 	"ChVII_QM6a",
// ];
// const chr_lengthList =  [6835803, 6234656, 5311445, 4556834, 4159965, 4000387, 3823438];

// const chr_nameList = [
// 	"Ch1_CBS1-1_unitig_0-RV_consensus",
// 	"Ch2_CBS1-1_unitig_1_consensus",
// 	"Ch3_CBS1-1_unitig_2-RV_consensus",
// 	"Ch4_CBS1-1_unitig_3-RV_consensus",
// 	"Ch5_CBS1-1_unitig_5_consensus",
// 	"Ch6_CBS1-1_unitig_6-RV_consensus",
// 	"Ch7_CBS1-1_unitig_4-RV_consensus"
// ];
// const chr_lengthList = [
// 	6822680,
// 	5559498,
// 	5258134,
// 	4872985,
// 	4096940,
// 	3741771,
// 	3967191
// ];


