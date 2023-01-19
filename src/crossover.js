// @ts-check

// 20210312: multi_align_to_segfile.js => remove SNP if trlomere, centromere, rDNA
// CO in rDNA-IGS ?: RAD51-WT+spo11-WT＃1 Ch12

// 20210408: add GCasso_snp, GCasso_indel

const NCO_MODE_A = "all";         // [SNV+SNP*]
const NCO_MODE_S = "snv";         // SNV+, SNV to SNV
const NCO_MODE_DS = "snv+snp*";   // SNP* to SNV+ to SNP* and value(SNV) == value(SNP)
const env_NCO_MODE = (() => {
	if (process.env.NCO_MODE &&
		[NCO_MODE_A, NCO_MODE_S, NCO_MODE_DS].includes(process.env.NCO_MODE)
	) {
		return process.env.NCO_MODE;
	}
	else {
		return null;
	}
})();
/** @type {NCO_MODE_A|NCO_MODE_S|NCO_MODE_DS} */
const NCO_MODE = /*20200721*/NCO_MODE_A;//env_NCO_MODE || NCO_MODE_S;

// // 20200720_v3
// const IGNORE_INDEL_DUPLICATION_POSITION = true;
// const REMOVE_4n0_markers = true;//is_4n0_markers
// const splice_23232 = false;

// Last
const IGNORE_INDEL_DUPLICATION_POSITION = true;
const IGNORE_INDEL_DUPLICATION_POSITION_V2 = false;
const REMOVE_4n0_markers = true;//is_4n0_markers
const splice_23232 = true;


const output_text_table = false;

console.log({
	NCO_MODE,
	IGNORE_INDEL_DUPLICATION_POSITION,
	IGNORE_INDEL_DUPLICATION_POSITION_V2,
	splice_23232,
});

if (typeof String.prototype.matchAll != "function") {
	const matchAll = require("string.prototype.matchall");
	matchAll.shim();
}

const fs = require("fs");


const { argv_parse } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { loadSegFile, SegRow } = require("./SegFile.js");
const { readFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const dataset_path = String(argv["-dataset"] || "");
const argv_segfile = String(argv["--segfile"] || "");//mafft_QM6a_WT_1_Ch1_segfile.txt
const argv_output_prefix = String(argv["--output-prefix"] || "");

const dataset = Dataset.loadFromFile(dataset_path);
dataset.load_GC_content();
dataset.load_rDNA_info();

const genome_info_list = dataset.loadGenomeInfoList();

class CrossoverAnalyser {
	constructor() {
		this.chrMinLength = Infinity;
		this.chrMaxLength = 0;

		//meta
		this.chrLength = {};

		/** @type {{[nChr:number]:SegRow[]}} */
		this.segFile = {};
		this.featureCounts = {};
		this.has_spo11_oligo_reads = false;

		//output
		this.co_list = [];
		this.nco_list = [];
		// this.tab_co_list = [];
		// this.tab_nco_list = [];
	}
	clear() {
		this.chrMinLength = Infinity;
		this.chrMaxLength = 0;

		//metaf
		this.chrLength = {};

		//data
		this.segFile = {};
		this.featureCounts = {};
	}
	clearOutput() {
		//output
		this.myBook = null;
		this.co_list = [];
		this.nco_list = [];
		// this.tab_co_list = [];
		// this.tab_nco_list = [];
	}
	
	_saveTableFile() {
		throw new Error("Not implement");

		console.log("write to file");

		// allow output empty array
		// if (this.co_list.length || this.nco_list.length) {
			{
				const file_co = `${dataset.output_path}/${argv_output_prefix}_co.json`;

				// const props = [
				// 	"chr", "chr_len", "pos", "co_type", "GCasso_tract_len", "GCasso_marker", "snp_start_out", "snp_start_in", "snp_end_in", "snp_end_out", "type", "why_remove", "before", "after", "before_ref1_snp", "before_ref2_snp", "after_ref1_snp", "after_ref2_snp",
				// ];
				// const list = this.co_list.map(row => props.reduce((obj, key, idx) => {
				// 	obj[key] = row[idx];
				// 	return obj;
				// }, {}));//.filter(co => !co.why_remove);
				// fs.writeFileSync(file_co, JSON.stringify(list, null, "\t"));

				fs.writeFileSync(file_co, JSON.stringify(this.co_list, null, "\t"));

				console.log("file_co", file_co);
			}
			// if (output_text_table) {
			// 	const file_co = `${dataset.output_path}/${argv_output_prefix}_co.txt`;
			//
			// 	//no tetrad
			//
			// 	fs.writeFileSync(file_co, "chr#	chr_len	pos(bp)	co_type	GCasso_tract_len	GCasso_marker#	snp_start_out	snp_start_in	snp_end_in	snp_end_out	type	why_remove	before	after	before_ref1_snp	before_ref2_snp	after_ref1_snp	after_ref2_snp\r\n");
			//
			// 	this.tab_co_list.forEach(row => {
			// 		fs.writeFileSync(file_co, row.join("\t") + "\r\n", { flag: "a" });
			// 	});
			//
			// 	console.log("file_co", file_co);
			// }

			{
				const file_nco = `${dataset.output_path}/${argv_output_prefix}_nco.json`;

				// const props = [
				// 	"chr", "chr_len", "pos", "GCasso_tract_len", "GCasso_marker",
				// 	"snp_start_out", "snp_start_in", "snp_end_in", "snp_end_out",
				// 	"type", "why_remove",
				// ];
				// const list = this.nco_list.map(row => props.reduce((obj, key, idx) => {
				// 	obj[key] = row[idx];
				// 	return obj;
				// }, {}));
				// fs.writeFileSync(file_nco, JSON.stringify(list, null, "\t"));

				fs.writeFileSync(file_nco, JSON.stringify(this.nco_list, null, "\t"));

				console.log("file_nco", file_nco);
			}
			// if (output_text_table) {
			// 	const file_nco = `${dataset.output_path}/${argv_output_prefix}_nco.txt`;
			// 	fs.writeFileSync(file_nco, "chr#	chr_len	pos(bp)	GCasso_tract_len	GCasso_marker#	snp_start_out	snp_start_in	snp_end_in	snp_end_out	type	why_remove\r\n");
			// 	this.tab_nco_list.forEach(row => {
			// 		fs.writeFileSync(file_nco, row.join("\t") + "\r\n", { flag: "a" });
			// 	});
			// 	console.log("file_nco", file_nco);
			// }
		// }
	}
	
	saveFile() {
		console.log("write to file");
	
		const file_co = `${dataset.output_path}/${argv_output_prefix}_co.json`;
		fs.writeFileSync(file_co, JSON.stringify(this.co_list, null, "\t"));
		console.log("file_co", file_co);
		
		const file_nco = `${dataset.output_path}/${argv_output_prefix}_nco.json`;
		fs.writeFileSync(file_nco, JSON.stringify(this.nco_list, null, "\t"));
		console.log("file_nco", file_nco);
	}

	load() {
		{
			const group = loadSegFile(argv_segfile);

			//group
			//make index
			Object.keys(group).forEach(chr => {
				const nChr = Number(chr);
				const _list = group[nChr].sort((a, b) => a.pos - b.pos);

				/** @type {SegRow[]} */
				let list;
				if (IGNORE_INDEL_DUPLICATION_POSITION_V2) {
					let p1s = new Set();
					let p2s = new Set();
					list = [];
					for (let i = 0; i < _list.length; ++i) {
						const seg = _list[i];
						if (!p1s.has(seg.pos) && !p2s.has(seg.ref2_pos)) {
							p1s.add(seg.pos);
							p2s.add(seg.ref2_pos);
							list.push(seg);
						}
					}
				}
				else if (IGNORE_INDEL_DUPLICATION_POSITION) {
					list = [];

					/** @type {Map<number, SegRow[]>} */
					let pos_map = new Map();
					_list.forEach(data => {
						data.merge_del = Number(data.merge_del);
						let ll = pos_map.get(data.pos);
						if (!ll) {
							ll = [];
							pos_map.set(data.pos, ll);
						}
						ll.push(data);
					});

					pos_map.forEach(ll => {
						// if (ll.find(a => 260115 <= a.pos && a.pos <= 260118)) {
						// 	debugger
						// }
						if (ll.length == 1) {
							list.push(ll[0]);
						}
						else if (ll.length > 0) {
							if (false) {
								let nco_list = ll.filter(data => data.type == "type=3" || data.type == "type=4");
								let snv_list = ll.filter(data => !data.type || (!data.isInDel() && data.type == "type=2"));

								if (nco_list.length) {
									list.push(nco_list.pop());//nco
								}
								else if (snv_list.length) {
									list.push(snv_list.pop());//match or 2:2
								}
							}
							else {
								const last = ll[ll.length - 1];
								last.merge_del = ll.length;
								list.push(last);
							}
						}
					});
					list.sort((a, b) => a.pos - b.pos);
				}
				else {
					list = _list;
				}

				list.forEach((data, index) => {
					data.$index = index;
				});

				this.segFile[nChr] = list;
			});
		}

		{
			const chrList = genome_info_list[0].chr_list;

			this.chrMinLength = Infinity;
			this.chrMaxLength = 0;

			chrList.forEach(row => {
				this.chrMinLength = Math.min(this.chrMinLength, row.length);
				this.chrMaxLength = Math.max(this.chrMaxLength, row.length);

				this.chrLength[row.index] = row;
			});
		}
	}

	/**
	 * @param {number} max_dist
	 */
	find_nco_close_co(max_dist) {
		const a_nco = this.nco_list;//.filter(a => !a.why_remove);
		const a_co = this.co_list;//.filter(a => !a.why_remove);
		
		a_nco.forEach(nco => {
			if (a_co.some(co => Math.abs(co.pos - nco.pos) <= max_dist)) {
				nco.nco_type = "NCO2";
			}
			else {
				nco.nco_type = "NCO1";
			}
		});
	}

	find_co_nco_simple() {
		console.log("step1: find all co / nco");

		function is_snp_equal(snp_1, snp_2) {
			return (
				(snp_1[1] == 2 || snp_2[1] == 2 || snp_1[1] == snp_2[1]) &&
				(snp_1[2] == 2 || snp_2[2] == 2 || snp_1[2] == snp_2[2]) &&
				(snp_1[3] == 2 || snp_2[3] == 2 || snp_1[3] == snp_2[3]) &&
				(snp_1[4] == 2 || snp_2[4] == 2 || snp_1[4] == snp_2[4])
			);
		}

		let results = Object.keys(this.chrLength).map(snChr => {
			const nChr = Number(snChr);

			const seg = this.segFile[nChr];
			if (!seg || seg.length <= 0) {
				//console.error("chr:", nChr, "if (!seg || seg.length <= 0) {", this.chrLength, Object.keys(this.segFile));
				return;
			}

			// const seq_list = (function () {
			// 	debugger;
			// 	const input_path = `${dataset.output_path}/mafft_ch${nChr}.fa`;
			// 	const fa = readFasta(input_path);
			// 	const keys = genome_info_list.map(a => a.chr_list[nChr - 1].chr);
			// 	return keys.map(k => fa[k]);
			// })();

			// remove indel 1:3/3:1 if in sss region
			// Strain-specific sequences

			class SSS_region {
				/**
				 * @param {number} s
				 * @param {number} e
				 */
				constructor(s, e) {
					/** ref pos start */
					this.start = s;
					/** ref pos end */
					this.end = e;
					/** number of snp */
					this.snp = 0;
					/** number of snp22 */
					this.snp22 = 0;
					/** SNP[] */
					this.state = [0, 0, 0, 0];
				}
			}

			// /** @type {SSS_region[]} */
			// const sss_region = [];
			{
				// /** @type {SSS_region} */
				// let current_sss = null;

				seg.forEach((data, index) => {
					let nco = [0, 0];//["red", "blue"]
					for (let row = 1; row <= 4; ++row) {
						let is_red = data[row];
						++nco[is_red];
					}

					data.$type = `${nco[1]}:${4 - nco[1]}`;
					//data.$snp = (data[4] << 3 | data[3] << 2 | data[2] << 1 | data[1] << 0);

					data.n_type = (nco[0] == 3 || nco[1] == 3) ? 3 : ((nco[0] == 4 || nco[1] == 4) ? 4 : 2);

					// // if (data.pos >= 834326 && data.pos <= 837512) {
					// // 	debugger;
					// // }

					// if (834329 <= data.pos && data.pos <= 837512) {
					// 	console.log(data)
					// 	// console.log(index, data.$pos)
					// 	// debugger
					// 	if (current_sss) {
					// 		console.log(data.pos - current_sss.end)
					// 	}
					// }
					// if (index == 30164) {
					// 	debugger
					// }

					// // if (data.indel == 'indel=2|2') {
					// if (data.isInDel() && (!current_sss || (data.pos - current_sss.end) == 1)) {
					// 	if (!current_sss) {
					// 		current_sss = new SSS_region(data.pos, data.pos);
					// 	}
					// 	current_sss.end = data.pos;
					// 	// if (data.pos == 606797) {
					// 	// 	debugger;
					// 	// }
					// 	if (data.n_type == 2) {
					// 		current_sss.state[0] += Number(data[1]);
					// 		current_sss.state[1] += Number(data[2]);
					// 		current_sss.state[2] += Number(data[3]);
					// 		current_sss.state[3] += Number(data[4]);
					// 		current_sss.snp22 += 1;
					// 	}
					// 	current_sss.snp += 1;
					// }
					// else if (current_sss) {
					// 	if (current_sss.snp22 == 1) {// remove small indel (1bp)
					// 		current_sss = null;
					// 	}
					// 	else {
					// 		current_sss.state[0] /= current_sss.snp22;
					// 		current_sss.state[1] /= current_sss.snp22;
					// 		current_sss.state[2] /= current_sss.snp22;
					// 		current_sss.state[3] /= current_sss.snp22;
					// 		// move
					// 		sss_region.push(current_sss);
					// 		current_sss = null;
					// 	}
					// }

					// if (data.merge_del > 0) {// ref1 del / ref2 ins
					// 	if (current_sss) {
					// 		sss_region.push(current_sss);
					// 		current_sss = null;
					// 	}
					// 	else {
					// 		sss_region.push(new SSS_region(data.pos, data.pos));
					// 	}
					// }
				});
			}

			// console.log(sss_region.find(a => a.start <= 1082685 && 1082685 <= a.end));
			// // sss_region.map(a => [a.start, a.end, a.snp].join("\t")).join("\n")
			// debugger;

			function snp_distance(snp_1, snp_2) {
				return Math.abs(snp_2.pos - snp_1.pos);
			}

			const closeCOsMinDistance = dataset.crossover.closeCOsMinDistance;

			/** @type {{[start:number]:any}} */
			let snp_block_map = closeCOsMinDistance != null ? (() => {
				const snp_block_map = {};

				let prev_snp;
				for (let i = 0; i < seg.length; ++i) {
					let snp_1 = seg[i];
					if (snp_1.n_type == 2) {
						prev_snp = snp_1;
						break;
					}
				}
				for (let i = 0; i < seg.length; ++i) {
					let snp_1 = seg[i];
					if (snp_1.n_type == 2) {
						if (!is_snp_equal(snp_1, prev_snp)) {
							snp_block_map[prev_snp.pos] = {
								start: prev_snp,
								end: snp_1,
								length: snp_1.pos - prev_snp.pos,
							};
							prev_snp = snp_1;
						}
					}
				}

				return snp_block_map;
			})() : null;

			//find simple CO

			// let paired_co_list = [];
			let simple_co_list = [];
			{
				// var erereerr = seg.filter(b => b.pos >= 611186 && b.pos <= 612505);

				for (let i = 0; i < seg.length - 1; ++i) {
					let snp_1 = seg[i];
					let snp_2 = seg[i + 1];

					if (snp_1.$pos == 157346 || snp_1.$pos == 157345 || snp_1.$pos == 157347) {
						debugger
					}

					if (snp_1.n_type == 2 && snp_2.n_type == 2) {
						if (is_snp_equal(snp_1, snp_2)) {
						}
						else {
							// let co_pair = [snp_1, snp_2];
							// paired_co_list.push(co_pair);
							simple_co_list.push([snp_2]);
						}
					}
				}
				// paired_co_list.forEach(([snp_1, snp_2]) => {
				// 	snp_1.$type_name = "CO";
				// 	snp_2.$type_name = "CO";
				// });
				simple_co_list.forEach(([snp_2]) => {
					snp_2.$type_name = "CO";
				});
			}

			//find CO with GC

			// [3569ac] => SNP22
			// [012478bdef] => SNP13 | SNP31

			let str = seg.map(raw => raw[4] << 3 | raw[3] << 2 | raw[2] << 1 | raw[1] << 0).map(a => a.toString(16)).join("");//make string
			// let _all_232 = [...str.matchAll(/((?:3|5|6|9|a|c)(?:0|1|2|4|7|8|b|d|e|f)+(?:3|5|6|9|a|c))/g)];//find 2:2 1:3|3:1|0:4|4:0 2:2
			// let _all_conco = _all_232.filter(a => a[1][0] != a[1][a[1].length - 1]).map(a => seg.slice(a.index + 1, a.index + a[1].length - 1));//.map(a => [a.index, a.index + a[1].length - 1]);//co(nco) -> remove 2:2, save GC
			let all_conco = [];

			if (!splice_23232) {
				// let all_232 = [...str.matchAll(/([3569ac][012478bdef]+[3569ac])/g)];
				
				let all_232 = [...str.matchAll(/([3569ac](?:[012478bdef]+[3569ac])+)/g)];

				all_conco = all_232.filter(a => a[1][0] != a[1][a[1].length - 1]).map(a => seg.slice(a.index + 1, a.index + a[1].length - 1));//.map(a => [a.index, a.index + a[1].length - 1]);//co(nco) -> remove 2:2, save GC
			}
			// split ⑵3⑵3② // ⑵3⑵ => NCO; ⑵3② => CO(NCO)
			else {
				let all_23232 = [...str.matchAll(/([3569ac](?:[012478bdef]+[3569ac])+)/g)];
				// let all_2conco = [];
				all_23232.filter(a => a[1][0] != a[1][a[1].length - 1]).map(aaa => {
					//         2 3 3 2 3 2
					// snp2:2  a     b   c
					// if a == b then a,b is NCO
					// if b != c then b,c is CO(NCO)
					// find all snp 2:2
					const ttt = [...aaa[1].matchAll(/[3569ac]/g)].map(a => a);

					// find diff snp
					for (let i = 0; i < ttt.length - 1; ++i) {
						const current = ttt[i];
						const next = ttt[i + 1];
						const c = current[0];
						const n = next[0];
						if (c != n) {
							const co = seg.slice(aaa.index + current.index + 1, aaa.index + next.index);
							all_conco.push(co);
						}
					}
				});
			}

			// let all_nco = [...str.matchAll(/((?:0|1|2|4|7|8|b|d|e|f)+)/g)].map(a => seg.slice(a.index, a.index + a[1].length));
			//
			// let all_nco = [...str.matchAll(/[012478bdef]+/g)].map(a => seg.slice(a.index, a.index + a[0].length));//20200527
			let all_nco = [...str.matchAll(/(0+|1+|2+|4+|7+|8+|b+|d+|e+|f+)/g)].map(a => seg.slice(a.index, a.index + a[1].length));//20200527
			console.log({
				"all_nco.length": all_nco.length,
				// snp_block_map[end_out.pos].length < closeCOsMinDistance
			});

			// var erereerr = all_nco.filter(a => a.some(b => b.pos >= 576021 && b.pos <= 578951));
			debugger

			all_conco.forEach(snp_list => {
				snp_list.forEach(snp => {
					snp.$type_name = "CO(NCO)";
				});
			});
			let paired_conco_list = all_conco;
			let merge_co_list = [
				...paired_conco_list,
				...simple_co_list
			].filter(snps => snps && snps.length).sort((a_snps, b_snps) => a_snps[0].pos - b_snps[0].pos);//sort co and co(nco)

			const directly_remove = true;

			if (closeCOsMinDistance == null) {
			}
			else {
				let remove_snps = {};
				{//remove 2:2 blocks if too short
					merge_co_list.forEach((list, snps_idx) => {
						if (list[0].$type_name == "CO") {
							let head = list[0];
							let foot = list[list.length - 1];
							let start_out = head;
							let end_out = foot;

							if (!snp_block_map[end_out.pos]) {
								//console.log("end_out.pos", end_out.pos);
							}
							else if (snp_block_map[end_out.pos].length < closeCOsMinDistance) {
								list.forEach(snp => {
									if (directly_remove) {
										delete snp.$type_name;
									}
									else {
										snp.why_remove = "closeCOsMinDistance";
									}
								});
								remove_snps[snps_idx] = list;
								//console.log("remove co 2:2 blocks if too short:", head);
							}
						}
						else {
							let head = list[0];
							let foot = list[list.length - 1];
							let start_out = seg[head.$index - 1];
							let end_out = seg[foot.$index + 1];

							if (end_out) {
								if (!snp_block_map[end_out.pos]) {
									//console.log("end_out.pos", end_out.pos);
								}
								else if (snp_block_map[end_out.pos].length < closeCOsMinDistance) {
									list.forEach(snp => {
										if (directly_remove) {
											delete snp.$type_name;
										}
										else {
											snp.why_remove = "closeCOsMinDistance";
										}
									});
									remove_snps[snps_idx] = list;
									//console.log("remove co(nco) 2:2 blocks if too short:", head);
								}
							}
						}
					});
				}
				console.log("merge_co_list.length", merge_co_list.length);
				console.log("remove_snps.length", Object.keys(remove_snps).length);
				{//remove co if too close
					let first_snaps = merge_co_list[0];
					for (let snps_idx = 1; snps_idx < merge_co_list.length; ++snps_idx) {
						const current_snps = merge_co_list[snps_idx];
						let foot;

						// if (current_snps[0].pos == 1196305) {
						// 	debugger;
						// }

						if (current_snps[0].$type_name == "CO") {
							foot = current_snps[current_snps.length - 1];
						}
						else {
							foot = seg[current_snps[current_snps.length - 1].$index + 1];
						}

						let first_head;

						if (first_snaps[0].$type_name == "CO") {
							first_head = first_snaps[0];
						}
						else {
							first_head = seg[first_snaps[0].$index - 1];
						}

						if (!first_head || !foot) {
							continue;
						}

						let prev_co_dist = snp_distance(first_head, foot);
						if (prev_co_dist < closeCOsMinDistance) {
							const prev_snps_idx = snps_idx - 1;
							const prev_snaps = merge_co_list[prev_snps_idx];

							remove_snps[prev_snps_idx] = prev_snaps;
							prev_snaps.forEach(snp => {
								if (directly_remove) {
									delete snp.$type_name;
								}
								else {
									snp.why_remove = "closeCOsMinDistance";
								}
							});
							//console.log("remove L:", prev_snaps[0], [first_head, foot]);

							if (is_snp_equal(first_head, foot)) {
								remove_snps[snps_idx] = current_snps;
								current_snps.forEach(snp => {
									if (directly_remove) {
										delete snp.$type_name;
									}
									else {
										snp.why_remove = "closeCOsMinDistance";
									}
								});
								//console.log("remove R:", current_snps[0]);
							}
						}
						else {
							first_snaps = current_snps;
						}
					}
				}
			}
			merge_co_list = merge_co_list.filter(list => list[0].$type_name);

			// find NCO, 2NCO

			/** @type {SegRow[][]} */
			let paired_nco_list = [];
			/** @type {SegRow[][]} */
			let paired_2nco_list = [];
			all_nco.forEach(_snp_list => {
				if (_snp_list.every(a => a.rip)) {
					return;
				}
				const snp_list = REMOVE_4n0_markers ? _snp_list.filter(a => !a.is_4n0_markers()) : _snp_list;// 20210127//.filter(a => !a.is_4n0_markers());// remove 4n:0 markers // why remove 4n:0 markers
				// if (!snp_list.length) {
				// 	fs.writeFileSync("ch" + nChr + "_removed_4n0.txt", JSON.stringify(_snp_list) + "\n", { flag: "a" });
				// 	return;
				// }//200603
				if (snp_list.length != _snp_list.length) {
					// const unknow_warn = JSON.stringify({ tetrad: dataset.name, nChr, _snp_list }) + "\n";
					const { pos: p11, ref2_pos: p12 } = _snp_list[0];
					const { pos: p21, ref2_pos: p22 } =_snp_list[_snp_list.length - 1];
					const unknow_warn = [dataset.name, nChr, _snp_list.length, snp_list.length, p11, p21, p12, p22].join("\t") + "\n";
					fs.writeFileSync("./unknow_warn.txt", unknow_warn, { flag: "a" });
				}

				if (snp_list.some(a => a.$type_name == "CO(NCO)")) {
					if (snp_list.every(a => a.$type_name == "CO(NCO)")) {
						//skip CO(NCO)
					}
					else {
						console.error("?? CO(NCO)", snp_list);
					}
				}
				else if (snp_list.every(a => a.n_type == 4)) {
					if (NCO_MODE == NCO_MODE_A) {
						const list = snp_list;
						const _n_list = snp_list.filter(a => a.isInDel());
						list.forEach(snp => {
							// if (snp.pos == 671437) {
							// 	console.log(snp);
							// 	console.log("_n_list.length", _n_list.length);
							// }
							// snp.$type_name = "2NCO" + (_n_list.length ? ("=" + _n_list.length + "/" + list.length) : "");
							snp.$type_name = "2NCO" + (_n_list.length == list.length ? " (indel)" : "");
						});

						paired_2nco_list.push(list);
					}
					else {
						const snv_list = snp_list.filter(a => !a.isInDel()); //200203 // remove every indel
						if (snv_list.length) {
							if (NCO_MODE == NCO_MODE_DS) {
								snp_list.forEach(snp => {
									snp.$type_name = "2NCO" + (snp_list.length ? ("=" + snv_list.length + "/" + snp_list.length) : "");
								});

								paired_2nco_list.push(snp_list);
							}
							else {// else if (NCO_MODE == NCO_MODE_S) {
								snv_list.forEach(snp => {
									snp.$type_name = "2NCO";
								});
								paired_2nco_list.push(snv_list);
							}
						}
					}
				}
				else {
					if (NCO_MODE == NCO_MODE_A) {
						const list = snp_list;
						const _n_list = snp_list.filter(a => a.isInDel());
						list.forEach(snp => {
							// if (snp.pos == 671437) {
							// 	console.log(snp);
							// 	console.log("_n_list.length", _n_list.length);
							// }
							// snp.$type_name = "NCO" + (_n_list.length ? ("=" + _n_list.length + "/" + list.length) : "");
							snp.$type_name = "NCO" + (_n_list.length == list.length ? " (indel)" : "");
						});

						paired_nco_list.push(list);
					}
					else {
						// if (snp_list[0].pos == 2042652) {
						// 	console.log(2042652, 2042652, 2042652, 2042652);
						// 	console.log(snp_list.filter(a => !a.isInDel()));
						// }
						const snv_list = snp_list.filter(a => !a.isInDel()); //200203 // remove every indel
						if (snv_list.length) {
							if (NCO_MODE == NCO_MODE_DS) {
								snp_list.forEach(snp => {
									snp.$type_name = "NCO" + (snp_list.length ? ("=" + snv_list.length + "/" + snp_list.length) : "");
								});

								paired_nco_list.push(snp_list);
							}
							else {// else if (NCO_MODE == NCO_MODE_S) {
								snv_list.forEach(snp => {
									snp.$type_name = "NCO";
								});
								paired_nco_list.push(snv_list);
							}
						}
					}
				}
			});

			debugger

			// check Q, C, S
			if (merge_co_list.length) {
				let co_prev_snp = seg[0];
				for (let index = 0; index < merge_co_list.length - 1; ++index) {
					const snp = merge_co_list[index][0];
					const next = merge_co_list[index + 1][0];
					before_after_crossover(snp, co_prev_snp, next);
					co_prev_snp = snp;
				}
				before_after_crossover(merge_co_list[merge_co_list.length - 1][0], co_prev_snp, seg[seg.length - 1]);
			}
			function before_after_crossover(snp, prev, next) {//QC_before_crossover
				let { ref: before, ref_snp: ref_snp_before, why_remove: why_remove_before, remark_before } = snp_qc(prev, snp, snp);
				let { ref: after, ref_snp: ref_snp_after, why_remove: why_remove_after, remark_after } = snp_qc(snp, next, snp);

				if (why_remove_before || why_remove_after) {
					if (why_remove_before == why_remove_after) {
						snp.why_remove = why_remove_before;
					}
					else {
						snp.why_remove = [why_remove_before, why_remove_after].join(",");
					}
				}
				//snp.co_snp = [snp.pos, next.pos, aa[1], aa[2], aa[3], aa[4]].join(",");
				snp.before = [before[1], before[2], before[3], before[4]].join(",");
				snp.after = [after[1], after[2], after[3], after[4]].join(",");

				snp.before_ref1_snp = ref_snp_before.map(a => a[0]);
				snp.before_ref2_snp = ref_snp_before.map(a => a[1]);

				snp.after_ref1_snp = ref_snp_after.map(a => a[0]);
				snp.after_ref2_snp = ref_snp_after.map(a => a[1]);

				if (remark_before) {
					snp.remark_before = remark_before;
				}
				if (remark_after) {
					snp.remark_after = remark_after;
				}
			}

			/**
			 * @param {SegRow} from_snp
			 * @param {SegRow} to_snp
			 * @param {SegRow} snp_header
			 */
			function snp_qc(from_snp, to_snp, snp_header) {
				let snps = seg.filter(a => a.pos >= from_snp.pos && a.pos <= to_snp.pos)
				.filter(a => {
					return (
						a.$type_name == null ||
						(a.$type_name.indexOf("CO") < 0 && a.$type_name.indexOf("NCO") < 0)
					);
				});// remove 22_SNP_MARKER if in CO

				if (!snps.length) {
					return {
						ref_snp: [
							[0, 0],
							[0, 0],
							[0, 0],
							[0, 0]
						],
						ref: [
							snp_header[1],
							snp_header[2],
							snp_header[3],
							snp_header[4],
						],
						// why_remove: `not enough snp ${from_snp.pos}~${to_snp.pos}`,
					};
				}

				let num_22 = 0;
				let num_not22 = 0;

				// number of ref2_snp
				let aa = {
					1: 0,
					2: 0,
					3: 0,
					4: 0
				};
				snps.forEach(_snp => {
					// if (_snp.$type_name) {
					// 	return;
					// }
					let [a, b, c, d] = [_snp[1] | 0, _snp[2] | 0, _snp[3] | 0, _snp[4] | 0];
					let n = a + b + c + d;
					if (n == 2) {// two ref2
						aa[1] += a;
						aa[2] += b;
						aa[3] += c;
						aa[4] += d;
						++num_22;
					}
					else {
						// if (!_snp.is_4n0_markers()) {// why remove 4n:0 markers
							++num_not22;
						// }
					}
				});

				if (num_not22) {
					console.warn({ num_not22 });
					// debugger;
				}

				// let hl = snps.length / 2;
				// let bb = {};
				// Object.keys(aa).forEach(k => bb[k] = aa[k] >= hl ? 1 : 0);
				// let ref_snp = Object.keys(aa).map(k => [snps.length - aa[k], aa[k]]);

				const [ref2_1, ref2_2, ref1_1, ref1_2] = Object.entries(aa).sort((a, b) => b[1] - a[1]).map(a => a[0]);
				const bb = {};
				bb[ref1_1] = 0;
				bb[ref1_2] = 0;
				bb[ref2_1] = 1;
				bb[ref2_2] = 1;

				/** @type {number[][]} [tetrad_idx][num_of_snp] */
				const ref_snp = Object.keys(aa).map(k => [snps.length - aa[k], aa[k]]);

				const hl = num_22 / 2;
				if (Object.keys(aa).every(k => bb[k] == (aa[k] >= hl ? 1 : 0))) {
				}
				else {
					console.log({
						ref1_start: from_snp.pos,
						ref1_end: to_snp.pos,
						ref_snp: ref_snp,
						ref: bb,
						num_not22,
					});
					// TDOD:
					// output log files
					// throw new Error("if (bb[ref2_1] == 1 && bb[ref2_2] == 1 && bb[ref1_1] == 0 && bb[ref1_2] == 0) {")
					return {
						ref_snp: ref_snp,
						ref: bb,
						remark: `not SNP 2:2 region ${from_snp.pos}~${to_snp.pos}`,
					};
				}
				return {
					ref_snp: ref_snp,
					ref: bb,
				};
			}

			// end check Q, C, S

			/** @type {SegRow[][]} */
			const merge_list = [
				...paired_nco_list,
				...paired_2nco_list,
				...merge_co_list
			].filter(snps => snps && snps.length);

			// merge_list.forEach(snps => {
			// 	const first = snps[0].pos;
			// 	const last = snps[snps.length - 1].pos;
			// 	let in_sss = sss_region.find(sss => {
			// 		// return first <= sss.end && last >= sss.start;// ??? 20210125
			// 		return (
			// 			(sss.start <= first && last <= sss.end) // NCO in indel
			// 			// (sss.start <= first && first <= sss.end) || (sss.start <= last && last <= sss.end) // snp in indel
			// 		);// ??? 20210125
			// 	});
			// 	if (in_sss != null) {
			// 		snps.forEach(snp => {
			// 			const sss_len = in_sss.end - in_sss.start + 1;
			// 			const range = `${in_sss.start}~${in_sss.end}`;
			// 			if (range == "606797~612995") {
			// 				debugger;
			// 			}
			// 			snp.why_remove = `Strain-specific sequences (${range})(${sss_len})`;
			// 		});
			// 	}
			// });

			// let g_region_id = 0;
			let ret_list = [];
			// function make_header(list) {
			// 	++g_region_id;
			// 	let head = list[0];
			// 	let foot = list[list.length - 1];
			// 	if (head.$type_name == "CO") {
			// 		_make_header(head, head, foot, foot);
			// 	}
			// 	else {
			// 		_make_header(seg[head.$index - 1], head, foot, seg[foot.$index + 1]);
			// 	}
			// }
			// function _make_header(start_out, start_in, end_in, end_out) {
			// 	start_out = start_out || start_in;
			// 	end_out = end_out || end_in;

			// 	head_comment.is_comment = true;

			// 	head_comment.g_region_id = start_in.region_id;
			// 	head_comment.start = head.pos;
			// 	head_comment.end = foot.pos;
			// 	head_comment.target = head;

			// 	head_comment.idx_start_out = idx_start_out;
			// 	head_comment.idx_start_in = idx_start_in;
			// 	head_comment.idx_end_in = idx_end_in;
			// 	head_comment.idx_end_out = idx_end_out;

			// 	head_comment["snp_start_out"] = start_out.pos;
			// 	head_comment["snp_start_in"] = head.pos;
			// 	head_comment["snp_end_in"] = foot.pos;
			// 	head_comment["snp_end_out"] = end_out.pos;
			// 	head_comment["pos(bp)"] = Math.round((mid2 + mid1) * 0.5);//191024
			// 	head_comment["tract_len"] = Math.round(Math.abs(mid2 - mid1));//191024
			// 	head_comment["tract_min_len"] = Math.abs(end_in.pos - start_in.pos);
			// 	head_comment["tract_max_len"] = Math.abs(end_out.pos - start_out.pos);
			// }
			merge_list.forEach((list, list_idx) => {
				let region_id = list_idx + 1;
				let head = list[0];
				let foot = list[list.length - 1];

				let idx_start_out, idx_start_in, idx_end_in, idx_end_out;
				if (head.$type_name == "CO") {
					idx_start_out = head.$index;
					idx_start_in = head.$index;
					idx_end_in = foot.$index;
					idx_end_out = foot.$index;
				}
				else {
					idx_start_out = head.$index - 1;
					idx_start_in = head.$index;
					idx_end_in = foot.$index;
					idx_end_out = foot.$index + 1;
				}

				let start_in = seg[idx_start_in];
				let end_in = seg[idx_end_in];
				let start_out = seg[idx_start_out] || start_in;
				let end_out = seg[idx_end_out] || end_in;

				let mid1 = (start_in.pos + start_out.pos) * 0.5;
				let mid2 = (end_in.pos + end_out.pos) * 0.5;

				list.forEach(a => {
					a.region_id = region_id;
					a.$type_name = head.$type_name;
				});


				let head_comment;
				if (head.$type_name == "CO(NCO)") {
					head_comment = {
						...start_out,
						$type_name: head.$type_name,
						co_type: 1,
						"GCasso_marker#": list.length,//191021
						"GCasso_snp": list.filter(a => !a.isInDel()).length,//20210408
						"GCasso_indel": list.filter(a => a.isInDel()).length,//20210408
					};//clone
				}
				else {
					head_comment = {
						...head,
						co_type: 0,
						// always zero
						"GCasso_marker#": head.$type_name == "CO" ? 0 : list.length,//191021
						"GCasso_snp": head.$type_name == "CO" ? 0 : list.filter(a => !a.isInDel()).length,//20210408
						"GCasso_indel": head.$type_name == "CO" ? 0 : list.filter(a => a.isInDel()).length,//20210408
					};//clone
				}

				{
					if (head.$type_name == "CO") {
						//head_comment.co_snp = start_out.co_snp || start_in.co_snp;
						head_comment.before = start_out.before || start_in.before;
						head_comment.after = start_out.after || start_in.after;

						head_comment.before_ref1_snp = start_out.before_ref1_snp || start_in.before_ref1_snp;
						head_comment.before_ref2_snp = start_out.before_ref2_snp || start_in.before_ref2_snp;
						head_comment.after_ref1_snp = start_out.after_ref1_snp || start_in.after_ref1_snp;
						head_comment.after_ref2_snp = start_out.after_ref2_snp || start_in.after_ref2_snp;
					}
					if (head.$type_name == "CO(NCO)") {
						//head_comment.co_snp = start_out.co_snp || start_in.co_snp;
						head_comment.before = start_out.before || start_in.before;
						head_comment.after = start_out.after || start_in.after;

						head_comment.before_ref1_snp = start_out.before_ref1_snp || start_in.before_ref1_snp;
						head_comment.before_ref2_snp = start_out.before_ref2_snp || start_in.before_ref2_snp;
						head_comment.after_ref1_snp = start_out.after_ref1_snp || start_in.after_ref1_snp;
						head_comment.after_ref2_snp = start_out.after_ref2_snp || start_in.after_ref2_snp;
					}
				}

				head_comment.is_comment = true;
				delete head_comment.is_start;
				delete head_comment.is_end;
				delete head_comment.is_single;

				head_comment.region_id = region_id;
				head_comment.start = head.pos;
				head_comment.end = foot.pos;
				head_comment.target = head;

				head_comment.idx_start_out = idx_start_out;
				head_comment.idx_start_in = idx_start_in;
				head_comment.idx_end_in = idx_end_in;
				head_comment.idx_end_out = idx_end_out;

				head_comment["snp_start_out"] = start_out.pos;
				head_comment["snp_start_in"] = head.pos;
				head_comment["snp_end_in"] = foot.pos;
				head_comment["snp_end_out"] = end_out.pos;
				//
				head_comment["ref2_snp_start_out"] = start_out.ref2_pos;
				head_comment["ref2_snp_start_in"] = head.ref2_pos;
				head_comment["ref2_snp_end_in"] = foot.ref2_pos;
				head_comment["ref2_snp_end_out"] = end_out.ref2_pos;
				//
				head_comment["$snp_start_out"] = start_out.$pos;
				head_comment["$snp_start_in"] = head.$pos;
				head_comment["$snp_end_in"] = foot.$pos;
				head_comment["$snp_end_out"] = end_out.$pos;
				//
				// if (head.$type_name == "CO") {
				// 	head_comment._snp_start_out = seg[head_comment.$index - 2].pos;
				// 	head_comment._snp_end_out = seg[head_comment.$index + 1].pos;
				// }
				// if (head.$type_name == "CO(NCO)") {
				// 	head_comment._snp_start_out = seg[head_comment.$index - 1].pos;
				// 	head_comment._snp_end_out = seg[head_comment.$index + 1].pos;
				// }

				head_comment["pos(bp)"] = Math.round((mid2 + mid1) * 0.5);//191024
				head_comment["tract_len"] = Math.round(Math.abs(mid2 - mid1));//191024
				head_comment["tract_min_len"] = Math.abs(end_in.pos - start_in.pos);
				head_comment["tract_max_len"] = Math.abs(end_out.pos - start_out.pos);

				// 20210312
				// @see {@link multi_align_to_segfile.js:REMOVE_SNP_IF_NOT_TARGET}

				/** @type {{ range: number[]|[number, number]; type: string; pos_type?: "tseta"|"ref1"|"ref2" }[] } */
				let remove_region_list = [];
				if (
					(dataset.rDNA != null && head_comment.chr == dataset.rDNA.nChr) || 
					(dataset.rDNA_info != null && head_comment.chr == dataset.rDNA_info.chr)
				) {
					if (dataset.rDNA.region != null) {
						remove_region_list.push({
							range: dataset.rDNA.region,
							type: "rDNA",
						});
					}
					if (dataset.rDNA.region_ref2 != null) {
						remove_region_list.push({
							range: dataset.rDNA.region_ref2,
							type: "rDNA",
							pos_type: "ref2",
						});
					}
				}
				// if (dataset.rDNA_info != null && head_comment.chr == (dataset.rDNA.nChr || dataset.rDNA_info.chr)) {
				// 	remove_region_list.push({
				// 		// range: dataset.rDNA_info.data[0].region,
				// 		range: [dataset.rDNA_info.region_start, dataset.rDNA_info.region_end],
				// 		// dataset.rDNA_info.data[0].region.map(a => ref1_pos_uint32array[a]).forEach(a => viewerState.pushPositionMarker(a))
				// 		type: "rDNA",
				// 		pos_type: "tseta",
				// 	});
				// }
				if (dataset.centromere && dataset.centromere[head_comment.chr]) {
					remove_region_list.push({
						range: dataset.centromere[head_comment.chr],
						type: "centromere",
					});
				}
				try {
					dataset.all_telomere.slice(0, 2).forEach((telomere, genome_idx) => {
						if (telomere && telomere[head_comment.chr]) {
							let [l_tele, r_tele] = telomere[head_comment.chr];
							remove_region_list.push({
								range: l_tele,
								type: "telomere",// + "," + ["R", genome_idx, head_comment.chr].join(","),
								pos_type: ["ref1", "ref2"][genome_idx],
							}, {
								range: r_tele,
								type: "telomere",// + "," + ["R", genome_idx, head_comment.chr].join(","),
								pos_type: ["ref1", "ref2"][genome_idx],
							});
						}
					});
				}
				catch (ex) {
					console.error(ex);
					console.log(dataset.telomere);
					console.log(head_comment.chr);
				}
				// if (!remove_region_list.length) {
				// 	if (
				// 		dataset.centromere && dataset.centromere[head_comment.chr] &&
				// 		dataset.telomere && dataset.telomere[head_comment.chr]
				// 	) {
				// 		console.log(dataset.telomere);
				// 		console.log(dataset.telomere[head_comment.chr]);
				// 		throw 12;
				// 	}
				// }
				remove_region_list.forEach(region => {
					let [r_start, r_end] = region.range;
					let [start, end] = (function () {
						if (region.pos_type == "tseta") {
							let $start = head_comment["$snp_start_out"];
							let $end = head_comment["$snp_end_out"];
							return [$start, $end];
						}
						else if (region.pos_type == "ref2") {
							let ref2_start = head_comment["ref2_snp_start_out"];
							let ref2_end = head_comment["ref2_snp_end_out"];
							return [ref2_start, ref2_end];
						}
						else {
							let _start = head_comment["snp_start_out"];
							let _end = head_comment["snp_end_out"];
							return [_start, _end];
						}
					})();
					if (start <= r_end && end >= r_start) {
						head_comment.why_remove = region.type;// + "|" + [r_start, r_end, start, end].join(",");
						//console.log(head_comment, region.type);
					}
				});

				ret_list.push(head_comment);

				if (head == foot && head.$type_name != "CO(NCO)" && head.$type_name != "CO") {//only nco 2nco//191024
					delete head.is_start;
					head_comment.is_start = true;
					head.is_single = true;
					foot.is_end = true;
					//head.region_id = region_id;

					ret_list.push(head);
				}
				else {
					if (head.$type_name == "CO(NCO)") {
						delete head.is_start;
						delete foot.is_end;
						//
						start_out.$type_name = head.$type_name;
						start_out.region_id = region_id;
						start_out.is_start = true;
						//
						end_out.$type_name = head.$type_name;
						end_out.region_id = region_id;
						end_out.is_end = true;
						//
						ret_list.push(start_out);
						//ret_list.push(head);
						//ret_list.push(foot);
						ret_list.push(end_out);
					}
					else {
						delete head.is_start;
						//delete foot.is_end;
						//
						head_comment.is_start = true;
						//
						//head.region_id = region_id;
						//head.is_start = true;
						//foot.region_id = region_id;
						foot.is_end = true;
						//
						ret_list.push(head);
						ret_list.push(foot);
					}
				}

				if (head_comment.merge_del)
					console.log(head_comment.merge_del);
			});

			return ret_list;
		});

		results = results.filter(a => a);
		if (results.length <= 0) {
			throw new Error("no result: results.length -> 0");
		}

		results.forEach(byChr => {
			const _sorted_filtered = byChr.filter(data => data.is_comment).sort((a, b) => a.snp_start_out - b.snp_start_out);

			// TODO: remove error co
			const sorted_filtered = _sorted_filtered;//.filter();

			sorted_filtered.forEach(data => {
				if (data.$type_name == "CO" || data.$type_name == "CO(NCO)") {
					// if (output_text_table) {
					// 	let row = [
					// 		//tetrad,//200217
					// 		data.chr,//                        "chr",
					// 		this.chrLength[data.chr].length,//    "chr_len", "pos", "co_type", "GCasso_tract_len", "GCasso_marker", "snp_start_out", "snp_start_in", "snp_end_in", "snp_end_out", "type", "why_remove", "before", "after", "before_ref1_snp", "before_ref2_snp", "after_ref1_snp", "after_ref2_snp",
					// 		data["pos(bp)"],
					// 		data.co_type,
					// 		data.tract_len,//GCasso_tract_len,
					// 		data["GCasso_marker#"],
					// 		data.snp_start_out, data.snp_start_in, data.snp_end_in, data.snp_end_out,
					// 		data.$type_name,
					// 		data.why_remove,
					// 		//data.co_snp,
					// 		data.before,
					// 		data.after,
					// 		data.before_ref1_snp,
					// 		data.before_ref2_snp,
					// 		data.after_ref1_snp,
					// 		data.after_ref2_snp,
					// 	];
					// 	this.tab_co_list.push(row);
					// }
					data.chr              = data.chr,
					data.chr_len          = this.chrLength[data.chr].length,
					data.pos              = data["pos(bp)"],
					data.co_type          = data.co_type,
					data.GCasso_tract_len = data.tract_len,
					data.GCasso_marker    = data["GCasso_marker#"],
					data.GCasso_snp       = data["GCasso_snp"],
					data.GCasso_indel     = data["GCasso_indel"],
					data.snp_start_out    = data.snp_start_out,
					data.snp_start_in     = data.snp_start_in,
					data.snp_end_in       = data.snp_end_in,
					data.snp_end_out      = data.snp_end_out,
					data.type             = data.$type_name,
					data.why_remove       = data.why_remove,
					data.before           = data.before,
					data.after            = data.after,
					data.before_ref1_snp  = data.before_ref1_snp,
					data.before_ref2_snp  = data.before_ref2_snp,
					data.after_ref1_snp   = data.after_ref1_snp,
					data.after_ref2_snp   = data.after_ref2_snp,
					this.co_list.push(data);
				}
				else if (
					(NCO_MODE == NCO_MODE_A  && (data.$type_name.indexOf("NCO") >= 0)) ||
					(NCO_MODE == NCO_MODE_DS && (data.$type_name.indexOf("NCO") >= 0)) ||
					(NCO_MODE != NCO_MODE_A  && (data.$type_name == "NCO" || data.$type_name == "2NCO"))
				) {
					// if (output_text_table) {
					// 	let row = [
					// 		//tetrad,
					// 		data.chr,
					// 		this.chrLength[data.chr].length,
					// 		data["pos(bp)"],
					// 		data.tract_len,//GCasso_tract_len,
					// 		data["GCasso_marker#"],
					// 		data.snp_start_out, data.snp_start_in, data.snp_end_in, data.snp_end_out,
					// 		data.$type_name,
					// 		data.why_remove,
					// 	];
					// 	this.tab_nco_list.push(row);
					// }
					data.chr              = data.chr;
					data.chr_len          = this.chrLength[data.chr].length;
					data.pos              = data["pos(bp)"];
					data.GCasso_tract_len = data.tract_len;
					data.GCasso_marker    = data["GCasso_marker#"];
					data.snp_start_out    = data.snp_start_out;
					data.snp_start_in     = data.snp_start_in;
					data.snp_end_in       = data.snp_end_in;
					data.snp_end_out      =  data.snp_end_out;
					data.type             = data.$type_name;
					data.why_remove       = data.why_remove;
					this.nco_list.push(data);
				}
				else {
					//console.log("not co/nco", data);
					throw new Error("ls CO/NCO");
				}
			});
		});
	}

	/**
	 * 
	 * @param {(Crossover|NonCrossover)[]} list
	 * @param {RegExp} start_regexp
	 * @param {RegExp} end_regexp
	 * @param {number} max_dist
	 */
	telomere_s(list, start_regexp, end_regexp, max_dist = 100) {
		if (!list.length) {
			return;
		}
		const nChr = Number(list[0].chr);
		const chr_idx = Number(list[0].chr) - 1;
		
		const ref_chr_list = genome_info_list[0].chr_list[chr_idx];

		// TODO: load ?? fasta
		// NOTICE !!!
		const seq_list = Object.values(readFasta(ref_chr_list.path));

		// TODO: user input telomere seq

		let left_telo;
		{
			const sss = seq_list[0].slice(0).replace(/-/g, "");
			const aa = [...sss.matchAll(RegExp(start_regexp))].map(m => m.index);
			const aaa = aa.map((v, i, a) => a[i + 1] - v)
			
			const left_last = aaa.findIndex(a => a > max_dist);
			left_telo = aa[left_last];

			list.filter((_, co_idx) => {
				const co = list[co_idx];
				if (!co.why_remove) {
					if (co.snp_end_in <= left_telo) {
						co.why_remove = "telomere";
					}
					return co.why_remove;
				}
			});
		}
		
		let right_telo;
		{
			const sss = seq_list[0].slice(0).replace(/-/g, "");
			const aa = [...sss.matchAll(end_regexp)].map(m => m.index).reverse();
			const aaa = aa.map((v, i, a) => -(a[i + 1] - v));
			const right_first = aaa.findIndex(a => a > max_dist);
			right_telo = aa[right_first];

			list.filter((_, co_idx) => {
				const co = list[co_idx];
				if (!co.why_remove) {
					if (co.snp_start_in >= right_telo) {
						co.why_remove = "telomere";
					}
					return co.why_remove;
				}
			});
		}

		fs.writeFileSync(`${dataset.tmp_path}/telomere_Ch${nChr}.json`, JSON.stringify([left_telo, right_telo]));
	}
}

let conco = new CrossoverAnalyser();

conco.load();
conco.find_co_nco_simple();

// dataset.NCO2.nco_close_co_dist
conco.find_nco_close_co(2000);


conco.saveFile();

module.exports = CrossoverAnalyser;
