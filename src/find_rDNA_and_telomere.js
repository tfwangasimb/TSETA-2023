// @ts-check

const fs = require("fs");
const Path = require("path");
const child_process = require("child_process");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { parse_blastn_results, BlastnCoord } = require("./blastn_util.js");
const { readFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}

const genome_info_list = dataset.loadGenomeInfoList();

main();

async function main() {
	const argv_genome_ref1 = String(argv["--parental-1"] || "");
	const argv_genome_ref2 = String(argv["--parental-2"] || "");
	
	const argv_telomere = String(argv["--parental-telomere"] || "");
	const argv_ref_rDNA = String(argv["--parental-rDNA"] || "");
	const argv_nChr_rDNA = Number(argv["--rDNA-ch"] || "");

	if (
		Number.isNaN(argv_nChr_rDNA) == true ||
		Number.isSafeInteger(argv_nChr_rDNA) == false ||
		argv_nChr_rDNA <= 0
	) {
		throw new Error("--rDNA-ch");
	}

	/** @typedef {{ [nChr:string]: [[left_start: number, left_end: number], [right_start: number, right_end: number]]; }} TelomereDef */
	/** @type {TelomereDef[]} TelomereDef[genome_idx] */
	const telomere_data = JSON.parse(fs.readFileSync(argv_telomere).toString()).slice(0, 2);// only 2 ref
	
	/** @typedef {[start: number, end: number]} rDNA_Def */
	/** @type {rDNA_Def[]} TelomereDef[genome_idx] */
	const rDNA_data = JSON.parse(fs.readFileSync(argv_ref_rDNA).toString());

	const ref1_ref2_seq_list = [
		Object.values(readFasta(argv_genome_ref1)),
		Object.values(readFasta(argv_genome_ref2)),
	];

	{
		// console.log("find rDNA");

		const rDNA_chrIdx = argv_nChr_rDNA - 1;

		const out = [];
		const task_genome = genome_info_list.map(async (genome_info, genome_idx) => {
			const chr_info = genome_info.chr_list[rDNA_chrIdx];
			const chr_fasta = chr_info.path;

			out[genome_idx] = [];

			const task_rDNA_ref = rDNA_data.map(async (ref_rDNA_region, ref_idx) => {
				const ref_seq_list = ref1_ref2_seq_list[ref_idx];
				const ref_chr_seq = ref_seq_list[rDNA_chrIdx];
				
				out[genome_idx][ref_idx] = [];
				
				const aln_list = merge_overlap_alignment(await blastn(chr_fasta, {
					query_seq: ref_chr_seq.slice(ref_rDNA_region[0], ref_rDNA_region[1]),
				}));
				aln_list.sort((a, b) => b.score - a.score);

				const aln = aln_list[0];

				out[genome_idx][ref_idx] = [genome_info.name, [argv_genome_ref1, argv_genome_ref2][ref_idx], aln.start, aln.end, aln_list.length, aln.score].join("\t");

				return aln;
			});

			const ref_result = await Promise.all(task_rDNA_ref);

			return [
				Math.min(...ref_result.map(a => a.start)),
				Math.max(...ref_result.map(a => a.end)),
			];
		});
		const genome_result = await Promise.all(task_genome);

		dataset.rDNA.nChr = argv_nChr_rDNA;

		// @ts-ignore
		dataset.rDNA.region = genome_result[0];
		// @ts-ignore
		dataset.rDNA.region_ref2 = genome_result[1];
		// @ts-ignore
		dataset.rDNA.all_region = genome_result;

		// console.log(dataset.rDNA.region);
		// console.log(out.flat(1).join("\n"));
	}

	{
		// console.log("find telomere");

		const table = [
			[
				"ch",
				...genome_info_list.map(info => ["left start", "left end", "right start", "right end"].map(a => `${info.name} ${a}`)).flat(1),
			],
		];

		const task_genome = genome_info_list.map(async (genome_info, genome_idx) => {
			const task_chr = Object.keys(telomere_data[0]).map(async (nChr) => {
				const chr_idx = Number(nChr) - 1;
				/** @type {any[]} */
				const table_row = table[nChr] = table[nChr] || [
					nChr,
				];
			
				const chr_info = genome_info.chr_list[chr_idx];
				const chr_fasta = chr_info.path;
				
				const sub_task = telomere_data.map(async (ref_telomere, ref_idx) => {
					const telomere = ref_telomere[nChr];
					const arm_task = telomere.map(async (arm_telomere, arm_idx) => {
						const ref_fa = ref1_ref2_seq_list[ref_idx];
						const ref_chr_seq = ref_fa[chr_idx];

						// console.log(ref_idx, chr_idx + 1, chr_fasta);
						
						const aln_list = merge_overlap_alignment((await blastn(chr_fasta, {
							query_seq: ref_chr_seq.slice(arm_telomere[0], arm_telomere[1]),
							subject_loc: (function () {
								const middle = Math.trunc(chr_info.length / 2);
								if (arm_idx == 0) {//left arm
									return [1, middle].join("-");
								}
								else if (arm_idx == 1) {//right arm
									return [middle, chr_info.length].join("-");
								}
								throw new Error("arm_telomere");
							})(),
						})).filter(a => a.align >= 1000 || Math.trunc(a.identity) > 90));
						if (aln_list.length == 0) {
							return null;
						}
						
						aln_list.sort((a, b) => b.score - a.score);;

						// if (arm_idx == 1 && Number(nChr) == 7) {
						// 	console.log(arm_telomere);
						// 	aln_list.forEach(a => console.log(a));
						// }

						if (arm_idx == 0) {
							aln_list.sort((a, b) => {
								return (b.score * (chr_info.length - b.start)) - (a.score * (chr_info.length - a.start));
							});
						}
						else if (arm_idx == 1) {
							aln_list.sort((a, b) => (b.score * b.end) - (a.score * a.end));
						}
						else {
							throw new Error("arm_telomere");
						}

						// // return {
						// // 	loc: sub_telomere,
						// // 	aln: aln_list[0],
						// // };

						return aln_list[0];
					});
					return await Promise.all(arm_task);
				});
				const [_ref1_telo, _ref2_telo] = await Promise.all(sub_task);
				const [left_arm, right_arm] = [
					[
						_ref1_telo[0],
						_ref2_telo[0],
					].filter(a => a),
					[
						_ref1_telo[1],
						_ref2_telo[1],
					].filter(a => a),
				];

				// return {
				// 	nChr: nChr,
				// 	left_arm: left_arm,
				// 	right_arm: right_arm,
				// };

				/** @type {number} */
				let left_end;
				/** @type {number} */
				let right_start;
				
				if (left_arm.length > 1 && left_arm[0].start <= left_arm[1].end && left_arm[0].end >= left_arm[1].start) {
					left_end = Math.max(...left_arm.map(a => a.end));
				}
				else if (left_arm.length == 1) {
					left_end = ss_sort(left_arm, 0)[0].end;
					console.log(genome_idx, nChr, "left", Math.max(...left_arm.map(a => a.end)), left_end);
				}
				else {
					left_end = 1;
					console.log("not found left arm telomere");
				}

				if (right_arm.length > 1 && right_arm[0].start <= right_arm[1].end && right_arm[0].end >= right_arm[1].start) {
					right_start = Math.min(...right_arm.map(a => a.start));
				}
				else if (right_arm.length == 1) {
					right_start = ss_sort(right_arm, 0)[0].start;
					console.log(genome_idx, nChr, "right", Math.min(...right_arm.map(a => a.start)), right_start);
				}
				else {
					right_start = chr_info.length;
					console.log("not found right arm telomere");
				}

				/**
				 * @param {_ref1_telo} aln_list
				 * @param {0|1} arm_idx
				 */
				function ss_sort(aln_list, arm_idx) {
					if (arm_idx == 0) {
						aln_list.sort((a, b) => {
							return (b.score * (chr_info.length - b.start)) - (a.score * (chr_info.length - a.start));
						});
					}
					else if (arm_idx == 1) {
						aln_list.sort((a, b) => (b.score * b.end) - (a.score * a.end));
					}
					return aln_list;
				}
				
				table_row.push(1, left_end);
				table_row.push(right_start, chr_info.length);

				return [
					// [genome_idx, nChr],
					[1, left_end],
					[right_start, chr_info.length],
				];
			});

			return await Promise.all(task_chr);
		});
		const final_telomere = await Promise.all(task_genome);

		// console.log(final_telomere);

		/** @type {{ [nChr:string]: number[][]; }[]} */
		const output_telomere = final_telomere.map((genome_telomere, genome_idx) => {
			/** @type {{ [nChr:string]: number[][]; }} */
			const chr_telomere = {};
			
			genome_telomere.forEach((telomere, chr_idx) => {
				chr_telomere[String(chr_idx + 1)] = telomere;
			});
			
			return chr_telomere;
		});

		// console.log(JSON.stringify(output_telomere, null, "\t"));
		
		// output_telomere.forEach((genome, genome_idx) => {
		// 	Object.keys(genome).forEach(nChr => {
		// 		const [left, right] = genome[nChr];
		// 		console.log([genome_idx, nChr, left, right].flat(1).join("\t"));
		// 	});
		// });

		try {
			const text_log_table = output_telomere.map((genome, genome_idx) => {
				return Object.keys(genome).map(nChr => {
					const [
						[, left_end],
						[right_start, ],
					] = genome[nChr];

					const [
						[
							r1_left,
							r1_right,
						],
						[
							r2_left,
							r2_right,
						]
					] = telomere_data.map(ref_genome => {
						const [
							[, ref_left_end],
							[ref_right_start, ],
						] = ref_genome[nChr];

						return [
							ref_left_end - left_end,
							ref_right_start - right_start,
						];
					});

					const d_left = [r1_left, r2_left].sort((a, b) => Math.abs(a) - Math.abs(b))[0];
					const d_right = [r1_right, r2_right].sort((a, b) => Math.abs(a) - Math.abs(b))[0];

					return [
						genome_idx, nChr, d_left, d_right
					].join("\t");
				}).join("\n");
			}).join("\n");

			fs.writeFileSync(`${dataset.tmp_path}/log_d_telomere.txt`, text_log_table);
		}
		catch (ex) {
		}

		dataset.telomere = output_telomere[0];
		
		dataset.all_telomere = output_telomere;
	}

	dataset.save(argv_dataset_path);
}

async function find_telomere_v1() {
	const argv_ref_genome = String(argv["--ref-genome"] || "");
	const argv_ref_telomere = String(argv["--ref-telomere"] || "");

	const telomere_data = JSON.parse(fs.readFileSync(argv_ref_telomere).toString());

	const ref_seq_list = Object.values(readFasta(argv_ref_genome));

	// const table = [
	// 	[
	// 		"ch",
	// 		...genome_info_list.map(info => ["left start", "left end", "right start", "right end"].map(a => `${info.name} ${a}`)).flat(1),
	// 	],
	// ];

	const task_genome = genome_info_list.map(async (genome_info, genome_idx) => {
		const task_chr = Object.keys(telomere_data).map(async (nChr) => {
			const chr_idx = Number(nChr) - 1;
			// /** @type {any[]} */
			// const table_row = table[chr_idx] || [
			// 	nChr,
			// ];
		
			const chr_info = genome_info.chr_list[chr_idx];
			const chr_fasta = chr_info.path;
			/** @type {[[number, number], [number, number]]} */
			const telomere = telomere_data[nChr];
	
			const sub_task = telomere.map(async sub_telomere => {
				const telomere_seq = ref_seq_list[chr_idx].slice(...sub_telomere);
				
				const aln_list = merge_overlap_alignment(await blastn(chr_fasta, {
					query_seq: telomere_seq,
				}));
				// return {
				// 	loc: sub_telomere,
				// 	aln: aln_list[0],
				// };
				return aln_list[0];
			});
			const [left_arm, right_arm] = await Promise.all(sub_task);
			
			// table_row.push(left_arm.aln.sstart, left_arm.aln.send);
			// table_row.push(right_arm.aln.sstart, right_arm.aln.send);

			// return {
			// 	nChr: nChr,
			// 	left_arm: left_arm,
			// 	right_arm: right_arm,
			// };
			return [
				// [genome_idx, nChr],
				[1, left_arm.end],
				[right_arm.start, chr_info.length],
			];
		});
		
		// table[chr_idx] = table_row;

		return await Promise.all(task_chr);
	});
	const final_telomere = await Promise.all(task_genome);

	// console.log(final_telomere);

	/** @type {{ [nChr:string]: number[][]; }[]} */
	const output_telomere = final_telomere.map((genome_telomere, genome_idx) => {
		/** @type {{ [nChr:string]: number[][]; }} */
		const chr_telomere = {};
		
		genome_telomere.forEach((telomere, chr_idx) => {
			chr_telomere[String(chr_idx + 1)] = telomere;
		});
		
		return chr_telomere;
	});

	console.log(JSON.stringify(output_telomere, null, "\t"));

	dataset.telomere = output_telomere[0];
	// @ts-ignore
	dataset.all_telomere = output_telomere;

	dataset.save(argv_dataset_path);
}

/**
 * @param {string} subject
 * @param {{ query?: string; query_loc?: string; query_seq?: string; subject_loc?: string; }} options
 */
async function blastn(subject, options) {
	const params = [
		"-subject", subject,
		"-evalue", "1e-20",// https://biology.stackexchange.com/questions/19338/what-e-value-blast-cut-off-is-best#:~:text=To%20find%20extremely,seeing%20your%20results.
		"-outfmt", "6",
		// "-subject_loc", subject_loc,
	];
	
	if (options.query) {
		params.push("-query", options.query);
		if (options.query_loc) {
			params.push("-query_loc", options.query_loc);
		}
	}

	if (options.subject_loc) {
		params.push("-subject_loc", options.subject_loc);
	}

	const blastn_proc = child_process.spawn("blastn", params);
	
	if (options.query_seq) {
		blastn_proc.stdin.end(options.query_seq);
	}

	const [
		blastn_list,
	] = await Promise.all([
		parse_blastn_results(blastn_proc.stdout),
		new Promise((resolve, rejects) => {
			blastn_proc.on("exit", function (code, signal) {
				console.error(blastn_proc.spawnargs.join(" "), code, signal);
				if (code || signal) {
					rejects({
						code, signal,
					});
				}
				else {
					resolve({
						code, signal,
					});
				}
			});
		})
	]);

	// blastn_list.sort((a, b) => b.align - a.align);

	return blastn_list;
}

/**
 * @param {BlastnCoord[]} list
 * @param {number} [tolerance]
 */
function merge_overlap_alignment(list, tolerance = 1) {
	/** @type {{ start: number; end: number; score: number; }[]} */
	const region_list = [];
	
	[...list].sort((a, b) => a.sstart - b.sstart).filter(aln => {
		const [start, end] = [aln.sstart, aln.send].sort((a, b) => a - b);

		const region = region_list.find(region => {
			return (
				(start - tolerance) <= region.end &&
				(end + tolerance) >= region.start
			);
		});
		if (region) {
			region.start = Math.min(region.start, start);
			region.end = Math.max(region.end, end);
			region.score = region.score + aln.score;
		}
		else {
			region_list.push({
				start: start,
				end: end,
				score: aln.score,
			});
		}
	});

	return region_list.sort((a, b) => b.score - a.score);
}





