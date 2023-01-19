
// @ts-check

const fs = require("fs");

const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { Dataset, GenomeInfo } = require("./dataset.js");
const { GC_Content_Data } = require("./GC_content_util.js");
const { readFasta } = require("./fasta_util.js");

const argv = argv_parse(process.argv);

if (process.argv.includes("--help")) {
	console.log([
		"--window, -w",
		"--minus, -m",
		"--no-telo",
		"--no-centro",
		"--verbose",
	].join("\n"));
}

const argv_window_size = Number(argv["--window"] || argv["-w"]) | 0;
const argv_AT_island_GC_minus = (Number(argv["--minus"] || argv["-m"]) | 0);//6 // AT_island => gc(window) <= (gc(all) - minus)

const argv_no_telo = !!argv["--no-telo"];
const argv_no_centro = !!argv["--no-centro"];

const VERBOSE = !!argv["--verbose"];

/**
A = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "A").length;
T = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "T").length;
G = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "G").length;
C = [...seq_list[0]].filter(a => a != "-").slice(0, 5000).filter(a => a == "C").length;
(G + C) / (A + T + G + C) * 100
 */

class AT_island_Data {
	constructor() {
		this.chr = 0;
		this.start = 0;
		this.end = 0;
		this.length = 0;
	}
}

main();

function main() {
	if (VERBOSE) {
		console.log({
			argv_AT_island_GC_minus, argv_window_size,
		});
	}

	const argv_dataset_path = String(argv["-dataset"] || "");
	const dataset = Dataset.loadFromFile(argv_dataset_path);
	if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
		throw new Error("dataset name no match file name");
	}

	program_log(`${dataset.name}.log.txt`, "start");

	if (argv_no_centro) {
		if (dataset.centromere) {
			dataset._old_centromere = dataset.centromere;
			delete dataset.centromere;
		}
	}
	if (argv_no_telo) {
		if (dataset.telomere) {
			dataset._old_telomere = dataset.telomere;
			delete dataset.telomere;
		}
	}

	const genome_info_list = dataset.loadGenomeInfoList();

	const ref1_name = dataset.ref;
	const ref2_name = dataset.parental_list[1];
	
	const gc_content_filepath = `${dataset.output_path}/${[ref1_name, ref2_name].join("_")}_GC_content.txt`;
	const at_island_filepath = `${dataset.output_path}/${ref1_name}_AT_island.txt`;

	const ref1_chr_list = genome_info_list[0];
	
	dataset.GC_Content_value = genome_info_list.map(_ => 0);
	dataset.min_gc_content = genome_info_list.map(_ => 0);

	let {
		all_gc_table: ref1_all_gc_table,
		all_island_table: ref1_all_island_table,
		all_gc_content: ref1_all_gc_content,
		at_island_maxmum_gc: ref1_at_island_maxmum_gc,
	} = make_table(ref1_chr_list);
	//
	dataset.GC_Content_value[0] = ref1_all_gc_content;
	dataset.min_gc_content[0] = ref1_at_island_maxmum_gc;
	
	let ref2_all_gc_table = null;
	if (ref2_name) {
		const ref2_chr_list = genome_info_list[1];
		let {
			all_gc_table: _ref2_all_gc_table,
			all_island_table: _ref2_all_island_table,
			all_gc_content: ref2_all_gc_content,
			at_island_maxmum_gc: ref2_at_island_maxmum_gc,
		} = make_table(ref2_chr_list);
		//
		dataset.GC_Content_value[1] = ref2_all_gc_content;
		dataset.min_gc_content[1] = ref2_at_island_maxmum_gc;
		
		_ref2_all_gc_table.forEach(group => {
			group.forEach(gc => {
				const ref2_chr_info = ref2_chr_list.chr_map[gc.chr];
				gc.name = ref2_name;
				gc.chr = ref2_chr_info.index;
				gc.start = gc.start + 1;
				gc.end = Math.min(gc.end + 1, ref2_chr_info.length);
			});
		});

		ref2_all_gc_table = _ref2_all_gc_table;
	}

	{
		ref1_all_gc_table.forEach(group => {
			group.forEach(gc => {
				const ref1_chr_info = ref1_chr_list.chr_map[gc.chr];
				gc.name = ref1_name;
				gc.chr = ref1_chr_info.index;
				gc.start = gc.start + 1;
				gc.end = Math.min(gc.end + 1, ref1_chr_info.length);
			});
		});
		
		const output_gc_table_header = [
			"name", "chr", "start", "end", "gc"
		];
		let all_gc_list = ref2_all_gc_table ? [].concat(...ref1_all_gc_table, ...ref2_all_gc_table) : [].concat(...ref1_all_gc_table);
		let text_all_gc_list = all_gc_list.map(row => output_gc_table_header.map(key => row[key]).join("\t")).join("\n");

		fs.writeFileSync(gc_content_filepath, text_all_gc_list);
	}

	{
		ref1_all_island_table.forEach(group => {
			group.forEach(island => {
				island.chr = ref1_chr_list.chr_map[island.chr].index;
				island.start = island.start + 1;
				island.end = island.end + 1;
			});
		});

		const output_at_table_header = [
			"chr", "start", "end", "length"
		];
		let all_island_list = [].concat(...ref1_all_island_table);
		let text_all_island_list = all_island_list.map(row => output_at_table_header.map(key => row[key]).join("\t")).join("\n");
		
		fs.writeFileSync(at_island_filepath, text_all_island_list);
	}

	dataset.all_centromere = genome_info_list.map(_ => []);
	dataset.all_telomere = genome_info_list.map(_ => []);

	const {
		centro: ref1_cen,
		telo: ref1_tel
	} = find_centro_telo(genome_info_list, ref1_all_island_table, ref1_chr_list);
	//
	if (!argv_no_centro) {
		ref1_cen.forEach(a => {
			dataset.centromere[a.chr] = [a.start, a.end];
			dataset.all_centromere[0][a.chr] = [a.start, a.end];
		});
	}
	if (!argv_no_telo) {
		ref1_tel.forEach(([a, b]) => {
			dataset.telomere[a.chr] = [
				[a.start, a.end],
				[b.start, b.end],
			];
			dataset.all_telomere[0][a.chr] = [
				[a.start, a.end],
				[b.start, b.end],
			];
		});
		console.warn("argv_no_telo:", argv_no_telo);
	}

	{
		genome_info_list.slice(1).forEach((chr_list, _ssIdx) => {
			const sIdx = _ssIdx + 1;
			
			const {
				all_island_table,
				all_gc_content: _all_gc_content,
				at_island_maxmum_gc: _at_island_maxmum_gc,
			} = make_table(chr_list);
			//
			dataset.GC_Content_value[sIdx] = _all_gc_content;
			dataset.min_gc_content[sIdx] = _at_island_maxmum_gc;

			all_island_table.forEach(group => {
				group.forEach(island => {
					island.chr = chr_list.chr_map[island.chr].index;
					island.start = island.start + 1;
					island.end = island.end + 1;
				});
			});
			const { centro, telo } = find_centro_telo(genome_info_list, all_island_table, chr_list);
			//
			if (!argv_no_centro) {
				centro.forEach(a => {
					dataset.all_centromere[sIdx][a.chr] = [a.start, a.end];
				});
			}
			if (!argv_no_telo) {
				telo.forEach(([a, b]) => {
					dataset.all_telomere[sIdx][a.chr] = [
						[a.start, a.end],
						[b.start, b.end],
					];
				});
			}
		});
	}

	dataset.GC_Content_filePath = gc_content_filepath;
	dataset.GC_Content_window = argv_window_size;
	dataset.GC_Content_m = argv_AT_island_GC_minus;

	const output_dataset = JSON.stringify(dataset, null, "\t");

	fs.writeFileSync(argv_dataset_path, output_dataset);
	
	console.log("save:", argv_dataset_path);
	
	console.log("next step:", "Sequential (5' to 3') slicing of the query chromosome into smaller fragments (~10 kb)");
	console.log("command:", `node ${__dirname}/slice_all.js -dataset ${argv_dataset_path}`);
	
	program_log(`${dataset.name}.log.txt`, "exit");
}

function find_centro_telo(genome_info_list, all_island_table, chr_list) {
	let centro = [];
	let telo = [];

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		const chrIdx = nChr - 1;
		const AT_desc_list = [...all_island_table[chrIdx]].sort((a, b) => b.length - a.length);

		if (!argv_no_centro) {
			const cen_range = AT_desc_list[0];
			try {
				let [at1, at2] = [AT_desc_list[0], AT_desc_list[1]].sort((a, b) => a.start - b.start);
				//console.log("ch", nChr, "at1, 500bp, at2", at1.start, at1.end, at2.start, at2.end);
				// 合併 2 個鄰近的 AT-island，2 個 AT-island 間最多能有 1 個 window (QM6a ChIV: 1482500-1559500,1560000-1659000)
				if ((at2.start - at1.end) <= Math.abs(at1.end - at1.start)) {
					cen_range.start = at1.start;
					cen_range.end = at2.end;
					cen_range.length = at2.end - at1.start;
				}
			}
			catch (ex) {
				console.error({
					ex,
					AT_desc_list
				});
			}
			//console.log("cen", cen_range.start, cen_range.end, cen_range.length);
			if (cen_range) {
				centro.push(cen_range);
			}
			else {
				centro.push(new AT_island_Data());
			}
		}

		if (!argv_no_telo) {
			try {
				let tel1 = Object.assign({}, all_island_table[chrIdx][0], {
					start: 1,
				});

				let tel2 = Object.assign({}, all_island_table[chrIdx][all_island_table[chrIdx].length - 1], {
					end: chr_list.chr_list[all_island_table[chrIdx][all_island_table[chrIdx].length - 1].chr - 1].length,
				});

				telo.push([tel1, tel2]);
			}
			catch (ex) {
				console.error(`nChr: ${nChr}`, ex);
				telo.push([0, 0]);
			}
		}
	}
	
	return {
		centro, telo,
	};
}

/**
 * 
 * @param {GenomeInfo} genome_info
 */
function make_table(genome_info) {
	console.log("input:", genome_info.name);

	const input_seq = genome_info.loadFasta();
	
	//console.log("input fa", Object.keys(input_seq));
	
	const all_gc_content = get_all_gc_content(input_seq);
	const at_island_maxmum_gc = all_gc_content - argv_AT_island_GC_minus;
	
	console.log("GC%", all_gc_content);
	console.log("AT-island maxmum GC%", at_island_maxmum_gc);

	/** @type {GC_Content_Data[][]} */
	let all_gc_table = [];
	/** @type {AT_island_Data[][]} */
	let all_island_table = [];

	Object.keys(input_seq).forEach(function (name) {
		/** @type {GC_Content_Data[]} */
		let gc_table = [];
		/** @type {AT_island_Data[]} */
		let island_table = [];

		let seq = input_seq[name];
		let counter = {
			A: 0,
			T: 0,
			G: 0,
			C: 0,
		};
		
		let i = 0, length = 0, start = 0;
		for (; i < seq.length; ++i) {
			const v = seq[i].toUpperCase();
			
			++counter[v];

			++length;
			if (length >= argv_window_size) {
				let gc = counter.G + counter.C;
				let gc_content = 100 * gc / (counter.A + counter.T + gc);

				gc_table.push({
					chr: name,
					start: start,
					end: i,
					gc: gc_content,
				});

				if (gc_content <= at_island_maxmum_gc) {
					island_table.push({
						chr: name,
						start: start,
						end: i,
						length: length,
					});
				}

				if ((i - start + 1) > argv_window_size) {
					console.error(start, i, argv_window_size);
					throw new Error("if ((start + 1) > argv_window_size) {");
				}

				length = 0;
				start = i + 1;//next
				counter = {
					A: 0,
					T: 0,
					G: 0,
					C: 0,
				};//clear
			}
		}
		if (length < argv_window_size) {
			let gc = counter.G + counter.C;
			let gc_content = 100 * gc / (counter.A + counter.T + gc);

			gc_table.push({
				chr: name,
				start: start,
				end: i,
				gc: gc_content,
			});
		}
		all_gc_table.push(gc_table);
		all_island_table.push(island_table);

		// console.log(name, "gc_table", gc_table.length);
		// console.log(name, "island_table", island_table.length);
	});

	// {
	// 	const output_gc_table_header = [
	// 		"chr", "start", "end", "gc"
	// 	];
	// 	let all_gc_list = [].concat(...all_gc_table);
	// 	let text_all_gc_list = all_gc_list.map(row => output_gc_table_header.map(key => row[key]).join("\t")).join("\n");
	// 	fs.writeFileSync(`${dataset.output_path}/${output_prifix}_GC_content.txt`, text_all_gc_list);
	// }

	function merge_island() {
		return all_island_table.map(island_list => {
			let megered_list = [];
			let list = island_list.slice(1);

			let current_island = island_list[0];
			list.forEach(next_island => {
				if ((current_island.end + 1) == next_island.start) {
					current_island.end = next_island.end;
					current_island.length = current_island.length + next_island.length;
				}
				else {
					megered_list.push(current_island);

					current_island = next_island;
				}
			});

			///.........

			return megered_list;
		});
	}

	//console.log(all_island_table[0].length);
	all_island_table = merge_island();
	
	// {
	// 	//console.log(all_island_table[0].length);
	// 	const output_at_table_header = [
	// 		"chr", "start", "end", "length"
	// 	];
	// 	let all_island_list = [].concat(...all_island_table);
	// 	let text_all_island_list = all_island_list.map(row => output_at_table_header.map(key => row[key]).join("\t")).join("\n");
	// 	//fs.writeFileSync(`${dataset.output_path}/${output_prifix}_AT_island.txt`, text_all_island_list);
	// }
	
	return {
		all_gc_table, all_island_table,

		all_gc_content, at_island_maxmum_gc,
	};
}

function get_all_gc_content(input_seq) {
	let counter = {
		A: 0,
		T: 0,
		G: 0,
		C: 0,
	};

	Object.keys(input_seq).forEach(function (name) {
		let seq = input_seq[name];
		for (let i = 0; i < seq.length; ++i) {
			const v = seq[i].toUpperCase();
			++counter[v];
		}
	});

	//console.log(counter);

	let all_gc = counter.G + counter.C;
	let all_atgc = counter.A + counter.T + all_gc;
	let all_gc_content = 100 * all_gc / all_atgc;

	// console.log("GC content", {
	// 	all_gc,
	// 	all_atgc,
	// 	all_gc_content,
	// 	"%": all_gc_content.toFixed(2),
	// });

	return all_gc_content;
}


