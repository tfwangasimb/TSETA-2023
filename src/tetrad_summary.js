// @ts-check

const argv_NCO_MODE = (() => {
	const regexp = /--nco-mode=(.*)/;
	const argv_ = process.argv.map(s => s.match(regexp)).filter(a => a).slice(-1)[0];
	if (argv_) {
		const argv = argv_[1];//(typeof process == "object" && process != null) ? 123 : null;
		process.env.NCO_MODE = argv;
		console.log("nco-mode =", argv);
	}
})();

const fs = require("fs");
const child_process = require("child_process");

const { romanize } = require("./tools/romanize.js");
const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { argv_parse } = require("./util.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");

const { FinalTableRow } = require("./final_table_format.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const argv_closeCOsMinDistance = Number(argv["-min-co"]);
const argv_output_summary_only = !!argv["--output-summary-only"];

const argv_snp_count_window = Number(argv["--snp-count-window"] | 0);

const argv_skip_exist = !!argv["--skip-exist"];
const VERBOSE = !!argv["--verbose"];

if (argv_skip_exist) {
	console.warn("argv_skip_exist");
}

const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}
dataset.load_rDNA_info();

if (Number.isSafeInteger(argv_closeCOsMinDistance)) {
	if (!dataset.crossover) {
		dataset.crossover = {};
	}
	dataset.crossover.closeCOsMinDistance = argv_closeCOsMinDistance;
	dataset.save(argv_dataset_path);
}
else {
	console.error("Required: -min-co (number >= 0)");
	console.error("-min-co:", "The closest distance between two adjacent crossovers");
}
dataset.load_GC_content();

const genome_info_list = dataset.loadGenomeInfoList();

const output_prefix = (dataset.name);

let cmds = [];

main();

async function main() {
	if (dataset.mode != "tetrad") {
		console.warn("switch tetrad mode");
		dataset.mode = "tetrad";
	}

	try {
		if (!argv_output_summary_only) {
			let fin_list = [];
			
			//for each chromosome
			let tasks = genome_info_list[0].chr_list.map(async function (chrName, chrIdx) {
				const nChr = chrIdx + 1;
				const seq = `${dataset.output_path}/mafft_ch${nChr}.fa`;

				const path_to_segfile = `${dataset.output_path}/${output_prefix}_seg_ch${nChr}.txt`;
				const path_to_co_list = `${dataset.output_path}/${output_prefix}_co_ch${nChr}_co.json`;
				const path_to_nco_list = `${dataset.output_path}/${output_prefix}_co_ch${nChr}_nco.json`;

				fin_list[chrIdx] = ["init"];

				if (!argv_skip_exist ||
					(argv_skip_exist && fs.existsSync(path_to_segfile))
				) {
					let multi_align_to_segfile_args =  [
						"--max-old-space-size=4096",
						`${__dirname}/multi_align_to_segfile.js`,
						"-dataset", argv_dataset_path,
						"-i", seq,
						"-chr", nChr,
						"--output-prefix", `${output_prefix}_seg_ch${nChr}`,
					];
					try {
						await spawnNodeAsync(multi_align_to_segfile_args);

						if (VERBOSE) {
							console.log("nChr", nChr, "segfile");
						}

						fin_list[chrIdx].push("segfile");
					}
					catch (ex) {
						console.error(ex);
						console.error("error", nChr, ["node", ...multi_align_to_segfile_args].join(" "));
						return;
					}
				}
				else {
					console.log("found:", path_to_segfile);
				}

				if (!argv_skip_exist ||
					argv_skip_exist && (fs.existsSync(path_to_co_list) && fs.existsSync(path_to_nco_list))
				) {
					let crossover_ars = [
						"--max-old-space-size=4096",
						`${__dirname}/crossover.js`,
						"-dataset", argv_dataset_path,
						"--segfile", path_to_segfile,
						"--output-prefix", `${output_prefix}_co_ch${nChr}`
					];
					try {
						await spawnNodeAsync(crossover_ars);

						if (VERBOSE) {
							console.log("nChr", nChr, "crossover");
						}

						fin_list[chrIdx].push("crossover");
					}
					catch (ex) {
						console.error(ex);
						console.error("error", nChr, ["node", ...crossover_ars].join(" "));
						return;
					}
				}
				else {
					console.log("found:", [
						path_to_co_list,
						path_to_nco_list,
					]);
				}

				for (let optional of ["", "--point_filter"]) {
					let tetrad_chr_summary_args = [
						"--max-old-space-size=8192",
						`${__dirname}/tetrad_chr_summary.js`,
						"-dataset", argv_dataset_path,
						"-chr", nChr,
						"--seq", seq,
						"--co-list", path_to_co_list,
						"--nco-list", path_to_nco_list,
						"--snp-count-window", `${argv_snp_count_window}`,
						"--output-prefix", `${output_prefix}_ch${nChr}${optional == "--point_filter" ? "_filtered" : ""}`,
						optional
					];
					try {
						await spawnNodeAsync(tetrad_chr_summary_args);

						if (VERBOSE) {
							console.log("nChr", nChr, "summary");
						}

						fin_list[chrIdx].push("summary");
					}
					catch (ex) {
						console.error(ex);
						console.error("error", nChr, ["node", ...tetrad_chr_summary_args].join(" "));
						return;
					}
				}
			});
			await Promise.all(tasks);

			if (VERBOSE) {
				console.log(fin_list);
			
				console.log(cmds.join("\n"));
			}
		}
	}
	catch (ex) {
		console.log(cmds.join("\n"));
		throw ex;
	}

	{
		let template_html = fs.readFileSync(`${__dirname}/template_viewer.html`).toString();
		let analyser_js = fs.readFileSync(`${__dirname}/analyser.js`).toString();
		let web_ui_js = fs.readFileSync(`${__dirname}/web_ui.js`).toString();
		let html2canvas_js = fs.readFileSync(`${__dirname}/html2canvas.js`).toString();
		let gff_js = fs.readFileSync(`${__dirname}/gff.js`).toString();

		const all_co = load_all_co("co");
		const all_nco = load_all_co("nco");
		{
			if (all_co) {
				const all_co_list = all_co.reduce((all, list) => all.concat(list), []);
				const header_map = {
					chr: "Chromosome number",
					chr_len: "Chromosome Length (bp)",
					pos: "Position (bp)",
					GCasso_tract_len: "GC tract length (bp)",
					GCasso_marker: "GC tract marker",
					snp_start_out: "snp_start_out",
					snp_start_in: "snp_start_in",
					snp_end_in: "snp_end_in",
					snp_end_out: "snp_end_out",
					type: "CO type",
					why_remove: "why_remove",
				};
				make_co_table(all_co_list, header_map, "co");
			}
			if (all_nco) {
				const all_nco_list = all_nco.reduce((all, list) => all.concat(list), []);
				const header_map = {
					chr: "Chromosome number",
					chr_len: "Chromosome Length (bp)",
					pos: "Position (bp)",
					GCasso_tract_len: "GC tract length (bp)",
					GCasso_marker: "GC tract marker",
					snp_start_out: "snp_start_out",
					snp_start_in: "snp_start_in",
					snp_end_in: "snp_end_in",
					snp_end_out: "snp_end_out",
					type: "NCO type",//type(NCO|2NCO)=InDels/SNPs
					why_remove: "why_remove",
					nco_type: "NCO type",
				};
				make_co_table(all_nco_list, header_map, "nco");
			}
		}

		dataset.genome_info_list = genome_info_list;
		dataset.crossover_list = all_co;
		dataset.non_crossover_list = all_nco;
		
		fs.writeFile(`${dataset.output_path}/${output_prefix}_co.json`, JSON.stringify(all_co, null, "\t"), function (err) {
			console.error(err);
		});
		fs.writeFile(`${dataset.output_path}/${output_prefix}_nco.json`, JSON.stringify(all_nco, null, "\t"), function (err) {
			console.error(err);
		});

		{
			dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
				const nChr = chrIdx + 1;
				return `mafft_ch${nChr}.fa`;
			});
			const output_path = `${dataset.output_path}/debug_${output_prefix}.html`;

			output_html(template_html, output_path, false);
		
			if (VERBOSE) {
				console.log("output debug viewer:", output_path);
			}
		}

		//single html
		{
			let all_seq = load_all_seq();
			dataset.results = genome_info_list[0].chr_list.map(function (chrName, chrIdx) {
				return all_seq[chrIdx];
			});
			const output_path = `${dataset.output_path}/${output_prefix}.html`;

			output_html(template_html, output_path, true);
		
			console.log("output viewer:", output_path);
		}
		
		function output_html(input_html, output_path, inline_script) {
			let str_json = JSON.stringify(dataset);// new copy

			let output_html = input_html.replace(`<script id="dataset.json" type="application/json"></script>`, `<script id="dataset.json" type="application/json">${str_json}</script>`);

			if (inline_script) {
				output_html = output_html.replace(`<script src="html2canvas.js"></script>`, `<script>${html2canvas_js}</script>`);
				output_html = output_html.replace(`<script src="gff.js"></script>`, `<script>${gff_js}</script>`);

				output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script>${analyser_js}</script>`);
				output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script>${web_ui_js}</script>`);
				fs.writeFile(`${dataset.output_path}/output.json`, str_json, function (err) {
					if (err) {
						throw err;
					}
				});
			}
			else {
				output_html = output_html.replace(`<script src="html2canvas.js"></script>`, `<script src="../src/html2canvas.js"></script>`);
				output_html = output_html.replace(`<script src="gff.js"></script>`, `<script src="../src/gff.js"></script>`);
				output_html = output_html.replace(`<script src="analyser.js"></script>`, `<script src="../src/analyser.js"></script>`);
				output_html = output_html.replace(`<script src="web_ui.js"></script>`, `<script src="../src/web_ui.js"></script>`);
			}

			fs.writeFile(output_path, output_html, function (err) {
				if (err) {
					throw err;
				}
			});
		}
	}

	output_final_table();
}

function load_all_seq() {
	let all_seq = [];
	genome_info_list[0].chr_list.forEach(function (chrName, chrIdx) {
		const nChr = chrIdx + 1;
		const seq = `${dataset.output_path}/mafft_ch${nChr}.fa`;
		
		const input_fasta = readFasta(seq);

		all_seq[chrIdx] = input_fasta;
	});
	return all_seq;
}

/**
 * @param {"co"|"nco"} type
 */
function load_all_co(type) {
	const chr_list = genome_info_list[0].chr_list;

	return chr_list.map(chr => {
		const nChr = chr.index;
		const co_list = _load_co_list(`${dataset.output_path}/${output_prefix}_co_ch${nChr}_${type}.json`);
		return co_list;
	});
	function _load_co_list(co_list_filename) {
		if (!fs.existsSync(co_list_filename)) {
			throw new Error(`Not found: ${co_list_filename}`);
		}

		const stat = fs.statSync(co_list_filename);
		if (stat.size <= 0) {
			throw new Error(`Empty file: ${co_list_filename}`);
		}

		try {
			const text = fs.readFileSync(co_list_filename).toString();
			const list = JSON.parse(text);
			list.forEach(co => {
				if (co.before) {
					co.before = co.before.split(",").map(n => Number(n));
				}
				if (co.after) {
					co.after = co.after.split(",").map(n => Number(n));
				}
			});
			return list;
		}
		catch (ex) {
			console.error("error:", "read or parse file:", co_list_filename);
			console.error(ex);
			return [];
		}
	}
}

function output_final_table() {
	const ordered_nChr_list = [];
	// const ref1_nChr_len_list = genome_info_list[0].chr_list.map(a => a.length);
	// const ref2_nChr_len_list = genome_info_list[1].chr_list.map(a => a.length);
	{
		// Warn: chr.index => chr.order
		const aaa = genome_info_list.map(cc => cc.chr_list.sort((a, b) => b.length - a.length).map(cc => Number(cc.index)));
		const bbb = aaa.map(a => a.join("\t"));
		const all_same_order = bbb.every(a => a == bbb[0]);

		if (!all_same_order) {
			console.warn(`order by ${dataset.ref} chr length ??`);
		}
		else {
			// chromosome 順序等於 chromosome length 長到短
			console.warn(`order by ${dataset.ref} chr length`);
		}

		ordered_nChr_list.push(...aaa[0]);
	}
	const final_table_all = _sub_final_table_1(ordered_nChr_list, "");
	const final_table_filtered = _sub_final_table_1(ordered_nChr_list, "_filtered");
	
	const output_head_list = (new FinalTableRow()).outputTableHeader(dataset.parental_list[0], dataset.parental_list[1]);

	const final_table = [
		final_table_all[0],// Chromosome
		final_table_all[1],// dad Chromosome len
		final_table_all[2],// mom Chromosome len
	];	
	for (let i = 3; i < output_head_list.length; ++i) {
		final_table.push(final_table_all[i]);
		final_table_filtered[i][0] += " (all - rDNA - telomere)";
		final_table.push(final_table_filtered[i]);
	}

	let final_text = final_table.map(row => row.join("\t")).join("\n");

	// let final_text = "";
	// final_text += output_head.join("\t") + "\n";
	// final_text += final_table.map(row => output_head.map(key => row[key]).join("\t")).join("\n");

	let output_path = `${dataset.output_path}/${output_prefix}_final_table.txt`;

	fs.writeFileSync(output_path, final_text);
	
	console.log("output final table:", output_path);
}

function _sub_final_table_1(ordered_nChr_list, sub_type) {
	const total = new FinalTableRow();
	const output_head = total.outputTableHeader(dataset.parental_list[0], dataset.parental_list[1]);

	/** @type {string[][]} */
	let final_table = [
		output_head,
	];

	ordered_nChr_list.map(function (nChr) {
		const tetrad_chr_summary = `${dataset.tmp_path}/table/${output_prefix}_ch${nChr}${sub_type}_summary.txt`;

		const text = fs.readFileSync(tetrad_chr_summary).toString();

		// @ts-ignore
		// let table = table_to_object_list(tsv_parse(text), output_head, { start_row: 1 });
		// output_head.forEach(k => {
		// 	table.forEach(row => {
		// 		total[k] += Number(row[k]);
		// 	});
		// });
		// final_table.push(...table);
		let objList = table_to_object_list(tsv_parse(text), 0); // read header




		// if (chrIdx == 0) {
		// 	console.table(objList);
		// 	console.table(output_head);
		// }
		total._keys().forEach((inKey, i) => {
			const outKey = output_head[i];
			objList.forEach(row => {
				total[inKey] += Number(row[outKey]);
			});
		});
		const chrPropName = output_head[0];
		objList.forEach(row => {
			row[chrPropName] = romanize(Number(row[chrPropName]));
		});
		const table = objList.map(row => output_head.map(k => row[k]));
		final_table.push(...table);
	});
	{
		const a = Number(total.simple_CO);
		const b = Number(total.CO_NCO);
		const c = Number(total.NCO_2p);
		total.CO_div_nCO = ((a + b) / (a + b + c)).toFixed(2);
	}

	total.Chromosome = "total";
	// @ts-ignore
	final_table.push(Object.values(total.outputTableRow(dataset.parental_list[0], dataset.parental_list[1])));

	// rotate -90deg;
	final_table = final_table[0].map((_, i) => final_table.map((__, j) => final_table[j][i]));
	return final_table;
}

/**
 * @param {*} _co_list
 * @param {*} header_map
 * @param {"co"|"nco"} type
 */
function make_co_table(_co_list, header_map, type) {
	const ordered_size_list = (function () {
		const li = [...genome_info_list[0].chr_list];
		// order by chr.length
		return [NaN, ...li.sort((a, b) => b.length - a.length).map(a => a.length)];
	})();

	const co_list = [..._co_list].sort((a, b) => ordered_size_list[b.nChr] - ordered_size_list[a.nChr]);
	
	const headers = Object.values(header_map);
	const co_table = [];
	co_table.push(headers);
	co_list.map(row => {
		const entries = Object.keys(header_map).map(from_key => {
			// const to_key = header_map[from_key];
			const value = row[from_key];
			// const entry = [to_key, value];
			// return entry;
			return value;
		});
		co_table.push(entries);
	});
	const text_co_table = co_table.map(row => row.join("\t")).join("\n");

	const output_path = `${dataset.output_path}/${output_prefix}_final_${type}_table.txt`;
	fs.writeFile(output_path, text_co_table, function (err) {
		if (err) {
			console.error(err);
		}
		console.log("output final table:", output_path);
	});
}

function spawnNodeAsync(args) {
	cmds.push(["node", ...args].join(" "));

	if (VERBOSE) {
		console.log("node", ...args);
	}

	return new Promise(function (resolve, reject) {
		let proc = child_process.spawn("node", args);
		proc.on("exit", function (code, signal) {
			if (code || signal) {
				/** @type {Buffer[]} */
				const stdout = [];
				/** @type {Buffer[]} */
				const stderr = [];

				proc.stdout.on("data", a => stdout.push(a));
				proc.stderr.on("data", a => stderr.push(a));

				proc.stdout.on("close", _ => {
					console.error("node " + args.join(" "));
					console.error(Buffer.concat(stdout).toString());
				});
				proc.stderr.on("close", _ => {
					console.error("node " + args.join(" "));
					console.error(Buffer.concat(stderr).toString());
				});

				reject({
					code, signal,
				});
			}
			else {
				resolve();
			}
		});
	});
}

