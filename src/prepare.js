// @ts-check

const fs = require("fs");
const Path = require("path");

const { inputFile, inputDirectory, inputNumber, inputText, inputSelect } = require("./interact-util.js").userInput;
const { argv_parse, array_groupBy, program_log } = require("./util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { GenomeDataSet, Dataset, GenomeInfo, MyDataInfo, MyFileInfo } = require("./dataset.js");
const { table_to_object_list, tsv_parse } = require("./tsv_parser.js");

const argv = argv_parse(process.argv);

let has_argv_dataset_path = false;


if (typeof argv["--init-table"] == "string") {
	const pp = argv["--init-table"];
	if (!fs.existsSync(pp)) {
		const text = [
			`taskName`, `tetrad_1`,
			`mode`,     `tetrad`,
			`ref`,      `dad`,
			`parental`, `dad`,  `dad.fa`,
			`parental`, `mom`,  `mom.fa`,
			`progeny`,  `F1-A`, `F1-A.fa`,
			`progeny`,  `F1-B`, `F1-B.fa`,
			`progeny`,  `F1-C`, `F1-C.fa`,
			`progeny`,  `F1-D`, `F1-D.fa`,
		].join("\n");
		fs.writeFileSync(pp, text);
	}
}
// else if (argv["--init-dataset"]) {
// }
else {
	main();
}

async function main() {
	const dataset = await load_or_input();

	program_log(`${dataset.name}.log.txt`, "start");

	if (!fs.existsSync(`${dataset.output_path}/`)) {//make output directory
		fs.mkdirSync(`${dataset.output_path}/`);
	}
	const loaded_data = make_genome_info(dataset);

	explode_genome(dataset, loaded_data);

	{
		dataset.info = dataset.info || new MyDataInfo();

		dataset.info.parental = dataset.info.parental || {};
		Object.keys(dataset.parental).forEach(k => {
			dataset.info.parental[k] = MyFileInfo.getInfo(dataset.parental[k]);
		});

		dataset.info.progeny = dataset.info.progeny || {};
		Object.keys(dataset.progeny).forEach(k => {
			dataset.info.progeny[k] = MyFileInfo.getInfo(dataset.progeny[k]);
		});
	}
	
	const output_filename = `${(dataset.name)}.json`;
	fs.writeFileSync(output_filename, JSON.stringify(dataset, null, "\t"));
	console.log("save:", output_filename);
	
	console.log("skip GC content %:", "Sequential (5' to 3') slicing of the query chromosome into smaller fragments (~10 kb)");
	console.log("command:", `node ${__dirname}/slice_all.js -dataset ${output_filename}`);
	
	console.log("next step:", "detect GC content %");
	console.log("command:", `node ${__dirname}/make_GC_AT_table.js -dataset ${output_filename} -w 500 -m 6 --no-telo --no-centro`);

	console.log("next step:", "detect GC content %, AT-rich blocks, centromere, telomeres");
	console.log("command:", `node ${__dirname}/make_GC_AT_table.js -dataset ${output_filename} -w 500 -m 6`);
	
	program_log(`${dataset.name}.log.txt`, "exit");
}

async function load_or_input() {
	let genome_dataset;
	try {
		if (argv["-dataset"]) {
			const argv_dataset_path = String(argv["-dataset"] || "");
			if (fs.existsSync(argv_dataset_path)) {
				genome_dataset = Dataset.loadFromFile(argv_dataset_path);
				has_argv_dataset_path = true;

				genome_dataset.parental_list ??= Object.keys(genome_dataset.parental);
				genome_dataset.progeny_list ??= Object.keys(genome_dataset.progeny);
			}
		}
		else if (argv["-table"]) {
			const path_to_table = String(argv["-table"] || "");
			if (fs.existsSync(path_to_table)) {
				const tab = table_to_object_list(tsv_parse(fs.readFileSync(path_to_table).toString()), ["property", "value_1", "value_2"]);
				// console.table(tab);
				genome_dataset = new GenomeDataSet();
				genome_dataset.name = tab.find(a => a.property == "taskName").value_1;
				genome_dataset.mode = tab.find(a => a.property == "mode").value_1;
				genome_dataset.ref = tab.find(a => a.property == "ref").value_1;

				genome_dataset.parental = Object.fromEntries(tab.filter(a => a.property == "parental").map(a => [a.value_1, a.value_2]));;
				genome_dataset.progeny = Object.fromEntries(tab.filter(a => a.property == "progeny").map(a => [a.value_1, a.value_2]));;
				
				genome_dataset.parental_list = Object.keys(genome_dataset.parental);
				genome_dataset.progeny_list = Object.keys(genome_dataset.progeny);
			}
			else {
				console.log("404:", path_to_table);
			}
		}
		genome_dataset.ref ??= genome_dataset.parental_list[0];
	}
	catch (ex) {
		genome_dataset = null;
		console.error(ex);
	}

	if (!genome_dataset) {
		genome_dataset = await userInput_genomeDataset();
	}

	return genome_dataset;
}

/**
 * @returns {Promise<GenomeDataSet>}
 */
async function userInput_genomeDataset() {
	let set_genomeName = new Set();

	let taskName = await inputText("Task name");

	/** @type {"tetrad"|"SNP"} */
	let mode = await inputSelect("mode", ["tetrad", "SNP"], "tetrad");
	
	// current version
	// genome \ mode | tetrad |  SNP
	// n parental    |    = 2 |  = 1
	// n progeny     |    = 4 | >= 1
	
	// future version
	// genome \ mode | tetrad |  SNP
	// n parental    |    = 2 | >= 1
	// n progeny     |   >= 4 | >= 1
	
	//let numParental = mode == "tetrad" ? 2 : await inputNumber("number of parental", { min: 1, max: 2, default: 1 });
	const numParental = await (async function () {
		if (mode == "tetrad") {
			return 2;
		}
		else {
			await inputNumber("number of parental", { min: 1, max: 2, default: 1, });
		}
	})();
	/** @type {string[]} */
	let parental_list = [];
	/** @type {{[genomeName:string]:string}} */
	let parental_map = {};
	let refName = await inputText(`${mode == "tetrad" ? "ref (parental 1)" : "ref"} name`);
	let refFasta = await inputFile(`${mode == "tetrad" ? "ref (parental 1)" : "ref"} genome fasta file`, `${refName}.genome.fa`);
	parental_list.push(refName);
	set_genomeName.add(refName);
	parental_map[refName] = refFasta;

	for (let  parentalIdx = 1; parentalIdx < numParental; ++parentalIdx) {
		let pName = await inputText(`parental ${(parentalIdx + 1)} name`);
		if (!set_genomeName.has(pName)) {
			parental_list.push(pName);
			set_genomeName.add(pName);
			
			let filepath = await inputFile(`parental ${(parentalIdx + 1)} genome fasta file`, `${pName}.genome.fa`);
			parental_map[pName] = filepath;
		}
		else {
			console.warn(`exist ${pName}`);
		}
	}

	//let numProgeny = await inputNumber("number of progeny", { min: mode == "tetrad" ? 4 : (numParental > 1 ? 0 : 1), max: 16 });
	let numProgeny = mode == "tetrad" ? 4 : await inputNumber("number of subject", { min: 1, max: null });
	/** @type {string[]} */
	let progeny_list = [];
	/** @type {{[genomeName:string]:string}} */
	let progeny_map = {};

	for (let  progenyIdx = 0; progenyIdx < numProgeny; ++progenyIdx) {
		let pName = await inputText(`${mode == "tetrad" ? "progeny" : "subject"} ${(progenyIdx + 1)} name`);
		if (!set_genomeName.has(pName)) {
			progeny_list.push(pName);
			set_genomeName.add(pName);
			
			let filepath = await inputFile(`${mode == "tetrad" ? "progeny" : "subject"} ${(progenyIdx + 1)} fasta file`, `${pName}.genome.fa`);
			progeny_map[pName] = filepath;
		}
		else {
			console.warn(`exist ${pName}`);
		}
	}

	let genome_set = new GenomeDataSet();
	genome_set.name = taskName;
	genome_set.mode = mode;
	genome_set.ref = refName;
	genome_set.parental_list = parental_list;
	genome_set.progeny_list = progeny_list;
	genome_set.parental = parental_map;
	genome_set.progeny = progeny_map;

	console.log(genome_set);
	
	process.stdin.pause();//stop input

	return genome_set;
}

/**
 * @param {GenomeDataSet} genome_dataset
 * @returns {{[genomeName:string]:{genome_name:string,chr_rawName_list:string[],fasta:{[chrName:string]:string}}}}
 */
function make_genome_info(genome_dataset) {
	const header = ["Index", "Chromosome", "Length", "raw_chr_name", "symbol"].join("\t") + "\n";

	const genomeNameList = [...genome_dataset.parental_list, ...genome_dataset.progeny_list];
	const genomeFileMap = {
		...genome_dataset.parental,
		...genome_dataset.progeny,
	};

	/** @type {{[genomeName:string]:{genome_name:string,chr_rawName_list:string[],fasta:{[chrName:string]:string}}}} */
	const loaded_data = {};
	
	let found_dup_name = false;
	const seq_name_set = new Set();

	genomeNameList.forEach(genome_name => {
		const fa = (function () {
			try {
				return readFasta(genomeFileMap[genome_name]);
			}
			catch (ex) {
				console.log("error: fasta:", genome_name, genomeFileMap[genome_name]);
			}
		})();
		const chr_rawName_list = Object.keys(fa);
		const chr_outputName_list = chr_rawName_list.map(seq_name => GenomeInfo.transformChrName(genome_name, seq_name));

		let inChr_found_dup_name = false;
		chr_outputName_list.forEach(name => {
			// console.log(genome_name, name);
			
			if (seq_name_set.has(name)) {
				found_dup_name = true;
				inChr_found_dup_name = true;
				console.log("found duplicate sequence name:", name, "in", genome_name);
			}
			else {
				seq_name_set.add(name);
			}
		});
		if (inChr_found_dup_name) {
			console.error({
				chr_outputName_list,
				seq_name_set,
			});
		}

		loaded_data[genome_name] = {
			genome_name: genome_name,
			fasta: fa,
			chr_rawName_list: chr_rawName_list,
		};

		let output_fname = `${genome_dataset.output_path}/${genome_name}.length.txt`;
		let out_text = "";
		out_text += header;
		out_text += chr_rawName_list.map((seq_name, idx) => [
			idx + 1,                     // number of chr
			chr_outputName_list[idx],    // chr name
			fa[seq_name].length,         // chr length
			seq_name,                    // raw chr name
			seq_name,                    // default: raw chr name
		].join("\t")).join("\n");
		fs.writeFileSync(output_fname, out_text);
		console.log("output:", output_fname);
	});

	if (found_dup_name ) {
		throw new Error("found duplicate sequence name");
	}

	return loaded_data;
}

/**
 * @param {GenomeDataSet} genome_dataset
 * @param {{[genomeName:string]:{genome_name:string,chr_rawName_list:string[],fasta:{[chrName:string]:string}}}} loaded_data
 */
function explode_genome(genome_dataset, loaded_data) {
	if (!fs.existsSync(`${genome_dataset.tmp_path}/`)) {
		fs.mkdirSync(`${genome_dataset.tmp_path}/`);
	}
	if (!fs.existsSync(`${genome_dataset.tmp_path}/fasta`)) {
		fs.mkdirSync(`${genome_dataset.tmp_path}/fasta`);
	}

	Object.keys(loaded_data).forEach(gName => {
		loaded_data[gName].chr_rawName_list.forEach((chr_rawName, idx) => {
			const nChr = idx + 1;

			const fname_fasta = GenomeInfo.makeChrFilePath(genome_dataset, gName, chr_rawName);
			
			if (!has_argv_dataset_path && fs.existsSync(fname_fasta)) {
				console.log("found:", fname_fasta);
			}
			else {
				const output_chr_name = GenomeInfo.transformChrName(gName, chr_rawName);
				console.log("output:", gName, nChr, `'${chr_rawName}' -> '${output_chr_name}'`, fname_fasta);
				saveFasta(fname_fasta, {
					[output_chr_name]: loaded_data[gName].fasta[chr_rawName],
				});
			}
		});
	});
}

