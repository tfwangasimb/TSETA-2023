//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList, MyCoord } = require("./load_frag_list.js");
const { validation_chr } = require("./validation_seq.js");
const { join_chr_frag } = require("./join_chr_frag.js");

const argv = argv_parse(process.argv);

const VERBOSE = !!argv["--verbose"];

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);
if (argv_dataset_path.endsWith(`${dataset.name}.json`) == false) {
	throw new Error("dataset name no match file name");
}

const start_frag = argv["--start-frag"] ? String(argv["--start-frag"]) : null;
const end_frag = argv["--end-frag"] ? String(argv["--end-frag"]) : null;

let mafft_output_directory = String(argv["-i"] || `${dataset.tmp_path}/mafft_seq_frag`);

const genome_info_list = dataset.loadGenomeInfoList();

if (fs.realpathSync(process.argv[1]) == __filename) {
	main();
}
else {
	debugger;
}

function main() {
	let genome_has_error = false;
	
	merge_chr_all_fa(start_frag, end_frag);
	
	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		try {
			let chr_has_error = validation_chr(nChr, dataset.tmp_path, false);
			genome_has_error = genome_has_error || chr_has_error;
		}
		catch (ex) {
			console.error(ex);
		}
	}

	if (genome_has_error) {
		console.log("has error, no output");
		return;
	}
	else {
		//clone all results to output path
		for (let i = 1; i <= genome_info_list[0].chr_list.length; ++i) {
			let input_path = `${dataset.tmp_path}/mafft_ch${i}.fa`;
			let output_path = `${dataset.output_path}/mafft_ch${i}.fa`;
			fs.createReadStream(input_path).pipe(fs.createWriteStream(output_path));//copy file
			console.log("output:", output_path);
		}
	}
	
	console.log("next step:", "define and align the 50S rDNA loci");
	console.log("command:", `node ${__dirname}/re_align_rDNA.js -dataset ${argv_dataset_path} -rDNA rDNA.fa -chr 6`);

	console.log("or next step:", "analysis");
	console.log("snp calling command:", `node ${__dirname}/snp_summary.js -dataset ${argv_dataset_path}`);
	console.log("tetrad analysis command:", `node ${__dirname}/tetrad_summary.js -dataset ${argv_dataset_path} -min-co 5000`);
}

/**
 * @param {string} start_frag
 * @param {string} end_frag
 */
function merge_chr_all_fa(start_frag, end_frag) {
	const all_chr_frag_list = loadFragIdList(dataset);//load_frag_id_list();

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		const frag_list = all_chr_frag_list[nChr];
		if (frag_list) {
			const id_list = frag_list.map(a => a.id);
			
			// const first_id = id_list.indexOf(start_frag);
			const last_id = end_frag ? id_list.indexOf(end_frag) : -1;
			const filtered_list = last_id >= 0 ? frag_list.splice(0, last_id) : frag_list;

			console.log({
				"id_list[0]": id_list[0],
				start_frag, end_frag,
				"id_list.length": id_list.length,
				"filtered_list.length": filtered_list.length,
			});

			join_chr_frag(dataset, nChr, filtered_list, `${dataset.tmp_path}/mafft_ch${nChr}.fa`);
		}
		else {
			console.log("skip", "ch", nChr);
		}
	}
}

