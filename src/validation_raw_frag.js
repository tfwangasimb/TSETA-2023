//@ts-check

const fs = require("fs");
const Path = require("path");

const { argv_parse, array_groupBy } = require("./util.js");
const { tsv_parse, table_to_object_list } = require("./tsv_parser.js");
const { Dataset } = require("./dataset.js");
const { BlastnCoord, execAsync, exec_blastn, parseBlastnResults, blastn_coord, isCollide, groupByOverlap } = require("./blastn_util.js");
const { readFasta, saveFasta } = require("./fasta_util.js");
const { loadFragIdList } = require("./load_frag_list.js");

const argv = argv_parse(process.argv);

const argv_dataset_path = String(argv["-dataset"] || "");
const dataset = Dataset.loadFromFile(argv_dataset_path);

const genome_info_list = dataset.loadGenomeInfoList();

if (fs.realpathSync(process.argv[1]) == __filename) {
	main();
}

function main() {
	validation_frag({ raw: true, spread_centromere: true });//raw frag
	validation_frag({ raw: false, spread_centromere: true });//mafft frag
}

function validation_frag(flags = { raw: true, spread_centromere: true }) {
	const all_chr_frag_list = loadFragIdList(dataset);

	const input_path_dir = flags.raw ? `${dataset.tmp_path}/seq_frag` : `${dataset.tmp_path}/mafft_seq_frag`;
	
	/** @type {{ [chrName: string]: string }[]} */
	let srcRaw_genome_list = [];
	genome_info_list.forEach((genomeInfo, i) => {
		srcRaw_genome_list[i] = genomeInfo.loadFasta();
	});

	for (let nChr = 1; nChr <= genome_info_list[0].chr_list.length; ++nChr) {
		/** @type {{ [chrName: string]: string }} */
		const output_src_fa = {};

		const chr_name_list = [];
		/** @type {{ [chrName: string]: string }} */
		const src_raw_seq = {};
		srcRaw_genome_list.forEach((src_raw_genome, i) => {
			const idxChr = nChr - 1;
			const chr_name = genome_info_list[i].chr_list[idxChr].chr;
			src_raw_seq[chr_name] = src_raw_genome[chr_name];
			chr_name_list.push(chr_name);
			output_src_fa[chr_name] = "";//init
		});
		
		if (all_chr_frag_list[nChr]) {
			try {
				const final_list = [];
				const good_list = [];

				for (let coord of all_chr_frag_list[nChr]) {
					const fragId = coord.id;
					const in_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}.fa`;
			
					/** @type {{ [chrName: string]: string }} */
					let in_fa;

					if ((!coord.centromere || flags.spread_centromere) && fs.existsSync(in_filename)) {
						if (fs.statSync(in_filename).size) {
							in_fa = readFasta(in_filename);
							final_list.push(in_filename);
						}
						else {
							console.log("invalid fasta", "ch", nChr, "frag", fragId, Path.resolve(in_filename));
							throw new Error("invalid fasta");
						}
						//console.log("ch", nChr, "frag", i);
					}
					else {
						const in_ref1_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}_ref1.fa`;
						const in_ref2_filename = `${input_path_dir}/${flags.raw ? "" : "mafft_"}ch${nChr}_${fragId}_ref2.fa`;
			
						if (fs.existsSync(in_ref1_filename) && fs.existsSync(in_ref2_filename)) {
							if (!fs.statSync(in_ref1_filename).size) {
								console.log("invalid fasta", "ch", nChr, "frag", fragId, "ref1", Path.resolve(in_ref1_filename));
								throw new Error("invalid fasta ref1");
							}
							if (!fs.statSync(in_ref2_filename).size) {
								console.log("invalid fasta", "ch", nChr, "frag", fragId, "ref2", Path.resolve(in_ref2_filename));
								throw new Error("invalid fasta ref2");
							}
							const in_ref1_fa = readFasta(in_ref1_filename);
							const in_ref2_fa = readFasta(in_ref2_filename);
							
							final_list.push([in_ref1_filename, in_ref2_filename]);
							
							//console.log("ch", nChr, "frag", i, "ref1 ref2");

							in_fa = Object.assign({}, in_ref1_fa, in_ref2_fa);
						}
						else {
							console.log({
								in_filename, in_ref1_filename, in_ref2_filename,
							});
							throw new Error("ref1 ?? ref2 ??");
						}
					}

					const all_good = chr_name_list.map(chr_seq_name => {
						/** @type {string} */
						const _seq = in_fa[chr_seq_name];
						if (_seq) {
							const seq = flags.raw ? _seq : _seq.replace(/-/g, "");

							const raw_start = output_src_fa[chr_seq_name].length;
							const raw_end = raw_start + seq.length;
							const raw_frag = src_raw_seq[chr_seq_name].slice(raw_start, raw_end);
							for (let ii = 0; ii < seq.length; ++ii) {
								if (seq[ii] != raw_frag[ii]) {
									console.error("error at", in_filename, ii);
									console.error("src", raw_frag.slice(ii - 10, ii + 10));
									console.error("res", seq.slice(ii - 10, ii + 10));
									console.error(">>>", "----------^");
									console.error({
										raw_start,
										raw_end,
									});
									throw new Error("invalid frag");
								}
							}
							
							output_src_fa[chr_seq_name] += seq;
							
							return {
								name: chr_seq_name,
								fragId: fragId,
								start: raw_start,
								end: raw_end,
								file: in_filename,
							};
						}
					});
					good_list[fragId] = all_good;
					if (all_good.every(a => a)) {
						console.log("ch", nChr, "-", fragId, "good");
					}
					else {
						console.log("ch", nChr, "-", fragId, "bad");
					}
				}

				let va_result = chr_name_list.every(chr_seq_name => {
					for (let ii = 0; ii < src_raw_seq[chr_seq_name].length; ++ii) {
						if (output_src_fa[chr_seq_name][ii] != src_raw_seq[chr_seq_name][ii]) {
							console.error("error at", ii);
							console.error("src", src_raw_seq[chr_seq_name].slice(ii - 10, ii + 10));
							console.error("res", output_src_fa[chr_seq_name].slice(ii - 10, ii + 10));
							console.error(">>>", "----------^");
							throw new Error("invalid full length");
						}
					}
				});
				
				fs.writeFileSync(`${dataset.tmp_path}/coord_from_frag_ch${nChr}.json`, JSON.stringify(good_list, null, "\t"));

				fs.writeFileSync(`${dataset.tmp_path}/mafft_ch${nChr}.json`, JSON.stringify(final_list, null, "\t"));

				console.log("ch", nChr, "ok");
			}
			catch (ex) {
				console.error(ex);
			}
		}//if (all_chr_frag_list[nChr]) {
		else {
			console.log("skip", "ch", nChr);
		}
	}
}


module.exports.validation_frag = validation_frag;

