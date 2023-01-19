// @ts-check

const fs = require("fs");
const Path = require("path");

const { tsv_parse, _table_to_object_list, table_to_object_list } = require("./tsv_parser.js");
const { parse_GC_Content_table, GC_Content_Data } = require("./GC_content_util.js");
const { readFasta } = require("./fasta_util.js");
const { loadSetting } = require("./setting.js");

const VERBOSE = process.argv.indexOf("--verbose") >= 0;

class ChromosomeData {
	constructor() {
		/** @type {number} - number of chromosome */
		this.index = null;
		/** @type {string} - chromosome name */
		this.chr = null;
		/** @type {number} - chromosome length*/
		this.length = null;
		/** @type {string} - chromosome fasta file path */
		this.path = null;
		
		/** @type {string} - raw chromosome name */
		this.raw_chr_name = null;

		/** @type {string} - symbol name (mix different data) */
		this.symbol = null;
	}

	loadSeq() {
		return readFasta(this.path)[this.chr];
	}
}
class ChromosomeInfo {
	/**
	 * @param {ChromosomeData[]} list
	 * @param {{ [chrName:string]: ChromosomeData }} map
	 */
	constructor(list, map) {
		/** @type {ChromosomeData[]} */
		this.list = list;
		
		/** @type {{ [chrName:string]: ChromosomeData }} */
		this.map = map;
	}
}

class GenomeInfo {
	/**
	 * @param {Dataset} dataset
	 * @param {string} genome_name
	 * @param {"parental"|"progeny"} type
	 */
	constructor(dataset, genome_name, type) {
		let info = GenomeInfo.loadChrLength(dataset, genome_name);
		
		this.name = genome_name;
		this.chr_list = info.list;
		this.chr_map = info.map;

		/** @type {{ [sChr:string]: string }} */
		this.fasta = {};

		/** @type {"parental"|"progeny"} */
		this.type = type;
	}
	
	loadFasta() {
		return this.chr_list.reduce((obj, chrInfo) => Object.assign(obj, readFasta(chrInfo.path)), {});
	}
	
	/**
	 * @param {Dataset} dataset
	 * @param {string} genome_name
	 * @returns {ChromosomeInfo}
	 */
	static loadChrLength(dataset, genome_name) {
		let text_tab = fs.readFileSync(`${dataset.output_path}/${genome_name}.length.txt`).toString();
		let tab = tsv_parse(text_tab);
		let rows = table_to_object_list(tab, ["index", "chr", "length", "raw_chr_name", "symbol"], { start_row: 1 });
		
		const _chrList = rows.map(row => {
			let data = new ChromosomeData();
			Object.assign(data, {
				index: Number(row.index),
				chr: String(row.chr),//seq name in fasta
				length: Number(row.length),
				path: GenomeInfo.makeChrFilePath(dataset, genome_name, String(row.raw_chr_name)),
				raw_chr_name: row.raw_chr_name,
				symbol: row.symbol,
			});
			return data;
		});
		const chrList = [..._chrList].sort((a, b) => a.index - b.index);

		/** @type {{ [chrName:string]: ChromosomeData }} */
		const chrMap = {};
		chrList.forEach(chrInfo => {
			chrMap[chrInfo.chr] = chrInfo;
		});
		
		return {
			map: chrMap,
			list: chrList,
		};
	}

	/**
	 * @param {string} genome_name
	 * @param {string} raw_chr_name
	 * @returns {string}
	 */
	static transformChrName(genome_name, raw_chr_name) {
		return `${genome_name}_${encodeURIComponent(raw_chr_name)}`;
	}

	/**
	 * @param {GenomeDataSet} dataset
	 * @param {string} genome_name
	 * @param {string} raw_chr_name
	 * @returns {string}
	 */
	static makeChrFilePath(dataset, genome_name, raw_chr_name) {
		return `${dataset.tmp_path}/fasta/${GenomeInfo.transformChrName(genome_name, raw_chr_name)}.fa`;
	}
}

class RibosomalDNA_Data {
	constructor() {
		/** @type {number} */
		this.nChr = null;
		/** @type {string} - rDNA fasta file path */
		this.sequence = null;

		/**
		 * ref1 pos
		 * @type {[number, number]}
		 */
		this.region = [0, 0];

		/**
		 * @type {[number, number]}
		 */
		this.region_ref2 = [0, 0];
	}
}

class RibosomalDNA_Position_info {
	constructor() {
		/** @type {number} */
		this.chr = null;

		/**
		 * @type {number}
		 */
		this.region_start = null;
		
		/**
		 * @type {number}
		 */
		this.region_end = null;
	}
	
	/**
	 * @deprecated
	 * @type {number}
	 */
	get alignment_start() {
		return this.region_start;
	}
	
	/**
	 * @deprecated
	 * @type {number}
	 */
	get alignment_end() {
		return this.region_end;
	}
}

class MyFileInfo {
	constructor() {
		this.path = "";
		this.size = 0;
		this.atime = "";
		this.mtime = "";
		this.ctime = "";
		this.birthtime = "";
	}

	/**
	 * @param {fs.PathLike} file
	 */
	static getInfo(file) {
		const stat = fs.statSync(file);
		return {
			path: file,
			size: stat.size,
			atime: (stat.atime).toString(),
			mtime: (stat.mtime).toString(),
			ctime: (stat.ctime).toString(),
			birthtime: (stat.birthtime).toString(),
		};
	}
}

class MyDataInfo {
	constructor() {
		/** @type {{[genomeName:string]:MyFileInfo}} - parental fasta file apth */
		this.parental = {};
		
		/** @type {{[genomeName:string]:MyFileInfo}} - progeny fasta file path */
		this.progeny = {};

		/** @type {{[genomeName:string]:MyFileInfo}} - progeny fasta file path */
		this.raw_progeny = {};
	}
}

class GenomeDataSet {
	constructor() {
		this.name = "";

		/** @type {string|"tetrad"|"SNP"} */
		this.mode = "tetrad";

		this.ref = "";

		/** @type {string[]} */
		this.parental_list = [];

		/** @type {{[genomeName:string]:string}} - parental fasta file apth */
		this.parental = {};

		/** @type {string[]} */
		this.progeny_list = [];

		/** @type {{[genomeName:string]:string}} - progeny fasta file path */
		this.progeny = {};

		this.info = new MyDataInfo();
	}
	
	/** @type {string} */
	get output_path() {
		return (this.name);
	}
	
	/** @type {string} */
	get tmp_path() {
		return `${this.output_path}/tmp`;
	}
}

class MafftOptions {
	constructor() {
		this.algorithm = "localpair";
		this.default_algorithm = "";
		this.maxIterate = 1000;
		this.thread = 20;
	}
}

class TSETA_Options {
	constructor() {
		this.RIP = true;
	}
}

/**
 * @type {Map<string, Dataset>}
 */
const loaded_Dataset = new Map();

class Dataset extends GenomeDataSet {
	constructor() {
		super();

		/** @type {number} */
		this.version = null;

		//user input data

		/** @type {TSETA_Options} */
		this.options = new TSETA_Options();
		
		/** @type {MafftOptions} */
		this.mafft = null;

		/** @type {{ closeCOsMinDistance: number }} */
		this.crossover = null;

		/** @type {{[parentalName:string]:{[nChr:number]:GC_Content_Data[]}}} */
		this.gc_content = null;
		// Object.defineProperty(this, "gc_content", {
		// 	get: function () {
		// 		throw new TypeError("gc");
		// 	},
		// });

		/** @type {number} */
		this.GC_Content_window = null;

		/** @type {number} */
		this.GC_Content_m = null;

		/** @type {number[]} */
		this.GC_Content_value = null;

		/**
		 * AT-island min gc content %
		 * @type {number[]}
		 */
		this.min_gc_content = null;

		/** @type {string} file path */
		this.GC_Content_filePath = null;

		/** @type {{[nChr:number]:number[]}} */
		this.centromere = {};

		/** @type {{[nChr:number]:number[][]}} */
		this.telomere = {};

		/** @type {{[nChr:number]:number[][]}[]} */
		this.all_telomere = [];

		/**
		 * CO/NCO ignore rDNA.region
		 * @type {RibosomalDNA_Data}
		 */
		this.rDNA = new RibosomalDNA_Data();

		// auto generate

		/**
		 * display data
		 * @deprecated use dataset.rDNA.region and dataset.rDNA.region_ref2 instead
		 * @type {RibosomalDNA_Position_info}
		 */
		this.rDNA_info = null;

		// internal property

		/** @type {string[]} */
		this.genomeNameList = [];
	}

	/**
	 * @returns {"tetrad"|"SNP"}
	 */
	auto_detect_mode() {
		if (this.parental_list.length == 2) {
			if (this.progeny_list.length % 4 == 0) {
				//Tetrad analysis
				return "tetrad";
			}
		}
		return "SNP";
		// if (this.parental_list.length == 2) {
		// 	if (this.progeny_list.length == 4) {
		// 		return "SNP CO InDel";
		// 	}
		// 	else if (this.progeny_list.length > 0) {
		// 		return "2 parental +progeny SNP";
		// 	}
		// 	else {
		// 		return "2 parental SNP";
		// 	}
		// }
		// else if (this.parental_list.length == 1 && this.progeny_list.length > 0) {
		// 	return "1 parental +progeny SNP";
		// }
		// else {
		// 	throw new Error("");
		// }
	}

	/**
	 * @param {number} nChr
	 * @param {number} ref1_pos
	 */
	isIn_rDNA_ref1(nChr, ref1_pos) {
		if (!this.rDNA) {
			throw new Error(`this.rDNA: ${this.name}`);
		}
		const rDNA_nChr = Number(this.rDNA.nChr);
		const ref1_start = Number(this.rDNA.region[0]);
		const ref1_end = Number(this.rDNA.region[1]);
		if (nChr == rDNA_nChr) {
			if (ref1_pos >= ref1_start && ref1_pos <= ref1_end) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * @param {number} nChr
	 * @param {number} pos
	 */
	isIn_rDNA(nChr, pos) {
		if (!this.rDNA_info) {
			throw new Error(`this.rDNA_info: ${this.name}`);
		}
		const rDNA_nChr = Number(this.rDNA_info.chr);
		const rDNA_ma_start = Number(this.rDNA_info.alignment_start);
		const rDNA_ma_end = Number(this.rDNA_info.alignment_end);
		if (nChr == rDNA_nChr) {
			if (pos >= rDNA_ma_start && pos <= rDNA_ma_end) {
				return true;
			}
		}
		return false;
	}

	/**
	 * @param {string} ref_name
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	getGCByPos(ref_name, nChr, ref1_pos) {
		if (this.gc_content != null &&
			this.gc_content[ref_name] != null &&
			this.gc_content[ref_name][nChr] != null
		) {
			let row = this.gc_content[ref_name][nChr].find(a => ref1_pos >= a.start && ref1_pos <= a.end);
			this.gc_content[ref_name][nChr][ref1_pos / this.GC_Content_window];
			if (row) {
				return row.gc;
			}
			else {
				console.error({
					GC_Content_filePath: this.GC_Content_filePath,
					ref_name, nChr, ref1_pos,
					"s names": Object.keys(this.gc_content),
					"chr names": Object.keys(this.gc_content[ref_name]),
					"rows.length": Object.keys(this.gc_content[ref_name][nChr].length),
					"chr.rows.length": Object.keys(this.gc_content[ref_name]).map(chrName => this.gc_content[ref_name][chrName].length).join(","),
					"max": this.gc_content[ref_name][nChr].sort((a, b) => b.end - a.end)[0],
				});
				//throw new Error("getGCByPos(ref_name, nChr, ref1_pos) {");
			}
		}
	}
	
	/**
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	isInTelomere(nChr, ref1_pos) {
		if (this.telomere[nChr]) {
			let [[start1, end1], [start2, end2]] = this.telomere[nChr];
			return (ref1_pos >= start1 && ref1_pos <= end1) || (ref1_pos >= start2 && ref1_pos <= end2);
		}
	}	
	/**
	 * @param {number} nChr
	 * @param {number} ref1_pos pos in bp
	 */
	isInCentromere(nChr, ref1_pos) {
		if (this.centromere[nChr]) {
			let [start, end] = this.centromere[nChr];
			return ref1_pos >= start && ref1_pos <= end;
		}
	}

	loadGenomeInfoList() {
		return [
			...this.parental_list.map(gName => new GenomeInfo(this, gName, "parental")),
			...this.progeny_list.map(gName => new GenomeInfo(this, gName, "progeny")),
		];
	}
	
	loadGenomeInfoMap() {
		/** @type {{[GenomeName:string]: GenomeInfo}} */
		let map = {};
		this.parental_list.map(gName => map[gName] = new GenomeInfo(this, gName, "parental"));
		this.progeny_list.map(gName => map[gName] = new GenomeInfo(this, gName, "progeny"));
		return map;
	}

	/**
	 * @param {string|number} nChr ch1 -> 1
	 */
	loadFasta(nChr) {
		/** @type {string} */
		const file_name = this["results"]?.[Number(nChr) - 1] ?? `mafft_ch${nChr}.fa`;
		const path_to_fasta = Path.resolve(this.output_path, file_name);
		if (fs.existsSync(path_to_fasta)) {
			return readFasta(path_to_fasta);
		}
		else {
			throw new Error(`404 Dataset::loadFasta: ${path_to_fasta}`);
		}
	}
	
	/** @type {string} */
	get output_path() {
		return this._get_output_path();
	}
	_get_output_path() {
		if (this.$path) {
			const dirname = Path.dirname(this.$path);
			return Path.join(dirname, (this.name));
		}
		else {
			return (this.name);
		}
	}

	/**
	 * @param {fs.PathLike} path
	 */
	save(path) {
		// this.dataset_path
		fs.writeFileSync(path, JSON.stringify(this, null, "\t"));
	}

	/**
	 * @param {string} dataset_path
	 * @param {boolean} [reload]
	 * @returns {Dataset}
	 */
	static loadFromFile(dataset_path, reload = false) {
		try {
			if (!reload && loaded_Dataset.get(dataset_path) instanceof Dataset) {
				if (VERBOSE) {
					console.warn("from cache");
				}
				return loaded_Dataset.get(dataset_path);
			}
			else {
				if (!fs.existsSync(dataset_path)) {
					console.error({
						error: "No such file or directory",
						path: dataset_path,
						absolute: Path.resolve(dataset_path),
					});
					throw new Error("No such file or directory");
				}
				let obj = JSON.parse(fs.readFileSync(dataset_path).toString());
				let dataset = Dataset.__fromObject(obj);

				/** @type {Dataset} */
				// @ts-ignore
				let loaded = VERBOSE ? (new Proxy(dataset, {
					get: function (target, propertyKey, receiver) {
						let value = Reflect.get(target, propertyKey, receiver);
						if (value === null || value === undefined) {
							console.log(`dataset["${propertyKey.toString()}"] =>`, value);
							debugger;
						}
						return value;
					},
				})) : dataset;
				
				//loaded.gc_content = Dataset.__load_GC_content(dataset.GC_Content_filePath);

				//loaded.rDNA_info = Dataset.__load_rDNA_info_fromDataset(dataset_path);

				loaded.genomeNameList = [].concat(loaded.parental_list, loaded.progeny_list);

				loaded._apply_default_options();

				Object.defineProperty(loaded, "$path", {
					enumerable: false,
					writable: false,
					configurable: false,
					value: dataset_path,
				});

				if (VERBOSE) {
					// @ts-ignore
					if (loaded_Dataset.get(dataset_path) instanceof Dataset) {
						// @ts-ignore
						console.warn("dataset loaded:", loaded_Dataset.get(dataset_path).$path);
					}
				}
				loaded_Dataset.set(dataset_path, loaded);

				global.dataset = loaded;

				return loaded;
			}
		}
		catch (ex) {
			console.error(ex);
			console.log("error:", dataset_path);
			console.log("cwd:", process.cwd());
			console.log("json:", fs.readFileSync(dataset_path).toString().length);
			throw ex;
		}
	}

	_apply_default_options() {
		const global_default = loadSetting();
		if (global_default && global_default.version != null) {
			this.version = global_default.version;
		}

		if (this.options == null) { 
			this.options = new TSETA_Options();
		}
		else {
			Object.assign(this.options, new TSETA_Options(), this.options);
		}
	}
	
	load_GC_content() {
		// try {
			const GC_Content_filePath = Path.resolve(this.output_path, "../", this.GC_Content_filePath);
			if (fs.existsSync(GC_Content_filePath)) {
				const text = fs.readFileSync(GC_Content_filePath).toString();
				// @ts-ignore
				let table = table_to_object_list(tsv_parse(text), ["name", "chr", "start", "end", "gc"]);
				this.gc_content = parse_GC_Content_table(table);
			}
			else {
				console.error("load_GC_content", "Not found:", GC_Content_filePath);
			}
		// }
		// catch (ex) {
		// 	console.error(ex);
		// }
	}

	/**
	 * @deprecated
	 */
	load_rDNA_info() {
		console.warn("load_rDNA_info was deprecated; use dataset.rDNA.region and dataset.rDNA.region_ref2 instead");
		if (fs.existsSync(`${this.output_path}/rDNA_info.json`)) {
			this.rDNA_info = JSON.parse(fs.readFileSync(`${this.output_path}/rDNA_info.json`).toString());
		}
	}

	/**
	 * @param {any} obj
	 * @returns {Dataset}
	 */
	static __fromObject(obj) {
		let ds = new Dataset();
		Object.assign(ds, obj);
		return ds;
	}
}

module.exports.MafftOptions = MafftOptions;

module.exports.GenomeDataSet = GenomeDataSet;
module.exports.Dataset = Dataset;

module.exports.GenomeInfo = GenomeInfo;

module.exports.RibosomalDNA_Data = RibosomalDNA_Data;
module.exports.RibosomalDNA_Position_info = RibosomalDNA_Position_info;

module.exports.MyDataInfo = MyDataInfo;
module.exports.MyFileInfo = MyFileInfo;
