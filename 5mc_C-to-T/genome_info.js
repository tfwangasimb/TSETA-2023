// @ts-check

/**
 * @interface T { [key: string]: string; }
 * @typedef StringMap
 * @property {string} [T]
 */

class ChrInfo {
	constructor() {
		this.name = "";
		this.length = 0;
	}
}

class GenomeInfo {
	/**
	 * @param {string} name - genome ID
	 * @param {string} file - fasta file path
	 * @param {ChrInfo[]} [chrInfo] - chrInfo
	 */
	constructor(name, file = "", chrInfo = null) {
		/**
		 * genome ID
		 * @type {string}
		 */
		this.name = name;
		
		/**
		 * fasta file path
		 * @type {string}
		 */
		this.file = file;
		
		/**
		 * chrInfo
		 * @type {ChrInfo[]}
		 */
		this.chrInfo = chrInfo;

		/**
		 * chr seq map
		 * @type {{ [chrName:string]:string }}}
		 */
		this.seq_map = null;
	}

	/**
	 * @param {string} name
	 * @param {string} file
	 */
	init(name, file) {
		this.name = name;
		this.file = file;
	}

	load() {
		const { readFasta } = require("../fasta_util.js");
		let fa = readFasta(this.file);
		this.seq_map = fa;
		this.chrInfo = Object.keys(fa).map(chrName => {
			let info = new ChrInfo();
			info.name = chrName;
			info.length = fa[chrName].length;
			return info;
		});
	}
}


module.exports.GenomeInfo = GenomeInfo;

