// @ts-check

const fs = require("fs");


/**
 * @template T
 */
class FinalTableRow {
	constructor() {
		/**************************************************
		 * Genomes
		 **************************************************/

		/** @type {number|T} Chromosome */
		this.Chromosome = 0;
		
		/** @type {number|T} Chromosome length */
		this.ref1_len = 0;

		/** @type {number|T} Chromosome length */
		this.ref2_len = 0;

		/**************************************************
		 * Dad, Mom, Individual nucleotides
		 **************************************************/
		
		/** @type {number|T} Q/C SNP (Individual nucleotides) */
		this.QC_SNV = 0;
		/** @type {number|T} Q/C InDel (Individual nucleotides) */
		this.QC_InDel = 0;
		/** @type {number|T} Q/C APD (Individual nucleotides) */
		this.QC_SNP = 0;

		/**************************************************
		 * TSETA: Individual recombination events (dad, mom, Fa, Fb, Fc and Fd)
		 **************************************************/

		/** @type {number|T} simple_CO */
		this.simple_CO = 0;
		
		/** @type {number|T} CO + CO(NCO) */
		this.CO_NCO = 0;

		/** @type {number|T} NCO */
		this.NCO = 0;
		/** @type {number|T} NCO 2 snp */
		this.NCO_2p = 0;

		/** @type {number|T} CO_div_nCO */
		this.CO_div_nCO = 0;

		// NCO = NCO1 + NCO2

		/** @type {number|T} NCO1 */
		this.NCO1 = 0;
		/** @type {number|T} NCO2 */
		this.NCO2 = 0;
		
		/** @type {number|T} NCO1 2 snp */
		this.NCO1_2p = 0;
		/** @type {number|T} NCO2 2 snp */
		this.NCO2_2p = 0;

		/**************************************************
		 * TSETA: Individual recombination events (dad, mom, Fa, Fb, Fc and Fd)
		 **************************************************/
		
		/** @type {number|T} SNP 2:2 */
		this.SNV_22 = 0;
		
		/** @type {number|T} NCO 3:1 */
		this.NCO_31 = 0;
		/** @type {number|T} NCO 4:0 */
		this.NCO_40 = 0;
		
		/** @type {number|T} SNP RIP (Q) */
		this.RIP_Q = 0;
		/** @type {number|T} SNP RIP (C) */
		this.RIP_C = 0;
		/** @type {number|T} SNP RIP (both Q and C) */
		this.RIP_QC = 0;
		
		/** @type {number|T} InDel RIP (Q) */
		this.RIP_2_Q = 0;
		/** @type {number|T} InDel RIP (C) */
		this.RIP_2_C = 0;
		
		/** @type {number|T} InDel 2:2, Strain-specific sequences 2:2 */
		this.sss22 = 0;
		//
		/** @type {number|T} InDel 1:3 */
		this.sss13 = 0;
		/** @type {number|T} InDel 0:4 */
		this.sss04 = 0;
		//
		/** @type {number|T} InDel 3:1, Strain-specific sequences 3:1 */
		this.sss31 = 0;
		/** @type {number|T} InDel 4:0, Strain-specific sequences 4:0 */
		this.sss40 = 0;
		
		/** @type {number|T} 1n:3 */
		this.del_1n3 = 0;
		/** @type {number|T} 2n:2 */
		this.del_2n2 = 0;
		/** @type {number|T} 3n:1 */
		this.del_3n1 = 0;
		/** @type {number|T} 4n:0 */
		this.del_4n0 = 0;

		/** @type {number|T} illegitimate mutation */
		this.IM1 = 0;
		// /** @type {number|T} illegitimate mutation (not 3:1) */
		// this.IM2 = 0;
		// /** @type {number|T} illegitimate mutation (not 4:0) */
		// this.IM3 = 0;
		// /** @type {number|T} illegitimate mutation deletion */
		// this.IM4 = 0;
		/** @type {number|T} illegitimate mutation InDel */
		this.IM2 = 0;
		/** @type {number|T} illegitimate mutation deletion */
		this.IM3 = 0;
	}

	/**
	 * @param {string} ref1_name
	 * @param {string} ref2_name
	 */
	__tableHeaderMap(ref1_name, ref2_name) {
		if (!ref1_name || !ref2_name) {
			throw new TypeError("ref1_name, ref2_name");
		}
		/** @type {FinalTableRow<string>} */
		const header_map = new FinalTableRow();// Make ordered map

		header_map.Chromosome = "Chromosome";
		header_map.ref1_len = `${ref1_name} (chromosome length)`;
		header_map.ref2_len = `${ref2_name} (chromosome length)`;
		
		header_map.QC_SNV = "SNP";
		header_map.QC_InDel = "InDel";
		header_map.QC_SNP = "APD";

		//

		header_map.simple_CO = "simple CO";// sCO
		header_map.CO_NCO = "CO(NCO)";// COn
		header_map.NCO = "NCO";// NCO1 + NCO2
		header_map.NCO_2p = "NCO(SNP≧2)";// NCO1 + NCO2 // ⩾
		header_map.CO_div_nCO = "CO/(CO+NCO(SNP≧2))"; // (sCO + COn) / (sCO + COn + NCO)
		
		header_map.NCO1 = "NCO1";
		header_map.NCO2 = "NCO2";
		header_map.NCO1_2p = "NCO1(SNP≧2)";
		header_map.NCO2_2p = "NCO2(SNP≧2)";
		
		header_map.SNV_22 = "SNP 2:2";
		header_map.NCO_31 = "SNP 1:3 or 3:1 (NCO)";
		header_map.NCO_40 = "SNP 4:0 (2NCO)";
		
		header_map.RIP_Q = `SNP RIP or C-to-T (${ref1_name})`;
		header_map.RIP_C = `SNP RIP or C-to-T (${ref2_name})`;
		header_map.RIP_QC = `SNP RIP or C-to-T (both ${ref1_name} and ${ref2_name})`;
		
		header_map.RIP_2_Q = `InDel RIP or C-to-T (${ref1_name})`;
		header_map.RIP_2_C = `InDel RIP or C-to-T (${ref2_name})`;
		
		header_map.sss22 = "InDel 2:2";// Strain-specific sequences 2:2 sss
		
		header_map.sss13 = "InDel 1:3";// Strain-specific sequences 1:3 sss
		header_map.sss04 = "InDel 0:4";// Strain-specific sequences 0:4 sss
		header_map.sss31 = "InDel 3:1";// Strain-specific sequences 3:1 sss
		header_map.sss40 = "InDel 4:0";// Strain-specific sequences 4:0 sss
		
		// header_map.del_1n3 = "1n:3";
		// header_map.del_2n2 = "2n:2";
		// header_map.del_3n1 = "3n:1";
		// header_map.del_4n0 = "4n:0";
		header_map.del_1n3 = "ID 1n:3";
		header_map.del_2n2 = "ID 2n:2";
		header_map.del_3n1 = "ID 3n:1";
		header_map.del_4n0 = "ID 4n:0";

		header_map.IM1 = "IM-1";
		header_map.IM2 = "IM-2";
		header_map.IM3 = "IM-3";
		// header_map.IM4 = "IM-4";
		
		Object.entries(header_map).forEach(([key, value]) => {
			if (!value) {
				// console.error();
				throw new Error(`final table header { "${key}": ${value} }`);
			}
		});

		return header_map;
	}

	/**
	 * @param {string} ref1_name
	 * @param {string} ref2_name
	 * @returns {{ [col_name:string]: string }}
	 */
	outputTableRow(ref1_name, ref2_name) {
		const header_map = this.__tableHeaderMap(ref1_name, ref2_name);
		// const col_name_list = this.colNameList();
		/** @type {{ [col_name:string]: string }} */
		const results = {};

		Object.keys(this).forEach(prop => {
			const name = header_map[prop];
			results[name] = this[prop];
		});

		// return col_name_list.map(k => results[k]).join("");

		return results;
	}
	/**
	 * @param {string} ref1_name
	 * @param {string} ref2_name
	 * @returns {string[]}
	 */
	outputTableHeader(ref1_name, ref2_name) {
		const header_map = this.__tableHeaderMap(ref1_name, ref2_name);
		return Object.keys(this).map(prop => {
			const name = header_map[prop];
			return name;
		});
	}
	
	/**
	 * @returns {string[]}
	 */
	_keys() {
		return Object.keys(this);
	}

	/**
	 * @param {string} ref1_name
	 * @param {string} ref2_name
	 * @returns {string[]}
	 */
	colNameList(ref1_name, ref2_name) {
		return Object.values(this.__tableHeaderMap(ref1_name, ref2_name));
	}
}

if (fs.realpathSync(process.argv[1]) == __filename) {
	const ttt = new FinalTableRow();
	const bbb = Object.keys(ttt.__tableHeaderMap());
	console.log({
		"FinalTableRow#_keys().sort()": ttt._keys().sort(),
		"Obkect.keys(FinalTableRow#__tableHeaderMap()).sort()": bbb.sort(),
	});
	if (ttt._keys().sort().toString() != bbb.sort().toString()) {
		throw new TypeError("FinalTableRow#__tableHeaderMap()");
	}
}

module.exports.FinalTableRow = FinalTableRow;

