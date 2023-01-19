// @ts-check

/**
 * line:1793 add_2NCO_mut_ref: skip all 04-40-nco
 * line:1179 push_all_22_IM_marker: skip all 04-40-nco
 */

// if (typeof window == "object") {
// 	window.version_2 = true;
// }

class RIP_data {
	constructor() {
		/** @type {number} pos index */
		this.pos = 0;
		
		this.ref1_pos = 0;
		
		/** @type {number} 0: no rip, 1: ref1, 2: ref2 */
		this.mut_ref = 0;

		this.ref1 = "";
		this.ref2 = "";
		this.a = "";
		this.b = "";
		this.c = "";
		this.d = "";

		/**
		 * @type {-1|0|1} gc_in_out - -1: out, 1: in, 0: none
		 */
		this.gc_in_out = 0;
	}
}

let align_start_index;
let align_end_index;

let seg_snp = [];
// /** @type {RIP_data[]} */
// let rip_list = [];

///** @type {Uint8Array} */
//let parental_cmp_uint8array = null;
let region_rect = [];
/** @type {Uint32Array} ref1 pos map to multialign */
let ref1_pos_uint32array = null;
/** @type {Uint32Array} multialign pos map to ref1 */
let pos_ref1_uint32array = null;
/** @type {Uint32Array} ref2 pos map to multialign */
let ref2_pos_uint32array = null;
/** @type {Uint32Array} multialign pos map to ref2 */
let pos_ref2_uint32array = null;
// /** @type {Uint8Array} */
// let rip_map_uint8array = null;
/** @type {Uint32Array} ref1 ref2 score */
let ref1_ref2_score_uint32array = null;

class MarkerValue {
	constructor() {
		/** @type {number} - ma pos */
		this.pos = 0;
		/** @type {number} */
		this.value = 0;
		/** @type {number} */
		this.ref1_pos = 0;
		/** @type {string} */
		this.type = "";
		/** @type {string} */
		this.name = "";

		/** @type {number} */
		this.order = 0; // === MarkerData # order
		
		/** @type {string[]} */
		this.seg = [];

		/** @type {string} */
		this.color = undefined;

		/** @type {number} mutated referece */
		this.mut_ref = undefined;
	}
}

class MarkerData {
	/**
	 * @param {string} name - display name
	 * @param {string} property - property name
	 * @param {number} order - display order
	 */
	constructor(name, property, order) {
		/** @type {string} - display name */
		this.name = name;

		/** @type {string} - property name */
		this.property = property;

		/** @type {number} - display order */
		this.order = order;

		/** @type {MarkerValue[]} */
		this.values = [];
	}
}

class MarkerList {
	constructor() {
		/** @type {MarkerData[]} */
		this.list = [];

		/** @type {{[propertyName:string]:MarkerData}} */
		this.map = {};
	}

	/**
	 * @param {string} name - display name
	 * @param {string} property - property name
	 * @param {any} [_]
	 */
	add_marker_type(name, property, _) {
		const marker_data = new MarkerData(name, property, this.list.length)
		this.list.push(marker_data);
		this.map[property] = marker_data;
		return marker_data;
	}

	/**
	 * @param {string} markerName
	 * @param {MarkerValue} value
	 */
	add_marker_value(markerName, value) {
		this.map[markerName].values.push(value);
	}
}

/** @type {MarkerList} */
let allMarker = null;

const ColorID = get_colorID_Data();
function get_colorID_Data() {
	let colorID = {
	};
	colorID.mask = 0b1111;
	colorID.dad = 0;
	colorID.mom = 1;
	colorID.identical = 2;
	colorID.dad_rip = 3;
	colorID.mom_rip = 4;
	colorID.diff = 5;
	colorID.none = 8;
	colorID.indel_bit = 0b10000000;// 0b1000_0000;
	colorID.indel_mask = 0b10000000;// 0b1000_0000;
	return colorID;
}
class AnalysisOptions {
	constructor() {
		/** @type {number} */
		this.version = null;
		
		/** @type {"tetrad"|"SNP"} */
		this.mode = "tetrad";

		/** @type {number} nChr */
		this.nChr = 0;

		/** @type {rDNA_info} */
		this.rDNA_info = null;

		/** @type {crossover_list[0]} - this chr co list*/
		this.co_list = null;

		/** @type {non_crossover_list[0]} - this chr co list*/
		this.nco_list = null;

		this.fill_prev_color = true;

		/** @type {function(string):void} */
		this._onProgressReport = null;

		/** @type {boolean} */
		this.show_rDNA_snp = false;
		/** @type {boolean} */
		this.show_rDNA_indel_only = false;
		/** @type {boolean} */
		this.show_rDNA_non_InDel = true;
	}
}

/**
 * @param {Partial<AnalysisOptions>} options
 */
function initAnalyser(options) {
	allMarker = new MarkerList();

	if (options.mode == "tetrad") {
		if (typeof window == "object" && typeof window.version_2 == "boolean" && window.version_2) {
			allMarker.add_marker_type("SNV 2:2", "22", "MARKER_22");
			allMarker.add_marker_type("SNV 3:1/1:3", "31", "MARKER_31");
			allMarker.add_marker_type("SNV 4:0", "40","MARKER_40");

			allMarker.add_marker_type("SNV RIP (Q)", "rip_Q", "MARKER_RIP");
			allMarker.add_marker_type("SNV RIP (C)", "rip_C", "MARKER_RIP");
			allMarker.add_marker_type("SNV RIP (Q and C)", "rip_QC", "MARKER_RIP");
		}
		else {
			allMarker.add_marker_type("SNP 2:2", "22", "MARKER_22");
			allMarker.add_marker_type("SNP 3:1/1:3", "31", "MARKER_31");
			allMarker.add_marker_type("SNP 4:0", "40","MARKER_40");

			allMarker.add_marker_type("SNP RIP (Q)", "rip_Q", "MARKER_RIP");
			allMarker.add_marker_type("SNP RIP (C)", "rip_C", "MARKER_RIP");
			allMarker.add_marker_type("SNP RIP (Q and C)", "rip_QC", "MARKER_RIP");
		}
		
		allMarker.add_marker_type("InDel RIP (Q)", "rip_2_Q", "MARKER_SSS_RIP");
		allMarker.add_marker_type("InDel RIP (C)", "rip_2_C", "MARKER_SSS_RIP");

		allMarker.add_marker_type("InDel 2:2", "sss_22", "MARKER_sss_22");// Strain-specific sequences
		
		allMarker.add_marker_type("InDel 0:4", "sss_04", "MARKER_sss_04");// Strain-specific sequences
		allMarker.add_marker_type("InDel 1:3", "sss_13", "MARKER_sss_13");// Strain-specific sequences
		allMarker.add_marker_type("InDel 3:1", "sss_31", "MARKER_sss_31");// Strain-specific sequences
		allMarker.add_marker_type("InDel 4:0", "sss_40", "MARKER_sss_40");// Strain-specific sequences
		
		allMarker.add_marker_type("illegitimate deletions 1n:3", "1n3", "MARKER_1n3");
		allMarker.add_marker_type("illegitimate deletions 2n:2", "2n2", "MARKER_2n2");
		allMarker.add_marker_type("illegitimate deletions 3n:1", "3n1", "MARKER_3n1");
		allMarker.add_marker_type("illegitimate deletions 4n:0", "4n0", "MARKER_4n0");

		if (typeof window == "object" && typeof window.version_2 == "boolean" && window.version_2) {
			allMarker.add_marker_type("illegitimate mutation", "illegitimate_mutation", "MARKER_ill");
		}
		else {
			allMarker.add_marker_type("IM-1 markers", "illegitimate_mutation", "MARKER_ill");
		}
		// allMarker.add_marker_type("IM-2 markers", "illegitimate_mutation_31", "MARKER_ill_31");
		// allMarker.add_marker_type("IM-3 markers", "illegitimate_mutation_40", "MARKER_ill_40");
		// allMarker.add_marker_type("IM-4 markers", "illegitimate_mutation_deletion", "MARKER_ill_del");	allMarker.add_marker_type("IM-2 markers", "illegitimate_mutation_31", "MARKER_ill_31");
		allMarker.add_marker_type("IM-2 markers", "illegitimate_mutation_indel", "MARKER_ill_indel");
		allMarker.add_marker_type("IM-3 markers", "illegitimate_mutation_deletion", "MARKER_ill_del");

		const IM_SNP_22 = allMarker.add_marker_type(`IM SNP 2:2`, `IM_SNP_22`, "");
		IM_SNP_22.color = "hsla(32, 89%, 43%, 1)";
		const ml = [
			"A",
			"T",
			"C",
			"G",
			"-"
		];
		ml.forEach(from_r => {
			ml.forEach(to_p => {
				if (from_r != to_p) {
					if (from_r == "C" && to_p == "T") {
						// skip RIP
					}
					else if (from_r == "G" && to_p == "A") {
						// skip RIP
					}
					else {
						const md = allMarker.add_marker_type(`${from_r}🢂${to_p}`, `IM_${from_r}_${to_p}`, "");
						md.color = "hsla(32, 89%, 43%, 1)";
					}
				}
			});
		});
	}
	else {
		allMarker.add_marker_type("(dummy markers)", "dummy", "MARKER_dummy");
		allMarker.add_marker_type("indel", "snp");
		allMarker.add_marker_type("snv", "snv");
	}
}

/**
 * @param {seq_list} seq_list
 * @param {Partial<AnalysisOptions>} options
 */
function loadCmpData(seq_list, options) {
	const result_data = calc_seg_reg(seq_list, options);

	const results = {
		allMarker,

		// argv_fill_prev_color: options.fill_prev_color,
		// argv_nChr: options.nChr,

		align_start_index,
		align_end_index,

		// rip_list,
		seg_snp,

		region_rect,
		
		/** @type {Uint32Array} ref1 pos map to multialign */
		ref1_pos_uint32array: ref1_pos_uint32array,
		/** @type {Uint32Array} multialign pos map to ref1 */
		pos_ref1_uint32array: pos_ref1_uint32array,
		/** @type {Uint32Array} ref2 pos map to multialign */
		ref2_pos_uint32array: ref2_pos_uint32array,
		/** @type {Uint32Array} multialign pos map to ref2 */
		pos_ref2_uint32array: pos_ref2_uint32array,
		
		spore_cmp: {
			s22: allMarker.map["22"],
			s31: allMarker.map["31"],
			s40: allMarker.map["40"],
			
			rip_Q: allMarker.map["rip_Q"],
			rip_C: allMarker.map["rip_C"],
			rip_QC: allMarker.map["rip_QC"],
			rip_2_Q: allMarker.map["rip_2_Q"],
			rip_2_C: allMarker.map["rip_2_C"],
			
			s1n3: allMarker.map["1n3"],
			s2n2: allMarker.map["2n2"],
			s3n1: allMarker.map["3n1"],
			s4n0: allMarker.map["4n0"],
			
			// Strain-specific sequences
			sss_22: allMarker.map["sss_22"],
			sss_13: allMarker.map["sss_13"],
			sss_04: allMarker.map["sss_04"],
			sss_31: allMarker.map["sss_31"],
			sss_40: allMarker.map["sss_40"],
			
			illegitimate_mutation_list: allMarker.map["illegitimate_mutation"],
			// illegitimate_mutation_31: allMarker.map["illegitimate_mutation_31"],
			// illegitimate_mutation_40: allMarker.map["illegitimate_mutation_40"],
			// illegitimate_mutation_deletion: allMarker.map["illegitimate_mutation_deletion"],
			illegitimate_mutation_indel: allMarker.map["illegitimate_mutation_indel"],
			illegitimate_mutation_deletion: allMarker.map["illegitimate_mutation_deletion"],
		},
		
		ColorID: ColorID,

		ref_snp_list: result_data.ref_snp_list,
	};

	const nChr = options.nChr;
	
	// /** @type {number} */
	// const gff_sIdx = (options.gff_sIdx);

	// /** @type {string} */
	// const gff_ref = dataset.parental_list[options.gff_sIdx];

	// if (options.find_rip) {
	// 	results.rip_seq_tuple = find_rip_out_AT_island(gff_ref, gff_sIdx, nChr, seq_list, 30, 30);
	// }

	if (options.strain_specific/* && typeof require == "function"*/) {
		// const child_process = require("child_process");

		const g1 = dataset.parental[dataset.parental_list[0]];
	
		const g2 = dataset.parental[dataset.parental_list[1]];

		const in_sss = options.strain_specific.profiles.map(sss_len_key => {
			/** @type {{ start:number; end:number; }[]} */
			const qqq_sss = options.strain_specific[dataset.parental_list[0]][sss_len_key];
			const ccc_sss = options.strain_specific[dataset.parental_list[1]][sss_len_key];
			
			options.strain_specific[dataset.parental_list[0]][sss_len_key] = qqq_sss.map(m => {
				const start = m.index;
				const end = m.index + m[0].length;

				// const q_start = pos_ref1_uint32array[start];
				// const q_end = pos_ref1_uint32array[end];
				// const cmd_q_c = `blastn -query ${g1} -subject ${g2} -query_loc ${q_start}-${q_end} -perc_identity 65 -outfmt 6`;
				// const q_c = child_process.execSync(cmd_q_c).toString();
				// const count = q_c.split("\t").length;
				// if (count) {
				// 	console.log("QM6a strain specific", nChr, q_start, q_end);
				// }

				return {
					start,
					end,
					// ref_start: q_start,
					// ref_end: q_end,
					// specific: q_c ? count : "strain",
				};
			});
			
			options.strain_specific[dataset.parental_list[1]][sss_len_key] = ccc_sss.map(m => {
				const start = m.index;
				const end = m.index + m[0].length;

				// const q_start = pos_ref2_uint32array[start];
				// const q_end = pos_ref2_uint32array[end];
				// const cmd_c_q = `blastn -query ${g2} -subject ${g1} -query_loc ${q_start}-${q_end} -perc_identity 65 -outfmt 6`;
				// const c_q = child_process.execSync(cmd_c_q).toString();
				// const count = c_q.split("\t").length;
				// if (count) {
				// 	console.log("CBS1-1 strain specific", nChr, q_start, q_end);
				// }

				return {
					start,
					end,
					// ref_start: q_start,
					// ref_end: q_end,
					// specific: c_q ? count : "strain",
				};
			});
		});
	}

	if (options.rip_repeat_methyl_gff && options.gff) {
		console.log("options.rip_repeat_methyl_gff && options.gff");

		results.IM_gff = {};

		// ["🢂", "RIP"].forEach(marker_class => {
		dataset.parental_list.forEach((gff_ref, gff_sIdx) => {
			/** @type {{ inout_gene: string; struct_type: string; gene_ID: any; gc_value: any; AT_island: string; snv_indel: string; s:string; pos:number; peek:string; mutation_type:string; mut_ref:number; repeat: number; sss_50bp: "-"|"+"; sss_100bp: "-"|"+"; }[]} */
			const IM_gff = [];
			results.IM_gff[gff_ref] = IM_gff;
			
			const gene_list = options.gff[gff_ref];
			
			const chr_length = options.chrInfo_list[gff_sIdx];

			const ref_pos_map = [
				ref1_pos_uint32array,
				ref2_pos_uint32array,
			][gff_sIdx];
			
			const pos_ref_map = [
				pos_ref1_uint32array,
				pos_ref2_uint32array,
			][gff_sIdx];

			gene_list.forEach(gene => {
				// gene.$length = ref_pos_map[gene.end - 1] - gene.$start + 1;
				gene.$start = ref_pos_map[gene.start - 1];
				gene.$end = ref_pos_map[gene.end - 1];
				gene.$length = gene.$end - gene.$start + 1;
			});

			// if (0) {
			// 	const cds = gene_list.find(a => a.type == "CDS");
			// 	const exon = gene_list.find(a => a.type == "exon");
			// 	const mRNA = gene_list.find(a => a.type == "mRNA");
			// 	console.error({
			// 		cds,
			// 		exon,
			// 		mRNA
			// 	});
			// }

			const sort_matrix = {
				"CDS": {
					"CDS": 0,
					"exon": -1,
					"mRNA": -1,
					"gene": -1,
				},
				"exon": {
					"CDS": 1,
					"exon": 0,
					"mRNA": -1,
					"gene": -1,
				},
				"mRNA": {
					"CDS": 1,
					"exon": 1,
					"mRNA": 0,
					"gene": -1,
				},
				"gene": {
					"CDS": 1,
					"exon": 1,
					"mRNA": 1,
					"gene": 0,
				},
			};

			// const gene_type = options.rip_gff == true ? null : options.rip_gff;

			if (options.rip_gff) {
				results.rip_gff = find_rip_gff(gff_sIdx, nChr, seq_list, gene_list, results, options);
			}

			const peek_seq_range = 10;

			/** @type {MarkerValue[]} */
			let all_IM_marker = [];

			// concat markers and sort
			const IM_markers_list = allMarker.list.filter(a => {
				return a.name.indexOf("🢂") >= 0 || a.name.indexOf("RIP") >= 0;
			});
			
			IM_markers_list.push(allMarker.map[40]);
			IM_markers_list.push(allMarker.map["sss_04"]);
			IM_markers_list.push(allMarker.map["sss_40"]);

			IM_markers_list.forEach(markerGroup => {
				const values = markerGroup.values.filter(marker => {
					return gff_sIdx == (marker.mut_ref - 1);
				});
				all_IM_marker = all_IM_marker.concat(values);
			});
			all_IM_marker.sort((a, b) => a.pos - b.pos);

			all_IM_marker.forEach(marker => {
				const gc_wnd = findGCwndByMarker(gff_ref, nChr, pos_ref_map, ref_pos_map, marker);
				
				const found_gene = findAllGeneByMarker(gene_list, pos_ref_map, marker);
				
				// TODO: repeat_segment, methyl_ratio

				const found_repeat = options.repeat_segment ? find_repeat_segment_by_marker(options.repeat_segment[gff_ref][0], marker, ref_pos_map, pos_ref_map) : [];

				const methyl_left = options.methyl_left;
				const methyl_right = options.methyl_right;
				const o_methyl_ratio = {};
				if (options.methratio && options.methratio[gff_ref]) {
					// const mix_methyl_sample = true;
					// if (mix_methyl_sample) {
						const array_methyl_ratio = find_methyl_ratio_by_marker(options.methratio[gff_ref], marker, methyl_left, methyl_right, pos_ref_map, chr_length);
						for (let i = methyl_left; i > 0; --i) {
							const wt_d4_d8 = `"${array_methyl_ratio[methyl_left - i].join("\n")}"`;
							o_methyl_ratio[`methyl_ratio_l${i}`] = wt_d4_d8;
						}
						for (let i = 0; i <= methyl_right; ++i) {
							const wt_d4_d8 = `"${array_methyl_ratio[methyl_left + i].join("\n")}"`;
							o_methyl_ratio[`methyl_ratio_r${i}`] = wt_d4_d8;
						}
					// }
					// else {
					// 	const methyl_samples = ["", "4", "8"];
					// 	methyl_samples.forEach((sample, sampleIdx) => {
					// 		const methyl_ratio = options.methratio[gff_ref][sampleIdx];
							
					// 		const array_methyl_ratio = find_methyl_ratio_single_by_marker(methyl_ratio, marker, methyl_left, methyl_right, pos_ref_map, chr_length);
					// 		for (let i = methyl_left; i > 0; --i) {
					// 			const wt_d4_d8 = array_methyl_ratio[methyl_left - i];
					// 			o_methyl_ratio[`methyl_ratio_${sample}l${i}`] = wt_d4_d8;
					// 		}
					// 		for (let i = 0; i <= methyl_right; ++i) {
					// 			const wt_d4_d8 = array_methyl_ratio[methyl_left + i];
					// 			o_methyl_ratio[`methyl_ratio_${sample}r${i}`] = wt_d4_d8;
					// 		}
					// 	});
					// }
				}

				/** @type {GFF_ROW[]} */
				const sorted_found_gene = found_gene.sort((a, b) => sort_matrix[a.type][b.type]); // sort exon mRNA intron

				/** @type {GFF_ROW} */
				const gene = sorted_found_gene[0];

				const inout_gene = gene ? "intragenic" : "extragenic";
				const struct_type = gene ? (gene.type == "mRNA" ? "intron" : gene.type) : "";
				// if (gene && !gene.attributes) {
				// 	console.log(gene);
				// 	throw gene.attributes;
				// }
				/** gene ID T */
				const gene_ID = gene ? gene.attributes.ID : "";
				
				const tpm = options.RNA_TPM ? (function (gene) {
					if (gene) {
						const geneID = gene_ID.match(/^.*_\d+/)[0];
						return options.RNA_TPM[gff_ref].get(geneID);
					}
				})(gene) : "";
				
				// if (gene_ID && gene_ID.indexOf("TRC1_009613") >= 0) {
				// 	debugger;
				// }
				marker.color = gene ? options.gene_rip_color[gene.type] : null;
				// if (!marker.gene_ID) {
				// 	marker.gene_ID = "";
				// }
				// marker.gene_ID += [marker.color, gene_ID].join("_") + ",";
				
				const gc_value = gc_wnd.gc;
				const AT_island = gc_wnd.gc <= dataset.min_gc_content[gff_sIdx] ? "+" : "-";
				const snv_indel = marker.seg.includes("-") ? "InDel" : "SNV";

				const in_sss = (function () {
					if (options.strain_specific) {
						return options.strain_specific.profiles.map(sss_len_key => {
							/** @type {{ start:number; end:number; }[]} */
							const strain_specific = options.strain_specific[gff_ref][sss_len_key];
							return strain_specific.find(({ start, end }) => {
								return marker.pos >= start && marker.pos <= end;
							});
						});
					}
				})();

				const peek_name = [
					dataset.parental_list[0][0],
					dataset.parental_list[1][0],
					"A",
					"B",
					"C",
					"D",
				].map(a => `${a}: `);
				const peek = seq_list.map((ss, si) => peek_name[si] + ss.slice(marker.pos - peek_seq_range, marker.pos + peek_seq_range + 1)).join("\r\n") + "\r\n" + "   " + " ".repeat(peek_seq_range) + "^";

				// const deletetion_sign = "Del";
				
				IM_gff.push({
					inout_gene,
					gene_ID,
					struct_type,
					gc_value,
					AT_island,
					snv_indel,

					mut_ref: marker.mut_ref,

					// // `"` + marker.name.replace(/-🢂/g, "-->").replace(/🢂-/g, "->-").replace(/🢂/g, "->") + `"`,
					// mutation_type: `"` + marker.name.replace(/-🢂/g, `${deletetion_sign}->`).replace(/🢂-/g, `->${deletetion_sign}`).replace(/🢂/g, "->") + `"`,
					mutation_type: marker.name,

					// repeat: found_repeat.filter(a => a.identity >= 65).length,
					// repeat: found_repeat.filter(a => a.identity >= 70).length,
					// repeat_70: found_repeat.filter(a => a.identity >= 65 && a.identity < 70).length,
					// repeat_65: found_repeat.filter(a => a.identity < 65).length,
					repeat: found_repeat.filter(a => a.identity >= 65).length,
					repeat_60: found_repeat.filter(a => a.identity < 65 && a.identity >= 60).length,
					// sss_50bp: in_sss[0] ? "+" : "-",
					// sss_100bp: in_sss[1] ? "+" : "-",
					sss_60bp: in_sss ? (in_sss[0] ? "+" : "-") : NaN,
					// sss_50bp: in_sss[0] ? in_sss[0].specific : "strain",
					// sss_100bp: in_sss[1] ? in_sss[0].specific : "strain",

					...o_methyl_ratio,
					...tpm,

					pos: marker.pos,
					s: marker.seg.join(","),
					peek: peek,//
				});
				
				// const out = {
				// 	"GC": gc_value,
				// 	"AT": AT_island,
				// 	indel: is_indel,
				// 	type: gene.type,
				// 	gc_wnd: gc_wnd,
				// 	gene: gene,
				// };
			});
		});
	}

	if (
		typeof window == "object" && typeof window.version_2 == "boolean" && window.version_2 &&
		options.mode != "SNP"
	) {
		// Q RIP
		if (allMarker.list[3])
			allMarker.list[3].values.push(...allMarker.list[6].values.splice(0));

		// C RIP
		if (allMarker.list[4])
			allMarker.list[4].values.push(...allMarker.list[7].values.splice(0));

		// illegitimate mutation <= IM-1, IM-2, IM-3
		if (allMarker.list[17])
			allMarker.list[17].values.push(...allMarker.list[18].values.splice(0), ...allMarker.list[19].values.splice(0));
	}

	// console.log({
	// 	sss_22: allMarker.map["sss_22"].values.length,
	// 	sss_31: allMarker.map["sss_31"].values.length,
	// 	sss_40: allMarker.map["sss_40"].values.length,
	// });
	return results;
}

/**
 * @param {seq_list} seq_list
 * @param {Partial<AnalysisOptions>} options
 */
function calc_seg_reg(seq_list, options) {
	if (options.mode == null) {
		throw new TypeError("options.mode");
	}
	console.info({
		"show_rDNA_snp": options.show_rDNA_snp,
		"show_rDNA_indel_only": options.show_rDNA_indel_only,
		"show_rDNA_non_InDel": options.show_rDNA_non_InDel,
	});

	init_viewModel(seq_list);

	const result_data = {
	};

	let co_detail_list = options.co_list ? get_co_detail(options.co_list) : [];
	console.log({
		co_detail_list,
	});

	let nco_detail_list = options.nco_list ? options.nco_list/*.map((nco, i) => {
		nco.start = nco.snp_start_out;
		nco.end = nco.snp_end_out;
		nco.raw_idx = i;//co_idx
		return nco;
	})*/ : [];
	nco_detail_list.forEach(nco => {
		delete nco.rip_snv22_markers;
		delete nco.maybe_rip_snv40_markers;
	});

	/** @type {co_detail_list} */
	let work_co_detail_queue = co_detail_list ? [...co_detail_list] : null;
	/** @type {co_detail_list[0]} */
	let current_co_detail = co_detail_list ? work_co_detail_queue.shift() : null;
	// console.log(current_co_detail);
	
	/** @type {nco_detail_list} */
	let work_nco_detail_queue = nco_detail_list ? [...nco_detail_list] : null;
	/** @type {nco_detail_list[0]} */
	let current_nco_detail = nco_detail_list ? work_nco_detail_queue.shift() : null;

	seg_snp = new Array(seq_list[0].length);
	
	// check flanking => remove flanking => telomere
	for (let i = 0; i < seq_list[0].length; ++i) {
		if (seq_list.every(aa => aa[i] != "-")) {
			align_start_index = i;
			break;
		}
	}
	for (let i = seq_list[0].length - 1; i >= 0 ; --i) {
		if (seq_list.every(aa => aa[i] != "-")) {
			align_end_index = i;
			break;
		}
	}

	console.log({
		align_start_index,
		align_end_index,
	});
	if (align_start_index == -1 || align_end_index == -1 || align_start_index == align_end_index) {
		throw new Error("if (align_start_index == -1 || align_end_index == -1 || align_start_index == align_end_index) {");
	}

	// begin make_seg

	let ref1_pos = 1;
	let ref2_pos = 1;

	let reverse_fill_left_flank = options.reverse_fill_left_flank;

	//args_colors
	let prev_22_pos = 0;
	/** prev snp color */
	let prev_22_color = [
		0, 1,
		0, 0, 0, 0
	];
	/**
	 * fill RIP next color, use prev color
	 * @type {[number, number, number, number, number, number]}
	 */
	let all_prev_color = [
		0, 1,
		0, 0, 0, 0
	];
	let prev_has_rip = 0;

	/**
	 * @type {MarkerData[][][]}
	 */
	const snp_marker_mat = [
		[ [], [], [], [], [] ],//snv
		[ [], [], [], [], [] ],//snv del
		[ [], [], [], [], [] ],//snp
		[ [], [], [], [], [] ],//IM-1 (parental != InDel, sss 2:2, 1 <= x < 3)
		[ [], [], [], [], [] ],//IM-2 (SSS 2:2 x >= 3)
	];

	if (options.mode == "tetrad") {
		snp_marker_mat[0][0][4] = allMarker.map["40"];// 4:0 2NCO (strict)
		snp_marker_mat[0][1][3] = allMarker.map["31"];// 3:1  NCO (strict)
		snp_marker_mat[0][2][2] = allMarker.map["22"];// 2:2 SNP  (strict)
		snp_marker_mat[0][3][1] = allMarker.map["31"];// 3:1  NCO (strict)
		snp_marker_mat[0][4][0] = allMarker.map["40"];// 4:0 2NCO  (strict)
	}
	const snp_marker_order_map = snp_marker_mat[0].map(a => a.map(b => b.order));

	if (options.mode == "tetrad") {
		snp_marker_mat[1][1][3] = allMarker.map["1n3"];// 1n:3 (p diff, s 1 del)
		snp_marker_mat[1][2][2] = allMarker.map["2n2"];// 2n:2 (p diff, s 2 del)
		snp_marker_mat[1][3][1] = allMarker.map["3n1"];// 3n:1 (p diff, s 3 del)
		snp_marker_mat[1][4][0] = allMarker.map["4n0"];// 4n:0 (p diff, s 4 del)

		// sss
		snp_marker_mat[2][0][4] = allMarker.map["sss_04"];// 4:0    (p 1 del, s 4 del)
		//
		snp_marker_mat[2][1][3] = allMarker.map["sss_13"];// 1:3    (p 1 del, s 3 del)
		snp_marker_mat[2][2][2] = allMarker.map["sss_22"];// 2:2    (p 1 del, s 2 del)
		snp_marker_mat[2][3][1] = allMarker.map["sss_31"];// 3:1    (p 1 del, s 1 del)
		snp_marker_mat[2][4][0] = allMarker.map["sss_40"];// 4:0    (p 1 del, s 0 del)
		
		// IM-1 SNV x:y
		snp_marker_mat[3][1][2] = allMarker.map["illegitimate_mutation"];// IM-1 x = 3
		snp_marker_mat[3][0][2] = allMarker.map["illegitimate_mutation"];// IM-1 x = 4 IM-snp-22
		snp_marker_mat[3][2][0] = allMarker.map["illegitimate_mutation"];// IM-1 x = 4 IM-snp-22
		snp_marker_mat[3][2][1] = allMarker.map["illegitimate_mutation"];// IM-1 x = 3

		// illegitimate_mutation_indel
		
		// snp_marker_mat[3][0][4] = allMarker.map["illegitimate_mutation_40"];// 4:0      (p 1 del, s 4 del)
		// snp_marker_mat[3][1][3] = allMarker.map["illegitimate_mutation_31"];// 3:1(1:3) (p 1 del, s 3 del)
		// snp_marker_mat[3][3][1] = allMarker.map["illegitimate_mutation_31"];// 3:1
		// snp_marker_mat[3][4][0] = allMarker.map["illegitimate_mutation_40"];// 4:0
	}

	/** @type {MarkerValue} */
	let prev_marker = null;

	//compare seq
	let ref1_seq = seq_list[0];
	let ref2_seq = seq_list[1];
	for (/** @type {number} bp index */let pos = 0; pos < ref1_seq.length; ++pos) {
		// if (pos == 1142910) {
		// 	debugger
		// }
		// if (pos == 1142914) {
		// 	debugger
		// }
		// if (pos == 1645832) {
		// 	debugger;
		// }
		const count_of_indel = seq_list.filter(ss => ss[pos] == "-").length;
		const ref1 = ref1_seq[pos];
		const ref2 = dataset.parental_list.length == 2 ? ref2_seq[pos] : null;
		
		// AnalysisOptions
		if (options.mode == "tetrad") {
			push_tetrad_data(
				Math.max(1, ref1 != "-" ? ref1_pos : (ref1_pos - 1)),
				Math.max(1, ref2 != "-" ? ref2_pos : (ref2_pos - 1))
			);
		}
		else {// if (options.mode == "SNP") {
			push_SNP_data(
				Math.max(1, ref1 != "-" ? ref1_pos : (ref1_pos - 1)),
				Math.max(1, ref2 != "-" ? ref2_pos : (ref2_pos - 1))
			);
		}

		/**
		 * @param {number} ref1_pos
		 * @param {number} ref2_pos
		 */
		function push_SNP_data(ref1_pos, ref2_pos) {
			let row = {
				chr: options.nChr,
				pos: pos,
			};

			// if (window.debugger_break_pos && pos == window.debugger_break_pos) {
			// 	debugger;// window.debugger_break_pos = 6666923
			// }

			let has_snp = false;
			let has_indel = false;

			seq_list.slice(1).forEach((seq, idx) => {
				const sid = Number(idx) + 1;
				const va = seq[pos];
				
				if (ref2 && ref1 == ref2) {
					// // if (dataset.parental_list.length > 1) {
					// 	row[sid] = seg_snp[pos -1]?.[sid] ?? ColorID.diff;
					// // }
					// // else {
					// // 	if (va == ref1 || va == ref2) {
					// // 		row[sid] = ColorID.identical;
					// // 		// row[sid] = ColorID.dad;
					// // 	}
					// // 	else {
					// // 		row[sid] = ColorID.diff;
					// // 		has_snp = true;
					// // 	}
					// // }
					row[sid] = prev_22_color[sid];
				}
				else {
					if (va == ref1) {
						row[sid] = ColorID.dad;
					}
					else if (ref2 && va == ref2) {
						row[sid] = ColorID.mom;
						has_snp = true;
					}
					else {
						row[sid] = (ref2 || ref1 == "-") ? ColorID.diff : ColorID.mom;
						has_snp = true;
					}
					if (ref1 != "-") {
						prev_22_color[sid] = row[sid];
					}
				}

				let is_indel = va == "-";
				has_indel = is_indel || has_indel;

				if (is_indel) {
					if (seq[pos - 1] == "-" || seq[pos + 1] == "-") {
						row[sid] = ColorID.none;
					}
				}
				
				row[sid] = row[sid] | (is_indel ? ColorID.indel_bit : 0);
			});
			let ref1_is_indel = ref1 == "-";

			if ("normal") {
				row[0] = ColorID.dad | (ref1_is_indel ? ColorID.indel_bit : 0);
			}
			else if ("SNP rebuild genome") {
				row[0] = has_snp ? ColorID.dad : ColorID.none;
			}

			seg_snp[pos] = row;

			// reverse fill
			if (reverse_fill_left_flank) {
				if (has_snp) {
					console.warn("reverse fill left-flank");
					for (let prev = pos - 1; prev > prev_22_pos; --prev) {
						for (let idx = 1; idx < seq_list.length; ++idx) {
							const sid = Number(idx) + 1;
							if ((seg_snp[prev][sid] & ColorID.mask) == ColorID.identical) {
								if (row[sid] & ColorID.indel_bit) {
									seg_snp[prev][sid] = row[sid] | ColorID.indel_bit;
								}
								else {
									seg_snp[prev][sid] = (row[sid] & ColorID.mask) | (seg_snp[prev][sid] & ColorID.indel_bit);
								}
							}
							else {
								break;
							}
						}
					}
					prev_22_pos = pos;

					reverse_fill_left_flank = false;
				}
			}
			
			if (has_snp) {
				if (ref1_is_indel || has_indel) {
					let marker = {
						type: "SNP",
						name: allMarker.map["snp"].name,
						pos: pos,
						value: 0,
						ref1_pos: ref1_pos,
					};
					allMarker.map["snp"].values.push(marker);
				}
				else {
					let marker = {
						type: "SNV",
						name: allMarker.map["snv"].name,
						pos: pos,
						value: 1,
						ref1_pos: ref1_pos,
					};
					allMarker.map["snv"].values.push(marker);
				}
			}

			if (!ref1_pos_uint32array[ref1_pos]) {
				ref1_pos_uint32array[ref1_pos] = pos + 1;
			}
			pos_ref1_uint32array[pos] = ref1_pos;
			
			if (ref2) {
				if (!ref2_pos_uint32array[ref2_pos]) {
					ref2_pos_uint32array[ref2_pos] = pos + 1;
				}
				pos_ref2_uint32array[pos] = ref2_pos;
			}
		}

		/**
		 * @param {number} ref1_pos
		 * @param {number} ref2_pos
		 */
		function push_tetrad_data(ref1_pos, ref2_pos) {
			let pos_in_nco = false;
			if (co_detail_list && current_co_detail) {
				if (ref1_pos < current_co_detail.end) {
					if (ref1_pos > current_co_detail.start) {
					}
					else {
						// console.log({
						// 	what: "CO(NCO) inner ??",
						// 	ref1_pos,
						// 	current_co_detail,
						// });
						// debugger;
					}
				}
				else {
					current_co_detail = work_co_detail_queue.shift();
					// console.log(current_co_detail);
				}
			}
			if (nco_detail_list && current_nco_detail) {
				if (ref1_pos < current_nco_detail.snp_end_out) {
					if (ref1_pos > current_nco_detail.snp_start_out) {
					}
					else {
						// console.log({
						// 	what: "CO(NCO) inner ??",
						// 	ref1_pos,
						// 	current_nco_detail,
						// });
						// debugger;
					}
				}
				if (ref1_pos > current_nco_detail.snp_end_out) {
					current_nco_detail = work_nco_detail_queue.shift();
				}
			}
			if (current_nco_detail) {
				// if (ref1_pos > current_nco_detail.snp_start_out &&	// 20200917
				// 	ref1_pos < current_nco_detail.snp_end_out	// 20200917
				if (ref1_pos >= current_nco_detail.snp_start_out &&	// 20200917
					ref1_pos <= current_nco_detail.snp_end_out	// 20200917
				) {
					pos_in_nco = true;
				}
			}
			
			let a = seq_list[2][pos];
			let b = seq_list[3][pos];
			let c = seq_list[4][pos];
			let d = seq_list[5][pos];
			/** @type {[string, string, string, string]} */
			const spores = [a, b, c, d];
			const spores_indel = spores.filter(a => a == "-");
			const spores_has_indel = spores_indel.length > 0;
			
			const ref_has_indel = ref1 == "-" || ref2 == "-";
			const cmp_ref1_ref2 = ref1 == ref2;

			if (ref1 == "-" || ref2 == "-") {
				let prev_score = ref1_ref2_score_uint32array[pos - 1] | 0;
				ref1_ref2_score_uint32array[pos] = Math.max(0, prev_score - 1);
			}
			else {
				let prev_score = ref1_ref2_score_uint32array[pos - 1] | 0;
				ref1_ref2_score_uint32array[pos] = Math.max(0, prev_score + (cmp_ref1_ref2 ? 1 : -2));
			}
			
			let row = {
				chr: options.nChr,
				pos: pos,
			};

			let is_rDNA = (() => {
				if (options.show_rDNA_snp) {//ignore_rDNA
					return false;
				}
				else if (options.rDNA_info && options.nChr == options.rDNA_info.chr) {
					const start = options.rDNA_info.region_start || options.rDNA_info.alignment_start;
					const end = options.rDNA_info.region_end || options.rDNA_info.alignment_end;
					if (pos >= start && pos <= end) {
						return true;
					}
				}
			})();
			if (is_rDNA && current_co_detail && !current_co_detail.is_CO) {
				if (options.fill_rDNA) {
					row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
					row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
					row[2] = current_co_detail.state[0] | (a == "-" ? ColorID.indel_bit : 0);
					row[3] = current_co_detail.state[1] | (b == "-" ? ColorID.indel_bit : 0);
					row[4] = current_co_detail.state[2] | (c == "-" ? ColorID.indel_bit : 0);
					row[5] = current_co_detail.state[3] | (d == "-" ? ColorID.indel_bit : 0);
				}
				else if (options.show_rDNA_non_InDel) {
					push_not_rDNA();

					// if (prev_marker.order != allMarker.map["22"].order ||
					// 	prev_marker.order != allMarker.map["sss_22"].order
					// ) {
					// 	// row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
					// 	// row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
					// 	// row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
					// 	// row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
					// }

					if (ref1 == "-" || ref2 == "-") {
						if (spores_has_indel) {
							// row[0] = ref1 == "-" ? ColorID.diff : row[0];
							// row[1] = ref2 == "-" ? ColorID.diff : row[1];
							row[2] = a == "-" ? (ColorID.diff | (a == "-" ? ColorID.indel_bit : 0)) : row[2];
							row[3] = b == "-" ? (ColorID.diff | (b == "-" ? ColorID.indel_bit : 0)) : row[3];
							row[4] = c == "-" ? (ColorID.diff | (c == "-" ? ColorID.indel_bit : 0)) : row[4];
							row[5] = d == "-" ? (ColorID.diff | (d == "-" ? ColorID.indel_bit : 0)) : row[5];
						}
						else {
							// row[0] = ref1 == "-" ? ColorID.diff : row[0];
							// row[1] = ref2 == "-" ? ColorID.diff : row[1];
							row[2] = ColorID.diff;
							row[3] = ColorID.diff;
							row[4] = ColorID.diff;
							row[5] = ColorID.diff;
						}
					}

					// // override seg
					// if (ref_has_indel || spores_has_indel) {
					// 	row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
					// 	row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
					// 	row[2] = current_co_detail.state[0] | (a == "-" ? ColorID.indel_bit : 0);
					// 	row[3] = current_co_detail.state[1] | (b == "-" ? ColorID.indel_bit : 0);
					// 	row[4] = current_co_detail.state[2] | (c == "-" ? ColorID.indel_bit : 0);
					// 	row[5] = current_co_detail.state[3] | (d == "-" ? ColorID.indel_bit : 0);
					// }
				}
				// if (ref_has_indel) {
				// 	row[2] = prev_22_color[2] | (a == "-" ? ColorID.indel_bit : 0);
				// 	row[3] = prev_22_color[3] | (b == "-" ? ColorID.indel_bit : 0);
				// 	row[4] = prev_22_color[4] | (c == "-" ? ColorID.indel_bit : 0);
				// 	row[5] = prev_22_color[5] | (d == "-" ? ColorID.indel_bit : 0);
				// }
				// else {
				// 	row[2] = a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff) | (a == "-" ? ColorID.indel_bit : 0);
				// 	row[3] = b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff) | (b == "-" ? ColorID.indel_bit : 0);
				// 	row[4] = c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff) | (c == "-" ? ColorID.indel_bit : 0);
				// 	row[5] = d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff) | (d == "-" ? ColorID.indel_bit : 0);
				// 	if ((row[2] + row[3] + row[4] + row[5])) {
				// 		prev_22_color[0] = 0;
				// 		prev_22_color[1] = 1;
				// 		prev_22_color[2] = row[2];
				// 		prev_22_color[3] = row[3];
				// 		prev_22_color[4] = row[4];
				// 		prev_22_color[5] = row[5];
				// 		all_prev_color[0] = 0;
				// 		all_prev_color[1] = 1;
				// 		all_prev_color[2] = row[2];
				// 		all_prev_color[3] = row[3];
				// 		all_prev_color[4] = row[4];
				// 		all_prev_color[5] = row[5];
				// 	}
				// }

				if (options.show_rDNA_indel_only && spores_indel.length) {
					const m_info = snp_marker_mat[1][spores_indel.length][4 - spores_indel.length];//xn:y
					const col_idx = m_info.order;
					const marker = Object.assign(make_marker(spores_has_indel, col_idx), {
						type: "",
						name: m_info.name,
						is_rDNA: true,
					});
					allMarker.list[col_idx].values.push(marker);
				}
			}
			else {
				push_not_rDNA();
			}//if (!is_rDNA) {

			seg_snp[pos] = row;

			if (!ref1_pos_uint32array[ref1_pos]) {
				ref1_pos_uint32array[ref1_pos] = pos + 1;
			}
			if (!ref2_pos_uint32array[ref2_pos]) {
				ref2_pos_uint32array[ref2_pos] = pos + 1;
			}
			pos_ref1_uint32array[pos] = ref1_pos;
			pos_ref2_uint32array[pos] = ref2_pos;
			
			/**
			 * @param {boolean} is_indel
			 * @param {number} markerOrder
			 * @returns {MarkerValue}
			 */
			function make_marker(is_indel, markerOrder) {
				let label = new MarkerValue();
				label.pos = pos;
				label.ref1_pos = ref1_pos;
				label.value = (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (is_indel ? ColorID.indel_bit : 0);
				label.seg = seq_list.map(ss => ss[pos]);
				label.order = markerOrder;
				return label;
				// return {
				// 	pos: pos,
				// 	ref1_pos,
				// 	value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (is_indel ? ColorID.indel_bit : 0),
				// };
			}
			
			function push_not_rDNA() {
				// if (pos == 4881783) {
				// 	debugger
				// }
				if (ref1 == "-" && ref2 == "-") {
					push_parental_del();
				}
				else {
					push_parental_not_del();
					if (current_co_detail && !current_co_detail.is_CO) {
						push_all_22_IM_marker(current_co_detail);
					}
				}
			}

			/**
			 * @see {@link maybe_rip_snv40_markers}
			 * @param {current_co_detail} co_detail
			 */
			function push_all_22_IM_marker(co_detail) {
				// if (pos == 1547) {
				// 	debugger;
				// }

				/**
				 * @typedef co_val
				 * @type {0|1}
				 * 0: 像 QM6a
				 * 1: 像 CBS1-1
				 */

				/**
				 * Crossover = [Q, C, Fa, Fb, Fc, Fd]
				 * @type {[co_val, co_val, co_val, co_val, co_val, co_val]}
				 */
				const co_state = (co_detail.state);
				
				// if (co_state.every(a => a == co_state[0])) {// nco 4:0
				// 	return;
				// }
				
				/** 使用 crossover 找出像 QM6a 的子代 */
				const ref1_sp = spores.filter((v, si) => co_state[si] == 0);
				/** 使用 crossover 找出像 CBS1-1 的子代 */
				const ref2_sp = spores.filter((v, si) => co_state[si] == 1);

				/** 像 QM6a 的子代，序列像 QM6a ？ */
				const c1 = ref1_sp.filter(v => v == ref1);
				/** 像 CBS1-1 的子代，序列像 CBS1-1 ？ */
				const c2 = ref2_sp.filter(v => v == ref2);
				
				/** 2NCO 4:0 => QCQQQQ, length == 4 */
				const nco40_1 = spores.filter(v => v == ref1);
				/** 2NCO 4:0 => QCCCCC, length == 4 */
				const nco40_2 = spores.filter(v => v == ref2);

				let is_im_22;
				let col_idx;
				let im_22_type;
				let mut_ref = 0;

				if ((
						c1.length == 2 &&	/** 像 QM6a 的子代，序列都像 QM6a */
						c2.length == 0	/** 像 CBS1-1 的子代，序列都不像 CBS1-1 */
					) &&
					(nco40_1.length == 2 || nco40_2.length == 2) /** 不是 2NCO */
				) {
					const p = ref2;
					const s = ref2_sp[0] || spores.filter(a => a != p)[0];// ref2_sp[0] || in 2NCO
					/** 像 CBS1-1 的子代，序列都不像 CBS1-1，子代序列都一樣？ */
					if (ref2_sp.slice(1).every(a => a == s)) {
						if ((p == "C" && s == "T") ||	// is RIP C->T ?
							(p == "G" && s == "A")	// is RIP G->A ?
						) {
							// skip RIP
						}
						else {
							im_22_type = `IM_${p}_${s}`;
							if (!allMarker.map[im_22_type]) {
								console.log({
									im_22_type,
									p, s,
									spores,
									co_state,
									ref1, ref2,

									ref1_sp,
									ref2_sp,
									c1,
									c2,
									nco40_1,
									nco40_2,

									nChr: options.nChr,
								});
								debugger;
							}
							col_idx = allMarker.map[im_22_type].order;
							is_im_22 = true;
							mut_ref = 2;//1 + 1;// dataset.parental_list[1];
						}
					}
					else {
						// global.err_IM = (global.err_IM | 0) + 1;
						// if (global.err_IM == 1) {
						// 	process.on("beforeExit", () => {
						// 		process.exitCode = global.err_IM;
						// 	});
						// }
						console.error("ref2 not IM 2:2", pos, ref1, ref2, ...spores);
						// throw new Error("ref2 not IM 2:2");
					}
				}
				
				if ((
						c1.length == 0 &&	/** 像 QM6a 的子代，序列都不像 QM6a */
						c2.length == 2	/** 像 CBS1-1 的子代，序列都像 CBS1-1 */
					) &&
					(nco40_1.length == 2 || nco40_2.length == 2)
				) {
					const p = ref1;
					const s = ref1_sp[0] || spores.filter(a => a != p)[0];// ref2_sp[0] || in 2NCO
					/** 像 QM6a 的子代，序列都不像 QM6a，子代序列都一樣？ */
					if (ref1_sp.slice(1).every(a => a == s)) {
						if ((p == "C" && s == "T") ||	// is RIP C->T ?
							(p == "G" && s == "A")	// is RIP G->A ?
						) {
							// skip RIP
						}
						else {
							im_22_type = `IM_${p}_${s}`;
							col_idx = allMarker.map[im_22_type].order;
							is_im_22 = true;
							mut_ref = 1;//0 + 1;// dataset.parental_list[0];
						}
					}
					else {
						// global.err_IM = (global.err_IM | 0) + 1;
						// if (global.err_IM == 1) {
						// 	process.on("beforeExit", () => {
						// 		process.exitCode = global.err_IM;
						// 	});
						// }
						console.error("ref1 not IM 2:2", pos, ref1, ref2, ...spores);
						// throw new Error("ref1 not IM 2:2");
					}
				}

				if (is_im_22) {
					const marker = Object.assign(make_marker(spores_has_indel, col_idx), {
						type: im_22_type,
						name: allMarker.list[col_idx].name,
						current_co_detail: co_detail,
						mut_ref: mut_ref,
					});

					allMarker.list[col_idx].values.push(marker);

					allMarker.map["IM_SNP_22"].values.push(marker);

					{
						if (options.nChr == 3 && marker.pos == 1008409) {
							console.log("? IM ?", marker);
						}
						if (options.nChr == 3 && marker.pos == 1008411) {
							console.log("? IM ?", marker);
						}
						if (options.nChr == 6 && marker.pos == 4355624) {
							console.log("? IM ?", marker);
						}
						if (options.nChr == 6 && marker.pos == 4355625) {
							console.log("? IM ?", marker);
						}
					}
				}
			}

			function push_parental_del() {
				const m_info = allMarker.map["illegitimate_mutation_deletion"];
				const col_idx = m_info.order;
				const marker = Object.assign(make_marker(spores_has_indel, col_idx), {
					type: "illegitimate_mutation_deletion",
					name: m_info.name,
				});
				allMarker.list[col_idx].values.push(marker);
				
				// console.info("IM-3", String(pos).padStart(8, " "), ref1, ref2, "|", ...spores);
				
				// if (options.fill_mode == "snp filter") {
					if (pos == 2671919) {
						debugger
					}
					if (0) {
						row[0] = ColorID.none | (ref1 == "-" ? ColorID.indel_bit : 0);
						row[1] = ColorID.none | (ref2 == "-" ? ColorID.indel_bit : 0);
						row[2] = ColorID.diff | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = ColorID.diff | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = ColorID.diff | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = ColorID.diff | (d == "-" ? ColorID.indel_bit : 0);
					}
					else {
						// row[0] = ColorID.none | ColorID.indel_bit;
						// row[1] = ColorID.none | ColorID.indel_bit;
						row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
						row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
						row[2] = (a == "-" ? (ColorID.none | ColorID.indel_bit) : ColorID.diff);
						row[3] = (b == "-" ? (ColorID.none | ColorID.indel_bit) : ColorID.diff);
						row[4] = (c == "-" ? (ColorID.none | ColorID.indel_bit) : ColorID.diff);
						row[5] = (d == "-" ? (ColorID.none | ColorID.indel_bit) : ColorID.diff);
					}
				// }
				// else {
				// 	push_seg_not_SNP();
				// }
			}

			function push_parental_not_del() {
				// if (pos == 2102452) {
				// 	debugger
				// }
				// if (pos == 1990741) {
				// 	debugger
				// }
				// if (pos == 2019326) {
				// 	debugger;
				// }
				// @type {0|1|2} - 0: no RIP; 1: ref 1 RIP; 2: ref 2 RIP
				// if (pos == 1115768) {
				// 	debugger
				// }
				// if (pos == 932744) {//909077
				// 	debugger;
				// }
				let mut_ref = co_detail_list && isRIP(current_co_detail, pos, a, b, c, d, spores, ref1, ref2, cmp_ref1_ref2, row, all_prev_color, ColorID);
				// if (mut_ref < 0) {
				// 	console.log({
				// 		pos,
				// 		ref1, ref2,
				// 		a, b, c, d,
				// 	});
				// }
				/** @type {-2|-1|0|1|2} -1: out; 1: in; 0: none; 1: ref 1 RIP; 2: ref 2 RIP */
				let gc_in_out = 0;
				// if (pos == 909075) {
				// 	debugger;
				// }
				
				//if (mut_ref && (!current_nco_detail || !(ref1_pos >= current_nco_detail.start && ref1_pos <= current_nco_detail.end)/*!current_nco_detail.every(snp => snp.is_RIP)*/)) {
				if (mut_ref > 0 && pos_in_nco && current_nco_detail && !current_nco_detail.is_rip) {
					// if (count_of_snp_marker == 1) {
					// 	// RIP ?
					// 	if (count_of_indel == 3) {
					// 		gc_in_out = 0;
					// 	}
					// 	else {
					// 		gc_in_out = -1;
					// 	}
					// }
					// else if (count_of_snp_marker > 1) {
					// 	// GC ?
					// 	gc_in_out = 1;
					// }
					
					// TODO: RIP 2NCO 誰優先，2NCO 的 marker 不能全部都是 GAAAAA or CTTTTT

					if (20200622) {
						// if (snv40)
						// 	current_nco_detail.maybe_rip_snv40_markers.push(marker);
						// else {
						// 	current_nco_detail.rip_snv22_markers.push(marker);
						// }
						gc_in_out |= 0x000000F0;
					}
					else {
						const count_of_snp_marker = current_nco_detail.GCasso_marker;
						if (count_of_snp_marker == 1) {//one 4:0 markers -> CTTTTT, GAAAAA is RIP
							// gc_in_out |= 0x0000000F;// 20200619 RIP or NCO ??
							// if (current_nco_detail.snp_start_in == 2042652) {
							// 	debugger
							// }
							current_nco_detail.is_rip = true;
						}
						else if (count_of_snp_marker > 1) {//4:0 next NCO -> [CAAAAA, CTTTTT]
							gc_in_out |= 0x000000F0;
						}
					}

					// gc_in_out = -1;// display GA-AAAXX, CT-TTTTT: red marker
					// gc_in_out = 1;// display GA-AAAAA, CT-TTTTT: greenyellow marker

					// if (ref1_pos >= current_nco_detail.snp_start_out && ref1_pos <= current_nco_detail.snp_end_out) {
						// // // console.log("RIP in GC out", {
						// // // 	"RIP in GC out": "RIP in GC out",
						// // // 	pos: pos,
						// // // 	ref1_pos: ref1_pos,
						// // // 	ref2_pos: ref2_pos,
						// // // 	mut_ref: mut_ref,
						// // // 	current_nco_detail: current_nco_detail,
						// // // });// 20200603
						// gc_in_out *= mut_ref; // is ref1 or ref2
					// }
					// if (ref1_pos >= current_nco_detail.snp_start_out && ref1_pos <= current_nco_detail.snp_end_out) {
					// 	console.log("RIP in GC out", {
					// 		"RIP in GC out": "RIP in GC out",
					// 		pos: pos,
					// 		ref1_pos: ref1_pos,
					// 		ref2_pos: ref2_pos,
					// 		mut_ref: mut_ref,
					// 		current_nco_detail: current_nco_detail,
					// 	});
					// 	gc_in_out = -mut_ref;
					// }
					// else if (ref1_pos >= current_nco_detail.snp_start_in && ref1_pos <= current_nco_detail.snp_end_in) {
					// 	console.log("RIP in GC", {
					// 		"RIP in GC": "RIP in GC",
					// 		pos: pos,
					// 		ref1_pos: ref1_pos,
					// 		ref2_pos: ref2_pos,
					// 		mut_ref: mut_ref,
					// 		current_nco_detail: current_nco_detail,
					// 	});
					// 	gc_in_out = mut_ref;
					// }
					// else {
					// 	gc_in_out = 1;
					// }
				}// end check rip in NCO
				
				// TODO: remove RIP in NCO
				// if (gc_in_out ??? ) {
				if (mut_ref > 0/* && !gc_in_out*/) {
					// if (mut_ref == -1) {
					// 	gc_in_out |= 0xF0000000;
					// }
					// else if (mut_ref == -2) {
					// 	gc_in_out |= 0x0F000000;
					// }

					push_seg_RIP(gc_in_out);// CCT*T*
					
					const ref_name = ["", "Q", "C", "QC"][mut_ref];
					const is_sss = spores_indel.length == 2;
					const col_idx = allMarker.map[(is_sss ? "rip_2" : "rip") + "_" + ref_name].order;//rip_2_Q, rip_2_C
					const marker = Object.assign(make_marker(spores_has_indel, col_idx), {
						type: allMarker.list[col_idx].property,//is_sss ? "RIP_2" : "RIP",
						name: allMarker.list[col_idx].name,
						mut_ref: mut_ref,
					});

					if (current_nco_detail) {
						if (current_nco_detail.pos == 3901276) {
							console.error(marker.pos);
							console.error(marker);
							debugger;
						}
						if (gc_in_out) {
							marker.gc_in_out = gc_in_out;
							marker.nco = current_nco_detail;
							if (spores.every(s => s == spores[0])) {
								if (!current_nco_detail.maybe_rip_snv40_markers) {
									current_nco_detail.maybe_rip_snv40_markers = [];
								}
								current_nco_detail.maybe_rip_snv40_markers.push(marker);
							}
							else {
								if (!current_nco_detail.rip_snv22_markers) {
									current_nco_detail.rip_snv22_markers = [];
								}
								current_nco_detail.rip_snv22_markers.push(marker);
							}
						}
					}
					allMarker.list[col_idx].values.push(marker);
				}//if (mut_ref) {
				else {
					// if (pos == 934591) {
					// 	debugger;
					// }
					// if (mut_ref && ref1_pos >= current_nco_detail.start && ref1_pos <= current_nco_detail.end) {
					// 	console.log("RIP in GC", {
					// 		"RIP in GC": "RIP in GC",
					// 		pos: pos,
					// 		ref1_pos: ref1_pos,
					// 		ref2_pos: ref2_pos,
					// 		mut_ref: mut_ref,
					// 		current_nco_detail: current_nco_detail,
					// 	});
					// }

					// if (pos == 2624209) {
					// 	console.log(prev_22_pos);
					// 	debugger
					// }

					/**
					 * TODO:
					 * A GGGGGGG A TCTCTC
					 * C --GGGGG C ----TC
					 * A ---GGGG A --TCTC
					 * C ----GGG C ----TC
					 * C -----GG C ----TC
					 * A ---GGGG A --TCTC
					 *   ^^^^^^    ^^^^
					 *   L = 1     L = 2
					 */

					const VERSION_20200720_v3 = options.version == 0 ? 0 : 3;
					 // KmerAnalysis (K = 1)
					if (/*prev_marker && *//*options.fill_mode == "poly snp filter"*/VERSION_20200720_v3) {// change seg color
						const label = push_spore_marker();
						if (label) {
							let { col_idx, marker } = label;
							if (!allMarker.list[col_idx]) {
								console.log({col_idx});
							}
							allMarker.list[col_idx].values.push(marker);
						}

						const change_20200727 = true;
						let n_indel_threshold = 1;
						const current_non_indel = seq_list.map(ss => ss[pos]).find(a => a != "-");
						const current_non_indel_list = seq_list.map(ss => ss[pos]).filter(a => a != "-" && a != current_non_indel);
						const next_non_indel_list = seq_list.map(ss => ss[pos + 1]).filter(a => a != "-" && a != current_non_indel);
						const prev_non_indel_list = seq_list.map(ss => ss[pos - 1]).filter(a => a != "-" && a != current_non_indel);

						if (options.rDNA_info && options.nChr == (options.rDNA_info.nChr || options.rDNA_info.chr)) {
							if (pos >= options.rDNA_info.region_start &&
								pos <= options.rDNA_info.region_end
							) {
								const rep_data = (options.rDNA_info.data);
								const [iii] = [...rep_data].map((v, i) => [i, v]).sort(([ia, a], [ib, b]) => a.repeats.length - b.repeats.length)[0];//.map(([ia, a]) => ia + ":" + a.repeats.length);
								const cp = [...options.rDNA_info.data[iii].alignment_repeats.slice(-1)[0]];
								const se = cp.sort((a, b) => b - a)[0];
								if (pos >= se) {
									n_indel_threshold = 2;
								}
							}
						}

						// check
						//-1 A-AAAA
						//+0 AAAAAA
						
						const _prev_s = prev_marker ? prev_marker.seg.filter(s => s != "-") : null;
						const prev_s = _prev_s ? _prev_s[0] : null;

						if (change_20200727 && 
							seq_list.filter(ss => ss[pos] == "-").length >= n_indel_threshold &&
							current_non_indel_list.length == 0 &&
							(next_non_indel_list.length == 0 || prev_non_indel_list.length == 0) &&
							label && (
								// any indel
								(label.col_idx == allMarker.map["sss_04"].order) ||
								(label.col_idx == allMarker.map["sss_13"].order) ||
								(label.col_idx == allMarker.map["sss_22"].order) || // !important
								(label.col_idx == allMarker.map["sss_31"].order) ||
								(label.col_idx == allMarker.map["sss_40"].order) ||
								
								(label.col_idx == allMarker.map["1n3"].order) ||
								(label.col_idx == allMarker.map["2n2"].order) || // !important
								(label.col_idx == allMarker.map["3n1"].order) ||
								(label.col_idx == allMarker.map["4n0"].order) ||
								
								(label.col_idx == allMarker.map["illegitimate_mutation"].order) ||
								(label.col_idx == allMarker.map["illegitimate_mutation_indel"].order) ||
								(label.col_idx == allMarker.map["illegitimate_mutation_deletion"].order)
							) &&
							seg_snp[pos - 1]
						) {
							// row[0] = ColorID.diff | (ref1 == "-" ? ColorID.indel_bit : 0);
							// row[1] = ColorID.diff | (ref2 == "-" ? ColorID.indel_bit : 0);
							// const temp_1 = ColorID.diff | (a == "-" ? ColorID.indel_bit : 0);
							// const temp_2 = ColorID.diff | (b == "-" ? ColorID.indel_bit : 0);
							// const temp_3 = ColorID.diff | (c == "-" ? ColorID.indel_bit : 0);
							// const temp_4 = ColorID.diff | (d == "-" ? ColorID.indel_bit : 0);
							row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
							row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
							
							const tmp2 = seg_snp[pos - 1][2];//(a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
							const tmp3 = seg_snp[pos - 1][3];//(b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
							const tmp4 = seg_snp[pos - 1][4];//(c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
							const tmp5 = seg_snp[pos - 1][5];//(d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);

							row[2] = a == "-" ? (ColorID.diff | ColorID.indel_bit) : tmp2;//(current_co_detail && !current_co_detail.is_CO ? current_co_detail.state[0] : tmp2);
							row[3] = b == "-" ? (ColorID.diff | ColorID.indel_bit) : tmp3;//(current_co_detail && !current_co_detail.is_CO ? current_co_detail.state[1] : tmp3);
							row[4] = c == "-" ? (ColorID.diff | ColorID.indel_bit) : tmp4;//(current_co_detail && !current_co_detail.is_CO ? current_co_detail.state[2] : tmp4);
							row[5] = d == "-" ? (ColorID.diff | ColorID.indel_bit) : tmp5;//(current_co_detail && !current_co_detail.is_CO ? current_co_detail.state[3] : tmp5);
							row.poly = current_non_indel;
							label.marker.poly = current_non_indel;
							// if (typeof window == "object" && window != null) {
							// 	console.log("A", pos, "poly", current_non_indel);
							// }

							// hide poly mer
							// allMarker.list.forEach(li => li.values.forEach(m => m.poly ? void(m.hide = 1) : void(0)));
						}
						else if (!change_20200727 &&
							prev_marker &&
							prev_marker.pos == (pos - 1) &&
							prev_s && prev_marker.seg.every(s => s == "-" || s == prev_s) &&
							seq_list.map(sss => sss[pos]).every(s => s == "-" || s == prev_s)
						) {
							if (
								(prev_marker.order == allMarker.map["sss_04"].order) ||
								(prev_marker.order == allMarker.map["sss_13"].order) ||
								(prev_marker.order == allMarker.map["sss_22"].order) ||
								(prev_marker.order == allMarker.map["sss_31"].order) ||
								(prev_marker.order == allMarker.map["sss_40"].order) ||
								
								(prev_marker.order == allMarker.map["1n3"].order) ||
								(prev_marker.order == allMarker.map["2n2"].order) ||
								(prev_marker.order == allMarker.map["3n1"].order) ||
								(prev_marker.order == allMarker.map["4n0"].order) ||
								
								(prev_marker.order == allMarker.map["illegitimate_mutation"].order) ||
								(prev_marker.order == allMarker.map["illegitimate_mutation_indel"].order) ||
								(prev_marker.order == allMarker.map["illegitimate_mutation_deletion"].order)
							) {
								// row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
								// row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
								// // row[2] = ColorID.diff | (a == "-" ? ColorID.indel_bit : 0);
								// // row[3] = ColorID.diff | (b == "-" ? ColorID.indel_bit : 0);
								// // row[4] = ColorID.diff | (c == "-" ? ColorID.indel_bit : 0);
								// // row[5] = ColorID.diff | (d == "-" ? ColorID.indel_bit : 0);
								
								// // row[0] = ColorID.identical;
								// // row[1] = ColorID.identical;
								// row[2] = prev_22_color[2]; //current_colorset["identical"]
								// row[3] = prev_22_color[3];
								// row[4] = prev_22_color[4];
								// row[5] = prev_22_color[5];

								const prev_seg = seg_snp[pos - 1];
								prev_seg[0] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.dad) | (ref1 == "-" ? ColorID.indel_bit : 0);
								prev_seg[1] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.mom) | (ref2 == "-" ? ColorID.indel_bit : 0);
								// prev_seg[2] = prev_22_color[2];
								// prev_seg[3] = prev_22_color[3];
								// prev_seg[4] = prev_22_color[4];
								// prev_seg[5] = prev_22_color[5];

								// row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
								// row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
								// row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
								// row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);

								// this scope need "prev snv22 color", not prev_22_color

								let dirty = false;
								if (prev_seg[2] != prev_22_color[2]) {
									prev_seg[2] = ColorID.diff | (a == "-" ? ColorID.indel_bit : 0);
									dirty = true;
								}
								if (prev_seg[3] != prev_22_color[3]) {
									prev_seg[3] = ColorID.diff | (b == "-" ? ColorID.indel_bit : 0);
									dirty = true;
								}
								if (prev_seg[4] != prev_22_color[4]) {
									prev_seg[4] = ColorID.diff | (c == "-" ? ColorID.indel_bit : 0);
									dirty = true;
								}
								if (prev_seg[5] != prev_22_color[5]) {
									prev_seg[5] = ColorID.diff | (d == "-" ? ColorID.indel_bit : 0);
									dirty = true;
								}
								if (dirty) {
									prev_marker.gc_in_out = 0x0F000000;
								}
							}
							// else {
								if (!cmp_ref1_ref2) {
									push_seg_SNP();
								}
								else {
									push_seg_not_SNP();
								}
							// }
						}
						else {
							if (!cmp_ref1_ref2) {
								push_seg_SNP();
							}
							else {
								push_seg_not_SNP();
							}
						}

						if (label) {
							prev_marker = label.marker;
						}
						else {
							prev_marker = null;
						}
					}
					// else if (options.fill_mode == "snp filter") {
					// 	const label = push_spore_marker();
					// 	if (label) {
					// 		let { col_idx, marker } = label;
					// 		if (!allMarker.list[col_idx]) {
					// 			console.log({col_idx});
					// 		}
					// 		allMarker.list[col_idx].values.push(marker);
							
					// 		if ((label.col_idx == allMarker.map["sss_04"].order) ||
					// 			// (label.col_idx == allMarker.map["sss_13"].order) ||
					// 			// (label.col_idx == allMarker.map["sss_22"].order) ||
					// 			// (label.col_idx == allMarker.map["sss_31"].order) ||
					// 			(label.col_idx == allMarker.map["sss_40"].order) ||
								
					// 			// (label.col_idx == allMarker.map["1n3"].order) ||
					// 			// (label.col_idx == allMarker.map["2n2"].order) ||
					// 			// (label.col_idx == allMarker.map["3n1"].order) ||
					// 			// (label.col_idx == allMarker.map["4n0"].order) ||
								
					// 			(label.col_idx == allMarker.map["illegitimate_mutation"].order) ||
					// 			(label.col_idx == allMarker.map["illegitimate_mutation_indel"].order) ||
					// 			(label.col_idx == allMarker.map["illegitimate_mutation_deletion"].order)
					// 		) {
					// 			row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
					// 			row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
					// 			row[2] = ColorID.diff | (a == "-" ? ColorID.indel_bit : 0);
					// 			row[3] = ColorID.diff | (b == "-" ? ColorID.indel_bit : 0);
					// 			row[4] = ColorID.diff | (c == "-" ? ColorID.indel_bit : 0);
					// 			row[5] = ColorID.diff | (d == "-" ? ColorID.indel_bit : 0);
					// 		}
					// 		else {
					// 			if (!cmp_ref1_ref2) {
					// 				push_seg_SNP();
					// 			}
					// 			else {
					// 				push_seg_not_SNP();
					// 			}
					// 		}
					// 	}
					// 	else {
					// 		if (!cmp_ref1_ref2) {
					// 			push_seg_SNP();
					// 		}
					// 		else {
					// 			push_seg_not_SNP();
					// 		}
					// 	}
					// 	prev_marker = label.marker;
					// }
					else {
						if (!cmp_ref1_ref2) {
							push_seg_SNP();
						}
						else {
							push_seg_not_SNP();
						}
						
						const label = push_spore_marker();
						if (label) {
							let { col_idx, marker } = label;
							if (!allMarker.list[col_idx]) {
								console.log({col_idx});
							}
							allMarker.list[col_idx].values.push(marker);
						}

						prev_marker = label ? label.marker : null;
					}
				}
				prev_has_rip = mut_ref;

				if (typeof window == "object" && window != null) {
					window.marker_mod = "mod_2";
				}
				function push_spore_marker() {
					const mk = push_spore_marker_mod_2();
					if (mk) {
						if (current_co_detail && !current_co_detail.is_CO) {
							const mut_ref = add_2NCO_mut_ref(current_co_detail);
							mk.marker.mut_ref = mut_ref;//mut_ref: mut_ref
						}
					}
					return mk;
				}
				
				/**
				 * add 2NCO mut_ref
				 * @see {@link push_all_22_IM_marker}
				 * @param {current_co_detail} co_detail
				 * @returns {0|1|2} 0 if no ref | 1 if ref1 | 2 if ref2
				 */
				function add_2NCO_mut_ref(co_detail) {
					const co_state = co_detail.state;
					// if (co_state.every(a => a == co_state[0])) {// nco 4:0
					// 	return;
					// }
					
					const ref1_sp = spores.filter((v, si) => co_state[si] == 0);
					const ref2_sp = spores.filter((v, si) => co_state[si] == 1);
	
					const c1 = ref1_sp.filter(v => v == ref1);
					const c2 = ref2_sp.filter(v => v == ref2);
					
					const nco40_1 = spores.filter(v => v == ref1);
					const nco40_2 = spores.filter(v => v == ref2);
	
					// let is_22;
					// let col_idx;
					// let type;
					/** @type {0|1|2} */
					let mut_ref = 0;
	
					if ((c1.length == 2 && c2.length == 0) &&
						(nco40_1.length == 2 || nco40_2.length == 2)
					) {
						const p = ref2;
						const s = ref2_sp[0];
						if ((p == "C" && s == "T") ||
							(p == "G" && s == "A")
						) {
							// skip RIP
						}
						else {
							// type = `IM_${p}_${s}`;

							// if (!allMarker.map[type]) {
							// 	console.error({
							// 		type,
									
							// 		co_state,
							// 		ref1, ref2,
							// 		ref1_sp, ref2_sp,
							// 		c1, c2,
							// 		nco40_1, nco40_2
							// 	});
							// 	debugger;
							// }

							// col_idx = allMarker.map[type].order;
							// is_22 = true;
							mut_ref = 2;//1 + 1;// dataset.parental_list[1];
						}
					}
					
					if ((c1.length == 0 && c2.length == 2) &&
						(nco40_1.length == 2 || nco40_2.length == 2)
					) {
						const p = ref1;
						const s = ref1_sp[0];
						if ((p == "C" && s == "T") ||
							(p == "G" && s == "A")
						) {
							// skip RIP
						}
						else {
							// type = `IM_${p}_${s}`;
							// col_idx = allMarker.map[type].order;
							// is_22 = true;
							mut_ref = 1;//0 + 1;// dataset.parental_list[0];
						}
					}

					if (ref1 != "-") {
						mut_ref = 1;
					}
					else if (ref2 != "-") {
						mut_ref = 2;
					}

					return mut_ref;
				}
				
				function push_spore_marker_mod_2() {
					// const is_snp = ref1 != ref2;
					// const is_snv = !spores_has_indel && is_snp;// parental is indel ?
					const s_1s = spores.filter(a => a == ref1);
					const s_2s = spores.filter(a => a == ref2);
					const s_diff = spores.filter(a => a != ref1 && a != ref2);
					const s_n_indel = s_diff.filter(a => a == "-");
					const n_1s = s_1s.length;
					const n_2s = s_2s.length;
					const n_12s = n_1s + n_2s;

					/** @type {number} */
					let col_idx = snp_marker_order_map[n_1s][n_2s];

					let marker;
					// if (pos == 1132 || pos == 1133 || pos == 1134) {
					// 	debugger;
					// }

					// if (cmp_ref1_ref2) {
					// 	col_idx = allMarker.map["4n0"].order;
					// 	if (col_idx != null) {
					// 		marker = Object.assign(make_marker(spores_has_indel, col_idx), {
					// 			type: "IM-1-3",//AAAAAx // InDel 1:3/3:1 IM
					// 			name: allMarker.list[col_idx].name,
					// 		});
					// 	}
					// }
					
					// ill3:1, ill4:0
					if (col_idx != null && cmp_ref1_ref2 &&
						s_diff.length > 0 && s_diff.length != s_n_indel.length &&
						n_12s <= 4 // ??
					) {
						const m_info = snp_marker_mat[3][2][0];
						col_idx = m_info.order;
						marker = Object.assign(make_marker(spores_has_indel, col_idx), {
							type: "IM-1-3",//AAAAAx
							name: allMarker.list[col_idx].name,
						});
						if (typeof window == "object" && window != null) {
							console.log(marker);
						}
					}
					else if (col_idx != null && ref_has_indel) {
						// sss2:2, sss3:1, sss4:0 // InDel
						{
							const m_info = snp_marker_mat[2][4 - spores_indel.length][spores_indel.length];//InDel x:y; x:match, y:del
							col_idx = m_info.order;
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								_type: "Strain-specific sequences",//m_info.property,
								type: "InDel",//InDel x:y; x:match, y:del
								name: m_info.name,
							});
						}

						// (sss2:2, sss3:1, sss4:0), (ill3:1, ill4:0)
						// if (s_indel.length <= 2) {// if (s_indel.length == 0 || s_indel.length == 1 || s_indel.length == 2) {
						// 	// if (pos == 1990741)
						// 	// 	debugger
						// 	// Strain-specific sequences
						// 	// sss2:2, sss3:1, sss4:0
						// 	const m_info = snp_marker_mat[2][s_indel.length][4 - s_indel.length];
						// 	col_idx = m_info.order;
						// 	marker = Object.assign(make_marker(spores_has_indel, col_idx), {
						// 		type: "Strain-specific sequences",//m_info.property,
						// 		name: m_info.name,
						// 	});
						// }
						// else { // if (s_indel.length >= 3) {
						// 	// ill3:1, ill4:0
						//
						// 	const m_info = snp_marker_mat[3][s_indel.length][4 - s_indel.length];
						// 	col_idx = m_info.order;
						// 	marker = Object.assign(make_marker(spores_has_indel, col_idx), {
						// 		type: "illegitimate_mutation",
						// 		name: m_info.name,
						// 	});
						// }
					}
					else if (col_idx != null && !cmp_ref1_ref2) {// else if (ref1 != ref2) {
						// snp 2:2, 3:1, 4:0
						if (typeof window == "object" && window != null) {
							if (cmp_ref1_ref2) { // ??
								console.log({
									"ref cmp": cmp_ref1_ref2,
									pos: pos,
									ref1, ref2,
									spores: spores.join(""),
								});
							}
						}
						// 2:2, 3:1, 4:0
						marker = Object.assign(make_marker(spores_has_indel, col_idx), {
							type: `SNV${n_1s}:${n_2s}`,
							name: allMarker.list[col_idx].name,
							n_1s,
							n_2s,
							ref1,
							ref2,
							spores,
						});
					}
					else {
						// (illegitimate), (1n:2, 2n:2, 3n:1, 4n:1)
						
						if (s_diff.length == s_n_indel.length) {// diff if del
							// 1n:2, 2n:2, 3n:1, 4n:1
							let n_del = spores_indel.length;
							if (n_del == 1) {
								col_idx = allMarker.map["1n3"].order;
							}
							else if (n_del == 2) {
								col_idx = allMarker.map["2n2"].order;
							}
							else if (n_del == 3) {
								col_idx = allMarker.map["3n1"].order;
							}
							else if (n_del == 4) {
								col_idx = allMarker.map["4n0"].order;
							}
							if (col_idx != null) {
								marker = Object.assign(make_marker(spores_has_indel, col_idx), {
									type: "progney del",// progney
									name: allMarker.list[col_idx].name,
								});
							}
						}
						else {// diff if diff
							if (ref_has_indel) {
								const m_info = snp_marker_mat[3][n_1s][n_2s];
								if (m_info) {
									col_idx = m_info.order;
									marker = Object.assign(make_marker(spores_has_indel, col_idx), {
										type: "IM-1-2",// illegitimate mutation // InDel 2:2 IM
										name: allMarker.list[col_idx].name,
									});

									// console.info("IM-1-2", String(pos).padStart(8, " "), ref1, ref2, "|", ...spores);
								}
								else {
									col_idx = allMarker.map["illegitimate_mutation_indel"].order;
									marker = Object.assign(make_marker(spores_has_indel, col_idx), {
										type: "IM-2",// illegitimate mutation
										name: allMarker.list[col_idx].name,
									});

									// console.info("IM-2", String(pos).padStart(8, " "), ref1, ref2, "|", ...spores);
								}
							}
							else {
								col_idx = allMarker.map["illegitimate_mutation"].order;
								marker = Object.assign(make_marker(spores_has_indel, col_idx), {
									type: "IM-1-1",// illegitimate mutation // SNV IM
									name: allMarker.list[col_idx].name,
								});
								
								// console.info("IM-1-1", String(pos).padStart(8, " "), ref1, ref2, "|", ...spores);
							}
						}
					}

					if (col_idx >= 0 && marker) {
						return {
							col_idx, marker,
						};
					}
					return;

					// let marker;

					if (col_idx != null) {
						// 2:2, 3:1, 4:0
						marker = Object.assign(make_marker(spores_has_indel, col_idx), {
							type: "SNV",
							name: allMarker.list[col_idx].name,
						});
					}
					else {
						let a_n1s = spores.filter(a => a != ref1);
						let a_n2s = spores.filter(a => a != ref2);

						let mut_n1s = a_n1s.filter(a => a != ref2 && ref2 != "-" && a != "-");	//mutation
						let mut_n2s = a_n2s.filter(a => a != ref1 && ref1 != "-" && a != "-");	//mutation
						let mut_arr = [
							mut_n1s, mut_n2s,
						];

						if (mut_arr.some(aa => aa.length > 0)) {
							//is rip
							col_idx = allMarker.map["illegitimate_mutation"].order;

							// illegitimate mutation
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								type: "illegitimate_mutation",
								name: "illegitimate_mutation",
							});
						}
						else {
							// progeny indel
							// 1n3 2n2 3n1 4n0
							let n_del = spores_indel.length;
							if (n_del == 1) {
								col_idx = allMarker.map["1n3"].order;
							}
							else if (n_del == 2) {
								col_idx = allMarker.map["2n2"].order;
							}
							else if (n_del == 3) {
								col_idx = allMarker.map["3n1"].order;
							}
							else if (n_del == 4) {
								col_idx = allMarker.map["4n0"].order;
							}
							if (col_idx != null) {
								marker = Object.assign(make_marker(spores_has_indel, col_idx), {
									type: "SNP",// mod_0
									name: allMarker.list[col_idx]. name,
								});
							}
						}
					}

					if (col_idx >= 0 && marker) {
						return {
							col_idx, marker,
						};
					}
				}

				function push_spore_marker_mod_0() {
					let s_indel = spores.filter(a => a == "-");
					let spores_has_indel = s_indel.length > 0;
					let is_snp = ref1 != ref2;
					let is_snv = !spores_has_indel && is_snp;// parental is indel ?
					let s_1s = spores.filter(a => a == ref1);
					let s_2s = spores.filter(a => a == ref2);
					let n_1s = s_1s.length;
					let n_2s = s_2s.length;

					const snp_marker_order_map = [ [], [], [], [], [] ];

					snp_marker_order_map[0][4] = allMarker.map["40"].order;// 4:0 NCO
					snp_marker_order_map[1][3] = allMarker.map["31"].order;// 3:1 NCO
					snp_marker_order_map[2][2] = allMarker.map["22"].order;
					snp_marker_order_map[3][1] = allMarker.map["31"].order;// 3:1 NCO
					snp_marker_order_map[4][0] = allMarker.map["40"].order;// 4:0 NCO

					/** @type {number} */
					let col_idx = snp_marker_order_map[n_1s][n_2s];

					let marker;

					if (col_idx != null) {
						// 2:2, 3:1, 4:0
						marker = Object.assign(make_marker(spores_has_indel, col_idx), {
							type: "SNV",
							name: allMarker.list[col_idx].name,
						});
					}
					else {
						let a_n1s = spores.filter(a => a != ref1);
						let a_n2s = spores.filter(a => a != ref2);

						let mut_n1s = a_n1s.filter(a => a != ref2 && ref2 != "-" && a != "-");	//mutation
						let mut_n2s = a_n2s.filter(a => a != ref1 && ref1 != "-" && a != "-");	//mutation
						let mut_arr = [
							mut_n1s, mut_n2s,
						];

						if (mut_arr.some(aa => aa.length > 0)) {
							//is rip
							col_idx = allMarker.map["illegitimate_mutation"].order;

							// illegitimate mutation
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								type: "illegitimate_mutation",
								name: "illegitimate_mutation",
							});
						}
						else {
							// progeny indel
							// 1n3 2n2 3n1 4n0
							let n_del = s_indel.length;
							if (n_del == 1) {
								col_idx = allMarker.map["1n3"].order;
							}
							else if (n_del == 2) {
								col_idx = allMarker.map["2n2"].order;
							}
							else if (n_del == 3) {
								col_idx = allMarker.map["3n1"].order;
							}
							else if (n_del == 4) {
								col_idx = allMarker.map["4n0"].order;
							}
							if (col_idx != null) {
								marker = Object.assign(make_marker(spores_has_indel, col_idx), {
									type: "SNP",// mod_0
									name: allMarker.list[col_idx]. name,
								});
							}
						}
					}

					if (col_idx >= 0 && marker) {
						return {
							col_idx, marker,
						};
					}
				}

				function push_spore_marker_mod_1() {
					let s_indel = spores.filter(a => a == "-");
					let spores_has_indel = s_indel.length > 0;
					let is_snp = ref1 != ref2;
					let is_snv = !spores_has_indel && is_snp;// parental is indel ?
					let s_1s = spores.filter(a => a == ref1);
					let s_2s = spores.filter(a => a == ref2);
					let n_1s = s_1s.length;
					let n_2s = s_2s.length;

					const snp_marker_order_map = [ [], [], [], [], [] ];

					snp_marker_order_map[0][4] = allMarker.map["40"].order;// 4:0 NCO
					snp_marker_order_map[1][3] = allMarker.map["31"].order;// 3:1 NCO
					snp_marker_order_map[2][2] = allMarker.map["22"].order;
					snp_marker_order_map[3][1] = allMarker.map["31"].order;// 3:1 NCO
					snp_marker_order_map[4][0] = allMarker.map["40"].order;// 4:0 NCO

					/** @type {number} */
					let col_idx = snp_marker_order_map[n_1s][n_2s];

					let marker;

					if (col_idx != null) {
						if (is_snv) {
							// 2:2, 3:1, 4:0
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								type: "SNV",
								name: allMarker.list[col_idx].name,
							});
						}
						else {
							// 1n3 2n2 3n1 4n0
							if (col_idx == allMarker.map["31"].order) {
								if (s_indel.length == 1) {
									col_idx = allMarker.map["1n3"].order;
								}
								else if (s_indel.length == 3) {
									col_idx = allMarker.map["3n1"].order;
								}
							}
							else if (col_idx == allMarker.map["40"].order) {
								col_idx = allMarker.map["4n0"].order;
							}
							else if (col_idx == allMarker.map["22"].order) {
								col_idx = allMarker.map["2n2"].order;
							}
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								type: "SNP",
								name: allMarker.list[col_idx].name,
							});
						}
					}
					else {
						let a_n1s = spores.filter(a => a != ref1);
						let a_n2s = spores.filter(a => a != ref2);

						let mut_n1s = a_n1s.filter(a => a != ref2 && ref2 != "-" && a != "-");	//mutation
						let mut_n2s = a_n2s.filter(a => a != ref1 && ref1 != "-" && a != "-");	//mutation
						let mut_arr = [
							mut_n1s, mut_n2s,
						];

						// let n1s_ins = a_n1s.filter(a => a != ref2 && ref2 == "-" && a != "-");	//progeny ins != r1, r2 del
						// let n2s_ins = a_n2s.filter(a => a != ref1 && ref1 == "-" && a != "-");	//progeny ins != r2, r1 del
						// let p_ins_arr = [
						// 	n1s_ins, n2s_ins,
						// ];

						if (mut_arr.some(aa => aa.length > 0)) {
							//is rip
							col_idx = allMarker.map["illegitimate_mutation"].order;

							// illegitimate mutation
							marker = Object.assign(make_marker(spores_has_indel, col_idx), {
								type: "illegitimate_mutation",
								name: "illegitimate_mutation",
							});
						}
						// else if (p_ins_arr.some(aa => aa.length > 0)) {
						// 	//progeny ins 1n3 2n2 3n1 
						// 	col_idx = markerTypes.map["PROGENY_INS"];
						// 	// progeny ins
						// 	marker = Object.assign(make_marker(is_indel, col_idx), {
						// 		type: "progeny ins",
						// 		name: "progeny ins",
						// 	});
						// }
						else {
							// progeny indel
							// 1n3 2n2 3n1 4n0
							let n_del = s_indel.length;
							if (n_del == 1) {
								col_idx = allMarker.map["1n3"].order;
							}
							else if (n_del == 2) {
								col_idx = allMarker.map["2n2"].order;
							}
							else if (n_del == 3) {
								col_idx = allMarker.map["3n1"].order;
							}
							else if (n_del == 4) {
								col_idx = allMarker.map["4n0"].order;
							}
							if (col_idx != null) {
								marker = Object.assign(make_marker(spores_has_indel, col_idx), {
									type: "progeny indel",
									name: allMarker.list[col_idx].name,
								});
							}
						}
					}

					if (col_idx >= 0 && marker) {
						return {
							col_idx, marker,
						};
					}

					// if ([0, 1, 2, 3, 4, 5].some(i => row[i] == ColorID.diff)) {
					// 	illegitimate_mutation.push({
					// 		pos: pos,//index
					// 		ref1_pos,
					// 	});
					// }
				}

				/**
				 * @param {-2|-1|0|1|2} gc_in_out - -1: out; 1: in; 0: none; 1: ref 1 RIP; 2: ref 2 RIP
				 */
				function push_seg_RIP(gc_in_out) {
					let tmp = [];
					tmp[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
					tmp[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
					tmp[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
					tmp[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
					for (let i = 2; i <= 5; ++i) {
						if (tmp[i] == prev_22_color[i]) {
							row[i] = prev_22_color[i];
						}
					}
					row.is_rip = mut_ref;
					// rip_list.push({
					// 	pos,//index
					// 	ref1_pos,
					// 	mut_ref: mut_ref,
					// 	ref1,
					// 	ref2,
					// 	a,
					// 	b,
					// 	c,
					// 	d,
					// 	gc_in_out,
					// });
				}

				function push_seg_SNP() {
					if (2020062410) {// 20200624_10
						row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
						row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
						row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
						prev_22_color[0] = row[0];
						prev_22_color[1] = row[1];
						prev_22_color[2] = row[2];
						prev_22_color[3] = row[3];
						prev_22_color[4] = row[4];
						prev_22_color[5] = row[5];
						all_prev_color[0] = row[0];
						all_prev_color[1] = row[1];
						all_prev_color[2] = row[2];
						all_prev_color[3] = row[3];
						all_prev_color[4] = row[4];
						all_prev_color[5] = row[5];
					}
					else {
						row[0] = ColorID.dad | (ref1 == "-" ? ColorID.indel_bit : 0);
						row[1] = ColorID.mom | (ref2 == "-" ? ColorID.indel_bit : 0);
						// if (options.fill_prev_color) {
						// 	row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
						// 	row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
						// 	row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
						// 	row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);

						// 	prev_22_color[0] = row[0];
						// 	prev_22_color[1] = row[1];
							
						// 	all_prev_color[0] = row[0];
						// 	all_prev_color[1] = row[1];

						// 	if (pos >= 1537153 && pos <= 1537741) {
						// 		debugger;
						// 	}

						const q2 = seq_list.slice(2).filter(sss => ref1 == sss[pos]).length;
						const c2 = seq_list.slice(2).filter(sss => ref2 == sss[pos]).length;
						const is_snp_22 = q2 == 2 && c2 == 2;

						// 	for (let iii = 2; iii <= 5; ++iii) {
						// 		if (row[iii] == (ColorID.diff & ColorID.indel_bit)) {
						// 			row[iii] = prev_22_color[iii];
						// 		}
						// 		else {
						// 			prev_22_color[iii] = row[iii];
						// 			all_prev_color[iii] = row[iii];
						// 		}
						// 	}
						// }
						// else {
							row[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff)) | (a == "-" ? ColorID.indel_bit : 0);
							row[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff)) | (b == "-" ? ColorID.indel_bit : 0);
							row[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff)) | (c == "-" ? ColorID.indel_bit : 0);
							row[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff)) | (d == "-" ? ColorID.indel_bit : 0);
							if (is_snp_22) {
								prev_22_pos = pos;
								prev_22_color[0] = row[0];
								prev_22_color[1] = row[1];
								prev_22_color[2] = row[2];
								prev_22_color[3] = row[3];
								prev_22_color[4] = row[4];
								prev_22_color[5] = row[5];
							}
							all_prev_color[0] = row[0];
							all_prev_color[1] = row[1];
							all_prev_color[2] = row[2];
							all_prev_color[3] = row[3];
							all_prev_color[4] = row[4];
							all_prev_color[5] = row[5];
						//}
						//seg.push(row);
					}

					//2:2, 3:1, 4:0
					//push_SNP_marker();//snp 1n:3 2n:2 3n:1 4n:0 marker

					function push_SNP_marker() {
						if ((row[2] & ColorID.mask) <= ColorID.mom &&
							(row[3] & ColorID.mask) <= ColorID.mom &&
							(row[4] & ColorID.mask) <= ColorID.mom &&
							(row[5] & ColorID.mask) <= ColorID.mom) {
							//let va = (row[2] & ColorID.mask) + (row[3] & ColorID.mask) + (row[4] & ColorID.mask) + (row[5] & ColorID.mask);
							let indel = ref1 == "-" || ref2 == "-" || a == "-" || b == "-" || c == "-" || d == "-";
							//if (indel || va != 2) {
							//const cmp_column_idx = [4, 2, 0, 2, 4];//2:2,3:1,4:0
							let sq = spores.filter(v => v == ref1);
							let sc = spores.filter(v => v == ref2);
							let nsq = sq.length;
							let nsc = sc.length;
							let col_idx = -1;
							if (ref1 != "-" && ref2 != "-") {
								if (nsq == 2 && nsc == 2) {
									//2:2
									col_idx = 0;
								}
								else if (nsq == 1 || nsc == 1) {//1:3 / 3:1
									//3:1
									col_idx = 1;
								}
								else if (nsq == 4 || nsc == 4) {//0:4 / 4:0
									//4:0
									col_idx = 2;
								}
							}
							else if ((nsq == 1 && ref1 == "-") || (nsc == 1 && ref2 == "-")) {// 1 indel
								//1n:3
								col_idx = 3;
							}
							else if ((nsq == 2 || nsc == 2) && indel) {
								//2n:2
								col_idx = 4;
							}
							else if ((nsq == 3 && ref1 == "-") || nsc == 3 && ref2 == "-") {
								//3n:1
								col_idx = 5;
							}
							else if ((nsq == 4 && ref1 == "-") || nsc == 4 && ref2 == "-") {
								//4n:0
								col_idx = 6;
							}
							// spore_cmp_array[cmp_column_idx[va] + (indel ? 1 : 0)].push({
							// 	pos: pos,
							// 	value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (indel ? ColorID.indel_bit : 0),
							// });
							if (col_idx >= 0 /* && pos >= align_start_index && pos <= align_end_index*/) {
								allMarker.list[col_idx].values.push(Object.assign(new MarkerValue(), {
									pos: pos,
									ref1_pos,
									value: (row[5] << 3 | row[4] << 2 | row[3] << 1 | row[2] << 0) | (indel ? ColorID.indel_bit : 0),
									type: "nco",
									name: "" + col_idx,
								}));
							}
							//}
						}//is marker
					}//push_marker
				}//push_SNP

				// reverse fill left-flank
				if (reverse_fill_left_flank) {
					if (!cmp_ref1_ref2) {
						console.warn("reverse fill left-flank");
						for (let prev = pos - 1; prev > prev_22_pos; --prev) {
							for (let idx = 1; idx < seq_list.length; ++idx) {
								const sid = Number(idx) + 1;
								if ((seg_snp[prev][sid] & ColorID.mask) == ColorID.identical) {
									if (seg_snp[prev][sid] & ColorID.indel_bit) {
										seg_snp[prev][sid] = row[sid] | ColorID.indel_bit;
									}
									else {
										seg_snp[prev][sid] = (row[sid] & ColorID.mask) | (seg_snp[prev][sid] & ColorID.indel_bit);
									}
								}
								else {
									break;
								}
							}
						}
						prev_22_pos = pos;

						reverse_fill_left_flank = false;
					}
				}
			}//notRDNA() push not del
			function push_seg_not_SNP() {
				if (options.fill_prev_color) {
					/*if ((
							(ref1 == "-" && ref2 != "-") ||
							(ref1 != "-" && ref2 == "-")
						) &&
						(
							spores.every(s => s == ref1) ||
							spores.every(s => s == ref2)
						)
					) {
						if ()
					}
					else */if (ref1 == "-" && ref2 == "-") {
						row[0] = ColorID.none;
						row[1] = ColorID.none;
						row[2] = a != "-" ? ColorID.diff : ColorID.none;
						row[3] = b != "-" ? ColorID.diff : ColorID.none;
						row[4] = c != "-" ? ColorID.diff : ColorID.none;
						row[5] = d != "-" ? ColorID.diff : ColorID.none;
						all_prev_color[0] = ColorID.none;
						all_prev_color[1] = ColorID.none;
						all_prev_color[2] = ColorID.none;
						all_prev_color[3] = ColorID.none;
						all_prev_color[4] = ColorID.none;
						all_prev_color[5] = ColorID.none;
					}
					else {
						/*let prev_is_indel_40 = false;
						// if (0) {
						// 	const prev_ss = seq_list.map(sss => sss[pos - 1]);
						// 	const prev_spores = prev_ss.slice(2);
						// 	// if prev_ss == `A-AAAA` or prev_ss == `-AAAAA`
						// 	if ((
						// 			(prev_ss[0] == "-" && prev_ss[1] != "-") ||
						// 			(prev_ss[0] != "-" && prev_ss[1] == "-")
						// 		) && (
						// 			prev_spores.every(s => s == ref1) ||
						// 			prev_spores.every(s => s == ref2)
						// 		)
						// 	) {
						// 		prev_is_indel_40 = true;
						// 	}
						// }
						// if (pos == 2328595) {
						// 	debugger;
						// }
						
						if (prev_is_indel_40) {
							if (change_202000730) {
							}
							else {
								row[2] = prev_has_rip ? prev_22_color[2] : ColorID.diff;
								row[3] = prev_has_rip ? prev_22_color[3] : ColorID.diff;
								row[4] = prev_has_rip ? prev_22_color[4] : ColorID.diff;
								row[5] = prev_has_rip ? prev_22_color[5] : ColorID.diff;
							}
						}
						else */{
							row[0] = ColorID.identical;
							row[1] = ColorID.identical;

							// if (pos == 29259) {
							// 	debugger;
							// }
							
							const change_202000730 = true;
							if (change_202000730) {
								// 20210125
								row[2] = (a == ref1 ? ColorID.identical : ColorID.diff) | (a == "-" ? ColorID.indel_bit : 0); //current_colorset["identical"]
								row[3] = (b == ref1 ? ColorID.identical : ColorID.diff) | (b == "-" ? ColorID.indel_bit : 0);
								row[4] = (c == ref1 ? ColorID.identical : ColorID.diff) | (c == "-" ? ColorID.indel_bit : 0);
								row[5] = (d == ref1 ? ColorID.identical : ColorID.diff) | (d == "-" ? ColorID.indel_bit : 0);
							}
							else {
								row[2] = prev_has_rip ? prev_22_color[2] : (a == ref1 ? ColorID.identical : ColorID.diff); //current_colorset["identical"]
								row[3] = prev_has_rip ? prev_22_color[3] : (b == ref1 ? ColorID.identical : ColorID.diff);
								row[4] = prev_has_rip ? prev_22_color[4] : (c == ref1 ? ColorID.identical : ColorID.diff);
								row[5] = prev_has_rip ? prev_22_color[5] : (d == ref1 ? ColorID.identical : ColorID.diff);
							}
							all_prev_color[0] = ColorID.identical;
							all_prev_color[1] = ColorID.identical;
							all_prev_color[2] = ColorID.identical;
							all_prev_color[3] = ColorID.identical;
							all_prev_color[4] = ColorID.identical;
							all_prev_color[5] = ColorID.identical;
						}
					}
				}
				else {
					if (ref1 == "-" && ref2 == "-") {
						row[0] = ColorID.none;
						row[1] = ColorID.none;
						row[2] = a != "-" ? ColorID.diff : ColorID.none;
						row[3] = b != "-" ? ColorID.diff : ColorID.none;
						row[4] = c != "-" ? ColorID.diff : ColorID.none;
						row[5] = d != "-" ? ColorID.diff : ColorID.none;
					}
					else {
						row[0] = ColorID.identical;
						row[1] = ColorID.identical;
						row[2] = ColorID.identical; //current_colorset["identical"]
						row[3] = ColorID.identical;
						row[4] = ColorID.identical;
						row[5] = ColorID.identical;
					}
				}
			}// push_seg_not_SNP()
		}

		if (ref1 != "-") {
			++ref1_pos;
		}
		if (ref2 != "-") {
			++ref2_pos;
		}
	}// for pos in seq_list[0]

	// TODO check maybe_rip_snv40_markers
	// maybe_rip_markers
	// throw new Error(maybe_rip_snv40_markers);

/*
2	1 rip* snv40	2NCO	
3	1 rip* snv40	2NCO	
4	1 rip* snv40	2NCO	
4	1 rip* snv40	2NCO	
6	1 rip* snv22	CTT-TT
6	1 rip* snv22		
6	1 rip* snv40		
6	1 rip* snv40		
6	5 rip* snv40	2NCO	2 * RIP SNV22
6	1 rip* snv22		
6	1 rip* snv22		
6	1 rip* snv22		
6	1 rip* snv40		
6	1 rip* snv40		
*/
	
	// convert 2NCO to RIP, RIP to 2NCO
	// check 2NCO
	nco_detail_list.filter(nco => nco.type == "2NCO").forEach(_2nco => {
		if (_2nco.maybe_rip_snv40_markers) {
			const count_of_snp_marker = _2nco.GCasso_marker;

			if (_2nco.maybe_rip_snv40_markers.length == count_of_snp_marker) {
				_2nco.is_rip = true;
				if (typeof window == "object" && window != null) {
					console.log("2NCO -> RIP ?", _2nco);
				}
			}
			else {
				_2nco.n_rip_snv40 = _2nco.maybe_rip_snv40_markers.length;

				_2nco.maybe_rip_snv40_markers.forEach(rip => {
					rip.is_2nco = true;
					const _2nco_marker = rip;

					const m_info = snp_marker_mat[0][4][0];// 2NCO
					
					const col_idx = m_info.order;
					const marker = Object.assign(new MarkerValue(), _2nco_marker, {
						type: "SNV4:0",
						name: allMarker.list[col_idx].name,
						not: _2nco_marker.name,
					});
					allMarker.list[col_idx].values.push(marker);

					if (typeof window == "object" && window != null) {
						console.log("RIP -> 2NCO ?", rip);
					}
				});
			}
		}
	});
	if (options.mode == "tetrad") {
		["rip_Q", "rip_C", "rip_QC"].forEach(rip_marker_name => {
			// replace elements
			const list = allMarker.map[rip_marker_name].values.splice(0);
			// const not_2nco_list = list.filter(rip => !rip.is_2nco);
			const not_2nco_list = list.filter(rip => {
				if (rip.is_2nco) {
					seg_snp[rip.pos][0] = seg_snp[rip.pos - 1][0];
					seg_snp[rip.pos][1] = seg_snp[rip.pos - 1][1];
					seg_snp[rip.pos][2] = seg_snp[rip.pos - 1][2];
					seg_snp[rip.pos][3] = seg_snp[rip.pos - 1][3];
					seg_snp[rip.pos][4] = seg_snp[rip.pos - 1][4];
					seg_snp[rip.pos][5] = seg_snp[rip.pos - 1][5];
					seg_snp[rip.pos].is_rip = undefined;
					return false;
				}
				else {
					return true;
				}
				// return !rip.is_2nco;
			});
			// allMarker.map[rip_marker_name].values.push(...not_2nco_list);
			allMarker.map[rip_marker_name].values = not_2nco_list;
		});
	}
 
	// begin merge_seg

	_merge_seg(seq_list, seg_snp);

	/**
	 * @param {number} pos
	 */
	 function is_rDNA(pos) {
		if (options.show_rDNA_snp) {//ignore_rDNA
			return false;
		}
		if (options.rDNA_info && options.nChr == options.rDNA_info.chr) {
			if (pos >= options.rDNA_info.alignment_start && pos <= options.rDNA_info.alignment_end) {
				return true;
			}
		}
	};
	/**
	 * @param {number} pos
	 */
	function is_telomere(pos) {
		if (options.telomere != null) {
			if (options.telomere[0] != null && options.telomere[1] != null) {
				const left_end = options.telomere[0][1];// left arm end (3')
				const right_start = options.telomere[1][0];// right arm start (5')
				const ref1_pos = pos_ref1_uint32array[pos - 1];
				return (
					// pos <= ref1_pos_uint32array[left_end - 1] || pos >= ref1_pos_uint32array[right_start - 1] ||
					ref1_pos <= left_end || ref1_pos >= right_start
				);
			}
		}
		return false;
	}
	/**
	 * @type {(pos: number) => boolean}
	 */
	const point_filter = (function () {
		if (options.point_filter) {
			return function point_filter(pos) {
				return !is_rDNA(pos) && !is_telomere(pos);
			};
		}
		else {
			return function () { return true; };
		}
	})();

	if (options.ref_snp_list) {
		result_data.ref_snp_list = [...seq_list[0]].map((q, i) => [i + 1, q, seq_list[1][i]]).filter(([pos, q, c]) => point_filter(Number(pos)) && q != c);
		if (result_data.ref_snp_list.length == 0) {
			console.error("result_data.ref_snp_list == 0");
		}
	}

	allMarker.list.forEach(marker_info => {
		const c1 = marker_info.values.length;
		marker_info.values = marker_info.values.filter(marker => point_filter(marker.pos));
		const c2 = marker_info.values.length;
		console.log(marker_info.name, c1, c2);
	});
	
	return result_data;
}

/** 
 * @param {any[]} co_list
 * @returns {{ start: number, end: number, state: number[], raw_idx: number, is_CO: boolean, }[]}
 */
function get_co_detail(co_list) {
	let co_detail_list = [];
	let prev_snp_r1_pos = 0;
	for (let co_idx = 0; co_idx < co_list.length; ++co_idx) {
		let co = co_list[co_idx];
		// before co
		co_detail_list.push({
			start: prev_snp_r1_pos,
			end: co.snp_start_in - 1,
			state: co.before,
			raw_idx: co_idx,//co_idx
		});
		
		if (co.snp_start_in != co.snp_end_in) {
			// in co
			co_detail_list.push({
				start: co.snp_start_in,
				end: co.snp_end_in,
				// state: co.before,
				is_CO: true,
				raw_idx: co_idx,//co_idx
			});
		}

		prev_snp_r1_pos = co.snp_end_in;
	}
	// after co
	if (co_list.length > 1) {
		let last_co = co_list[co_list.length - 1];
		co_detail_list.push({
			start: prev_snp_r1_pos,
			end: Infinity,
			state: last_co.after,
			raw_idx: co_list.length,//co_idx
		});
	}
	return co_detail_list;
}

/**
 * @param {{ start: number, end: number, state: number[], raw_idx: number, is_CO: boolean, }} co_detail
 * @param {number} pos
 * @param {string} a s1
 * @param {string} b s2
 * @param {string} c s3
 * @param {string} d s4
 * @param {[string, string, string, string]} spores [a, b, c, d]
 * @param {string} ref1 ref1
 * @param {string} ref2 ref2
 * @param {boolean} cmp_ref1_ref2 ref1 == ref2
 * @param {*} row seg row
 * @param {[number, number, number, number, number, number]} all_prev_color 
 * @param {{[id:string]:number}} ColorID
 * @returns {0|1|2|3} 1: ref1, 2: ref2, 3: ref2
 */
function isRIP(co_detail, pos, a, b, c, d, spores, ref1, ref2, cmp_ref1_ref2, row, all_prev_color, ColorID) {
	if (co_detail && co_detail.is_CO && co_detail.state) {
		const co_state = co_detail.state;

		let seg_ref1 = 0;
		let seg_ref2 = 0;
		let mut_ref1 = 0;
		let mut_ref2 = 0;

		// let mut_progeny = [];
		if (pos == 610720) {
			debugger;
		}
		
		let spores_rip = co_state.map((spore_val, spore_idx) => {
			let spore_s = spores[spore_idx];
			if (spore_val == 0) {
				if (ref1 == spore_s) {
					++seg_ref1;

					return {
						ref: -1,
						type: null
					};
				}
				else if (ref1 == "G" && spore_s == "A") {
					++mut_ref1;
					return {
						ref: 1,
						type: "A"
					};
				}
				else if (ref1 == "C" && spore_s == "T") {
					++mut_ref1;
					return {
						ref: 1,
						type: "T"
					};
				}
				return {
					//illeg
					ref: -1,
					type: null
				};
			}
			else {
				if (ref2 == spore_s) {
					++seg_ref2;
					return {
						ref: -2,
						type: null
					};
				}
				else if (ref2 == "G" && spore_s == "A") {
					++mut_ref2;
					return {
						ref: 2,
						type: "A"
					};
				}
				else if (ref2 == "C" && spore_s == "T") {
					++mut_ref2;
					return {
						ref: 2,
						type: "T"
					};
				}
				return {
					//illeg
					ref: -2,
					type: null
				};
			}
		});
		if (mut_ref1 >= 2 || mut_ref2 >= 2) {
			// if (ref1_pos == 1071253) {
			// 	debugger;
			// }

			// row[0] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.dad) | (ref1 == "-" ? ColorID.indel_bit : 0);
			// row[1] = (cmp_ref1_ref2 ? ColorID.identical : ColorID.mom) | (ref2 == "-" ? ColorID.indel_bit : 0);

			row[0] = cmp_ref1_ref2 ? ColorID.identical : (ref1 != "-" ? ColorID.dad : all_prev_color[0]);
			row[1] = cmp_ref1_ref2 ? ColorID.identical : (ref2 != "-" ? ColorID.mom : all_prev_color[1]);

			// indel or diff

			let rs_indel = 0;
			let rs_diff = 0;

			spores_rip.forEach((rip, spore_idx) => {
				let bp = spores[spore_idx];
				if (rip.ref == -1) {
					row[2 + spore_idx] = (cmp_ref1_ref2 ? ColorID.identical : (bp == ref1 ? ColorID.dad : (bp == ref2 ? ColorID.mom : ColorID.diff))) | (a == "-" ? ColorID.indel_bit : 0);
					if (bp != ref1 && bp != ref2) {
						rs_diff++;
					}
				}
				else if (rip.ref == -2) {
					row[2 + spore_idx] = (cmp_ref1_ref2 ? ColorID.identical : (bp == ref1 ? ColorID.dad : (bp == ref2 ? ColorID.mom : ColorID.diff))) | (a == "-" ? ColorID.indel_bit : 0);
					if (bp != ref1 && bp != ref2) {
						rs_diff++;
					}
				}
				else if (rip.ref == 1) {
					row[2 + spore_idx] = ColorID.dad_rip;
				}
				else if (rip.ref == 2) {
					row[2 + spore_idx] = ColorID.mom_rip;
				}

				rs_indel += bp == "-" ? 1 : 0;
			});

			if (typeof window == "object" && window != null) {
				if (rs_diff != 2 && rs_indel != 2) {
					console.log(
						"pos:" + pos,
						"ill m:" + rs_diff,
						ref1, ref2, ...spores
					);
				}
				if (rs_indel == 1) {
					console.log(
						"1n:1 rip",
						"pos:" + pos,
						"rs_indel:", rs_indel,
						ref1, ref2, ...spores
					);
				}
				else if (rs_indel == 2) {
					// console.log(
					// 	"sss 2 rip",
					// 	"pos:" + pos,
					// 	"rs_indel:", rs_indel,
					// 	ref1, ref2, ...spores
					// 	);
				}
				else {
					console.log(
						"ill n + rip",
						"pos:" + pos,
						"rs_indel:", rs_indel,
						ref1, ref2, ...spores
					);
				}
			}

			if (rs_indel <= 2 || !rs_diff) {
				if (mut_ref1 >= 2 && mut_ref2 >= 2) {
					return 3;
				}
				else {
					return (mut_ref1 >= 2 ? 1 : (mut_ref2 >= 2 ? 2 : 0));
				}
			}
			else {
				return 0;
			}
			// return (mut_ref1 >= 2 ? 1 : (mut_ref2 >= 2 ? 2 : 0)) * ((rs_indel == 1 || has_illegitimate) ? -1 : 1);
		}
		else {
			return 0;
		}
		// let num_spores_rip = spores_rip.filter(a => a).length;
		// if (num_spores_rip == 2) {
		// 	// 2:2 rip
		// 	//return true;
		// }
		// else if (num_spores_rip > 2) {
		// 	spores.filter(v => ref1 == v).length == num_spores_rip;
		// 	spores.filter(v => ref2 == v).length == num_spores_rip;
		// 	//return true;
		// }
		// else if (num_spores_rip == 1) {
		// 	console.log({
		// 		type: "illeg", ref1_pos, ref1, ref2, a, b, c, d, co_state,
		// 	});
		// 	return false;
		// }
		// else {
		// 	return false;//not rip
		// }
	}
	else {
		let is_rip = 0;
		let rip;
		rip = check_rip(ref1, ref2, ColorID.dad_rip, ColorID.diff);
		if (rip == ColorID.dad_rip) {
			row[0] = cmp_ref1_ref2 ? ColorID.dad : ColorID.identical;
			row[1] = ref2 != "-" ? ColorID.mom : all_prev_color[1];
			is_rip = 1;
		}
		else if (rip == ColorID.mom_rip) {
			row[0] = ref1 != "-" ? ColorID.dad : all_prev_color[0];
			row[1] = cmp_ref1_ref2 ? ColorID.mom : ColorID.identical;
			is_rip = 2;
		}
		else {
			rip = check_rip(ref2, ref1, ColorID.mom_rip, ColorID.diff);
			if (rip == ColorID.dad_rip) {
				row[0] = cmp_ref1_ref2 ? ColorID.dad : ColorID.identical;
				row[1] = ref2 != "-" ? ColorID.mom : all_prev_color[1];
				is_rip = 1;
			}
			else if (rip == ColorID.mom_rip) {
				row[0] = ref1 != "-" ? ColorID.dad : all_prev_color[0];
				row[1] = cmp_ref1_ref2 ? ColorID.mom : ColorID.identical;
				is_rip = 2;
			}
		}

		return is_rip;
	}

	/**
	 * @param {string} ref_1
	 * @param {string} ref_2
	 * @param {number} col
	 * @param {number} col2
	 */
	function check_rip(ref_1, ref_2, col, col2) {
		let rip = 0;
		if (ref_1 == "G") { //GA AAAA
			if (ref_2 == "A") {
				if (a == "A" && b == "A" && c == "A" && d == "A") {
					//G->A
					rip = col;
					//parental_cmp_uint8array[pos] = col;
					row[2] = col;
					row[3] = col;
					row[4] = col;
					row[5] = col;
				}
			}
			else { //G? AA??
				let count_a = spores.filter(s => s == "A").length;
				let count_x = spores.filter(s => s == ref_2).length;
				if (count_a >= 2) {//2:2, 3:1
					if (count_x == (4 - count_a)) {//2:2, 3:1
						//G->A
						rip = col;
						//parental_cmp_uint8array[pos] = col;
						row[2] = a == "A" ? col : (a == ref_1 ? ColorID.dad : (a == ref_2 ? ColorID.mom : col2)) | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = b == "A" ? col : (b == ref_1 ? ColorID.dad : (b == ref_2 ? ColorID.mom : col2)) | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = c == "A" ? col : (c == ref_1 ? ColorID.dad : (c == ref_2 ? ColorID.mom : col2)) | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = d == "A" ? col : (d == ref_1 ? ColorID.dad : (d == ref_2 ? ColorID.mom : col2)) | (d == "-" ? ColorID.indel_bit : 0);
					}
					else {
						//G??A
						rip = 0; //fill default color
					}
				}
			}
		}
		else if (ref_1 == "C") { //CT TTTT
			if (ref_2 == "T") {
				if (a == "T" && b == "T" && c == "T" && d == "T") {
					//C->T
					rip = col;
					//parental_cmp_uint8array[pos] = col;
					row[2] = col;
					row[3] = col;
					row[4] = col;
					row[5] = col;
				}
			}
			else { //C? TT??
				let count_t = spores.filter(s => s == "T").length;
				let count_x = spores.filter(s => s == ref_2).length;
				if (count_t >= 2) {//2:2, 3:1
					if (count_x == (4 - count_t)) {//2:2, 3:1
						//C->T
						rip = col;
						//parental_cmp_uint8array[pos] = col;
						row[2] = a == "T" ? col : (a == ref_1 ? ColorID.dad : (a == ref_2 ? ColorID.mom : col2)) | (a == "-" ? ColorID.indel_bit : 0);
						row[3] = b == "T" ? col : (b == ref_1 ? ColorID.dad : (b == ref_2 ? ColorID.mom : col2)) | (b == "-" ? ColorID.indel_bit : 0);
						row[4] = c == "T" ? col : (c == ref_1 ? ColorID.dad : (c == ref_2 ? ColorID.mom : col2)) | (c == "-" ? ColorID.indel_bit : 0);
						row[5] = d == "T" ? col : (d == ref_1 ? ColorID.dad : (d == ref_2 ? ColorID.mom : col2)) | (d == "-" ? ColorID.indel_bit : 0);
					}
					else {
						//C??T
						rip = 0; //fill default color
					}
				}
			}
		}
		return rip;
	}
}

function init_viewModel(seq_list) {
	region_rect = seq_list.map(_ => []);
	//parental_cmp_uint8array = new Uint8Array(seq_list[0].length);
	ref1_pos_uint32array = new Uint32Array(seq_list[0].length);
	pos_ref1_uint32array = new Uint32Array(seq_list[0].length);
	ref2_pos_uint32array = new Uint32Array(seq_list[0].length);
	pos_ref2_uint32array = new Uint32Array(seq_list[0].length);
	//rip_map_uint8array = new Uint8Array(seq_list[0].length);
	ref1_ref2_score_uint32array = new Uint32Array(seq_list[0].length);
	allMarker.list.forEach(a => a.values.splice(0));//clear all
	align_start_index = -1;
	align_end_index = -1;
}

function _merge_seg(seq_list, seg_snp, max_reg_len) {
	for (let sid = 0; sid < seq_list.length; ++sid) {
		let start = 0;
		let prev = seg_snp[start][sid];
		let i = 1;
		for (; i < seq_list[0].length; ++i) {
			let col = seg_snp[i][sid];
			if (prev != col) {
				region_rect[sid].push({
					start: start,
					end: i,
					col: prev,
				});
				start = i;
				prev = col;
			}
		}
		// for (; i < seq_list[0].length;) {
		// 	let col = seg_snp[i][sid];
		// 	if (prev != col) {
		// 		region_rect[sid].push({
		// 			start: start,
		// 			end: i,
		// 			col: prev,
		// 		});
		// 		start = i;
		// 		prev = col;
		// 		continue;
		// 	}
		// 	++i;
		// }
		if (region_rect[sid].length && i != region_rect[sid][region_rect[sid].length - 1].end) {
			region_rect[sid].push({
				start: start,
				end: i,
				col: prev,
			});
		}
	}
}

if (typeof window != "object") {
	module.exports.loadCmpData = loadCmpData;
	module.exports.initAnalyser = initAnalyser;
}

function check_4n0() {
	return [...seq_list[0]].map((a, i) => seq_list.slice(2).every(b => b == "-") ? i : -1).filter(i => i >= 0);
}
function check_CCTTTT() {
	return [...seq_list[0]].map((q, i) => [
		i,
		seq_list.slice(0, 2).every(ss => ss[i] == "C") &&
		seq_list.slice(2).every(ss => ss[i] == "T")
	]).filter(a => a[1]);
}
function check_GGAAAA() {
	return [...seq_list[0]].map((q, i) => [
		i,
		seq_list.slice(0, 2).every(ss => ss[i] == "G") &&
		seq_list.slice(2).every(ss => ss[i] == "A")
	]).filter(a => a[1]);
}


/**
 * @param {any} seq_list
 * @param {number} ma_start
 * @param {number} ma_end
 * @param {number} [pos]
 */
function extract_seq(seq_list, ma_start, ma_end, pos = Math.trunc((ma_start + ma_end) / 2)) {
	const seq_names = [...dataset.parental_list, ...dataset.progeny_list];

	const range_info = [
		["a", ma_start, ma_end].join("-"),
		["Q", pos_ref1_uint32array[ma_start], pos_ref1_uint32array[ma_end]].join("-"),
		["C", pos_ref2_uint32array[ma_start], pos_ref2_uint32array[ma_end]].join("-"),
	].join(" ");

	const fa = seq_list.map((ss, i) => {
		const head = ">" + seq_names[i] + " " + range_info + "\n";
		const body = ss.slice(ma_start - 1, ma_end).replace(/-/g, "");
		return head + body;
	}).join("\n");

	return fa;
}

/**
 * @param {string} ref
 * @param {number} sIdx
 * @param {number} nChr
 * @param {any} seq_list
 * @param {number} view_range
 * @param {number} extract_range
 */
function find_rip_out_AT_island(ref, sIdx, nChr, seq_list, view_range = 30, extract_range = 500) {
	const seq_names = [...dataset.parental_list, ...dataset.progeny_list];
	const ref1_name = seq_names[0];

	// const fm = {
	// 	"A": 0,
	// 	"T": 0,
	// 	"C": 0,
	// 	"G": 0,
	// };
	// [...seq_list[sIdx]].forEach(a => ++fm[a]);
	// // ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06)
	// const block_min_gc = ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06) * 100;

	const block_min_gc = dataset.min_gc_content[sIdx];
	
	/** @type {"SNV RIP"|"RIP"} */
	const rip_snp_type = "RIP";
	const rip_markers_list = allMarker.list.filter(a => a.name.indexOf(rip_snp_type) >= 0);

	const map_pos_ref1 = pos_ref1_uint32array;
	const map_ref1_pos = ref1_pos_uint32array;

	console.log({
		block_min_gc: block_min_gc,
		rip_markers_list: rip_markers_list.length,
	});

	const result = rip_markers_list.map(grp => {
		return grp.values.map(marker => {
			const gc_wnd = findGCwndByMarker(ref, nChr, map_pos_ref1, map_ref1_pos, marker);
			// console.log("s1", nChr, a.pos, gc_wnd.gc);
			return {
				pos: a.pos,
				// seq: seq_list.map(ss => ss.slice(a.pos - range, a.pos + range)).join("\n") + "\n" + " ".repeat(range) + "^",
				// p0s: seq_list.map(ss => ss[a.pos]).join(""),
				gc_wnd: gc_wnd,
				type: grp.name,
			};
		}).filter(a => a.gc_wnd.gc > block_min_gc);//.forEach(a => console.log(a.pos, a.gc) + console.log(a.seq))
	});

	let fin_output_text = result.map(grp => {
		return grp.map(mark => {
			const pos = mark.pos;
			const ref1_pos = map_pos_ref1[mark.pos];
			const view = seq_list.map(ss => ss.slice(mark.pos - view_range, mark.pos + view_range)).join("\n") + "\n" + " ".repeat(view_range) + "^";
			const gc_value = mark.gc_wnd.gc;
			const gc_ref1_start = mark.gc_wnd.start;
			const gc_ref1_end = mark.gc_wnd.end;
			const ript_type = mark.type;

			const out_seq = extract_seq(seq_list, pos - extract_range, pos + extract_range, ref1_pos);

			// console.log("s2", nChr, pos);

			const info = {
				["Ch"]: nChr,
				["a pos"]: pos,
				[`${ref1_name} pos`]: ref1_pos,
				["GC%"]: gc_value,
				[`GC window ${ref1_name} start`]: gc_ref1_start,
				[`GC window ${ref1_name} end`]: gc_ref1_end,
				["ript_type"]: ript_type,
			};
			let output_text = [
				JSON.stringify(info, null, "\t"),
				view,
			].join("\n");

			const table_row = Object.keys(info).map(k => info[k]);

			return [output_text, out_seq, String(pos), table_row];
		});
	});
	
	// console.log({
	// 	fin_output_text: fin_output_text.length,
	// 	fin_output_text_map: fin_output_text.map(a => a.length).join("-"),
	// });

	return fin_output_text;
}

/**
 * @param {number} gff_sIdx
 * @param {number} nChr
 * @param {seq_list} seq_list
 * @param {any} gff
 * @param {any} analysis_results
 * @param {Partial<AnalysisOptions>} options
 */
function find_rip_gff(gff_sIdx, nChr, seq_list, gff, analysis_results, options) {
	const gff_ref = dataset.parental_list[gff_sIdx];

	// const seq_names = [...dataset.parental_list, ...dataset.progeny_list];
	// const ref_seq_name = seq_names[sIdx];

	console.log(Object.keys(analysis_results).map(k => [k, typeof analysis_results[k]]));

	// const fm = {
	// 	"A": 0,
	// 	"T": 0,
	// 	"C": 0,
	// 	"G": 0,
	// };
	// [...seq_list[0]].forEach(a => ++fm[a]);
	// // ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06)
	// const block_min_gc = ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06) * 100;
	
	/** @type {"SNV RIP"|"RIP"} */
	const rip_snp_type = "RIP";
	const rip_markers_list = analysis_results.allMarker.list.filter(a => a.name.indexOf(rip_snp_type) >= 0);

	const map_pos_ref = [
		analysis_results.pos_ref1_uint32array,
		analysis_results.pos_ref2_uint32array,
	][gff_sIdx];
	const map_ref_pos = [
		analysis_results.ref1_pos_uint32array,
		analysis_results.ref2_pos_uint32array,
	][gff_sIdx];

	console.log({
		// block_min_gc: block_min_gc,
		rip_markers_list: rip_markers_list.length,
	});

	const results = [];
	
	rip_markers_list.map(grp => {
		return grp.values.map(marker => {
			const gene = findAllGeneByMarker(gff, map_pos_ref, marker)[0];
			
			const gc_wnd = findGCwndByMarker(gff_ref, nChr, map_pos_ref, map_ref_pos, marker);
			
			const AT_island = gc_wnd.gc <= dataset.min_gc_content[gff_sIdx] ? "+" : "-";

			marker.color = gene ? options.gene_rip_color[gene.type] : "";
			// if (marker.color) {
			// 	console.log(gene.type, marker.pos);
			// }
			
			return {
				pos: marker.pos,
				marker: marker,
				gene: gene,
				rip_type: grp.name,

				gc_value: gc_wnd.gc,
				AT_island,
			};
		}).filter(a => a && a.gene).forEach(a => {
			results.push(a);
		});
	});

	results.sort((a, b) => a.pos - b.pos);

	// results.forEach(pair => {
	// 	const {
	// 		gene,
	// 		marker,
	// 	} = pair;
	// });

	return results;
}

/**
 * @param {GFF_ROW[]} gene_list
 * @param {Uint32Array|number[]} map_pos_ref
 * @param {MarkerValue} marker
 */
function findAllGeneByMarker(gene_list, map_pos_ref, marker) {
	return gene_list.filter(gene => {
		if (gene.start <= map_pos_ref[marker.pos] && map_pos_ref[marker.pos] <= gene.end ||
			gene.$start <= marker.pos && marker.pos <= gene.$end) {
			return true;
		}
	});
}

/**
 * @param {string} ref
 * @param {number} nChr
 * @param {Uint32Array|number[]} map_pos_ref
 * @param {Uint32Array|number[]} map_ref_pos
 * @param {MarkerValue} marker
 */
function findGCwndByMarker(ref, nChr, map_pos_ref, map_ref_pos, marker) {
	const gc_wnd = dataset.gc_content[ref][nChr].find(b => {
		return (
			b.start <= map_pos_ref[marker.pos] && map_pos_ref[marker.pos] <= b.end ||
			map_ref_pos[b.start] <= marker.pos && marker.pos <= map_ref_pos[b.end]
		);
	});
	if (!gc_wnd) {
		console.error(marker);
		console.error({
			ref, nChr
		});
	}
	return gc_wnd;
}

/**
 * @param {RepeatSegment[]} repeat_segment
 * @param {MarkerValue} marker
 * @param {Uint32Array|number[]} ref_pos_map
 * @param {Uint32Array|number[]} pos_ref_map
 */
function find_repeat_segment_by_marker(repeat_segment, marker, ref_pos_map, pos_ref_map) {
	const m_ref_pos = pos_ref_map[marker.pos]

	return repeat_segment.filter(function (segment) {
		const start = ref_pos_map[segment.q_start];
		const end = ref_pos_map[segment.q_end];

		if ((m_ref_pos >= segment.q_start && m_ref_pos <= segment.q_len) ||
			(marker.pos >= start && marker.pos <= end)
		) {
			return true;
		}
	});
}

/**
 * @param {Float32Array[]} me_list
 * @param {MarkerValue} marker
 * @param {number} left_size
 * @param {number} right_size
 * @param {Uint32Array|number[]} pos_ref_map
 * @param {number} chr_length
 */
function find_methyl_ratio_by_marker(me_list, marker, left_size, right_size, pos_ref_map, chr_length) {
	// const ref_pos = pos_ref_map[marker.pos];
	// const results = [];

	// for (let i = ref_pos - 1; i > 0; i--) {
	// 	const val = methyl_ratio[i];
	// 	if (!Number.isNaN(val) && Number.isFinite(val)) {
	// 		results.unshift((val * 100).toFixed(0));
	// 	}
	// 	if (results.length >= left_size) {
	// 		break;
	// 	}
	// }
	// for (let i = left_size - results.length; i > 0; --i) {
	// 	results.unshift("");
	// }

	// const limit = left_size + 1 + right_size;

	// for (let i = ref_pos; i < chr_length; ++i) {
	// 	const val = methyl_ratio[i];
	// 	if (!Number.isNaN(val) && Number.isFinite(val)) {
	// 		results.push((val * 100).toFixed(0));
	// 	}
	// 	if (results.length >= limit) {
	// 		break;
	// 	}
	// }

	// if (results.length != 11) {
	// 	console.log(results);
	// }

	// return results;

	const ref_pos = pos_ref_map[marker.pos];
	
	// const methyl_ratio = me_list[0];
	// return [
	// 	...methyl_ratio.subarray(ref_pos - left_size, ref_pos + right_size + 1)
	// ].map(val => Number.isNaN(val) ? "" : (val * 100).toFixed(0));

	const results = [];
	const start = ref_pos - left_size;
	const end = ref_pos + right_size
	for (let pos = start; pos <= end; ++pos) {
		const aa = [];
		for (let j = 0; j < me_list.length; ++j) {
			const val = getVal(me_list[j], pos);
			aa.push(val);
		}
		results.push(aa);
	}
	/**
	 * @param {Float32Array} arr
	 * @param {number} i
	 */
	function getVal(arr, i) {
		if (i < 0 || i >= arr.length) {
			return "";
		}
		else {
			const val = arr[i];
			if (!Number.isNaN(val)) {
				return (val * 100).toFixed(0);
			}
			else {
				return "N";
			}
		}
	}
	return results;

	// const results = [];
	// const start = ref_pos - left_size;
	// const end = ref_pos + right_size
	// for (let i = start; i <= end; ++i) {
	// 	console.log(getVal(me_list[0], i));
	// }
	// function getVal(arr, i) {
	// const val = arr[i];
	// return val !== undefined ? (val * 100).toFixed(0) : "";
	// }
}
/**
 * @param {Float32Array} methyl_ratio
 * @param {MarkerValue} marker
 * @param {number} left_size
 * @param {number} right_size
 * @param {Uint32Array|number[]} pos_ref_map
 * @param {number} chr_length
 */
function find_methyl_ratio_single_by_marker(methyl_ratio, marker, left_size, right_size, pos_ref_map, chr_length) {
	const ref_pos = pos_ref_map[marker.pos];

	const results = [];
	const start = ref_pos - left_size;
	const end = ref_pos + right_size
	for (let pos = start; pos <= end; ++pos) {
		const val = getVal(methyl_ratio, pos);
		results.push(val);
	}
	/**
	 * @param {Float32Array} arr
	 * @param {number} i
	 */
	function getVal(arr, i) {
		if (i < 0 || i >= arr.length) {
			return "";
		}
		else {
			const val = arr[i];
			if (!Number.isNaN(val)) {
				return (val * 100).toFixed(0);
			}
			else {
				return "";
			}
		}
	}
	return results;
}

function test_find_methyl_ratio_by_marker() {
	const ref_pos = 6;
	const chr_length = 12;
	const methyl_ratio = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
	const results = [];
	
	for (let i = ref_pos - 1; i > 0; i--) {
		const val = methyl_ratio[i];
		if (!Number.isNaN(val) && Number.isFinite(val)) {
			results.unshift(val);
		}
		if (results.length >= 5) {
			break;
		}
	}
	for (let i = 5 - results.length; i > 0; --i) {
		results.unshift(NaN);
	}

	for (let i = ref_pos; i < chr_length; ++i) {
		const val = methyl_ratio[i];
		if (!Number.isNaN(val) && Number.isFinite(val)) {
			results.push(val);
		}
		if (results.length >= 11) {
			break;
		}
	}
	results
}

