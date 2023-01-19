// @ts-check
// 20230119: 2023/01/19 clean

// 5mC
// TSETA dataset: QM6a_CBS1-1_WT_dim2_rid1_dimrid
//  all_preset_20221104() // set viewport
//  _load_presset_20220916() // load 5mc, C-to-T
//  _load_presset_20221025_QM6a_RNA() // load RNA

// ssDNA
// TSETA dataset: Tr_ssDNA_parent
//  preset_ssDNA_venn_minimal()
// 	ssDNA_peak_venn()
// 	ssDNA_peak_Rvenn()
// 	ssDNA_20211014()
// 	ssDNA_AgeI_all_chr()

[...document.querySelectorAll("style")].forEach(a => a.innerHTML = a.innerHTML.replace("transition: 1s all;", ""));

let promise_load_task = Promise.resolve();// init, lock layout

// window.version_2 = true;//-20210111
let CNChuang_spore_1 = false;

const time_app_start = new Date();

let CONFIG_LOAD_METH = false;

const METHYL_MODE = true;

/**
 * @type {{ [refId:string]:{ [seqid:string]:GFF_ROW[] } }}
 */
let gff_data_map = null;

// /** @type {{ [refId:string]: { [seqid: string]: { [queryId: number]: GFF_ROW[]; } } }} - queryId = gene.start / 100 */
// let gff_data_queryMap = null;

const dataset = JSON.parse(document.getElementById("dataset.json").innerText);
dataset._gc_content = JSON.parse(JSON.stringify(dataset.gc_content));//deep clone plain data


const RNA_cmp_Pvalue = 0.01;
const RNA_cmp_FC = 10;//
const RNA_cmp_C = Math.log2(RNA_cmp_FC).toFixed(6);

// draw_GTF_gene_s
const blend_lncRNA_color = {
	[0b000]: { fill: "#AAAAAA7F", stroke: "#AAAAAAFF" },

	[0b001]: { fill: "#0000FF7F", stroke: "#0000FFFF" },// veg
	[0b010]: { fill: "#007F007F", stroke: "#007F00FF" },// glucose
	[0b011]: { fill: "#00FFFF7F", stroke: "#00FFFFFF" },
	[0b100]: { fill: "#FF00007F", stroke: "#FF0000FF" },//cellulose
	[0b101]: { fill: "#FF00FF7F", stroke: "#FF00FFFF" },
	[0b110]: { fill: "#FF7F007F", stroke: "#FF7F00FF" },

	// [0b001]: { fill: "#0000FF7F", stroke: "#0000FFFF" },
	// [0b010]: { fill: "#007F007F", stroke: "#007F00FF" },
	// [0b011]: { fill: "#007FFF7F", stroke: "#007FFFFF" },
	// [0b100]: { fill: "#FF00007F", stroke: "#FF0000FF" },
	// [0b101]: { fill: "#FF00FF7F", stroke: "#FF00FFFF" },
	// [0b110]: { fill: "#FF7F007F", stroke: "#FF7F00FF" },

	[0b111]: { fill: "#0000007F", stroke: "#000000FF" },
};

class CmpSubset {
	id = "";
	logFC = 0;
	logCPM = 0;
	PValue = 0;
	FDR = 0;

	/**
	 * @param {Partial<CmpSubset>} obj
	 */
	constructor(obj) {
		this.id = obj.id;
		this.logFC = obj.logFC;
		this.logCPM = obj.logCPM;
		this.PValue = obj.PValue;
		this.FDR = obj.FDR;
	}
}

/**
 * @type {{ [up: string]: { [down: string]: CmpSubset[]; }; }}
 */
const RNA_cmp = {
};


class GeneClusterInfo {
	/** @type {Map<string,GFF_ROW>} <cluster geneId, annotation> */
	data = new Map();

	name = "";

	/** @type {string[]} */
	urls = [];

	/**
	 * @param {string} name
	 */
	constructor(name) {
		this.name = name;
	}
}

/** @type {GeneClusterInfo[]} */
const gene_cluster = [
];

/**
 * @param {GTF_ROW[]} gff3
 */
async function load_gene_cluster_veg_glu_Cel(gff3) {
	// V, G, C
	// const cl_list = [
	// 	{ type: ["C", ">", "V", "=", "G"], subcluster: [], }, // --+
	// 	{ type: ["V", ">", "C", "=", "G"], subcluster: [], }, // +--
	// 	{ type: ["V", ">", "G", ">", "C"], subcluster: [], }, // +--
	// 	{ type: ["G", ">", "C", "=", "V"], subcluster: [], }, //
	// 	{ type: ["C", "=", "V", ">", "G"], subcluster: [], }, //
	// 	{ type: ["C", "=", "G", ">", "V"], subcluster: [], }, //
	// ];
	// const cl_list = [
	// 	{ type: ["C", "=", "V", ">", "G"], subcluster: [9, 8, 2], }, // C=V>G	++-	9	8	2
	// 	{ type: ["C", ">", "V", "=", "G"], subcluster: [7, 1, 4], }, // C>V=G	+--	7	1	4
	// 	{ type: ["V", ">", "C", "=", "G"], subcluster: [5, 11],   }, // V>C=G	-+-	5	11
	// 	{ type: ["C", "<", "V", "=", "G"], subcluster: [6],       }, // C<V=G	-++	6
	// 	{ type: ["G", ">", "V", "=", "C"], subcluster: [3],       }, // G>V=C	--+	3
	// 	{ type: ["C", "=", "G", ">", "V"], subcluster: [10],      }, // C=G>V	+-+	10
	// ];
	const cl_list = [
		{ type: ["C", "=", "V", ">", "G"], subcluster: [9, 10, ], }, // C=V>G	++-	9	8	2
		{ type: ["C", ">", "V", "=", "G"], subcluster: [1, 4, 6], }, // C>V=G	+--	7	1	4
		{ type: ["V", ">", "C", "=", "G"], subcluster: [3, 7,  ], }, // V>C=G	-+-	5	11
		{ type: ["C", "<", "V", "=", "G"], subcluster: [8,     ], }, // C<V=G	-++	6
		{ type: ["G", ">", "V", "=", "C"], subcluster: [2, 5,  ], }, // G>V=C	--+	3
		{ type: ["C", "=", "G", ">", "V"], subcluster: [11,    ], }, // C=G>V	+-+	10
	];

	const name_map = {
		"V": "veg",
		"C": "cellulose",
		"G": "glucose",
		">": "＞",
		"<": "＜",
		"=": "＝",
	};
	gene_cluster.splice(0);// clear

	// make ordered list
	cl_list.forEach(info => {
		const cluster_name = info.type.map(a => name_map[a]).join(" ");

		const subcluster = new GeneClusterInfo(cluster_name);

		gene_cluster.push(subcluster);
	});

	const tasks = cl_list.map(async (info, cluster_idx) => {
		const subcluster = gene_cluster[cluster_idx];

		const sub_task = info.subcluster.map(async subcluster_id => {
			const url = `${LNCRNA_ROOT_DIR}/lncRNA_QM6a/subcluster_${subcluster_id}_log2_medianCentered_fpkm.matrix`;

			subcluster.urls.push(url);

			/** @type {string} */
			const text = await fetchData(url, "text");
			text.split("\n").slice(1).forEach(line => {
				const cols = line.split("\t");
				const gene_name = cols[0].trim();

				const cluster_gene = gff3.find(a => a.type == "gene" && a.geneID == gene_name);
				if (cluster_gene) {
					subcluster.data.set(gene_name, cluster_gene);
					cluster_gene.cluster_idx = cluster_idx;
				}
			});
		});

		await Promise.all(sub_task);
	});

	await Promise.all(tasks);

	return gene_cluster;
}

/**
 * @param {GTF_ROW[]} gff3
 */
 async function load_gene_cluster_glu_Cel(gff3) {
	const cl_list = [
		{ type: ["C", ">", "G"], subcluster: [1, ], },
		{ type: ["G", ">", "C"], subcluster: [2, ], }, // C=G>V	+-+	10
	];

	const name_map = {
		"C": "cellulose",
		"G": "glucose",
		">": "＞",
		"<": "＜",
	};
	gene_cluster.splice(0);// clear

	// make ordered list
	cl_list.forEach(info => {
		const cluster_name = info.type.map(a => name_map[a]).join(" ");

		const subcluster = new GeneClusterInfo(cluster_name);

		gene_cluster.push(subcluster);
	});

	const tasks = cl_list.map(async (info, cluster_idx) => {
		const subcluster = gene_cluster[cluster_idx];

		const sub_task = info.subcluster.map(async subcluster_id => {
			const url = `20210322_glu_vs_cel/QM6a/subcluster_${subcluster_id}_log2_medianCentered_fpkm.matrix`;

			subcluster.urls.push(url);

			/** @type {string} */
			const text = await fetchData(url, "text");
			text.split("\n").slice(1).forEach(line => {
				const cols = line.split("\t");
				const gene_name = cols[0].trim();

				const cluster_gene = gff3.find(a => a.type == "gene" && a.geneID == gene_name);
				if (cluster_gene) {
					subcluster.data.set(gene_name, cluster_gene);
					cluster_gene.cluster_idx = cluster_idx;
				}
			});
		});

		await Promise.all(sub_task);
	});

	await Promise.all(tasks);

	return gene_cluster;
}

let LNCRNA_ROOT_DIR = "lncRNA";
class LncRNAData {
	/**
	 * @param {string} path
	 * @param {string} genome
	 * @param {string} $name
	 * @param {{ has_gtf: boolean; has_coverage: boolean; has_gff3: boolean; }}
	 */
	constructor(path, genome, $name, { has_gtf, has_coverage, has_gff3 }) {
		this.$name = $name;
		this.path = path;
		this.genome = `${LNCRNA_ROOT_DIR}/${path}/${genome}`;
		this.transcripts = `${LNCRNA_ROOT_DIR}/${path}/${$name}.tx.fa`;
		this.lncRNA_list = `${LNCRNA_ROOT_DIR}/${path}/${$name}_lncRNA_list.txt`;

		if ((has_gtf || has_gff3) && !this.lncRNA_list) {
			throw new TypeError(`this.lncRNA_list = '${this.lncRNA_list}'`);
		}

		if (has_gtf) {
			this.gtf = `${LNCRNA_ROOT_DIR}/${path}/${$name}.gtf`;
		}
		if (has_coverage) {
			this.coverage = `${LNCRNA_ROOT_DIR}/${path}/${$name}_coverage.uint32`;
		}
		// if (has_gff3) {
		// 	this.gff3 = `${$name}.gff3`;
		// }

		/**
		 * @type {GTF_ROW[]}
		 */
		this.gtf_data = null;

		/**
		 * @type {GTF_ROW[]}
		 */
		this.gff3_data = null;

		//

		this.tpm = `${LNCRNA_ROOT_DIR}/${path}/${$name}_TPM.txt`;

		/** @type {{ [geneID:string]:number; }} */
		this.tpm_data = null;

		/** @type {{ [geneID:string]:{ [sampleName:string]:number; }; }} */
		this.cmp_tpm_data = null;

		/** @type {number} */
		this.max_tpm_value = 0;

		//

		this.tx_tpm = `${LNCRNA_ROOT_DIR}/${path}/${$name}_tx_TPM.txt`;

		/** @type {{ [geneID:string]:number; }} */
		this.tx_tpm_data = null;

		/** @type {number} */
		this.tx_max_tpm_value = 0;
	}
}

let QM6a_transposome = {};

// const chr_new_name = {
// 	"ChI_QM6a": "QM6a_1",
// 	"ChII_QM6a": "QM6a_2",
// 	"ChIII_QM6a": "QM6a_3",
// 	"ChIV_QM6a": "QM6a_4",
// 	"ChV_QM6a": "QM6a_5",
// 	"ChVI_QM6a": "QM6a_6",
// 	"ChVII_QM6a": "QM6a_7",
// };
const chr_new_name = {
	"ChI_QM6a": "scaffold_1",
	"ChII_QM6a": "scaffold_2",
	"ChIII_QM6a": "scaffold_3",
	"ChIV_QM6a": "scaffold_4",
	"ChV_QM6a": "scaffold_5",
	"ChVI_QM6a": "scaffold_6",
	"ChVII_QM6a": "scaffold_7",
};

const name_mapto_html = {
	"QM6a": `QM6a (<i>MAT1-2</i>)`,
	"CBS1-1": `CBS1-1 (<i>MAT1-1</i>)`,
	// "WTH109": ``,
	// "WTH111": ``,
	// "WTH115": ``,
	// "WTH119": ``,
};

/**
 * rename_html_seq has been deprecated, use rename_seq_row_header instead
 * @deprecated
 * @param {boolean} insert_prefix
 */
function rename_html_seq(insert_prefix = false) {
	console.warn("rename_html_seq has been deprecated, use rename_seq_row_header instead");
	return rename_seq_row_header(insert_prefix = false);
}

function rename_seq_row_header(insert_prefix = false) {
	try {
		["Dad", "Mom"].forEach((parental, pid) => {
			document.querySelectorAll(`.${parental}-name`).forEach(el => {
				if (el.parentElement.classList.contains("gc-plot")) {
					el.innerHTML = "";
					return;
				}

				const old_name = dataset.parental_list[pid].replace("_", " ");
				const new_name = name_mapto_html[old_name];

				if (insert_prefix) {
					if (new_name != null) {
						el.innerHTML = `${parental} ${new_name}`;
					}
					else {
						el.innerHTML = `${parental} (${old_name})`;
					}
				}
				else {
					el.innerHTML = new_name != null ? new_name : old_name;
				}
			});
		});

		["F1_A", "F1_B", "F1_C", "F1_D"].forEach((progeny, pid) => {
			document.querySelectorAll(`.${progeny}-name`).forEach(el => {
				if (insert_prefix) {
					el.innerHTML = `${progeny.replace("_", " ")} (${dataset.progeny_list[pid].replace(/_/g, " ")})`;
				}
				else {
					el.innerHTML = dataset.progeny_list[pid].replace(/_/g, " ");
				}
			});
		});
	}
	catch (ex) {
		console.error("warn:", ex);
	}
}

// CONFIG_LOAD_METH = false;// dataset.mode == "tetrad";

const reverse_parental = false;
if (reverse_parental) {
	const new_parental_list = [...dataset.parental_list].reverse();
	dataset.ref = new_parental_list[0];
	dataset.parental_list = new_parental_list;

	// const old_parental_map = Object.assign({}, dataset.parental);
	// dataset.parental = new_parental_list.reduce((obj, key) => {
	// 	obj[key] = old_parental_map[key];
	// 	return obj;
	// }, {});

	// genomeNameList;
	// genomeFileMap;
}

let seq_list = [""];
// let methy_list = [];

class ChromosomeCoordinate {
	mouse_bp_pos = 1;
	_bp_start = 0;
	_bp_end = 0;

	get bp_start() { return this._bp_start; }
	get bp_end() { return this._bp_end; }

	set bp_start(val) {
		this._bp_start = val;
	}
	set bp_end(val) {
		this._bp_end = val;
	}
	
	b_draw_selected_range = false;
	mouse_bp_pos_0 = 0;

	/** @type {Uint32Array[]} */
	map_pos_to_ref = [];

	/** @type {Uint32Array[]} */
	map_ref_to_pos = [];

	/**
		* 20220912 set_bounding
	 * make_pos_ref_map_list
	 * make_ref_pos_map_list
	 * @param {string[]} seq_list
	 * @see {g_methylRenderer.setLoadingBounding}
	 */
	setSeqList(seq_list) {
		for (let ref_idx = 0; ref_idx < dataset.genome_info_list.length; ++ ref_idx) {
			const map_pos_to_ref = multialign_to_chrPos_posMap(seq_list[ref_idx]);
			const map_ref_to_pos = chrPos_to_multialign_posMap(seq_list[ref_idx]);

			this.map_pos_to_ref[ref_idx] = map_pos_to_ref;
			this.map_ref_to_pos[ref_idx] = map_ref_to_pos;
		}
	}

	/**
	 * @param {number} ref_idx
	 * @param {number} pos
	 * @returns {number} this.map_pos_to_ref[ref_idx][pos - 1];
	 */
	pos_to_ref(ref_idx, pos) {
		return this.map_pos_to_ref[ref_idx][pos - 1];
	}
	
	/**
	 * @param {number} ref_idx
	 * @param {number} ref_pos
	 * @returns {number} this.map_ref_to_pos[ref_idx][ref_pos];
	 */
	ref_to_pos(ref_idx, ref_pos) {
		return this.map_ref_to_pos[ref_idx][ref_pos];
	}
}
let g_chrCoord = new ChromosomeCoordinate();

let snp_ratio = null;

let max_gc_content_window = 3000;

class RepeatSegment {
	constructor() {
		/** TSETA pos */
		this.start = 0;

		/** TSETA pos */
		this.end = 0;

		this.q_start = 0;
		this.q_end = 0;

		this.identity = 0;

		this.alignment = 0;

		this.q_len = 0;

		this.s_len = 0;
	}
}
/**
 * @param {string} url
 * @param {Uint32Array|number[]} ref_pos_map
 * @param {number} chr_length
 */
async function _loadRepeatSegment(url, ref_pos_map, chr_length) {
	/** @type {string} */
	const text = await fetchData(url, "text");
	const rows = text.split("\n").map(line => {
		const row = line.trim().split("\t");
		if (row.length <= 11) {
			return;
		}
		// return [
		// 	// row[0], // q
		// 	// row[1], // s
		// 	// row[2], // identity
		// 	// row[3], // alignment length
		// 	// row[4], // mismatch
		// 	// row[5], // gap
		// 	row[6], // q.start
		// 	row[7], // q.end
		// 	// row[8], // s.start
		// 	// row[9], // s.end
		// 	// row[10], // evalue
		// 	// row[11], // score
		// ];

		const q_start = Number(row[6].trim());
		const q_end = Number(row[7].trim());
		const s_start = Number(row[8].trim());
		const s_end = Number(row[9].trim());

		const q_len = Math.abs(q_end - q_start) + 1;
		const s_len = Math.abs(s_end - s_start) + 1;

		const start = ref_pos_map[q_start - 1];
		const end = ref_pos_map[q_end - 1];
		const identity = Number(row[2].trim());
		const alignment = Number(row[3].trim());// alignment (end - start + 1) >= 100

		// if (identity >= 65) {
			if (alignment < chr_length) {
				return Object.assign(new RepeatSegment(), {
					start: start,
					end: end,
					q_start: q_start,
					q_end: q_end,
					identity: identity,
					alignment: alignment,
					q_len: q_len,
					s_len: s_len,
				});
			}
		// }
	}).filter(a => a);
	return rows;
}
/**
 * @param {string} ref
 * @param {number} nChr
 * @param {1|2} in_or_out
 * @param {Uint32Array|number[]} ref_pos_map
 */
async function LoadRepeatSegment(ref, nChr, in_or_out, ref_pos_map) {
	// appendMarkersTable(ref, "");

	// const url = `repeat_segment/${ref}_ch${nChr}_repeat_${in_or_out}.txt`;
	// const url = `repeat_segment/${ref}_ch${nChr}_repeat_i65_${["inout", "in", "out"][in_or_out]}.txt`;
	const url = `repeat_segment/${ref}_ch${nChr}_repeat_blastn.txt`;
	const chr_length = analysis_options.chrInfo_list[dataset.parental_list.indexOf(ref)].length;

	return await _loadRepeatSegment(url, ref_pos_map, chr_length);
}
// function appendMarkersTable(ref, title) {
// 	const el_markers_table = document.getElementById("markers_table");
// 	const el_div = document.createElement("div");
// 	el_markers_table.before(el_div);
//
// 	const innerText = `${ref} ${title}`;
// 	// el_div.outerHTML = `<div style="line-height: ${seg_row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;
// 	el_div.outerHTML = `<div class="markers"><span contenteditable="true" style="vertical-align: middle;">${innerText}</span></div>`;
// }
/**
 * @type {{ [refId:string]: RepeatSegment[][] }}
 */
let repeat_segment = {};

/** @type {Plot_Helper} */
let plot_Helper = null;

/**
 * @type {{ [refId:string]: { [sampleId:string]: Uint32Array } }}
 */
let w_rnaCoverage_ui32a_map = {};

/**
 * small size => chr.size / rna_window_size
 * @type {{ [refId:string]: { [sampleId:string]: Uint32Array } }}
 */
let m_rnaCoverage_ui32a_map = {};

/**
 * rnaCoverage_max_map[refId][sampleName];
 * @type {{ [refId:string]: { [sampleId:string]: number } }}
 */
let rnaCoverage_max_map = {};

/**
 * rnaCoverage_norm_map[refId][sampleName];
 * @type {{ [refId:string]: { [sampleId:string]: number } }}
 */
let rnaCoverage_norm_map = {};

/**
 * min/mean/max size
 * @type {number}
 */
let rna_window_size = 100;

/**
 * min/mean/max size
 * @type {"min"|"mean"|"max"}
 */
let rna_func_type = "max";

/**
 * @see {@link window.rna_reads_func}
 * @see {@link viewerState.rna_reads_max_display_value}
 * @see {@link viewerState.rna_reads_low_cutoff_display_value}
 * @see {@link viewerState.rna_reads_high_cutoff_display_value}
 */
let rna_display_log_func = Math.log10;// window.rna_reads_func
let rna_display_log_base = 10;
let rna_display_exp = 4;// rna_display_log_base ** rna_display_exp >= maximum value
let rna_display_middleValue_list = [
	1,
	2 / 3,
	// 0.5,
	1 / 3,
	0,
];

let need_merge_rna_reads = true;// init m_rnaCoverage_ui32a_map

const methyl_ref_id_map = {// display methyl name
	// "QM6a": "QM6a",
	// "CBS1-1": "CBS1-1",
	[dataset.parental_list[0]]: dataset.parental_list[0],
	[dataset.parental_list[1]]: dataset.parental_list[1],
};
// const methyl_ref_idx_map = {
// 	// "QM6a": 0,
// 	// "CBS1-1": 1,
// 	[dataset.parental_list[0]]: 0,
// 	[dataset.parental_list[1]]: 1,
// };
// const methyl_ref_idx_map = Object.fromEntries([
// 	...dataset.parental_list,
// 	...dataset.progeny_list,
// ].map((v, i) => {
// 	return [v, i];
// }));
/**
 * @param {string} genomeName
 */
function get_genome_index(genomeName) {
	return [
		...dataset.parental_list,
		...dataset.progeny_list,
	].findIndex(a => a == genomeName);
}
const sample_type_list = [
	"veg",
	"D4",
	"D8",
];
const MethyContextTypeList = [
	"CG", "CHG", "CHH",
];
/**
 * @typedef {{ pos: number; value: number }[]} MethyDiffList
 * @type {{ [refId:string]: { [sampleId:string]: { CG: MethyDiffList, CHG: MethyDiffList, CHH: MethyDiffList } } }}
 */
let methy_diff_map = null;
let methy_BS_nBS_diff_map = null;
let a_methy_nBS_diff_map = null;
// const methy_info = {
// 	"QM6a": {
// 		"rid-QM6a": "Q_qr.txt",
// 		"WT-D4": "Q_d4.txt",
// 		"rid-D4": "Q_d4r.txt",
// 		"WT-D8": "Q_d8.txt",
// 		"rid-D8": "Q_d8r.txt",
// 		"CBS": "Q_c.txt",
// 		"rid-CBS": "Q_cr.txt",
// 	},
// 	"CBS1-1": {
// 		"WT-QM6a": "C_q.txt",
// 		"rid-QM6a": "C_qr.txt",
// 		"WT-D4": "C_d4.txt",
// 		"rid-D4": "C_d4r.txt",
// 		"WT-D8": "C_d8.txt",
// 		"rid-D8": "C_d8r.txt",
// 		"rid-CBS": "C_cr.txt",
// 	},
// };

// BS_oBS-seq20201219
// M4393	QM6a-1	No Oxi	Do BS_NUGEN
// M4394	QM6a-2	Do Oxi	Do BS_NUGEN
// M4395	QM6a-3	No Oxi	Do BS_Zymo
// M4396	QM6a-4	Do Oxi	Do BS_Zymo

let methy_sample_mapto_html = {
	// "veg": "veg",
	"$Q_bs": "Q bs",// 20200806
	// "$Q_bs1": "Q bs 1",// BS_oBS-seq20201219
	// "$Q_bs2": "Q bs 2 o",// BS_oBS-seq20201219
	"$Q_bs3": "Q bs 3",// BS_oBS-seq20201219
	"$Q_bs4": "Q bs 4 o",// BS_oBS-seq20201219

	"rid_veg": "veg <i>rid1</i>△",
	"dim2_veg": "veg <i>dim2</i>△",// 20200806 // Q / bs_dim2

	"D4": "D4",
	"rid_D4": "D4 <i>rid1</i>△",

	"D8": "D8",
	"rid_D8": "D8 <i>rid1</i>△",
};

// filename = {strain}_{name}.{type}.{format}
// strain = Q | C
// type = methratio
// format = float32
let methy_sample_fileName_map = {
	// [methy_sample_mapto_html.veg]: "wt",
	[methy_sample_mapto_html.$Q_bs]: "bs",// 20200806 // QM6a_bs.methratio.float32
	// [methy_sample_mapto_html.$Q_bs1]: "bs-1",// QM6a_bs-1.methratio.float32
	// [methy_sample_mapto_html.$Q_bs2]: "bs-2",// QM6a_bs-2.methratio.float32
	[methy_sample_mapto_html.$Q_bs3]: "bs-3",// QM6a_bs-3.methratio.float32
	[methy_sample_mapto_html.$Q_bs4]: "bs-4",// QM6a_bs-4.methratio.float32

	[methy_sample_mapto_html.rid_veg]: "rid",
	[methy_sample_mapto_html.dim2_veg]: "bs_dim2",// 20200806 // QM6a_bs_dim2.methratio.float32

	[methy_sample_mapto_html.D4]: "d4",
	[methy_sample_mapto_html.rid_D4]: "d4r",

	[methy_sample_mapto_html.D8]: "d8",
	[methy_sample_mapto_html.rid_D8]: "d8r",

	// "rid(veg)_CG": "rid_CG",//Qrid CG
	// "Qdim2 NP30": "NP30",//NP
	// "Qrid NP32": "NP32",//NP
};
const disabled_methyl = "disabled";
let methy_refMap_sampleList = {//see line:4270
	"QM6a": [
		// methy_sample_mapto_html.veg,//QM6a
		methy_sample_mapto_html.$Q_bs,// 20200806
		// methy_sample_mapto_html.$Q_bs1,// QM6a_bs-1.methratio.float32
		// methy_sample_mapto_html.$Q_bs2,// QM6a_bs-2.methratio.float32
		methy_sample_mapto_html.$Q_bs3,// QM6a_bs-3.methratio.float32
		methy_sample_mapto_html.$Q_bs4, // QM6a_bs-4.methratio.float32

		disabled_methyl ?? methy_sample_mapto_html.rid_veg,//QM6a
		disabled_methyl ?? methy_sample_mapto_html.dim2_veg,// 20200806

		disabled_methyl ?? methy_sample_mapto_html.D4,
		disabled_methyl ?? methy_sample_mapto_html.rid_D4,
		disabled_methyl ?? methy_sample_mapto_html.D8,
		disabled_methyl ?? methy_sample_mapto_html.rid_D8,

		// "rid(veg)_CG",//Qrid CG
		// "Qdim2 NP30",//NP
		// "Qrid NP32",//NP
	],
	"CBS1-1": [
		disabled_methyl ?? methy_sample_mapto_html.veg,//CBS1-1
		disabled_methyl ?? methy_sample_mapto_html.rid_veg,//CBS1-1

		disabled_methyl ?? methy_sample_mapto_html.D4,
		disabled_methyl ?? methy_sample_mapto_html.rid_D4,
		disabled_methyl ?? methy_sample_mapto_html.D8,
		disabled_methyl ?? methy_sample_mapto_html.rid_D8,
	],
};
// const methy_refMap_sampleList = {
// 	"QM6a": [
// 		"WT(veg)",//QM6a
// 		// "Q bs",// 20200806
// 		"rid(veg)",//QM6a
// 		"Q bs dim2",// 20200806

// 		"WT(D4)",
// 		"rid(D4)",
// 		"WT(D8)",
// 		"rid(D8)",

// 		// "rid(veg)_CG",//Qrid CG
// 		// "Qdim2 NP30",//NP
// 		// "Qrid NP32",//NP
// 	],
// 	"CBS1-1": [
// 		"WT(veg)",//CBS1-1
// 		"rid(veg)",//CBS1-1

// 		"WT(D4)",
// 		"rid(D4)",
// 		"WT(D8)",
// 		"rid(D8)",
// 	],
// };
const display_methy_refMap_sampleList = {//see line:4270
	"veg": {
		"QM6a": [
			{ title: methy_sample_mapto_html.veg, display: false, },//QM6a
			{ title: methy_sample_mapto_html.$Q_bs, display: true, },// 20200806
			{ title: methy_sample_mapto_html.$Q_bs1, display: true, },// QM6a_bs-1.methratio.float32
			{ title: methy_sample_mapto_html.$Q_bs2, display: true, },// QM6a_bs-2.methratio.float32
			{ title: methy_sample_mapto_html.$Q_bs3, display: true, },// QM6a_bs-3.methratio.float32
			{ title: methy_sample_mapto_html.$Q_bs4, display: true, }, // QM6a_bs-4.methratio.float32
			{ title: methy_sample_mapto_html.rid_veg, display: false, },//QM6a
			{ title: methy_sample_mapto_html.dim2_veg, display: false, },// 20200806
		],
		"CBS1-1": [
			{ title: methy_sample_mapto_html.veg, display: false, },//CBS1-1
			{ title: methy_sample_mapto_html.rid_veg, display: false, },//CBS1-1
		],
	},
	"D4": {
		"QM6a": [
			{ title: methy_sample_mapto_html.D4, display: false, },
			{ title: methy_sample_mapto_html.rid_D4, display: false, },
		],
		"CBS1-1": [
			{ title: methy_sample_mapto_html.D4, display: false, },
			{ title: methy_sample_mapto_html.rid_D4, display: false, },
		],
	},
	"D8": {
		"QM6a": [
			{ title: methy_sample_mapto_html.D8, display: false, },
			{ title: methy_sample_mapto_html.rid_D8, display: false, },
		],
		"CBS1-1": [
			{ title: methy_sample_mapto_html.D8, display: false, },
			{ title: methy_sample_mapto_html.rid_D8, display: false, },
		],
	}
};
// const methy_refMap_rowNameList = {
// 	"QM6a": {
// 		ratio: [
// 			"WT(veg)",
// 			"rid(veg)",
// 			"WT(D4)",
// 			"rid(D4)",
// 			"WT(D8)",
// 			"rid(D8)",
// 		],
// 		diff: [
// 			"rid(veg)-WT(veg)",
// 			"WT(D4)-WT(veg)",
// 			"rid(D4)-WT(veg)",
// 			"WT(D8)-WT(veg)",
// 			"rid(D8)-WT(veg)",
// 		],
// 	},
// 	"CBS1-1": {
// 		ratio: [
// 			"WT-CBS",
// 			"rid-CBS",
// 			"WT-D4",
// 			"rid-D4",
// 			"WT-D8",
// 			"rid-D8",
// 		],
// 		diff: [
// 			"rid-CBS",
// 			"WT-D4",
// 			"rid-D4",
// 			"WT-D8",
// 			"rid-D8",
// 		],
// 	},
// };
/**
 * @type {{ [refId:string]: { [sampleId:string]: { pos: number; value: number }[] } }}
 */
let methy_ratio_map = null;
/**
 * @type {{ pos: number; value: number }[]}
 */
let NP32_NP30_diff = null;//new Float32Array(methy_ratio_map["QM6a"]["Qrid NP32"].length)

/**
 * @type {{ [refId:string]: { [sampleId:string]: { pos: number; value: number }[] } }}
 */
let nBS_methy_ratio_map = null;

/**
 * @type {{ [refId:string]: { [sampleId:string]: { pos: number; value: number }[] } }}
 */
let BS_nBS_match_methy_map = null;

// /** @type {boolean[]} */
// const display_sample = new Array(Object.values(methy_refMap_sampleList.Q).length);


class module_Methyl_sampleData {
	/** @type {""|"QM6a"|"CBS1-1"} */
	ref = "";

	sample = "";

	/**
	 * @type {null|""|"BS"|"oxBS"|"ctrl"}
	 * null <- ""
	 */
	seq_type = "";

	/** display name */
	name = "";

	/** @description data description */
	description = "";

	/** @type {string[]} user data */
	tags = [];
	
	static SN = 0;
	// order = module_Methyl_sampleData.SN++;
	// get order() {
	// 	throw Error("No use");
	// }
	order = 0;

	value_desc = "";

	/** @type {string} */
	html_value_desc = null;

	/** @type {HTMLCanvasElement} */
	canvas_value_desc = null;

	/** @type {string} file path */
	url = "";

	/** is region or loci */
	region = false;

	/** min display width (pixel) */
	strong = 0;// strong

	// /**
	//  * @param {number} value
	//  * @returns {boolean}
	//  */
	// value_filter = function (value) {
	// 	return true;
	// }

	/** @type {number} */
	fontSize_scale = null;

	/**
	 * @default 15
	 */
	row_height = 1.5;

	/** @type {string} */
	border_color = null;

	/** @type {boolean} */
	mid_line = true;

	/** @type {boolean} */
	density_to_opacity = true;

	/**
	 * @type {(view_length: number) => number}
	 * @returns {number} [0, 1]
	 */
	func_density_to_opacity = null;

	range_repeat_density = 0;

	/** @type {number} */
	max_display_value = 0;

	value_normalize = null;

	split_strand = false;

	__hide = false;
	get hide() {
		return this.__hide;
	}
	set hide(val) {
		this.__hide = val;
		if (val == false) {
			console.warn(val ? "hide" : "show", new Error(`${this.name}: ${this.url}`), this);
		}
		else {
			console.info(val ? "hide" : "show", new Error(`${this.name}: ${this.url}`), this);
		}
	}

	/** merge point mode */
	fast = false;
	
	/** merge point pre bin */
	fast_binSize = 10;

	/**
	 * { "#E00000", v => v > 0 },
	 * { "#0000E0", v => v < 0 },
	 * { "#000000", v => v == 0 },
	 * @type {module_Methyl_sampleData_RenderingCondition[]}
	 */
	rendering_condition = [];

	display_minus_value = false;//-1 ~ +1 if true; 0 ~ 1 if false
	remove_on_reload = false;

	// /**
	//  * strokeStyle: string | CanvasGradient | CanvasPattern
	//  * @type {{ value: number; strokeStyle: string | CanvasGradient | CanvasPattern; }[]}
	//  */
	// additional_horizontal_line = [];

	/** @type {((self: module_Methyl_sampleData, row_width: number, row_height: number) => void)[]} */
	beforeRenderList = [];

	/** @type {((self: module_Methyl_sampleData, row_width: number, row_height: number) => void)[]} */
	afterRenderList = [];

	getValueTransformer() {
		if (this.max_display_value > 0) {
			return v => {
				const sign = Math.sign(v);
				return sign * Math.min(this.max_display_value, Math.abs(v)) / this.max_display_value;
			};
		}
		else if (this.value_normalize > 0) {
			const value_normalize = this.value_normalize;//Math.max(1, Math.abs(this.value_normalize));
			return v => v / value_normalize;
		}
		else {
			return v => v;
		}
	}

	/**
	 * @example
	 * if (ref_to_pos_map == null) { global_pos = data.pos }
	 * else { global_pos = ref_to_pos_map[local_pos] }
	 * @type {Uint32Array|number[]}
	 */
	ref_to_pos_map = null;

	/** @type {Partial<module_Methyl_ratioData>[][]} */
	data = [];

	/**
	 * render time elapsed (millisecond)
	 */
	_elapsed_ms = 0;

	/**
	 * @param {Partial<module_Methyl_sampleData>=} clone_source
	 */
	constructor(clone_source) {
		if (clone_source) {
			Object.assign(this, clone_source);
		};

		this.__initiator = new Error();
	}

	/**
	 * @param {string} innerHTML
	 */
	async apply_valueDesc_HTML_to_canvas(innerHTML) {
		this.html_value_desc = innerHTML;
		return await this._make_value_desc_canvas();
	}

	/**
	 * @description html2canvas only one instance in global
	 */
	async _make_value_desc_canvas() {
		const el_con_html_value_desc = (function () {
			let el = document.getElementById("el_con_html_value_desc");
			if (!el) {
				el = document.createElement("div");
				el.id = "el_con_html_value_desc";
				document.body.append(el);
			}
			return el;
		})();

		const hidden = document.createElement("div");
		el_con_html_value_desc.append(hidden);

		hidden.dataset.id = this.name;
		hidden.innerHTML = this.html_value_desc;

		await delayFrame(); //?? html2canvas

		const el = hidden.querySelector("*");
		try {
			/** @type {HTMLCanvasElement} */
			const canvas = await html2canvas_async(el); // html2canvas use Promise polyfill
			await delayFrame();
			this.canvas_value_desc = canvas;

			return canvas;
		}
		catch (ex) {
			console.error(ex);
			debugger;
		}
		// finally {
		// 	el_con_html_value_desc.remove();
		// }
	}
	
	/**
	 * @type {undefined | null | (x1: number, x2: number, min_width: number, row_height: number, value: number, display_minus_value: boolean) => void}
	 */
	drawShape = null;

	// value_filter(value) {
	// 	throw new Error("Method not implemented.");
	// }

	// /**
	//  * @param {module_MethylData} data
	//  */
	// static from(data) {
	// 	const ret = new module_MethylData();
	// 	ret.ref = data.ref;
	// 	// ret.strand = data.strand;
	// 	ret.name = data.name;
	// 	ret.pos_map = data.pos_map;
	// 	ret.file = data.file;
	// 	ret.region = data.region;
	// 	return ret;
	// }
}

class module_Methyl_sampleData_RenderingCondition {
	color = "black";

	/** @type {(point: module_Methyl_ratioData) => boolean} */
	condition = null;

	/**
	 * low density: 1
	 * high density: 1/4
	 */
	min_width = 1;

	/**
	 * @param {Partial<module_Methyl_sampleData_RenderingCondition>} assign_source
	 */
	constructor(assign_source) {
		if (assign_source != null) {
			Object.assign(this, assign_source);
		}
	}
}

class module_Methyl_ratioData {
	// nChr = 0;

	start = 0;
	end = 0;

	/** ratio */
	value = 0;

	strand = 0;

	// /**
	//  * @param {Partial<module_MethRatioData>} data
	//  */
	// constructor(data) {
	// 	if (data) {
	// 		this.nChr = Number(data.nChr);
	// 		this.start = Number(data.start);
	// 		this.end = Number(data.end);
	// 		this.value = Number(data.value);
	// 	}
	// }
}

// 20221003: red line
// average to red line
class module_MethylRenderer_averageToLine {
	// /** @type {module_Methyl_sampleData} */
	// target = null;

	/**
	 * @param {module_Methyl_sampleData} target
	 */
	constructor(target) {
		this._target = target;
		this._average = 0;

		this.enable = false;
		/** @type {string | CanvasGradient | CanvasPattern} */
		this.strokeStyle = "red";
	}

	get target() { return this._target; };
	get average() { return this._average; };
	get value() { return this._average; };

	/**
	 * @param {module_Methyl_sampleData} config
	 * @param {number} row_width
	 * @param {number} row_height
	 */
	render(config, row_width, row_height) {
		const y = row_height - config.getValueTransformer()(this.value) * row_height;
		main_ctx.beginPath();
		main_ctx.moveTo(0, y);
		main_ctx.lineTo(row_width, y);
		main_ctx.strokeStyle = this.strokeStyle;
		main_ctx.stroke();

		// registerBeforeRender(...)
	}

	/**
	 * @param {number} [start_pos]
	 * @param {number} [end_pos]
	 * @see {draw_methy_ratio_row}
	 */
	_updateRange(start_pos, end_pos) {
		// g_chrCoord: ChromosomeCoordinate
		// g_methylRenderer._real_pos // ref_pos
		// methyl_dataset_list[1].data[5].reduce((acc, v) => acc + v.value, 0) / methyl_dataset_list[1].data[5].length;
		
		const filtered_chr_data = (() => {
			if (start_pos && end_pos && start_pos != end_pos) {
				const target_genomeIdx = dataset.genomeNameList.indexOf(this.target.ref);
				const target_start_pos = g_chrCoord.pos_to_ref(target_genomeIdx, start_pos);
				const target_end_pos = g_chrCoord.pos_to_ref(target_genomeIdx, end_pos)

				return this.target.data[viewerState.nChr - 1].filter(a => {
					return Number(a.end) >= target_start_pos && Number(a.start) <= target_end_pos;
				});
			}
			else {
				return this.target.data[viewerState.nChr - 1];
			}
		})();

		this._average = filtered_chr_data.reduce((acc, v) => acc + Number(v.value), 0) / filtered_chr_data.length;
		
		return this._average;
	}
	
	/**
	 * @param {module_Methyl_sampleData[]} sampleData_list
	 * @param {number} [start_pos]
	 * @param {number} [end_pos]
	 */
	updateAvgLine(sampleData_list, start_pos, end_pos) {
		if (this.enable && this.target) {
			this._updateRange(start_pos, end_pos);
			
			const render = (cfg, row_width, row_height) => this.render(cfg, row_width, row_height);
			sampleData_list.forEach(cfg => {
				cfg.afterRenderList.push(render);
			});
		}
	}
}

class module_MethylRenderer {
	debug = false;

	// /**
	//  * @type {Promise<void>}
	//  */
	// promise = null;

	/**
	 * @param {module_Methyl_sampleData[]} methyl_dataset_list
	 */
	constructor(methyl_dataset_list) {
		this.methyl_dataset_list = methyl_dataset_list;
		this.display = true;
		this.sort_layout = false;
	}

	/**
	 * @param {HTMLElement} el_appendRow
	 * @param {number} row_width
	 * @param {number} seg_row_height
	 * @param {number} seg_row_separate
	 */
	async draw_all_methyl_row(el_appendRow, row_width, seg_row_height, seg_row_separate) {
		if (!this.display) {
			return;
		}

		const ref_pos_map_list = make_ref_pos_map_list();

		const list = this.sort_layout == true ? 
			methyl_dataset_list.slice(0).sort((a, b) => a.order - b.order)
			:
			methyl_dataset_list;

		if (this.sort_layout == true && list.some(a => a.order == null)) {
			console.warn("undefined order");
		}

		// 20221003: red line
		// average to red line
		// redline { avg_target, avg_val }

		// methyl_dataset_list.reduce(async (prev_promise, config, config_idx) => {
		await list.reduce(async (prev_promise, config, config_idx) => {
			await prev_promise;

			const chr_idx = viewerState.nChr - 1;
			if (!config.data.length || !config.data[chr_idx] || config.hide) {
				return;
			}

			const ratio_rows = config.data[chr_idx];
			const row_height = Math.trunc(seg_row_height * (config.row_height ?? 1.5));

			const innerText = `${config.name}`;
			// const el_div = document.createElement("div");
			// el_appendRow.append(el_div);
			// const css_fontSize = config.fontSize_scale > 0 ? `font-size: ${config.fontSize_scale}em; ` : "";
			// /**
			//  * margin-bottom: ${seg_row_separate - 1}px;
			//  * if element has no content then line-height: 1px;
			//  * 0.5 is Math.round result
			//  */
			// el_div.outerHTML = `
			// 	<div title="methy ratio" style="${css_fontSize}line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;">
			// 		<span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span>
			// 	</div>`;
			const el = appendLeftTitle(innerText, row_height, seg_row_separate, 1, (config.fontSize_scale ?? 1) * 100);
			el.dataset.__initiator = config.__initiator.stack;
			el.title = `${config_idx}. ${config.sample}`;

			// always find map
			config.ref_to_pos_map = ref_pos_map_list[dataset.genome_info_list.findIndex(a => a.name == config.ref)];

			await draw_methy_ratio_row(config, ratio_rows, row_width, row_height, "#0123AB", $color_set_view.dad_bk);

			main_ctx.translate(0, 1 * (row_height + seg_row_separate));
		}, Promise.resolve());
	}

	/**
	 * @param {module_Methyl_sampleData[]} md_list
	 */
	make_methyl_Chm_row(md_list = methyl_dataset_list) {
		if (md_list.length == 0) {
			return;
		}
		md_list.push(make_methyl_Chm_row(
			md_list[8], // Q2+
			md_list[10], // Q2o+
			"QM6a-2 " + html_Chm + " +",
			make_methyl_value_desc("Chm", "Cm + Chm + C")
		));
		md_list.push(make_methyl_Chm_row(
			md_list[9], // Q2-
			md_list[11], // Q2o-
			"QM6a-2 " + html_Chm + " -",
			make_methyl_value_desc("Chm", "Cm + Chm + C")
		));

		md_list.push(make_methyl_C_row(
			md_list[8],
			"QM6a-2 C +",
			make_methyl_value_desc("C", "Cm + Chm + C")
		));
		md_list.push(make_methyl_C_row(
			md_list[9],
			"QM6a-2 C -",
			make_methyl_value_desc("C", "Cm + Chm + C")
		));

		function make_methyl_Chm_row(dataset_1, dataset_2, name, html_value_desc) {
			const result = new module_Methyl_sampleData({
				ref: "QM6a",
				sample: "NP43",
				name: name,
				split_strand: true,
				html_value_desc: html_value_desc,
				url: "",
				region: false,
				display_minus_value: true,
				remove_on_reload: true,
			});

			for (let chr_idx = 0; chr_idx < dataset.genome_info_list[0].chr_list.length; ++chr_idx) {
				/** @type {module_Methyl_ratioData[]} */
				const pm = [];

				dataset_1.data[chr_idx].forEach(a => {
					pm[a.start] = { ...a };// clone
				});

				dataset_2.data[chr_idx].forEach(a => {
					if (pm[a.start]) {
						pm[a.start].value = pm[a.start].value - a.value;
					}
					else {
						pm[a.start] = { ...a };// clone
						pm[a.start].value = 0 - a.value;
					}
				});

				result.data[chr_idx] = pm.filter(a => a);
			}

			return result;
		}

		function make_methyl_C_row(dataset_1, name, html_value_desc) {
			const result = new module_Methyl_sampleData({
				ref: "QM6a",
				sample: "NP43",
				name: name,
				split_strand: true,
				html_value_desc: html_value_desc,
				url: "",
				region: false,
				// display_minus_value: false,
				remove_on_reload: true,
			});

			for (let chr_idx = 0; chr_idx < dataset.genome_info_list[0].chr_list.length; ++chr_idx) {
				result.data[chr_idx] = dataset_1.data[chr_idx].map(a => {
					return {
						...a,
						value: 1 - Math.abs(a.value),
					};
				});
			}

			return result;
		}
	}

	/**
	 * @param {module_Methyl_sampleData[]} conifg_list
	 */
	split_strand(conifg_list) {
		if (conifg_list.length == 0) {
			return;
		}
		conifg_list.forEach(config => {
			if (!config.split_strand) {
				return;
			}
			console.log(new Error(20221111).stack, config);
			// alert(new Error(20221111));

			{
				const nc_plus = new module_Methyl_sampleData({
					...config,// clone
					name: config.name + " " + "+",
					split_strand: false,
					remove_on_reload: true,
					data: config.data.map(ds => {
						return ds.filter(a => a.strand > 0);
					}),
				});
				conifg_list.push(nc_plus);
			}

			{
				const nc_minus = new module_Methyl_sampleData({
					...config,// clone
					name: config.name + " " + "-",
					split_strand: false,
					remove_on_reload: true,
					data: config.data.map(ds => {
						return ds.filter(a => a.strand < 0);
					}),
				});
				conifg_list.push(nc_minus);
			}

			config.split_strand = false;
			config.hide = true;
		});
	}

	// http://localhost:8080/tseta/20200720_v3_QCt/debug_20200720_v3_QCt.html
	/**
	 * @param {module_Methyl_sampleData[]} conifg_list
	 * @param {boolean} [is_reload]
	 */
	async __html_to_canvas(conifg_list = methyl_dataset_list, is_reload = false) {
		// load data
		for (let config of conifg_list) {
			if ((is_reload || !config.canvas_value_desc) && config.html_value_desc) {
				await config._make_value_desc_canvas();// html2canvas only one instance in global
			}
		}
	}

	/**
	 * @param {module_Methyl_sampleData[]} conifg_list
	 * @param {boolean} [is_reload]
	 */
	async load(conifg_list = methyl_dataset_list, is_reload = false) {
		// load data
		const tasks = conifg_list.map(async (config, config_idx) => {
			if (config.data.length && !is_reload) {
				console.warn("skip load:", config_idx, config.name, config.data.length, new Error("skip load").stack);
				return;
			}

			try {
				if (config.url.endsWith(".tsv")) {
					config.data = await this.load_np_meth_tsv(config);
				}
				else if (config.url.endsWith(".float32")) {
					config.data = await this.load_meth_ratio_float32(config);
				}
				// else if (config.url.endsWith(".bed")) {
				// 	config.data = await this._load_tsv(config);
				// }
				else {
					throw new Error("file fmt");
				}
			}
			catch (ex) {
				console.error(config_idx, config.name);
				console.error(ex);
			}
		});

		tasks.push(this.__html_to_canvas(conifg_list, is_reload));

		await Promise.all(tasks);
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<Float32Array>}
	 * @see {@link draw_methy_ratio_row}
	 */
	async load_float32Array(config) {
		if (this.debug) {
			return new Float32Array(0);
		}

		if (config instanceof module_Methyl_sampleData) {
			const a = (await this.load_TypedArray(config, Float32Array));
			// return a instanceof Float32Array ? a : undefined;
			return a;
		}
		else {
			throw new TypeError();
		}
	}

	/**
	 * use NaN padding gap
	 * read Uint32Array to Float32Array
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<Float32Array>}
	 * @see {@link draw_methy_ratio_row}
	 */
	async load_uint32Array(config) {
		if (this.debug) {
			return new Float32Array(0);
		}

		if (config instanceof module_Methyl_sampleData) {
			const a = (await this.load_TypedArray(config, Uint32Array));
			return a || new Float32Array(0);
			// return a instanceof Uint32Array ? a : undefined;
		}
		else {
			throw new TypeError();
		}
	}

	/** @type {{ start: number; end: number; }[]} */
	_real_pos = [];

	isLimitLoad() {
		return this._real_pos != null && this._real_pos.length == dataset.genome_info_list.length;
	}

	/**
	 * 20220912 set_bounding
	 * @param {number} ma_start
	 * @param {number} ma_end
	 */
	setLoadingBounding(ma_start = g_chrCoord.bp_start, ma_end = g_chrCoord.bp_end) {
		// copy_visable_seq
		// const ma_seq_start_idx = ma_start - 1;
		// const ma_seq_end_idx = ma_end - 1;

		for (let ref_idx = 0; ref_idx < dataset.genome_info_list.length; ++ ref_idx) {
			// const pos_to_ref = multialign_to_chrPos_posMap(seq_list[ref_idx]);
			// const ref_to_pos = chrPos_to_multialign_posMap(seq_list[ref_idx]);

			const ref_start = g_chrCoord.pos_to_ref(ref_idx, ma_start);
			const ref_end = g_chrCoord.pos_to_ref(ref_idx, ma_end);

			this._real_pos[ref_idx] = {
				start: ref_start,
				end: ref_end,
			};
			// this._real_pos[ref_idx].start = ref_start;
			// this._real_pos[ref_idx].end = ref_end;
		}

		const _ref1_len = pos_ref1_uint32array[g_chrCoord.bp_end - 1] - pos_ref1_uint32array[g_chrCoord.bp_start - 1] + 1;
		const ref1_len = this._real_pos[0].end - this._real_pos[0].start + 1;

		if (ref1_len != _ref1_len) {
			console.warn({
				visible_length: _ref1_len,
				loading_length: ref1_len,
			});
		}

		return ref1_len == _ref1_len;
	}

	/**
	 * @template Ctor
	 * @param {{ sample: string; url: string; ref: string; }|module_Methyl_sampleData} config
	 * @param {Ctor extends Float32ArrayConstructor ? Float32ArrayConstructor : (Ctor extends Uint32ArrayConstructor ? Uint32ArrayConstructor : never} typedArray
	 * @see {@link draw_methy_ratio_row}
	 * @see {module_MethylRenderer.load_meth_ratio_float32}
	 */
	async load_TypedArray(config, typedArray) {
		if (this._real_pos.length != dataset.genome_info_list.length) {
			return await this._load_TypedArray_full(config, typedArray);
		}
		else {
			return await this._load_TypedArray_visible(config, typedArray);
		}
	}

	/**
	 * @template Ctor
	 * @param {{ sample: string; url: string; ref: string; }|module_Methyl_sampleData} config
	 * @param {Ctor extends Float32ArrayConstructor ? Float32ArrayConstructor : (Ctor extends Uint32ArrayConstructor ? Uint32ArrayConstructor : never} typedArray
	 * @see {@link draw_methy_ratio_row}
	 */
	async _load_TypedArray_full(config, typedArray) {
		// throw new Error("ref pos to tseta");
		console.info(timeElapsed(), "loading: methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		
		const ref = config.ref;
		if (ref == null) {
			throw new TypeError("required ref");
		}

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
		/** @type {number[]} */
		const chrLengthList = dataset.genome_info_list[ref_idx].chr_list.map(a => a.length);

		const chrOffsetInByte = [];
		chrOffsetInByte[0] = 0;
		const whole_genome_size = chrLengthList.reduce((pos, len, i) => {
			let offset = pos + (len * 4);
			chrOffsetInByte[i + 1] = offset;
			return offset;
		}, 0);

		const chr_idx = analysis_options.nChr - 1;
		const offsetStart = chrOffsetInByte[chr_idx];
		const chr_len = chrLengthList[chr_idx];
		// const f32a = new typedArray(buffer, offsetStart, chr_len);
		
		// "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_BS.methratio.float32"

		if (this.debug) {
			return new typedArray(0);
		}
		else {
			// get file full size
			const hd = (await fetch(config.url, {
				"method": "HEAD",
			}));
			const got_whole_genome_size = Number(hd.headers.get("content-length"));

			if (got_whole_genome_size != whole_genome_size) {
				debugger;
				throw new Error(`No match size, ${config.url}, got ${got_whole_genome_size}, expected ${whole_genome_size}`);
			}

			/**
			 * 206: Partial Content
			 * @type {ArrayBuffer}
			 */
			const buffer = await (await fetch(config.url, {
				headers: {
					"Range": `bytes=${offsetStart}-${offsetStart + (chr_len * typedArray.BYTES_PER_ELEMENT)}`,
				}
			})).arrayBuffer();
			if (buffer && buffer.byteLength == (chr_len * typedArray.BYTES_PER_ELEMENT + 1)) {
				document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;
				const f32a = new typedArray(buffer, 0, chr_len);
				
				return f32a;
			}
			else {
				return new typedArray(0);
			}
		}
	}

	/**
	 * 20220912 set_bounding
	 * @template Ctor
	 * @param {{ sample: string; url: string; ref: string; }|module_Methyl_sampleData} config
	 * @param {Ctor extends Float32ArrayConstructor ? Float32ArrayConstructor : (Ctor extends Uint32ArrayConstructor ? Uint32ArrayConstructor : never} typedArray
	 * @see {@link draw_methy_ratio_row}
	 */
	async _load_TypedArray_visible(config, typedArray) {
		// throw new Error("ref pos to tseta");
		console.info(timeElapsed(), "loading: methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);

		const ref = config.ref;
		if (ref == null) {
			throw new TypeError("required ref");
		}

		// if (this._real_pos.length != dataset.genome_info_list.length) {
		// 	this.setLoadingBounding(g_chrCoord.bp_start, g_chrCoord.bp_end);
		// }

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
		/** @type {number[]} */
		const chrLengthList = dataset.genome_info_list[ref_idx].chr_list.map(a => a.length);

		const chrOffsetInByte = [];
		chrOffsetInByte[0] = 0;
		const whole_genome_size = chrLengthList.reduce((pos, len, i) => {
			let offset = pos + (len * typedArray.BYTES_PER_ELEMENT);
			chrOffsetInByte[i + 1] = offset;
			return offset;
		}, 0);

		const chr_idx = analysis_options.nChr - 1;
		/** offset in bytes */
		const chrOffsetStart = chrOffsetInByte[chr_idx];
		/** length in elements */
		const chr_len = chrLengthList[chr_idx];
		// const f32a = new typedArray(buffer, offsetStart, chr_len);

		const local_start = this._real_pos[ref_idx].start - 1;
		const local_end = this._real_pos[ref_idx].end - 1;
		const local_len = Math.min(local_end - local_start + 1, chr_len);

		const local_offsetStart = chrOffsetInByte[chr_idx] + local_start * typedArray.BYTES_PER_ELEMENT;
		const local_lengthInBytes = local_len * typedArray.BYTES_PER_ELEMENT;
		const local_offsetEnd = local_offsetStart + local_lengthInBytes;
		
		const chrEndOffset = chrOffsetStart + (chr_len * typedArray.BYTES_PER_ELEMENT);
		if (local_offsetEnd > chrEndOffset) {
			debugger;
		}
		
		// "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_BS.methratio.float32"

		if (this.debug) {
			return new typedArray(0);
		}
		else {
			// get file full size
			const hd = (await fetch(config.url, {
				"method": "HEAD",
			}));
			const got_whole_genome_size = Number(hd.headers.get("content-length"));

			if (got_whole_genome_size != whole_genome_size) {
				debugger;
				throw new Error(`No match size, ${config.url}, got ${got_whole_genome_size}, expected ${whole_genome_size}`);
			}
			
			/** @type {ArrayBuffer} */
			const buffer = await (await fetch(config.url, {
				headers: {
					"Range": `bytes=${local_offsetStart}-${local_offsetEnd}`,
				}
			})).arrayBuffer();

			if (buffer && buffer.byteLength == (local_len * typedArray.BYTES_PER_ELEMENT + 1)) {
				document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;

				const local_f32a = new typedArray(buffer, 0, local_len);
				const f32a = new Float32Array(chr_len);// full length // use NaN and Infinity
				
				// if (config.sample.startsWith("ddr_D8")) {// ddr_D8 ctrl
				// 	const all_nan = local_f32a.every(a => Number.isNaN(a));
				// 	if (all_nan) {
				// 		debugger;//20220921 all NaN
				// 	}
				// }

				f32a.fill(NaN, 0, f32a.length);// fill NaN, ignore point
				if (!Number.isNaN(f32a[0])) { // NaN can't compare, no compare operator
					debugger;
				}
				f32a.set(local_f32a, local_start);

				return f32a;
			}
			else {
				debugger;
				return new typedArray(0);
			}
		}
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<module_Methyl_ratioData[][]>}
	 * @see {@link draw_methy_ratio_row}
	 * @see {module_MethylRenderer.load_TypedArray}
	 */
	async load_meth_ratio_float32(config) {
		if (this._real_pos.length != dataset.genome_info_list.length) {
			return await this._load_meth_ratio_float32_full(config);
		}
		else {
			return await this._load_meth_ratio_float32_visible(config);
		}
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<module_Methyl_ratioData[][]>}
	 * @see {@link draw_methy_ratio_row}
	 * @see {module_MethylRenderer.load_TypedArray}
	 */
	async _load_meth_ratio_float32_full(config) {
		// throw new Error("ref pos to tseta");
		console.info(timeElapsed(), "loading: methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		
		const ref = config.ref;
		if (ref == null) {
			throw new TypeError("required ref");
		}

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
		/** @type {number[]} */
		const chrLengthList = dataset.genome_info_list[ref_idx].chr_list.map(a => a.length);

		const chrOffsetInByte = [];
		chrOffsetInByte[0] = 0;
		chrLengthList.reduce((pos, len, i) => {
			let offset = pos + (len * 4);
			chrOffsetInByte[i + 1] = offset;
			return offset;
		}, 0);

		const chr_idx = analysis_options.nChr - 1;
		const offsetStart = chrOffsetInByte[chr_idx];
		const chr_len = chrLengthList[chr_idx];
		// const f32a = new Float32Array(buffer, offsetStart, chr_len);
		
		// "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_BS.methratio.float32"

		const data = chrLengthList.map(_ => []);

		if (this.debug) {
			// return new Float32Array(0);
			return data;
		}
		
		/** @type {ArrayBuffer} */
		const buffer = await (await fetch(config.url, {
			headers: {
				"Range": `bytes=${offsetStart}-${offsetStart + (chr_len * Float32Array.BYTES_PER_ELEMENT)}`,
			}
		})).arrayBuffer();
		if (buffer && buffer.byteLength == (chr_len * Float32Array.BYTES_PER_ELEMENT + 1)) {
			document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;
			
			if (buffer && buffer.byteLength != 0) {
				const f32a = new Float32Array(buffer, 0, chr_len);

				const fff = new Float32Array(1);
				const uuu = new Uint32Array(fff.buffer);
				
				/** @see {@link _load_meth_ratio_float32} */
				f32a.forEach((value, index) => {
					if (!Number.isNaN(value) && value != 0) {
						const pos = Number(index) + 1;// .methratio.txt, base position = 0
						const r = new module_Methyl_ratioData();
						r.start = pos;
						r.end = r.start;
						r.value = Math.abs(value);
						
						fff[0] = value;
						r.strand = (uuu[0] >>> 31) ? -1 : 1;

						data[chr_idx].push(r);
					}
				});
			}
		}
		// @ts-ignore
		console.info(timeElapsed(), "loading: loaded methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		document.getElementById("status").innerText = "loading: loaded methy ratio:" + config.sample;

		await delayFrame(); // unlock UI // globalThis.gc()
		return data;
	}

	/**
	 * 20220912 set_bounding
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<module_Methyl_ratioData[][]>}
	 * @see {@link draw_methy_ratio_row}
	 * @see {module_MethylRenderer.load_TypedArray}
	 */
	async _load_meth_ratio_float32_visible(config) {
		// throw new Error("ref pos to tseta");
		console.info(timeElapsed(), "loading: methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);

		const ref = config.ref;
		if (ref == null) {
			throw new TypeError("required ref");
		}

		// if (this._real_pos.length != dataset.genome_info_list.length) {
		// 	this.setLoadingBounding(g_chrCoord.bp_start, g_chrCoord.bp_end);
		// }

		if (config.sample == "QM6a WT_D8 ctrl") {
			debugger;
		}

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
		/** @type {number[]} */
		const chrLengthList = dataset.genome_info_list[ref_idx].chr_list.map(a => a.length);

		const chrOffsetInByte = [];
		chrOffsetInByte[0] = 0;
		chrLengthList.reduce((pos, len, i) => {
			let offset = pos + (len * Float32Array.BYTES_PER_ELEMENT);
			chrOffsetInByte[i + 1] = offset;
			return offset;
		}, 0);

		const chr_idx = analysis_options.nChr - 1;
		const chrOffsetStart = chrOffsetInByte[chr_idx];
		const chr_len = chrLengthList[chr_idx];
		// const f32a = new Float32Array(buffer, offsetStart, chr_len);
		
		const local_start = this._real_pos[ref_idx].start - 1;
		const local_end = this._real_pos[ref_idx].end - 1;
		const local_len = Math.min(local_end - local_start + 1, chr_len);

		const local_offsetStart = chrOffsetInByte[chr_idx] + local_start * Float32Array.BYTES_PER_ELEMENT;
		const local_lengthInBytes = local_len * Float32Array.BYTES_PER_ELEMENT;
		const local_offsetEnd = local_offsetStart + local_lengthInBytes;

		const chrEndOffset = chrOffsetStart + (chr_len * Float32Array.BYTES_PER_ELEMENT);
		if (local_offsetEnd > chrEndOffset) {
			debugger;
		}

		// "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_BS.methratio.float32"

		const data = chrLengthList.map(_ => []);

		if (this.debug) {
			// return new Float32Array(0);
			return data;
		}

		/** @type {ArrayBuffer} */
		const buffer = await (await fetch(config.url, {
			headers: {
				"Range": `bytes=${local_offsetStart}-${local_offsetEnd}`,
			}
		})).arrayBuffer();
		if (buffer && buffer.byteLength == (local_len * Float32Array.BYTES_PER_ELEMENT + 1)) {
			document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;
			const f32a = new Float32Array(buffer, 0, local_len);

			// const fff = new Float32Array(1);
			const uuu = new Uint32Array(buffer, 0, local_len);

			/** @see {@link _load_meth_ratio_float32} */
			f32a.forEach((value, index) => {
				if (!Number.isNaN(value) && value != 0) {
					const pos = Number(index) + 1;// .methratio.txt, base position = 0
					const r = new module_Methyl_ratioData();
					r.start = pos + local_start;
					r.end = r.start;
					r.value = Math.abs(value);

					// fff[0] = value;
					// r.strand = (uuu[0] >>> 31) ? -1 : 1;
					r.strand = (uuu[index] >>> 31) ? -1 : 1;

					data[chr_idx].push(r);
				}
			});
		}

		// @ts-ignore
		console.info(timeElapsed(), "loading: loaded methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		document.getElementById("status").innerText = "loading: loaded methy ratio:" + config.sample;
		await delayFrame(); // unlock UI // globalThis.gc()

		return data;
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 * @returns {Promise<module_Methyl_ratioData[][]>}
	 * @see {@link draw_methy_ratio_row}
	 */
	async old_load_meth_ratio_float32(config) {
		// throw new Error("ref pos to tseta");
		console.info(timeElapsed(), "loading: methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);

		/** @type {ArrayBuffer} */
		const buffer = await fetchData(config.url, "arraybuffer");
		document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;

		const ref = config.ref;
		if (ref == null) {
			throw new TypeError("required ref");
		}

		// uint32 uuu = 0;
		// float32* fff = &uuu;
		// *fff = value;
		// uint32 signbit = fff >> 31;
		const fff = new Float32Array(1);
		const uuu = new Uint32Array(fff.buffer);

		const data = _load_meth_ratio_float32(ref, buffer, (value, index) => {
			if (!Number.isNaN(value) && value != 0) { // && config.value_filter(value)
				// /** @type {number} */
				// const pos = pos_map[index];
				// methy_ratio_map[refId][sampleName].push({
				// 	// pos,
				// 	// value,
				// });
				const pos = index;
				const r = new module_Methyl_ratioData();
				// r.start = pos_map[pos];
				r.start = pos;
				r.end = r.start;
				r.value = value;

				fff[0] = value;
				r.strand = (uuu[0] >>> 31) ? -1 : 1;

				return r;
			}
		});

		// @ts-ignore
		console.info(timeElapsed(), "loading: loaded methy ratio:" + config.sample, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		document.getElementById("status").innerText = "loading: loaded methy ratio:" + config.sample;
		await delayFrame(); // unlock UI // globalThis.gc()

		return data;
	}

	/**
	 * load_nanopore_tsv
	 * @param {module_Methyl_sampleData} config
	 * @see {@link draw_methy_ratio_row}
	 */
	async load_np_meth_tsv(config) {
		throw new Error("ref pos to tseta");
		const text = await fetchData(config.url, "text");
		document.getElementById("status").innerText = "warn: ref1 pos:" + config.sample;

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == config.ref);
		// const pos_map = [
		// 	ref1_pos_uint32array,
		// 	ref2_pos_uint32array,
		// ][ref_idx];

		/** @type {module_Methyl_ratioData[][]} */
		const chr_meth_list = dataset.genome_info_list[ref_idx].chr_list.map(_ => []);

		if (config.region) {
			const _rows = parseTable(text, ["nChr", "start", "end", "value"]);
			_rows.forEach(row => {
				const r = new module_Methyl_ratioData();
				const nChr = Number(row.nChr);

				// r.start = pos_map[Number(row.start)];
				// r.end = pos_map[Number(row.end)];
				r.start = Number(row.start);
				r.end = Number(row.end);
				r.value = Number(row.value);

				chr_meth_list[nChr - 1].push(r);
			});
		}
		else {
			const _rows = parseTable(text, ["nChr", "start", "value"]);
			_rows.forEach(row => {
				const r = new module_Methyl_ratioData();
				const nChr = Number(row.nChr);

				// r.start = pos_map[Number(row.start)];
				r.start = Number(row.start);
				r.value = Number(row.value);
				r.end = r.start;

				chr_meth_list[nChr - 1].push(r);
			});
		}

		await delayFrame(); // unlock UI // globalThis.gc()

		return chr_meth_list;
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 */
	async load_bedgraph(config) {
		config.data = await this._load_tsv(config, ["sChr", "start", "end", "value"]);
		return config;
	}

	/**
	 * @param {module_Methyl_sampleData} config
	 */
	async load_bed(config) {
		config.data = await this._load_tsv(config, ["sChr", "start", "end"]);
		return config;
	}

	/**
	 * load bed
	 * @param {module_Methyl_sampleData} config
	 * @param {string[]|["sChr", "start", "end", "tagName", "value"]} headers
	 * param {(nChr: number, start: number, end: number, strValue: string) => module_Methyl_ratioData} callbackfn
	 * @param {(chrIdx: number, data: module_Methyl_ratioData, rawData: any) => boolean} [onLoadRow]
	 *
	 * @see {@link https://genome.ucsc.edu/goldenPath/help/bedgraph.html|BedGraph Track Format}
	 * @example load_bedgraph = config => _load_tsv(config, ["sChr", "start", "end", "value"])
	 *
	 * @see {@link https://www.biostars.org/p/10141/#10145|Bedgraph To Bed}
	 * @example load_bed = config => _load_tsv(config, ["sChr", "start", "end"])
	 */
	async _load_tsv(config, headers, onLoadRow) {
		const text = await fetchData(config.url, "text");

		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == config.ref);

		/** @type {{ [sChr: string]: number }} */
		const sChr_to_chrIdx = Object.fromEntries(dataset.genome_info_list[ref_idx].chr_list.map((v, i) => [v.symbol, i]));

		// const ref_pos_map = make_ref_pos_map_list()[ref_idx];

		/** @type {module_Methyl_ratioData[][]} */
		const chr_meth_list = dataset.genome_info_list[ref_idx].chr_list.map(_ => []);

		// ChIII_QM6a	98059	98060	QM6a_rad51_sae2_peak_1	35.487

		const _rows = parseTable(text, headers).filter(a => !a.sChr?.startsWith("#"));
		_rows.forEach(row => {
			const chrIdx = sChr_to_chrIdx[row.sChr];
			if (chrIdx != null) {// "sChr" "Ch"
				const r = new module_Methyl_ratioData();
				// r.start = ref_pos_map[Number(row.start) - 1];
				// r.end = ref_pos_map[Number(row.end) - 1];
				r.start = Number(row.start);
				r.end = Number(row.end);
				r.value = Number(row.value);

				// // nChr, start, end, strValue
				// const r = callbackfn(chrIdx + 1, Number(row.start), Number(row.end), row.value);

				if (onLoadRow == null || onLoadRow(chrIdx, r, row) == true) {
					chr_meth_list[chrIdx].push(r);
				}
			}
			else {
				console.warn({ "row.sChr": row.sChr, chrIdx, row });
			}
		});

		return chr_meth_list;
	}

	// _load_tsv_primitive() {
	// }
}

var html_Cm = `<ruby style="padding-right: 0.5em;">C<rt style="margin-right: -2.5em">m</rt></ruby>`;
var html_Chm = `<ruby style="padding-right: 0.5em;">C<rt style="margin-right: -2.5em">hm</rt></ruby>`;
var html_CmG = `<ruby style="padding-right: 0.5em;">C<rt style="margin-right: -2.5em">m</rt>G</ruby>`;
var html_CmhG = `<ruby style="padding-right: 0.5em;">C<rt style="margin-right: -2.5em">mh</rt>G</ruby>`;

/**
 * @param {string} ref
 * @param {ArrayBuffer} buffer
 * @param {(value: number, index: number) => any} func_transform
 */
function _load_meth_ratio_float32(ref, buffer, func_transform) {
	const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
	// const pos_map = [
	// 	ref1_pos_uint32array,
	// 	ref2_pos_uint32array,
	// ][ref_idx];
	/** @type {number[]} */
	const chrLengthList = dataset.genome_info_list[ref_idx].chr_list.map(a => a.length);

	const chrOffsetInByte = [];
	chrOffsetInByte[0] = 0;
	chrLengthList.reduce((pos, len, i) => {
		let offset = pos + (len * 4);
		chrOffsetInByte[i + 1] = offset;
		return offset;
	}, 0);

	const data = chrLengthList.map((chr_len, chrIdx) => {
		// const chrIdx = analysis_options.nChr - 1;
		const offsetStart = chrOffsetInByte[chrIdx];

		/** array view */
		const f32a = new Float32Array(buffer, offsetStart, chr_len);

		console.log(chrIdx, chr_len,{
			"f32a.length": f32a.length,
			offsetStart,
		});

		// methy_ratio_map[refId][sampleName] = [];
		/** @type {module_Methyl_ratioData[]} */
		const rows = [];

		f32a.forEach((value, index) => {
			const r = func_transform(value, index);
			if (r != null) {
				rows.push(r);
			}
		});

		return rows;
	});
	return data;
}

/**
 * @param {string} row1
 * @param {string} row2
 */
function make_methyl_value_desc(row1, row2) {
	function rl(s) {
		return s.replace(/Cm/g, html_Cm).replace(/Chm/g, html_Chm).replace(/CmG/g, html_CmG).replace(/CmhG/g, html_CmhG);
	}

	return (
`<div style="display: inline-block; text-align: center;">
	<div style="display: flex; flex-direction: column;">
		<div style="flex: 1; border-bottom: 1px solid;">
			${rl(row1)}
		</div>
		<div style="flex: 1;">
			${rl(row2)}
		</div>
	</div>
</div>`);
}

const inner_html_asdasdasdasd = `
<div>
	<style>
		rt {
		padding-left: 0em;
		margin-right: -2.5em;
		}
		ruby {
		margin-right: 0.5em;
		}
		table {
		border-spacing: 0px;
		border-collapse: collapse;
		line-height: 1em;
		}
	</style>
	<table style="height: 100%;">
		<tr>
			<td style="padding-right: 0.5em;">
				<table style="height: 100%;">
				<tr>
					<td style="vertical-align: baseline;">100</td>
				</tr>
				<tr>
					<td style="vertical-align: bottom;">0</td>
				</tr>
				</table>
			</td>
			<td style="padding-right: 0.5em;">
				<table style="text-align: center;">
				<tr style="border-bottom: 1px solid black;">
					<td><ruby>C<rt>m</rt></ruby> + <ruby>C<rt>hm</rt></ruby></td>
				</tr>
				<tr>
					<td><ruby>C<rt>m</rt></ruby> + <ruby>C<rt>hm</rt></ruby> + C</td>
				</tr>
				</table>
			</td>
			<td>%</td>
		</tr>
	</table>
</div>
`;

/**
 * @type {module_Methyl_sampleData[]}
 */
const methyl_dataset_list = [];

async function loadJSON(url) {
	return await (await fetch(url)).json();
}

/**
 * @param {"WT"|"dim2"|"rid1"|"dimrid"} type
 * @param {2|4|5|6|8} d
 * @param {"+"|"-"} strStrand
 */
function get_5hmC_by_BS_minus_oxBS(type, d, strStrand) {
	const chr_idx = viewerState.nChr - 1;

	const vec1Strand = strStrand == "+" ? 1 : (strStrand == "-" ? -1 : 0);

	const cfg_bs = methyl_dataset_list.find((v, i) => v.name == `${type} D${d} BS ${strStrand}`);
	const cfg_ox = methyl_dataset_list.find((v, i) => v.name == `${type} D${d} oxBS ${strStrand}`);
	const cfg_hmc = { ...cfg_bs };

	cfg_hmc.name = `${type} D${d} 5hmC ${strStrand}`;

	const hmc_arraybuffer = new Float32Array(dataset.genome_info_list[0].chr_list[0].length + 1);
	hmc_arraybuffer.fill(NaN);

	cfg_bs.data[chr_idx].forEach(bs_data => {
		hmc_arraybuffer[bs_data.start] = bs_data.value;
	});
	cfg_ox.data[chr_idx].forEach(oxbs_data => {
		if (oxbs_data.value > 0) {
			hmc_arraybuffer[oxbs_data.start] = 0;
		}
		// const bsVal = hmc_arraybuffer[oxbs_data.start];
		// hmc_arraybuffer[oxbs_data.start] = (Number.isNaN(bsVal) ? 0 : bsVal) - oxbs_data.value;
	});

	cfg_hmc.data = cfg_hmc.data.map(a_ => []);
	const hmc = cfg_hmc.data[chr_idx];

	hmc_arraybuffer.forEach((val, idx) => {
		if (!Number.isNaN(val) && val > 0) {
			hmc.push({
				start: idx,
				end: idx,
				strand: vec1Strand,
				value: val,
			});
		}
	});

	cfg_hmc.split_strand = false;
	// cfg_hmc.display_minus_value = true;

	methyl_dataset_list.push(cfg_hmc);
}

/**
 * @param {""|"QM6a"|"CBS1-1"} refId
 * @param {"WT"|"Pd"|"rid1"|"ddr"} type
 * @param {"BS"|"oxBS"|"ctrl"} seq_type
 * @param {number} threshold_RIP [0, 1]
 * @param {number[]} ignore_target_A_T
 * @see {@link draw_methy_ratio_row}
 */
async function load_methyl_group_by_RIP(refId = "QM6a", type, seq_type, threshold_RIP = 1, ignore_target_A_T) {
	const root = refId == "np59_Cdimrid" && type == "Pd" ? "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd" : {
		"QM6a": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a",
		"Qdim": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim",
		"Qrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid",
		"Crid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid",
		"np128_Qdimrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid",
		"np59_Cdimrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid",
	}[refId];

	await [
		"2", "4", "5", "6", "8"
	].reduce(async (prev_promise, dt) => {
		await prev_promise;

		const sample_group = `${type}_D${dt}`;
		const sample_id = `${sample_group}_${seq_type}`;

		const sample = await load_BS_nBS_methylratio_by_RIP(
			refId,
			sample_group,
			seq_type,
			`${root}/${refId}_${sample_id}.methratio.float32`,
			threshold_RIP,
			ignore_target_A_T
		);
		sample.name = `D × M ${sample.name}`;
	}, Promise.resolve());
}

/**
 * sample *ctrl
 * @param {""|"QM6a"|"CBS1-1"} refId
 * @param {string} sample_group
 * @param {string} seq_type
 * @param {string} url
 * @param {number} threshold_RIP [0, 1]
 * @param {number[]} ignore_target_A_T
 * @param {boolean} one_minus_value
 * @see {@link draw_methy_ratio_row}
 * @see {@link load_methratioFloat32_depthUint32_by_list} dname_mapping
 * @see {@link module_Methyl_sampleData}
 * @alias load_nBS_methylratio_by_RIP
 */
async function load_methyl_by_RIP(refId = "QM6a", sample_group, seq_type, url, threshold_RIP = 1, ignore_target_A_T, one_minus_value = false) {
	return await load_BS_nBS_methylratio_by_RIP(refId, sample_group, seq_type, url, threshold_RIP, ignore_target_A_T, one_minus_value);
}

/**
 * sample *ctrl
 * @param {""|"QM6a"|"CBS1-1"} refId
 * @param {string} sample_group
 * @param {string} seq_type
 * @param {string} url
 * @param {number} threshold_RIP [0, 1]
 * @param {number[]} ignore_target_A_T RIP target genome index (not include ref)
 * @param {boolean} one_minus_value
 * @see {@link draw_methy_ratio_row}
 * @see {@link load_methratioFloat32_depthUint32_by_list} dname_mapping
 * @see {@link module_Methyl_sampleData}
 * @alias load_methyl_by_RIP
 */
async function load_BS_nBS_methylratio_by_RIP(refId = "QM6a", sample_group, seq_type, url, threshold_RIP = 1, ignore_target_A_T, one_minus_value = false) {
	const dname_mapping = {
		"rid1": "rid1Δ",
		"Pd": "dim2Δ",
		"ddr": "dim2Δrid1Δ",
	};
	const display_name = sample_group.split("_").map(a => dname_mapping[a] ?? a).join(" ");

	const parental_role = {
		"QM6a": "D",
		"CBS1-1": "M",
		"Qdim": "D",
		"Qrid": "D",
		"Crid": "M",
		"np128_Qdimrid": "D",
		"np59_Cdimrid": "M",
	};

	const seqType_displayName = {
		"BS": "+BS",
		"ctrl": "-BS",
	}[seq_type];

	if (seqType_displayName == null) {
		debugger;
	}

	const sample = new module_Methyl_sampleData({
		ref: refId,
		sample: `${refId} ${sample_group} ${seq_type}`,
		name: `${display_name} ${seqType_displayName} >> ${parental_role[refId]}`,
		seq_type: seq_type,
		url: url,
		region: false,
		description: "dim2_rid1_dimrid/results",
		density_to_opacity: false,
		func_density_to_opacity: (view_length) => {
			if (g_methylRenderer.isLimitLoad()) {
				return 1;
			}
			else if (view_length > 3_000_000) {
				return 0.72;
			}
			else if (view_length > 1_000_000) {
				return 0.8;
			}
			else if (view_length > 500_000) {
				return 0.9;
			}
			else {
				return 1;
			}
		},
	});

	return await _load_BS_nBS_methylratio_by_RIP(sample, threshold_RIP, ignore_target_A_T, one_minus_value);
}

/**
 * sample *ctrl
 * @param {module_Methyl_sampleData} sample
 * @param {number} threshold_RIP [0, 1]
 * @param {number[]} ignore_target_A_T RIP target genome index (not include ref)
 * @param {boolean} one_minus_value
 * @see {@link draw_methy_ratio_row}
 * @see {@link load_methratioFloat32_depthUint32_by_list} dname_mapping
 * @see {@link module_Methyl_sampleData}
 * @alias load_methyl_by_RIP
 */
async function _load_BS_nBS_methylratio_by_RIP(sample, threshold_RIP = 1, ignore_target_A_T, one_minus_value = false) {
	const TAG_ZERO_READ = 0;
	const TAG_5mC = 1;
	const TAG_5hmC = 2;
	const TAG_RIP = 3;
	const TAG_REMOVED_RIP = 4;
	const TAG_C = 5;

	const refId = sample.ref;
	const sample_group = sample.sample;
	const seq_type = sample.seq_type;
	const url = sample.url;

	sample.rendering_condition = [
		new module_Methyl_sampleData_RenderingCondition({ color: "#FF9F", condition: v => v.value > 0 && v.t == TAG_ZERO_READ,   min_width: 1, desc: "TAG_ZERO_READ", }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#F00F", condition: v => v.value > 0 && v.t == TAG_5mC,         min_width: 1, desc: "TAG_5mC", }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#080F", condition: v => v.value > 0 && v.t == TAG_5hmC,        min_width: 1, desc: "TAG_5hmC", }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#000F", condition: v => v.value > 0 && v.t == TAG_RIP,         min_width: 1, desc: "TAG_RIP", }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#F0F0", condition: v => v.value > 0 && v.t == TAG_REMOVED_RIP, min_width: 1, desc: "TAG_REMOVED_RIP", }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#000F", condition: v => v.value > 0 && v.t == TAG_C,           min_width: 1, desc: "TAG_C", }),
	];

	// methyl_dataset_list.filter(c => c.rendering_condition[0].color == "#FF9").filter(c => {
	// 	c.rendering_condition.push(c.rendering_condition.shift());// move to back
	// 	return c;
	// });
	sample.data = dataset.genome_info_list[0].chr_list.map(_ => []);

	// if (sample.name == "dim2Δrid1Δ D8 -BS >> D") { // config.sample == "ddr_D8 ctrl"
	// 	debugger;// 20220921 all NaN
	// }

	const check_readCoverage = false;
	
	const [
		f32a_ctrl,
		ui32a_coverage,
	] = await Promise.all([
		load_float32(),
		check_readCoverage ? load_Uint32() : null,
	]);
	if (f32a_ctrl == null || (check_readCoverage && ui32a_coverage == null)) {
		return;
	}

	function is_valid(v) {
		return !Number.isNaN(v) && Number.isFinite(v) && v != 0;
	}

	const rip_pair = {
		"QM6a": { src: "QM6a", dst: "CBS1-1", },
		"CBS1-1": { src: "QM6a", dst: "CBS1-1", },
		"Qdim": { src: "Qdim", dst: "np59_Cdimrid", },
		"Qrid": { src: "Qrid", dst: "Crid", },
		"Crid": { src: "Qrid", dst: "Crid", },
		"np128_Qdimrid": { src: "np128_Qdimrid", dst: "np59_Cdimrid", },
		"np59_Cdimrid": { src: "np128_Qdimrid", dst: "np59_Cdimrid", },
	};
	const ripSrc_idx = dataset.genome_info_list.findIndex(a => a.name == rip_pair[sample.ref].src);
	const ripSrc_to_pos = make_ref_pos_map_list()[ripSrc_idx];

	// // 20220920
	// ignore_target_A_T = [
	// 	dataset.genome_info_list.findIndex(a => a.name == rip_pair[sample.ref].dst),
	// ];

	function in_visible_region(pos) {//2867589 // 2737467
		if (g_methylRenderer.isLimitLoad()) {
			// const idx = dataset.genome_info_list.findIndex(v => v.name == refId);
			// const region_start;
			// const region_end;
			return pos > g_chrCoord.bp_start && pos < g_chrCoord.bp_end;
		}
		else {
			return true;
		}
	}
	const is_RIP = (function () {
		if (ignore_target_A_T == null || ignore_target_A_T.length <= 0) {
			return function is_RIP_tr(ref_pos) {
				return false;
			};
		}
		else {
			return function is_RIP(ref_pos) {
				// if (ref_pos == 341676) {
				// 	debugger
				// }
				// if (ignore_target_A_T == null || ignore_target_A_T.length <= 0) {
				// 	return true;
				// }
				const pos = ripSrc_to_pos[ref_pos];
				if (in_visible_region(pos)) {
					const q_val = seq_list[ripSrc_idx][pos - 1];// seq[idx]
			
					return ignore_target_A_T.some(ripDst_idx => {
						const c_val = seq_list[ripDst_idx][pos - 1];// seq[idx]
						if (q_val == "C" && c_val == "T") {
							return true;
						}
						else if (q_val == "G" && c_val == "A") {
							return true;
						}
						return false;
					});
				}
				return false;
			};
		}
	})();

	const data_arr = sample.data[viewerState.nChr - 1];
	if (check_readCoverage && ui32a_coverage) {
		ui32a_coverage.forEach((depth, idx) => {
			if (depth == 0) {
				const ref_pos = Number(idx) + 1;
				const val = new module_Methyl_ratioData();
				val.start = ref_pos;
				val.end = ref_pos;
				delete val.strand;
				val.value = 1;
				val.t = TAG_ZERO_READ;
				data_arr.push(val);
			}
		});
	}
	f32a_ctrl.forEach((_, idx) => {
		const ctrl = Math.abs(f32a_ctrl[idx]);
		const ref_pos = Number(idx) + 1;

		if (is_valid(ctrl) && ctrl <= threshold_RIP) {
			const val = new module_Methyl_ratioData();
			val.start = ref_pos;
			val.end = ref_pos;
			
			delete val.strand;
			
			if (one_minus_value) {
				val.value = 1 - ctrl;
			}
			else {
				val.value = ctrl;
			}

			if (dataset.mode != "single" && is_RIP(ref_pos)) {
				val.t = TAG_REMOVED_RIP;
			}
			else {
				val.t = TAG_C;// TAG_RIP
			}
			
			data_arr.push(val);
		}
	});

	if (data_arr.length == 0) {
		console.log("if (data_arr.length == 0) {", sample.url);
		debugger;// 20220921
	}

	methyl_dataset_list.push(sample);

	async function load_float32() {
		const f32_sample = new module_Methyl_sampleData({
			// ...sample,
			ref: refId,
			sample: `${sample_group} ${seq_type}`,
			name: `${sample_group.replace(/_/g, " ")} ${seq_type}`,//seq_type <- seqType_displayName
			url: url,
		});
		const a = await g_methylRenderer.load_float32Array(f32_sample);
		sample.f32 = f32_sample;
		await delayFrame();
		return a;
	}
	async function load_Uint32() {
		const ui32_sample = new module_Methyl_sampleData({
			// ...sample,
			ref: refId,
			sample: `${sample_group} ${seq_type}`,
			name: `${sample_group.replace(/_/g, " ")} ${seq_type}`,//seq_type <- seqType_displayName
			url: url.replace(/\.methratio\.float32$/, ".depth.uint32"),
		});
		const a = await g_methylRenderer.load_uint32Array(ui32_sample);
		sample.ui32 = ui32_sample;
		await delayFrame();
		return a;
	}

	return sample;
}

/**
 * @param {number[]} ignore_are_not_C_G
 */
function remove_parental_are_not_C_G(ignore_are_not_C_G) {
	methyl_dataset_list.forEach(sample => {
		_remove_parental_are_not_C_G(sample, ignore_are_not_C_G);
	});
}

/**
 * @param {module_Methyl_sampleData} sample
 * @param {number[]} ignore_are_not_C_G
 */
function _remove_parental_are_not_C_G(sample, ignore_are_not_C_G) {
	const ref_idx = dataset.genome_info_list.findIndex(a => a.name == sample.ref);
	if (ref_idx < 0) {
		return;
	}
	const ref_to_pos = make_ref_pos_map_list()[ref_idx];

	function is_C_G(ref_pos) {
		if (ignore_are_not_C_G == null || ignore_are_not_C_G.length <= 0) {
			return true;
		}
		const pos = ref_to_pos[ref_pos];
		// const q_val = seq_list[ref_idx][pos - 1];// seq[idx]

		const is_c = ignore_are_not_C_G.every(target_idx => {
			const c_val = seq_list[target_idx][pos - 1];// seq[idx]
			return c_val == "C";
		});

		const is_g = ignore_are_not_C_G.every(target_idx => {
			const c_val = seq_list[target_idx][pos - 1];// seq[idx]
			return c_val == "G";
		});

		return (is_c ^ is_g);
	}

	const data_arr = sample.data[viewerState.nChr - 1];

	data_arr.forEach(point => {
		point.value = Math.abs(point.value);

		if (!is_C_G(point.start)) {
			point.value = -point.value;
		}
	});
}

/**
 * @param {string} name
 * @returns {string}
 */
function ref_to_display_name(name) {
	const display_name = dataset.display_name?.[name] ?? name;
	return display_name;
}

/**
 * @param {number} target_A_T target genome index
 * @param {boolean} insert_last
 * @example
 * add_RIP_marker(0, dataset.genome_info_list.findIndex(a => a.name == "np128_Qdimrid"))
 */
function add_RIP_marker(parental_index, target_A_T, insert_last = true) {
	const p_name = ref_to_display_name(dataset.genome_info_list[parental_index].name);
	const t_name = ref_to_display_name(dataset.genome_info_list[target_A_T].name);
	
	if (parental_index != 0) {
		console.log(`add_RIP_marker: Why ${p_name} vs ${t_name} ?`);
		debugger;
	}

	const sample = new module_Methyl_sampleData({
		ref: "",
		sample: `${p_name}-${t_name}`,
		name: `${p_name} vs ${t_name} C-to-T`,
		split_strand: false,
		value_desc: null,
		html_value_desc: null,
		url: null,
		region: false,
		description: null,
		density_to_opacity: false,
		func_density_to_opacity: (view_length) => {
			return 1;
		},
		rendering_condition: [
			new module_Methyl_sampleData_RenderingCondition({ color: "#000", condition: () => true, }),
		],
	});
	sample.data = dataset.genome_info_list[0].chr_list.map(_ => []);

	const chr_idx = viewerState.nChr - 1;

	const data = sample.data[chr_idx];
	for(let pos_idx = 1; pos_idx <= seq_list[parental_index].length; ++pos_idx) {
		const p = seq_list[parental_index][pos_idx];
		const d = seq_list[target_A_T][pos_idx];
		if (
			(p == "C" && d == "T") ||
			(p == "G" && d == "A")
		) {
			const a = new module_Methyl_ratioData();
			a.start = pos_idx + 1;
			a.end = a.start;
			a.value = 1;
			delete a.strand;
			data.push(a);
		}
	}

	if (insert_last) {
		methyl_dataset_list.push(sample);
	}

	return sample;
}

function func_density_to_opacity(view_length) {
	return 0.8
}

async function load_dim2_rid1_dimrid_methratio(refId = "QM6a", type = "WT", seq_type = "BS", is_reload = false) {
	if (type == "WT") {
		if (seq_type == "BS") {// 20200806 WT vs dim2
			const veg_bs = new module_Methyl_sampleData({
				ref: "QM6a",
				sample: "QM6a WT veg", // #2
				name: "QM6a WT veg +BS >> D", // #2
				url: "data/BS_QM6a/QM6a_methyl.methratio.float32",
				region: false,
				description: "20200806/BatMeth2/",
				density_to_opacity: false,
				func_density_to_opacity: func_density_to_opacity,
			});
			await load_methratioFloat32_depthUint32(seq_type, veg_bs);
		}
		if (seq_type == "BS") {
			if (false) {
				// CBS1-1 WT veg BS map to QM6a
				const veg_bs = new module_Methyl_sampleData({
					ref: "QM6a",
					sample: "QM6a CBS1-1 WT veg",
					name: "CBS1-1 WT veg +BS >> D",
					url: "data/BS_20200411/QM6a_CBS_WT_BS.methratio.float32",
					region: false,
					description: "BS-seq/20200411_BS_results/",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, veg_bs);
			}
			else {
				const veg_bs = new module_Methyl_sampleData({
					ref: "CBS1-1",
					sample: "CBS1-1 WT veg",
					name: "CBS1-1 WT veg +BS >> M",
					url: "data/BS_20200411/CBS1-1/20200902_C_wt.methratio.float32",
					region: false,
					description: "BS-seq/20200411_BS_results/",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, veg_bs);
			}
		}
		if (seq_type == "BS") {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"WT_D2_BS",
				"WT_D4_BS",
				"WT_D5_BS",
				"WT_D6_BS",
				"WT_D8_BS",
			]);
		}
		else if (seq_type[0] == "o" || seq_type == "ox" || seq_type == "oxBS") {
			{
				const veg_obs = new module_Methyl_sampleData({
					ref: "QM6a",
					sample: "NP43",
					name: "QM6a oxBS",
					// split_strand: true,
					// // value_desc: "M / (M + H + C)",
					// html_value_desc: make_methyl_value_desc("Cm", "Cm + Chm + C"),
					url: "data/BS_oBS/QM6a+oBS-Zymo.methratio.float32",
					region: false,
					description: "BS-seq/BS_oBS-merged_reads",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, veg_obs);
			}
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"WT_D2_oxBS",
				"WT_D4_oxBS",
				"WT_D5_oxBS",
				"WT_D6_oxBS",
				"WT_D8_oxBS",
			]);
		}
		else {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"WT_D2_ctrl",
				"WT_D4_ctrl",
				"WT_D5_ctrl",
				"WT_D6_ctrl",
				"WT_D8_ctrl",
			]);
		}
	}
	else if (type == "dim2") {
		if (seq_type == "BS") {
			if (refId == "QM6a") {
				const sample = new module_Methyl_sampleData({
					ref: "QM6a",
					sample: "QM6a dim2 veg",
					name: "QM6a dim2Δ veg +BS",
					// split_strand: true,
					// // value_desc: "(M + H) / (M + H + C)",
					// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
					url: `data/BS_QM6a/QM6a_dim2_methyl.methratio.float32`,
					region: false,
					description: "BS-seq/BS-seq_20200806/BatMeth2",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, sample);
			}
			else if (refId == "Qdim") {
				// alert("QM6a dim2Δ >> np30");//???????
				
				const sample = new module_Methyl_sampleData({
					ref: "Qdim",
					sample: "QM6a dim2 veg",
					name: "QM6a dim2Δ veg +BS",
					// split_strand: true,
					// // value_desc: "(M + H) / (M + H + C)",
					// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
					url: `data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_veg_BS.methratio.float32`,
					region: false,
					description: "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, sample);
			}

			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"Pd_D2_BS",
				"Pd_D4_BS",
				"Pd_D5_BS",
				"Pd_D6_BS",
				"Pd_D8_BS",
			]);
		}
		else if (seq_type[0] == "o" || seq_type == "ox" || seq_type == "oxBS") {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"Pd_D2_oxBS",
				"Pd_D4_oxBS",
				"Pd_D5_oxBS",
				"Pd_D6_oxBS",
				"Pd_D8_oxBS",
			]);
		}
		else {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"Pd_D2_ctrl",
				"Pd_D4_ctrl",
				"Pd_D5_ctrl",
				"Pd_D6_ctrl",
				"Pd_D8_ctrl",
			]);
		}
	}
	else if (type == "rid1") {
		if (seq_type == "BS") {
			/** @see {@link load_methly_preset_20210924} */
			await load_QM6a_CBS1_rid1_BS(refId, seq_type);
		}
		else if (seq_type[0] == "o" || seq_type == "ox" || seq_type == "oxBS") {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"rid1_D2_oxBS",
				"rid1_D4_oxBS",
				"rid1_D5_oxBS",
				"rid1_D6_oxBS",
				"rid1_D8_oxBS",
			]);
		}
		else {
			alert("deprecated, use load_methly_preset_20210924 instead.");

			{
				const sample = new module_Methyl_sampleData({
					ref: refId,
					sample: "QM6a CBS1-1 rid1 veg",
					name: "CBS1-1 rid1 veg -BS >> D",
					// split_strand: true,
					// // value_desc: "(M + H) / (M + H + C)",
					// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
					url: `data/BS_20200411/QM6a_Crid_ctrl.methratio.float32`,
					region: false,
					description: "/wcli/BatMeth2/",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, sample);
			}

			/**
			 * @example
			 * rid_D[1-8]_ctrl low density
			 */
			var x;
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"rid1_D2_ctrl",
				"rid1_D4_ctrl",
				"rid1_D5_ctrl",
				"rid1_D6_ctrl",
				"rid1_D8_ctrl",
			]);
		}
	}
	else if (type == "dimrid") {
		if (seq_type == "BS") {
			try {
				// 20220928 ddr parental
				const sample = new module_Methyl_sampleData({
					ref: "np128_Qdimrid",
					name: "dim2Δ ridΔ veg +BS >> D",// dim2 rid1 (1-2) veg +BS >> D
					sample: "dim2Δ ridΔ veg",
					url: `data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_veg_BS.methratio.float32`,
					region: false,
					description: "",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, sample);
			}
			catch {
				alert("dim2 rid1 (1-1) veg +BS >> M");
			}

			if ("ddr 1-1") {
				// 20220928 ddr parental
				const sample = new module_Methyl_sampleData({
					ref: "np59_Cdimrid",
					name: "dim2Δ ridΔ veg +BS >> M",// dim2 rid1 (1-2) veg +BS >> D
					sample: "dim2Δ ridΔ veg",
					url: `data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_veg_BS.methratio.float32`,
					region: false,
					description: "",
					density_to_opacity: false,
					func_density_to_opacity: func_density_to_opacity,
				});
				await load_methratioFloat32_depthUint32(seq_type, sample);
			}

			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"ddr_D2_BS",
				"ddr_D4_BS",
				"ddr_D5_BS",
				"ddr_D6_BS",
				"ddr_D8_BS",
			]);
		}
		else if (seq_type[0] == "o" || seq_type == "ox" || seq_type == "oxBS") {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"ddr_D2_oxBS",
				"ddr_D4_oxBS",
				"ddr_D5_oxBS",
				"ddr_D6_oxBS",
				"ddr_D8_oxBS",
			]);
		}
		else {
			await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
				"ddr_D2_ctrl",
				"ddr_D4_ctrl",
				"ddr_D5_ctrl",
				"ddr_D6_ctrl",
				"ddr_D8_ctrl",
			]);
		}
	}
	// await drawFrame();

	// viewerState.resizeCanvas();
	// await drawFrame();
}

/**
 * @param {string} refId
 * @param {"BS"} seq_type
 * @see {@link load_methly_preset_20210924}
 */
async function load_QM6a_CBS1_rid1_BS(refId, seq_type) {
	const parental_map_to_QM6a = false;
	if (parental_map_to_QM6a) {
		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: "QM6a rid1 veg",
			name: "QM6a rid1 veg +BS >> QM6a",
			// split_strand: true,
			// // value_desc: "(M + H) / (M + H + C)",
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: `data/BS_20200411/QM6a_rid1_BS.methratio.float32`,
			region: false,
			description: "/wcli/BatMeth2/",
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});
		await load_methratioFloat32_depthUint32(seq_type, sample);
	}
	if (parental_map_to_QM6a) {
		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: "QM6a CBS1-1 rid1 veg",
			name: "CBS1-1 rid1 veg +BS >> QM6a",
			// split_strand: true,
			// // value_desc: "(M + H) / (M + H + C)",
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: `data/BS_20200411/QM6a_Crid_BS.methratio.float32`,
			region: false,
			description: "/wcli/BatMeth2/",
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});
		await load_methratioFloat32_depthUint32(seq_type, sample);
	}
	
	{
		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: "QM6a rid1 veg",
			name: "QM6a rid1 veg +BS >> D",
			url: `data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_BS.methratio.float32`,// 20220930
			region: false,
			description: "BS-seq/20200411_BS_results/BatMeth2_parental/results/",
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});
		await load_methratioFloat32_depthUint32(seq_type, sample);
	}
	{
		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: "CBS1-1 rid1 veg",
			name: "CBS1-1 rid1 veg +BS >> M",
			url: `data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_BS.methratio.float32`,// 20220930
			region: false,
			description: "BS-seq/20200411_BS_results/BatMeth2_parental/results/",
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});
		await load_methratioFloat32_depthUint32(seq_type, sample);
	}

	await load_methratioFloat32_depthUint32_by_list(refId, seq_type, [
		"rid1_D2_BS",
		"rid1_D4_BS",
		"rid1_D5_BS",
		"rid1_D6_BS",
		"rid1_D8_BS",
	]);
}

/**
 * @param {string} refId
 * @param {"BS"|"oxBS"|"ctrl"} seq_type
 * @param {string[]} data_list
 * @see {@link load_BS_nBS_methylratio_by_RIP} dname_mapping
 * @see {@link load_methyl_group_by_RIP} np59_Cdimrid Pd
 * @see {@link module_Methyl_sampleData}
 */
async function load_methratioFloat32_depthUint32_by_list(refId, seq_type, data_list, is_reload = false) {
	const root = refId == "np59_Cdimrid" && data_list.every(a => a.indexOf("Pd") >= 0) ? "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd" : {
		"QM6a": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a",
		"Qdim": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim",
		"Qrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid",
		"Crid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid",
		"np128_Qdimrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid",
		"np59_Cdimrid": "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid",
	}[refId];
	const parental_role = {
		"QM6a": "D",
		"CBS1-1": "M",
		"Qdim": "D",
		"Qrid": "D",
		"Crid": "M",
		"np128_Qdimrid": "D",
		"np59_Cdimrid": "M",
	};

	for (let sample_id of data_list) {
		const dname_mapping = {
			"rid1": "rid1Δ",
			"Pd": "dim2Δ",
			"ddr": "dim2Δrid1Δ",
		}
		const display_name = sample_id.split("_").map(a => dname_mapping[a] ?? a).join(" ");

		const seqType_displayName = {
			"BS": "+BS",
			"ctrl": "-BS",
		}[seq_type];

		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: `${refId} ${sample_id} BS`,
			name: `D × M ${display_name} ${seqType_displayName} >> ${parental_role[refId]}`,
			// split_strand: true,
			// // value_desc: "(M + H) / (M + H + C)",
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: `${root}/${refId}_${sample_id}.methratio.float32`,
			region: false,
			description: "dim2_rid1_dimrid/results",
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});
		await load_methratioFloat32_depthUint32(seq_type, sample);
	}

	await g_methylRenderer.__html_to_canvas(methyl_dataset_list, is_reload);

	await delayFrame();
}


/**
 * @param {string} refId
 * @param {"QM6a"|"CBS1"|"Qdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid"|"D × M"} sampleGenome
 * @param {"BS"|"oxBS"|"ctrl"} seq_type
 * @param {"WT" | "Pd" | "rid1" | "ddr"} mut
 * @param {("veg" | "D2" | "D4" | "D5" | "D6" | "D8")[]} day_list
 * @param {?string} [_root_path]
 * @see {@link load_BS_nBS_methylratio_by_RIP} dname_mapping
 * @see {@link load_methyl_group_by_RIP} np59_Cdimrid Pd
 * @see {@link module_Methyl_sampleData}
 * @see {@link load_methratioFloat32_depthUint32_by_list} // ignore_target_A_T, rollback 20220920
 */
 async function __load_methratioFloat32_depthUint32_by_list(refId, sampleGenome, seq_type, mut, day_list, _root_path) {
	const root = _root_path || {
		"WT": "data/methyl_20221007/WT",
		"Pd": "data/methyl_20221007/dim2",
		"rid1": "data/methyl_20221007/rid1",
		"ddr": "data/methyl_20221007/dim2rid1",
	}[mut];

	const mut_display = {
		"WT": "WT",
		"Pd": "dim2Δ",
		"rid1": "rid1Δ",
		"ddr": "dim2Δrid1Δ",
	}[mut];

	const parental_role = {
		"QM6a": "D",
		"CBS1-1": "M",
		"Qdim": "D",
		"Qrid": "D",
		"Crid": "M",
		"np128_Qdimrid": "D",
		"np59_Cdimrid": "M",
	};

	const ref_display_name = {
		"np30": "QM6a dim2Δ",//QdimΔ // np30
		"np31": "CBS1-1 rid1Δ",//CridΔ // np31
		"np32": "QM6a rid1Δ",//QridΔ // np32
		"np128_Qdimrid": "dim2Δrid1Δ (1-2)", // np128_Qdimrid
		"np59_Cdimrid": "dim2Δrid1Δ (1-1)", // np59_CdimΔridΔ

		"QM6a": "QM6a",
		"CBS1-1": "CBS1-1",

		"Qdim": "QM6a dim2Δ",
		"Crid": "CBS1-1 rid1Δ",
		"Qrid": "QM6a rid1Δ",
	};
	
	/**
	 * @see {load_methratioFloat32_depthUint32}
	 * @see {__load_methratioFloat32_depthUint32_by_list}
	 */
	const rip_target = seq_type == "ctrl" ? [
		dataset.genomeNameList.indexOf({
			"QM6a": "CBS1-1",
			"CBS1-1": "QM6a",
			"Qdim": "np59_Cdimrid",
			"Qrid": "Crid",
			"Crid": "Qrid",
			"np128_Qdimrid": "CBS1-1",
			"np59_Cdimrid": "np128_Qdimrid",
		}[refId]),
	] : [];

	for (let day of day_list) {
		const display_name = [
			(ref_display_name[sampleGenome] ?? sampleGenome),
			day != "veg" ? mut_display : undefined,
			day,
			_root_path ? "#2" : undefined,
		].filter(a => a).join(" ");

		const seqType_displayName = {
			"BS": "+BS",
			"ctrl": "-BS",
			"oxBS": "oBS",
		}[seq_type];

		const sample = new module_Methyl_sampleData({
			ref: refId,
			sample: `${refId} ${mut} ${day} ${seqType_displayName}`,
			name: `${display_name} ${seqType_displayName} >> ${parental_role[refId]}`,
			seq_type: seq_type,
			// split_strand: true,
			// // value_desc: "(M + H) / (M + H + C)",
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: `${root}/${refId}_${mut}_${day}_${seq_type}.methratio.float32`,
			region: false,
			description: root,
			density_to_opacity: false,
			func_density_to_opacity: func_density_to_opacity,
		});

		// await load_methratioFloat32_depthUint32(seq_type, sample);
		// await load_BS_nBS_methylratio_by_RIP(sample.ref, sample.sample, seq_type, sample.url, 1.0, []);
		await _load_BS_nBS_methylratio_by_RIP(sample, 1.0, rip_target);
	}

	// await g_methylRenderer.__html_to_canvas(methyl_dataset_list, is_reload);

	await delayFrame();
}

/**
 * add row methratio.float32, add row depth.uint32
 * @warp load_BS_nBS_methylratio_by_RIP
 * @param {"BS"|"oxBS"|"ctrl"} seq_type
 * @param {module_Methyl_sampleData} sample
 * @see {@link load_BS_nBS_methylratio_by_RIP}
 */
async function load_methratioFloat32_depthUint32(seq_type, sample) {
	if (1) {// 20220916, 20220920
		// if (!rip_target || rip_target.length == 0) {
		// 	debugger
		// 	alert(`seq_type=${seq_type}; sample=${sample};`);
		// }
		const rip_target = [
			dataset.genomeNameList.indexOf({
				"QM6a": "CBS1-1",
				"CBS1-1": "QM6a",
				"Qdim": "np59_Cdimrid",
				"Qrid": "Crid",
				"Crid": "Qrid",
				"np128_Qdimrid": "CBS1-1",
				"np59_Cdimrid": "np128_Qdimrid",
			}[sample.ref]),
		];
		const nnn = await load_BS_nBS_methylratio_by_RIP(sample.ref, sample.sample, seq_type, sample.url, 1.0, rip_target);
		nnn.name = sample.name;
		return nnn;
	}

	const check_readCoverage = false;
	let ui32a_coverage;

	[
		sample.data,
		ui32a_coverage,
	] = await Promise.all([
		g_methylRenderer.load_meth_ratio_float32(sample),
		// load_Float32(),
		check_readCoverage ? load_Uint32() : null,
	]);

	const TAG_ZERO_READ = 0;
	const TAG_5mC = 1;
	const TAG_5hmC = 2;
	const TAG_RIP = 3;

	sample.rendering_condition = [
		new module_Methyl_sampleData_RenderingCondition({ color: "#FF9", condition: v => v.t == TAG_ZERO_READ, min_width: 1, }),
	];
	if (seq_type == "BS") {
		sample.rendering_condition.push(
			new module_Methyl_sampleData_RenderingCondition({ color: "#080", condition: v => v.t != TAG_ZERO_READ, min_width: 1 / 4, }),
		);
	}
	else if (seq_type[0] == "o" || seq_type == "ox" || seq_type == "oxBS") {
		sample.rendering_condition.push(
			new module_Methyl_sampleData_RenderingCondition({ color: "#F00", condition: v => v.t != TAG_ZERO_READ, min_width: 1 / 4, }),
		);
	}
	else {
		sample.rendering_condition.push(
			new module_Methyl_sampleData_RenderingCondition({ color: "#000", condition: v => v.t != TAG_ZERO_READ, min_width: 1 / 4, }),
		);
	}

	const data_arr = sample.data[viewerState.nChr - 1];
	if (check_readCoverage && ui32a_coverage) {
		ui32a_coverage.forEach((depth, idx) => {
			if (depth == 0) {
				const ref_pos = Number(idx) + 1;
				const val = new module_Methyl_ratioData();
				val.start = ref_pos;
				val.end = ref_pos;
				delete val.strand;
				val.value = 1;
				val.t = TAG_ZERO_READ;
				data_arr.push(val);
			}
		});
	}

	// function is_valid(v) {
	// 	return !Number.isNaN(v) && v != 0 && Number.isFinite(v);
	// }
	// f32a_ratio.forEach((_, idx) => {
	// 	const ctrl = Math.abs(f32a_ratio[idx]);
	// 	const ref_pos = Number(idx) + 1;

	// 	if (dataset.mode != "single" && !check_RIP(ref_pos)) {
	// 		return;
	// 	}

	// 	if (is_valid(ctrl) && ctrl <= threshold_RIP) {
	// 		const val = new module_Methyl_ratioData();
	// 		val.start = ref_pos;
	// 		val.end = ref_pos;
	// 		delete val.strand;
	// 		val.value = 1 - ctrl;
	// 		val.t = TAG_RIP;
	// 		data_arr.push(val);
	// 	}
	// });

	const list = [sample];
	// g_methylRenderer.split_strand(list);

	if (list.some(v => v.data[viewerState.nChr - 1].length == 0)) {
		debugger;// 20220921
	}

	await delayFrame();

	methyl_dataset_list.push(...list.filter(a => !a.hide));

	async function load_Uint32() {
		const a = await g_methylRenderer.load_uint32Array(new module_Methyl_sampleData({
			// ...sample,
			ref: sample.ref,
			sample: sample.name,
			name: sample.name.replace(/_/g, " "),
			url: sample.url.replace(/\.methratio\.float32$/, ".depth.uint32"),
		}));
		await delayFrame();
		return a;
	}
}

/**
 * @see {@link capture_all_5mC}
 */
async function load_5mc_debug(idx, nChr) {
	g_methylRenderer.debug = true;
	g_methylRenderer.display = false;

	const refType = [
		"WT",
		"Qdim",
		"Qrid",
		"np128_Qdimrid",
		"np59_Cdimrid",
	][idx];

	await load_methly_preset_20210924(refType, true);
	document.querySelector("#el_input_chr").value = nChr;
	document.querySelector("#el_input_chr").oninput(null);
	await promise_load_task;
	await load_methly_preset_20210924(refType, false);

	g_methylRenderer.display = true;
	drawFrame()
}

/**
 * @see {@link mark_methyl_point}
 * @param {"WT"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid"} ref_type
 * @param {boolean} init_dataset
 */
async function load_methly_preset_20210924(ref_type, init_dataset) {
	document.querySelector("#data-rows > div").style.width = "18em";

	dataset.display_name = {
		"QM6a": "QM6a",
		"CBS1-1": "CBS1-1",
		
		"np43_QM6a": "QM6a WT (D-np43)",
		"np42_CBS1-1": "CBS1-1 WT (M-np42)",

		"np30": "QM6a dim2Δ (D-np30)",//QdimΔ // np30
		"np31": "CBS1-1 rid1Δ (M-np31)",//CridΔ // np31
		"np32": "QM6a rid1Δ (D-np32)",//QridΔ // np32
		"np128_Qdimrid": "dim2Δrid1Δ (D-np128)", // np128_Qdimrid
		"np59_Cdimrid": "dim2Δrid1Δ (M-np59)", // np59_CdimΔridΔ

		"Qdim": "QM6a dim2Δ (D-np30)",
		"Crid": "CBS1-1 rid1Δ (M-np31)",
		"Qrid": "QM6a rid1Δ (D-np32)",
	};

	if (ref_type == "WT") {
		if (init_dataset) {
			RemoveGenomeHelper.rename("np30", "Qdim").rename("np31", "Crid").rename("np32", "Qrid");
			RemoveGenomeHelper.remove("Qdim").remove("Crid").remove("Qrid").remove("np128_Qdimrid").remove("np59_Cdimrid");
		}
		else {
			setPlotTitle("WT(D)", "WT(M)");

			const QM6a = dataset.genome_info_list.findIndex(a => a.name == "QM6a");
			const CBS1 = dataset.genome_info_list.findIndex(a => a.name == "CBS1-1");
			if (QM6a != 0) {
				alert("QM6a must is first parental");
			}

			add_RIP_marker(QM6a, CBS1);

			await load_dim2_rid1_dimrid_methratio("QM6a", "WT", "BS", true);
			await load_BS_nBS_methylratio_by_RIP("QM6a", "QM6a WT veg", "ctrl", "data/BS_QM6a/QM6a_control.methratio.float32", 1.0, [1]);
			await load_BS_nBS_methylratio_by_RIP("QM6a", "CBS1-1 WT veg", "ctrl", "data/nBS/QM6a_CBS_WT_ctrl.methratio.float32", 1.0, [1]);
			await load_methyl_group_by_RIP("QM6a", "WT", "ctrl", 1.0, [1]);
		}
	}
	else if (ref_type == "Qdim" || ref_type == "Cdim") {
		if (init_dataset) {
			RemoveGenomeHelper.rename("np30", "Qdim").rename("np31", "Crid").rename("np32", "Qrid");
			RemoveGenomeHelper.remove("Crid").remove("Qrid");
			RemoveGenomeHelper.remove("np128_Qdimrid").rebuild();
		}
		else {
			const Qdim = dataset.genome_info_list.findIndex(a => a.name == "Qdim");
			const np59_Cdimrid = dataset.genome_info_list.findIndex(a => a.name == "np59_Cdimrid");

			setPlotTitle(dataset.genome_info_list[Qdim].name, dataset.genome_info_list[np59_Cdimrid].name);

			add_RIP_marker(0, Qdim);
			add_RIP_marker(0, np59_Cdimrid);

			await [
				ref_type == "Cdim" ? "np59_Cdimrid" : (ref_type == "Qdim" ? "Qdim" : "never"),
				// "Qdim",
				// "Cdim",
			].reduce(async (prev_promise, ref_id) => {
				await prev_promise;

				// if (ref_id != "Qdim") {
				// 	throw Error("Qdim only");
				// }

				const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_id);
				if (ref_idx >= 0) {
					await load_dim2_rid1_dimrid_methratio(ref_id, "dim2", "BS", true);

					if (ref_id == "QM6a") {
						await load_BS_nBS_methylratio_by_RIP("QM6a", "QM6a dim2Δ veg", "ctrl", "data/BS_QM6a/QM6a_dim2_control.methratio.float32", 1.0, [Qdim, np59_Cdimrid]);
					}
					else if (ref_id == "Qdim") {
						// alert("QM6a dim2Δ >> np30");
						await load_BS_nBS_methylratio_by_RIP("Qdim", "QM6a dim2Δ veg", "ctrl", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_veg_ctrl.methratio.float32", 1.0, [Qdim, np59_Cdimrid]);
					}
					
					await load_methyl_group_by_RIP(ref_id, "Pd", "ctrl", 1.0, [Qdim, np59_Cdimrid]);
				}
			}, Promise.resolve());
		}
	}
	else if (ref_type == "Qrid" || ref_type == "Crid") {
		if (init_dataset) {
			RemoveGenomeHelper.rename("np30", "Qdim").rename("np31", "Crid").rename("np32", "Qrid");
			RemoveGenomeHelper.remove("Qdim").remove("np128_Qdimrid").remove("np59_Cdimrid");
			RemoveGenomeHelper.move_to_last("Crid").rebuild();
		}
		else {
			const Qrid = dataset.genome_info_list.findIndex(a => a.name == "Qrid");
			const Crid = dataset.genome_info_list.findIndex(a => a.name == "Crid");

			setPlotTitle(dataset.genome_info_list[Qrid].name, dataset.genome_info_list[Crid].name);

			add_RIP_marker(0, Qrid);
			add_RIP_marker(0, Crid);

			await [
				ref_type,
				// "Qrid",
				// "Crid",
			].reduce(async (prev_promise, ref_id) => {
				await prev_promise;

				const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_id);
				if (ref_idx >= 0) {
					/** @see {@link load_QM6a_CBS1_rid1_BS} */
					await load_dim2_rid1_dimrid_methratio(ref_id, "rid1", "BS", true);
					await load_BS_nBS_methylratio_by_RIP("Qrid", "QM6a rid1Δ veg", "ctrl", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_ctrl.methratio.float32", 1.0, [Qrid, Crid]); // 20220930
					await load_BS_nBS_methylratio_by_RIP("Crid", "CBS1-1 rid1Δ veg", "ctrl", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_ctrl.methratio.float32", 1.0, [Qrid, Crid]); // 20220930
					await load_methyl_group_by_RIP(ref_id, "rid1", "ctrl", 1.0, [Qrid, Crid]);
				}
			}, Promise.resolve());
		}
	}
	else if (ref_type == "np128_Qdimrid" || ref_type == "np59_Cdimrid") {
		if (init_dataset) {
			RemoveGenomeHelper.rename("np30", "Qdim").rename("np31", "Crid").rename("np32", "Qrid");
			RemoveGenomeHelper.remove("Qdim").remove("Crid").remove("Qrid");
			RemoveGenomeHelper.move_to_last("np59_Cdimrid").rebuild();
		}
		else {
			const np128_Qdimrid = dataset.genome_info_list.findIndex(a => a.name == "np128_Qdimrid");
			const np59_Cdimrid = dataset.genome_info_list.findIndex(a => a.name == "np59_Cdimrid");

			setPlotTitle(dataset.genome_info_list[np128_Qdimrid].name, dataset.genome_info_list[np59_Cdimrid].name);

			add_RIP_marker(0, np128_Qdimrid);
			add_RIP_marker(0, np59_Cdimrid);

			await [
				ref_type,
				// "np128_Qdimrid",
				// "np59_Cdimrid",
			].reduce(async (prev_promise, ref_id) => {
				await prev_promise;

				const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_id);
				if (ref_idx >= 0) {
					// add_RIP_marker(0, dataset.genome_info_list.findIndex(a => a.name == ref_id));
					
					// await load_dim2_rid1_dimrid_methratio(ref_id, "dimrid", "BS", true);

					// 20220928 ddr parental
					if (ref_id == "np128_Qdimrid") {
						await load_BS_nBS_methylratio_by_RIP(ref_id, "dim2Δ rid1Δ veg", "BS", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_veg_BS.methratio.float32", 1.0, [np128_Qdimrid, np59_Cdimrid]);
					}
					else if (ref_id == "np59_Cdimrid") {
						await load_BS_nBS_methylratio_by_RIP(ref_id, "dim2Δ rid1Δ veg", "BS", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_veg_BS.methratio.float32", 1.0, [np128_Qdimrid, np59_Cdimrid]);
					}

					await load_methratioFloat32_depthUint32_by_list(ref_id, "BS", [
						"ddr_D2_BS",
						"ddr_D4_BS",
						"ddr_D5_BS",
						"ddr_D6_BS",
						"ddr_D8_BS",
					]);// see load_dim2_rid1_dimrid_methratio

					// 20220928 ddr parental
					if (ref_id == "np128_Qdimrid") {
						await load_BS_nBS_methylratio_by_RIP(ref_id, "dim2Δ rid1Δ veg", "ctrl", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_veg_ctrl.methratio.float32", 1.0, [np128_Qdimrid, np59_Cdimrid]);
					}
					else if (ref_id == "np59_Cdimrid") {
						await load_BS_nBS_methylratio_by_RIP(ref_id, "dim2Δ rid1Δ veg", "ctrl", "data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_veg_ctrl.methratio.float32", 1.0, [np128_Qdimrid, np59_Cdimrid]);
					}

					await load_methyl_group_by_RIP(ref_id, "ddr", "ctrl", 1.0, [np128_Qdimrid, np59_Cdimrid]);
				}
			}, Promise.resolve());

			// remove_parental_are_not_C_G([np128_Qdimrid, np59_Cdimrid]);
		}
	}
	else {
		debugger;
	}

	/**
	 * @param {string} ref1_name
	 * @param {string} ref2_name
	 */
	function setPlotTitle(ref1_name, ref2_name) {
		const ref1_start = pos_ref1_uint32array[g_chrCoord.bp_start - 1];
		const ref1_end = pos_ref1_uint32array[g_chrCoord.bp_end - 1];
		const loc = `ch${viewerState.nChr}:${ref1_start}-${ref1_end}`;

		const d1 = ref_to_display_name(ref1_name);
		const d2 = ref_to_display_name(ref2_name);
		const title = `${d1} x ${d2} ${loc}`;
		viewerState.setPlotTitle(title);
		return title;
	}
 }

/**
 * @param {"WT"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid"} type
 * @see {@link load_5mc_debug}
 */
async function capture_all_5mC(type, fixedPlotHeight = 1266) {
	// viewerState.fixedPlotHeight = fixedPlotHeight;
	// viewerState.resizeCanvas();
	// await drawFrame();

	// document.querySelectorAll("canvas").forEach(a => a.height = 1266);

	if (["WT", "Qdim", "Cdim", "Qrid", "Crid", "np128_Qdimrid", "np59_Cdimrid"].includes(type)) {
		await load_methly_preset_20210924(type, true);

		await captureScreenAll(true, async () => {
			methyl_dataset_list.splice(0);// clear

			await load_methly_preset_20210924(type, false);
		});
	}
	else {
		throw ["WT", "Qdim", "Cdim", "Qrid", "Crid", "np128_Qdimrid", "np59_Cdimrid"];
	}
}

/**
 * TSETA:QM6a_CBS1-1_WT_dim2_rid1_dimrid
 * @param {"WT"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid"} type
 * @param {number} chr
 */
 async function capture_chr_5mC(type, chr) {
	if (["WT", "Qdim", "Cdim", "Qrid", "Crid", "np128_Qdimrid", "np59_Cdimrid"].includes(type)) {
		document.querySelectorAll("canvas").forEach(a => a.height = 1266);

		// await load_methly_preset_20210924(type, true);

		// await captureScreenAll(async () => {
		// 	methyl_dataset_list.splice(0);// clear

		// 	await load_methly_preset_20210924(type, false);
		// });
	}
	else {
		throw ["WT", "Qdim", "Cdim", "Qrid", "Crid", "np128_Qdimrid", "np59_Cdimrid"];
	}
}

/**
 * TSETA:Tr_ssDNA_parent
 */
async function capture_all_ssDNA(single = true, start_chr, end_chr, fork_feat, download) {
	if (single) {
		dataset.mode = "single";
		hide_all_marker();
		viewerState.$display_SNP_density = false;
	}
	else {
		await load_methly_preset_20210924("WT", true);
		viewerState.$display_SNP_density = false;

		dataset.progeny_list.forEach(a => {
			RemoveGenomeHelper.rename(a, a.replace(/_f$/, "♀"));
			RemoveGenomeHelper.rename(a, a.replace(/_m$/, "♂"));
		});
	}

	// 20211220
	viewerState.plot_list.push(new PlotInfo({
		title: "",
		row_height: 2,// px = row_height * viewerState.seg_row_height
		row_separate: viewerState.seg_row_separate,
		rowspan: 1,
		p_fontSize: 100,
		func: (info) => {
			const ctx = main_ctx;
			const height = info.row_height * viewerState.seg_row_height;

			ctx.save();
			
			try {
				const text_list = [
					"red = p-value ≦ 0.05",
					"blue = 0.05 < p-value ≦ 0.1",
				];
				const a = Math.max(...text_list.map(a => ctx.measureText(a)).map(a => a.actualBoundingBoxAscent + a.actualBoundingBoxDescent));
				
				ctx.translate(0, info.row_separate);

				// const fontSize = info.p_fontSize && info.p_fontSize != 100 ? info.p_fontSize + "%" : "";
				// ctx.font = `${fontSize} Arial`;
				ctx.font = `${Math.ceil(viewerState.seg_row_height * 0.8)}px Arial`;
				
				// main
				ctx.fillStyle = "red";
				ctx.fillText(text_list[0], 0, height * 1 / 4);

				ctx.fillStyle = "blue";
				ctx.fillText(text_list[1], 0, height * 3 / 4);
			}
			finally {
				ctx.restore();
			}

			// auto extend canvas height
			ctx.translate(0, height);
		},
	}));

	// methyl_dataset_list.splice(0);

	await ssDNA_20211014(true, { row_height: 3, fork_feat: fork_feat }, true);// clear and rename

	// document.querySelector("#el_input_chr").value = 7;
	// document.querySelector("#el_input_chr").oninput(null);

	const canvas_list = await captureScreenAll(true, async () => {
		await promise_load_task;

		await ssDNA_20211014(true, { row_height: 3, fork_feat: fork_feat }, false);// load chr data

		viewerState.setPlotTitle(`ssDNA ch${viewerState.nChr}`);

		console.log(viewerState.nChr);
	}, start_chr, end_chr, download, {
		rep: 1,
		// rep: 2,
		// beforeCapture: async function beforeCapture(rep, nChr) {
		// 	if (rep == 1) {
		// 		methyl_dataset_list.forEach(a => a.hide = !a.hide);
		// 	}
		// },
	});

	const key = `screenshot_${new Date().getTime()}`;
	window[key] = canvas_list;

	console.log("save screenshot", key);

	return canvas_list;
}

async function compare_ssDNA_peak_seq(reload) {
	const ref_idx = 0;
	const ref_name = dataset.genome_info_list[ref_idx].name;

	// ssDNA_20211014.lg = [4];
	// ssDNA_20211014.file_extname = "bed";
	// if (reload || methyl_dataset_list.length == 0) {
	// 	await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, false);
	// }
	
	await Promise.all(dataset.genome_info_list[ref_idx].chr_list.map(async (ch_info, chr_idx) => {
		await load_multi_alignment_results(chr_idx);
	}));
	
	const ref1_chr_seq = window.ch ?? await Promise.all(dataset.genome_info_list[ref_idx].chr_list.map(async info => parseFasta(await fetchData("../" + info.path, "text"))));
	window.ch = ref1_chr_seq;
		
	const row_list = methyl_dataset_list
		.filter(a => !a.hide)
		.filter(a => a.tags[0] == "rad51_vs_sae2" || a.tags[0] == "spo11rad51_vs_spo11sae2");

	row_list.forEach((cfg, data_row_idx) => {
		const out_fa = {};

		dataset.genome_info_list[ref_idx].chr_list.forEach((ch_info, chr_idx) => {
			const nChr = chr_idx + 1;
			const sChr = ch_info.chr;
			const seq = ref1_chr_seq[chr_idx][sChr];
			
			if (cfg.data[chr_idx] == null) {
				console.warn(chr_idx + 1);
				return;
			}

			cfg.data[chr_idx].map(peak => {
				const start = peak.start;
				const end = peak.end;
				const peak_faName = `${ref_name}_ch${nChr}_${start}_${end}`;
				out_fa[peak_faName] = seq.slice(start - 1, end);
			});
		});

		const tt = Object.entries(out_fa).map(a => ">" + a[0] + "\n" + a[1] + "\n").join("\n");
		downloadTextFile(`${ref_name}_ssDNA_${cfg.tags.join("_")}_peak.fa`, tt);
	});

	// $(`<div>${methyl_dataset_list[0].name}</div>`).text().replace(/[\r\n\t]/g, "")
}

/**
 * @param {number} genomeIdx
 * @param {string} url
 * @param {(aln: BlastnCoord, idx: number) => boolean} filter
 */
async function _load_view_blastnResult(genomeIdx, url, filter) {
	const results = BlastnCoord.load_tab(await fetchData(url, "text")).tab.filter(filter);

	// addViewBlastnResult(genomeIdx, results, 0, 0);

	const ref_pos = genomeIdx >= 0 ? make_ref_pos_map_list()[genomeIdx] : null;
	const aaa = new Uint8Array(seq_list[0].length);
	
	results.forEach(aln => {
		let ss, se;
		if (ref_pos) {
			ss = ref_pos[aln.qstart - 1];
			se = ref_pos[aln.qend - 1];
			[ss, se] = [ss, se].sort((a, b) => a - b);
			viewerState.pushRangeMarker(ss, se);
		}
		else {
			ss = aln.qstart;
			se = aln.qend;
		}
		
		results.forEach(a => {
			aaa.fill(1, ss - 1, se);
		});
	});

	let range = {
		start: -1,
		end: -1,
	};
	/** @type {range[]} */
	const lst = [];
	range = null;
	aaa.forEach((flag, idx) => {
		if (flag) {
			if (range == null) {
				range = {
					start: idx + 1,
					end: idx + 1,
				}
			}
			else {
				range.end = idx + 1;
			}
		}
		else if (range) {
			lst.push(range);
			range = null;
		}
	});

	return lst;
}

async function ssDNA_all_fa() {
	const cmp_list = [
		"rad51_vs_sae2",
		"spo11rad51_vs_sae2",
		// "spo11rad51_vs_spo11sae2",
		// "spo11sae2_vs_sae2",
	];

	await [
		"QM6a",
		"CBS1-1",
	].reduce(async (prev_promise, ref_id) => {
		await prev_promise;
		
		const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_id);
		const ref_seq_list = dataset.genome_info_list[ref_idx].chr_list.map((ch_info, chr_idx) => {
			const sChr = ch_info.chr;
			const seq = dataset.results[chr_idx][sChr];
			return seq.replace(/-/g, "");
		});

		await cmp_list.reduce(async (prev_promise, cmp, cmp_idx) => {
			await prev_promise;

			const cfg = methyl_dataset_list.find(a => a.hide && a.ref == ref_id && a.tags[0] == cmp);
			
			const out_fa = {};
			cfg.data.forEach((rr, chr_idx) => {
				const nChr = chr_idx + 1;
				const chr_seq = ref_seq_list[chr_idx];
				rr.forEach(point => {
					const start = point.start;
					const end = point.end;
					const point_faName = `${ref_id}_ch${nChr}_${start}_${end}`;
					const point_seq = chr_seq.slice(start, end);
					// if (point_seq < 10) {
					// 	console.log(cfg.name);
					// 	debugger;
					// }
					out_fa[point_faName] = point_seq;
				});
			});

			const tt = Object.entries(out_fa).map(a => ">" + a[0] + "\n" + a[1] + "\n").join("\n");
			downloadTextFile(`${ref_id}_${cmp}_peak_seq.fa`, tt);
			
			await timeout(1000);
		}, Promise.resolve());
		
		await timeout(1000);
	}, Promise.resolve());

	/** @param {number} ms */
	function timeout(ms) {
		return new Promise(resolve => setTimeout(resolve, ms));
	}
}

// QM6a overlap 2 CBS1-1
// ssDNA_peak_venn.output[0][0].filter(a => a.id == "PCG_ch1_178")

/**
 * ext0, ext50, ext0_all_spo11rad51
 * @param {boolean} all_spo11rad51
 * @param {string} prefix
 * @param {number} extends_length
 * @param {boolean} download
 * @param {boolean} debug
 */
async function ssDNA_type_fa(all_spo11rad51, prefix, extends_length, download = false, debug = true) {
	const ref_list = [
		"QM6a",
		"CBS1-1",
	];
	const type_list = [
		"PCG",
		"PCG+intergenic",
		"intergenic",
		"ATisland",
	];
	const cmp_list = [
		"rad51_vs_sae2",
		"spo11rad51_vs_sae2",
		// "spo11rad51_vs_spo11sae2",
		// "spo11sae2_vs_sae2",
	];
	
	const fa_peak_map = await _ssDNA_type_fa(all_spo11rad51, extends_length, ref_list, type_list, cmp_list, debug);

	const {
		out_map,
		peak_map,
	} = fa_peak_map;
	
	if (debug == false) {
		if (download) {
			cmp_list.forEach((cmp, cmp_idx) => {
				type_list.forEach(type => {
					// ref_list.forEach(ref_id => {
					// 	const ffname = `${ref_id}_${cmp}_${type}_peak_ext${extends_length}_seq.fa`;
					// 	const out_fa = out_map[cmp][type][ref_id];
					// 	saveFa(ffname, out_fa);
					// });
					
					const peak_sum = new Set([
						peak_map[cmp][type][ref_list[0]],
						peak_map[cmp][type][ref_list[1]]
					].map(b => b.map(a => [...a.peakList_map.values()].flat(1))).flat(2).map(a => a.$peak)).size;

					const ffa = ref_list.reduce((fa, ref_id) => {
						return Object.assign(fa,  out_map[cmp][type][ref_id]);
					}, {});
					const ffname = [
						prefix,
						`${cmp}_${type}_peak${peak_sum}_ext${extends_length}`,
						"seq.fa",
					].filter(a => a).join("_");

					saveFa(ffname, ffa);
				});
			});
		}
		return fa_peak_map;
	}
	else {
		console.log(fa_peak_map);
		window.__fa_peak_map = fa_peak_map;
	}

	/**
	 * @param {string} ffname
	 * @param {*} out_fa
	 */
	function saveFa(ffname, out_fa) {
		const tt = Object.entries(out_fa).map(a => ">" + a[0] + "\n" + a[1] + "\n").join("\n");
		downloadTextFile(ffname, tt);
	}
}
/**
 * @param {boolean} all_spo11rad51
 * @param {number} extends_length
 * @param {string[]} ref_list
 * @param {string[]} type_list
 * @param {string[]} cmp_list
 * @param {boolean} debug
 */
async function _ssDNA_type_fa(all_spo11rad51, extends_length, ref_list, type_list, cmp_list, debug = true) {
	console.log({
		extends_length
	});
	
	/**
	 * @type {{ [ref_cmp: string]: (groups: string[]) => any }}
	 */
	const filter_map = {
		"QM6a_rad51_vs_sae2": make_exclude_fn([
			  "QM6a_spo11rad51_vs_sae2",
			"CBS1-1_spo11rad51_vs_sae2",
		]),
		"CBS1-1_rad51_vs_sae2": make_exclude_fn([
			  "QM6a_spo11rad51_vs_sae2",
			"CBS1-1_spo11rad51_vs_sae2",
		]),
	};
	if (all_spo11rad51 == false) {
		filter_map["QM6a_spo11rad51_vs_sae2"] = make_exclude_fn([
			  "QM6a_rad51_vs_sae2",
			"CBS1-1_rad51_vs_sae2",
		]);
		filter_map["CBS1-1_spo11rad51_vs_sae2"] = make_exclude_fn([
			  "QM6a_rad51_vs_sae2",
			"CBS1-1_rad51_vs_sae2"
		]);
	}
	else if (all_spo11rad51) {
		filter_map["QM6a_spo11rad51_vs_sae2"] = make_some_fn([
			  "QM6a_spo11rad51_vs_sae2",
			"CBS1-1_spo11rad51_vs_sae2",
		]);
		filter_map["CBS1-1_spo11rad51_vs_sae2"] = make_some_fn([
			  "QM6a_spo11rad51_vs_sae2",
			"CBS1-1_spo11rad51_vs_sae2",
		]);
	}
	function make_some_fn(some_list) {
		/** @param {string[]} groups */
		return function fn(groups) {
			return some_list.some(some => groups.includes(some));
		};
	}
	function make_exclude_fn(exclude_list) {
		/** @param {string[]} groups */
		return function fn(groups) {
			return exclude_list.every(exclude => groups.includes(exclude) == false);
		};
	}

	/**
	 * @type {{ [refName: string]: string[]; }}
	 */
	const seq_map = {};
	if (!debug) {
		ref_list.forEach(ref_id => {
			seq_map[ref_id] = [];
			
			const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_id);
			dataset.genome_info_list[ref_idx].chr_list.map((ch_info, chr_idx) => {
				const sChr = ch_info.chr;
				const seq = dataset.results[chr_idx][sChr].replace(/-/g, "");
				seq_map[ref_id][chr_idx] = seq;
			});
		});
	}
	
	/** @type {{ [cmp: string]: { [type: string]: { [ref_id: string]: {}; }; }; }} */
	const out_map = {};
	/** @type {{ [cmp: string]: { [type: string]: { [ref_id: string]: any[]; }; }; }} */
	const peak_map = {};

	// init
	cmp_list.forEach((cmp, cmp_idx) => {
		out_map[cmp] = {};
		peak_map[cmp] = {};
		type_list.forEach(type => {
			out_map[cmp][type] = {};
			peak_map[cmp][type] = {};
			ref_list.forEach(ref_id => {
				out_map[cmp][type][ref_id] = null;
				peak_map[cmp][type][ref_id] = null;
			});
		});
	});

	cmp_list.forEach((cmp, cmp_idx) => {
		type_list.forEach(type => {
			ref_list.forEach(ref_id => {
				const cfg = methyl_dataset_list.find(a => !a.hide && a.ref == ref_id && a.tags[0] == cmp && a.tags[1] == type);

				const out_fa = {};
				const info_list = [];
				
				cfg.data.forEach((rr, chr_idx) => {
					const nChr = chr_idx + 1;
					const chr_seq = debug ? null : seq_map[ref_id][chr_idx];
					rr.forEach(point => {
						if (!point.$peak) {
							// console.log(cfg);
							// debugger
							return;
						}
						else if (point.$peak.id.startsWith(type)) {
							const cbfn_filter = filter_map[`${ref_id}_${cmp}`];
							if (filter_map == null ||
								// exclude_list.every(exclude => point.$peak.groups.includes(exclude) == false)
								cbfn_filter(point.$peak.groups)
							) {
								if (debug == false) {
									const start = Math.max(0, point.start - extends_length);
									const end = point.end + extends_length;
									// 20211109
									const point_faName = `${ref_id}_ch${nChr}_${start}_${end}_${cmp}_${type} ${point.$peak.id} extend(${extends_length})`;
									const point_seq = chr_seq.slice(start, end);
									// if (point_seq < 10) {
									// 	console.log(cfg.name);
									// 	debugger;
									// }
									out_fa[point_faName] = point_seq;
								}
								info_list.push(point.$peak);
							}
						}
					});
				});

				// const ffname = `${ref_id}_${cmp}_${type}_peak_ext${extends_length}_seq.fa`;
				// // if (ffname == "CBS1-1_rad51_vs_sae2_PCG_peak_ext0_seq.fa") {
				// // 	debugger
				// // }
				// if (debug == false) {
				// 	const tt = Object.entries(out_fa).map(a => ">" + a[0] + "\n" + a[1] + "\n").join("\n");
				// 	downloadTextFile(ffname, tt);
				// }
				// // else {
				// // 	console.log(ffname, out_fa.llist);
				// // 	window.__out_fa = out_fa.lilst;
				// // }

				out_map[cmp][type][ref_id] = out_fa;
				peak_map[cmp][type][ref_id] = info_list;
			});
		});
	});

	return {
		out_map,
		peak_map,
	};

	/** @param {number} ms */
	function timeout(ms) {
		return new Promise(resolve => setTimeout(resolve, ms));
	}
}
function ssDNA_cmp_type_seq_motifs() {
	const c1_list = [];
	const c2_list = [];
	[
		"QM6a",
		// "CBS1-1",
	].forEach(ref_id => {
		[
			"PCG",
			"intergenic",
			"PCG+intergenic",
			"ATisland",
		].forEach(type => {
			[
				"rad51_vs_sae2",
				"spo11rad51_vs_sae2",
				// "spo11rad51_vs_spo11sae2",
				// "spo11sae2_vs_sae2",
			].forEach(cmp => {
				type = type.replace(/&/g, "+");
				const pp_name = `${ref_id}_${cmp}_${type}_peak_seq`;
				const fname = `${pp_name}.fa`;
				const cmd1 = `findMotifs.pl ${fname} fasta ./homer_${pp_name}/`;
				c1_list.push(cmd1);
				const cmd2 = `meme ${fname} -dna -oc ./meme_${pp_name} -nostatus -time 14400 -mod zoops -nmotifs 100 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0`;
				c2_list.push(cmd2);
			});
		});
	});
	const tttt = c1_list.join(" &\n") + "\n" + c2_list.join(" &\n");
	document.querySelector("textarea").value = tttt;
	return tttt;
}

async function ssDNA_length_table() {
	const cmp_list = [
		"rad51_vs_sae2",
		"spo11rad51_vs_sae2",
		// "spo11rad51_vs_spo11sae2",
		// "spo11sae2_vs_sae2",
	];

	return await [
		"PCG",
		"intergenic",
		"PCG+intergenic",
		"ATisland",
	].reduce(async (prev_promise, type) => {
		await prev_promise;
		await [
			"QM6a",
			"CBS1-1",
		].reduce(async (prev_promise, ref_id) => {
			await prev_promise;
			await cmp_list.reduce(async (prev_promise, cmp, cmp_idx) => {
				await prev_promise;
				const cfg = methyl_dataset_list.find(a => !a.hide && a.ref == ref_id && a.tags[0] == cmp && a.tags[1] == type);
				cfg.data.forEach((peaks, ch_idx) => {
					peaks.forEach(peak => {
						peak.nChr = ch_idx + 1;
					});
				});
				//20211105
				download(`${type}_${ref_id}_${cmp}_peak_length.txt`, cfg.data.flat(1));
				await timeout(1000);
			}, Promise.resolve());
			await timeout(1000);
		}, Promise.resolve());
		await timeout(1000);
	}, Promise.resolve());

	/** @param {number} ms */
	function timeout(ms) {
		return new Promise(resolve => setTimeout(resolve, ms));
	}

	/**
	 * @param {string} file_name
	 * @param {Partial<module_Methyl_ratioData>[]} peaks
	 */
	function download(file_name, peaks) {
		const lines = peaks.map(peak => [
			peak.nChr,
			peak.start,
			peak.end,
			peak.end - peak.start,
		].join("\t")).join("\n");
	
		const header = [
			"ch",
			"start",
			"end",
			"length",
			"\n"
		].join("\t");
	
		downloadTextFile(file_name, header + lines);
	}
}

example_ssDNA_peak_venn.loaded = false;
/**
 * @param {number} nChr
 * @param {number} min_pvalue
 * @param {number} min_val
 */
async function example_ssDNA_peak_venn(nChr, min_pvalue, min_val) {
	if (nChr == null) {
		throw new Error("nChr");
	}
	if (min_val == null) {
		throw new Error("min_val");
	}
	if (min_pvalue == null) {
		throw new Error("min_pvalue");
	}
	if (example_ssDNA_peak_venn.loaded == false) {
		[...dataset.progeny_list].forEach(a => RemoveGenomeHelper.remove(a));
		RemoveGenomeHelper.rebuild();

		document.querySelector("#el_input_chr").value = nChr
		document.querySelector("#el_input_chr").oninput(null)
		await promise_load_task
		
		await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
			await load_multi_alignment_results(chr_idx);
		}));
		
		ssDNA_20211014.file_extname = "txt"
		ssDNA_20211014.lg = [2];
		ssDNA_20211014.show_pvalue = false;
		await ssDNA_20211014(true, { row_height: 3, fork_feat: true, load_QM6a: false, load_CBS1: false }, false);
		viewerState.resizeCanvas()
		await drawFrame()
	}
	
	return ssDNA_peak_venn(min_pvalue, min_val);
}

/**
 * 20211224
 * @param {number} min_pvalue ks_test_pvalue ≦ min_pvalue
 * @param {number} min_val (default = 2) coverage a - b
 * @param {boolean} bUseAvg
 * @param {(peak: Partial<module_Methyl_ratioData>) => boolean} cond
 * @param {number} min_cons_len
 * 
 * venn 的數量會少於原本的peak
 * p1 p2 會被當作同一個
 * Qr    **p1**      **p2**
 * Qsr      ***********
 */
function ssDNA_peak_venn(min_pvalue, min_val, bUseAvg, cond, min_cons_len = 40) {
	if (min_pvalue == null) {
		throw new Error("min_pvalue");
	}
	if (methyl_dataset_list.every(a => a.tags[1] == null)) {
		throw new Error("require ssDNA_20211014.fork_feat");
	}
	console.log({
		min_cons_len
	});

	/** @type {any[][]} */
	const output_info_table = [];
	
	const cmp_list = [
		"rad51_vs_sae2",
		"spo11rad51_vs_sae2",
		// "spo11rad51_vs_spo11sae2",
		// "spo11sae2_vs_sae2",
	];

	const ref1_id = "QM6a";
	const ref2_id = "CBS1-1";

	const ref_name_map = {
		[ref1_id]: "Q",
		[ref2_id]: "C",
	};

	const output = [
		"PCG",
		"PCG+intergenic",
		"intergenic",
		"ATisland",
	].map(type => {
		const ref1_idx = dataset.genome_info_list.findIndex(a => a.name == "QM6a");
		const ref2_idx = dataset.genome_info_list.findIndex(a => a.name == "CBS1-1");

		const group_name_list = [ref1_id, ref2_id].map(a => cmp_list.map(b => [a, b].join("_"))).flat(1);// 20211111

		class PeakSets_Base {
			id = "";
			start = 0;
			end = 0;
			
			/** @type {Map<string, Partial<module_Methyl_ratioData>[]>} */
			peakList_map = new Map([ref1_id, ref2_id].map(a => cmp_list.map(b => [[a, b].join("_"), []])).flat(1));
			
			/** @type {string[]} venn */
			groups = Array(group_name_list.length);// 20211111 // fill hole
		}

		/** @type {PeakSets_Base[][]} */
		const lalala = dataset.genome_info_list[0].chr_list.map((_, chr_idx) => {
			const nChr = chr_idx + 1;
			const ref_to_pos_list = _make_chr_ref_pos_map_list(nChr);// 20221111
			const ref1_to_pos_map = ref_to_pos_list[ref1_idx];
			const ref2_to_pos_map = ref_to_pos_list[ref2_idx];

			const _pos_to_ref_list = _make_chr_pos_ref_map_list(nChr);
			const _pos_to_ref1_map = _pos_to_ref_list[ref1_idx];
			const _pos_to_ref2_map = _pos_to_ref_list[ref2_idx];

			class PeakSets extends PeakSets_Base {
				id = `${type}_ch${nChr}_${PeakSets.list.length + 1}`;
				
				get ref1_start() {
					return _pos_to_ref1_map[this.start];
				}
				get ref1_end() {
					return _pos_to_ref1_map[this.end];
				}
				get ref1_length() {
					const [x1, x2] = [this.ref1_start, this.ref1_end];
					// if (x1 == x2 || x1 == 0 || x2 == 0) {
					// 	debugger;
					// }
					return x2 - x1 + 1;
				}

				get ref2_start() {
					return _pos_to_ref2_map[this.start];
				}
				get ref2_end() {
					return _pos_to_ref2_map[this.end];
				}
				get ref2_length() {
					const [y1, y2] = [this.ref2_start, this.ref2_end];
					// if (y1 == y2 || y1 == 0 || y2 == 0) {
					// 	debugger;
					// }
					return this.ref2_end - this.ref2_start + 1;
				}

				/** @type {PeakSets[]} */
				static list = [];

				/**
				 * @param {Partial<module_Methyl_ratioData>} point
				 * @param {string} group
				 * @param {number} group_idx
				 */
				static add_peak(point, group, group_idx) {
					// if (nChr == 7 && point.start <= 2996565 && 2996565 <= point.end) {
					// 	debugger
					// }
					if (
						// @see peak.ks_test_pvalue = 0.1;// fill
						point.ks_test_pvalue <= min_pvalue &&
						point.value >= min_val &&
						(cond == null || cond(point))
					) {
						const found = PeakSets.list.find(ps => point.$start >= ps.start && point.$end <= ps.end);
						if (found != null) {
							// if (found.ref1_length == 0 || found.ref2_length == 0) {
							// 	debugger
							// }
							point.$peak = found;

							found.start = Math.min(found.start, point.$start);
							found.end = Math.max(found.end, point.$end);
							found.groups[group_idx] = group;
							found.peakList_map.get(group).push(point);
						}
						else {
							const ps = new PeakSets();
							point.$peak = ps;

							ps.start = point.$start;
							ps.end = point.$end;
							ps.peakList_map.get(group).push(point);
							ps.groups[group_idx] = group;
							PeakSets.list.push(ps);
						}
					}
				}
			}

			/**
			 * enum all peak
			 * @type {Partial<module_Methyl_ratioData>[]}
			 */
			let all_peak_list = [];// not result

			cmp_list.forEach((cmp, cmp_idx) => {
				function is_use_avg(data_name) {
					return (
						(bUseAvg && data_name.indexOf("≥2 samples") < 0) ||
						(!bUseAvg && data_name.indexOf("≥2 samples") >= 0)
					);
				}

				// ?? remove if chr is empty ??

				const q_ = methyl_dataset_list.find(a => !a.hide && is_use_avg(a.name) && a.ref == ref1_id && a.tags[0] == cmp && a.tags[1] == type);
				const c_ = methyl_dataset_list.find(a => !a.hide && is_use_avg(a.name) && a.ref == ref2_id && a.tags[0] == cmp && a.tags[1] == type);

				// cache position
				const q_list = q_.data[chr_idx].filter(point => {
					const pos_start = ref1_to_pos_map[point.start];
					const pos_end = ref1_to_pos_map[point.end];
					const ref2_start = _pos_to_ref2_map[pos_start];
					const ref2_end = _pos_to_ref2_map[pos_end];
					const ref2_len = ref2_end - ref2_start;
					delete point.$peak;
					
					point.$_len = ref2_len;
					if (ref2_len >= min_cons_len) {
						point.$group_idx = cmp_idx;
						point.$group = [ref1_id, cmp].join("_");
						point.$start = pos_start;
						point.$end = pos_end;
						return true;
					}
				});
				const c_list = c_.data[chr_idx].filter(point => {
					const pos_start = ref2_to_pos_map[point.start];
					const pos_end = ref2_to_pos_map[point.end];
					const ref1_start = _pos_to_ref1_map[pos_start];
					const ref1_end = _pos_to_ref1_map[pos_end];
					const ref1_len = ref1_end - ref1_start;
					delete point.$peak;
					
					point.$_len = ref1_len;
					if (ref1_len >= min_cons_len) {
						point.$group_idx = cmp_list.length + cmp_idx;
						point.$group = [ref2_id, cmp].join("_");
						point.$start = pos_start;
						point.$end = pos_end;
						return true;
					}
				});
				
				all_peak_list = all_peak_list.concat(q_list, c_list);
			});
			
			// add peak to group
			all_peak_list.sort((a, b) => a.$start - b.$start).forEach(point => {
				PeakSets.add_peak(point, point.$group, point.$group_idx);
			});

			return PeakSets.list.filter(ps => {
				if (type == "PCG") {
					// peakSet 'inside PCG' or 'overlap PCG+intergenic'
					const aaa = [...ps.peakList_map.values()].flat(1);
					// always true
					const bbb = aaa.every(point => point.in_PCG > 0);
					if (!bbb) {
						console.log(aaa.length);
						debugger;
					}
					return bbb;
				}
				else {
					return true;
				}
			});
		});//each chr

		const download = true;
		if (!download) {
			// /** @type {Set<PeakSets_Base>} */
			// const _all_peak = new Set();

			// cmp_list.forEach(cmp => {
			// 	const q_ = methyl_dataset_list.find(a => !a.hide && a.ref == ref1_id && a.tags[0] == cmp && a.tags[1] == type);
			// 	const c_ = methyl_dataset_list.find(a => !a.hide && a.ref == ref2_id && a.tags[0] == cmp && a.tags[1] == type);
			// 	const qq = q_.data.map(mm => mm.map(pt => _all_peak.add(pt.$peak)));
			// 	const cc = c_.data.map(mm => mm.map(pt => _all_peak.add(pt.$peak)));
			// });

			// const all_peak = lalala.flat(1);
			// all_peak.every(peak => {
			// 	if (_all_peak.has(peak)) {
			// 		return true;
			// 	}
			// 	else {
			// 		console.error("lost peak", peak);
			// 		debugger;
			// 	}
			// });

			// downloadTable(`${type}_peaks.txt`, [...all_peak.values()]);
		}
		else {
			cmp_list.forEach(cmp => {
				const q_ = methyl_dataset_list.find(a => !a.hide && a.ref == ref1_id && a.tags[0] == cmp && a.tags[1] == type);
				const c_ = methyl_dataset_list.find(a => !a.hide && a.ref == ref2_id && a.tags[0] == cmp && a.tags[1] == type);

				const qq = q_.data.map(mm => distinct(mm.filter(pt => pt.$_len >= min_cons_len).map(pt => pt.$peak))).flat(1);
				const cc = c_.data.map(mm => distinct(mm.filter(pt => pt.$_len >= min_cons_len).map(pt => pt.$peak))).flat(1);

				/**
				 * @template T
				 * @param {Array<T>} list
				 */
				function distinct(list) {
					return Array.from(new Set(list));
				}

				downloadTable(`${type}_${ref1_id}_${cmp}.txt`, qq);
				downloadTable(`${type}_${ref2_id}_${cmp}.txt`, cc);
			});
		}

		/**
		 * @param {string} file_name
		 * @param {PeakSets_Base[]} peaks
		 */
		function downloadTable(file_name, peaks) {
			const filtered_peaks = peaks.filter(a => a != null);

			output_info_table.push([file_name, filtered_peaks]);

			const cmp_venn_strMap = [];
			for (let i = 0; i < (cmp_list.length - 1); ++i) {
				cmp_venn_strMap[i] = "";
				cmp_list.map((cmp, cmp_idx) => {
					if (i & (1 << cmp_idx)) {
						cmp_venn_strMap[i] = [cmp_venn_strMap[i], cmp].join(",");
					}
				});
			}
			
			const lines = filtered_peaks.map($peak => [
				$peak.id,
				//
				$peak.ref1_start,
				$peak.ref1_end,
				$peak.ref1_length,
				//
				$peak.ref2_start,
				$peak.ref2_end,
				$peak.ref2_length,
				//
				...Array.from($peak.groups),//20211111 // fill holes
			].join("\t")).join("\n");

			const header = [
				"id",
				//
				`${ref1_id} start`,
				`${ref1_id} end`,
				`${ref1_id} length`,
				//
				`${ref2_id} start`,
				`${ref2_id} end`,
				`${ref2_id} length`,
				//
				...group_name_list,//20211111 // R.lang read.table: read.table columns = column names
			].join("\t") + "\n";

			downloadTextFile(file_name, header + lines);
		}

		return lalala;
	});

	ssDNA_peak_venn.output = output;

	ssDNA_peak_venn.output_info_table = output_info_table;
	
	return output;
}

function ssDNA_peak_Rvenn() {
	/**
Rscript venn.R 2021112_1038.PCG \
 "Q rad51_vs_sae2" \
 "Q spo11rad51_vs_sae2" \
 "C rad51_vs_sae2" \
 "C spo11rad51_vs_sae2" \
 peaks/PCG_QM6a_rad51_vs_sae2.txt \
 peaks/PCG_QM6a_spo11rad51_vs_sae2.txt \
 peaks/PCG_CBS1-1_rad51_vs_sae2.txt \
 peaks/PCG_CBS1-1_spo11rad51_vs_sae2.txt
Rscript venn.R 2021112_1043.intergenic "Q rad51_vs_sae2" "Q spo11rad51_vs_sae2" "C rad51_vs_sae2" "C spo11rad51_vs_sae2" peaks/intergenic_QM6a_rad51_vs_sae2.txt peaks/intergenic_QM6a_spo11rad51_vs_sae2.txt peaks/intergenic_CBS1-1_rad51_vs_sae2.txt peaks/intergenic_CBS1-1_spo11rad51_vs_sae2.txt
Rscript venn.R 2021112_1043.ATisland "Q rad51_vs_sae2" "Q spo11rad51_vs_sae2" "C rad51_vs_sae2" "C spo11rad51_vs_sae2" peaks/ATisland_QM6a_rad51_vs_sae2.txt peaks/ATisland_QM6a_spo11rad51_vs_sae2.txt peaks/ATisland_CBS1-1_rad51_vs_sae2.txt peaks/ATisland_CBS1-1_spo11rad51_vs_sae2.txt
	 */
	const d = new Date();
	const out = [
		[d.getFullYear(), d.getMonth() + 1, d.getDate()].map(a => a.toString().padStart(2, "0")).join(""),
		[d.getHours(), d.getMinutes()].map(a => a.toString().padStart(2, "0")).join(""),
	].join("_");
	
	const cmpList = [
		"rad51_vs_sae2",
		"spo11rad51_vs_sae2",
		// "spo11rad51_vs_spo11sae2",
		// "spo11sae2_vs_sae2",
	];

	const cmd_list = [
		"PCG",
		"intergenic",
		"PCG+intergenic",
		"ATisland",
	].map(type => {
		const ref1_id = "QM6a";
		const ref2_id = "CBS1-1";
		const ref_list = [
			ref1_id,
			ref2_id,
		];

		const ref_name_map = {
			[ref1_id]: "Q",
			[ref2_id]: "C",
			// [ref1_id]: ref1_id,
			// [ref2_id]: ref2_id,
		};

		const names = ref_list.map(ref => {
			return cmpList.map(cmp => `"${ref_name_map[ref]} ${cmp}"`);
		}).flat(1).join(" ");

		const list = ref_list.map(ref => {
			return cmpList.map(cmp => `"peaks/${type}_${ref}_${cmp}.txt"`);
		}).flat(1).join(" ");

		const cmd = `Rscript venn.R "${out}.${type}" ${names} ${list}`;

		return cmd;
	});

	document.querySelector("textarea").value = cmd_list.join("\n");

	return cmd_list;
}

async function ssDNA_avg_20211224_venn_peak_sum() {
	{
		const nChr = 2;
		
		[...dataset.progeny_list].forEach(a => RemoveGenomeHelper.remove(a));
		RemoveGenomeHelper.rebuild();

		document.querySelector("#el_input_chr").value = nChr
		document.querySelector("#el_input_chr").oninput(null)
		await promise_load_task
		
		await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
			await load_multi_alignment_results(chr_idx);
		}));
		
		ssDNA_20211014.file_extname = "txt"
		ssDNA_20211014.lg = [2];
		ssDNA_20211014.show_pvalue = false;
		await ssDNA_20211014(true, { row_height: 3, fork_feat: true, load_QM6a: true, load_CBS1: true, load_dup3: false }, false);
		viewerState.resizeCanvas()
		await drawFrame()
	}
	// ssDNA_peak_venn(0.05, 0, true, null);

	ssDNA_peak_link_avg_and_2sample("QM6a");
	ssDNA_peak_link_avg_and_2sample("CBS1-1");
	
	ssDNA_peak_venn(0.05, 0, true, peak => peak.trip?.value >= 2);

	/**
	 * QM6a  rad51Δ vs sae2Δ scaled coverage (cutoff 2) average
	 * methyl_dataset_list.filter(a => a.name.indexOf("≥") < 0)[0].data[1].filter(a => a.value >= 0 && a.ks_test_pvalue <= 0.05).filter(a => a.$_len < 40 && !a.AT_island?.size).length == 95
	 */
	
	/**
	 * @see {@link ssDNA_peak_venn} venn 的數量會少於原本的peak
	 */
	
	// const bbb = ssDNA_peak_venn.output_info_table.map(a => [a[0], a[1].length].join("\t")).join("\n");
	// document.querySelector("textarea").value = aaa;

	// const aaa = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") < 0 && a.tags[1] != null).map(a => {
	// 	const line_head = [a.ref, ...a.tags].join(" ");
	// 	const cc = a.data.reduce((acc, peaks) => acc + peaks.reduce((acc2, peak) => acc2 + (peak.ks_test_pvalue <= 0.05 && peak.$_len >= 40 && peak.value >= 0 ? 1 : 0), 0), 0);
	// 	return [line_head, cc].join("\t");
	// }).join("\n");
	
	// document.querySelector("textarea").value = aaa + "\n\n\n\n" + bbb;
}

async function ssDNA_3_20211224_venn_peak_sum() {
	const nChr = 2;

	{
		[...dataset.progeny_list].forEach(a => RemoveGenomeHelper.remove(a));
		RemoveGenomeHelper.rebuild();

		document.querySelector("#el_input_chr").value = nChr
		document.querySelector("#el_input_chr").oninput(null)
		await promise_load_task
		
		await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
			await load_multi_alignment_results(chr_idx);
		}));
		
		ssDNA_20211014.file_extname = "txt"
		ssDNA_20211014.lg = [2];
		ssDNA_20211014.show_pvalue = false;
		await ssDNA_20211014(true, { row_height: 3, fork_feat: true, load_QM6a: false, load_CBS1: false, load_dup3: true }, false);
		viewerState.resizeCanvas()
		await drawFrame()
	}
	ssDNA_peak_venn(0.05, 0, false)
	
	/**
	 * @see {@link ssDNA_peak_venn} venn 的數量會少於原本的peak
	 */

	// const bbb = ssDNA_peak_venn.output_info_table.map(a => [a[0], a[1].length].join("\t")).join("\n");

	// const aaa = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") > 0 && a.tags[1] != null).map(a => {
	// 	const line_head = [a.ref, ...a.tags].join(" ");
	// 	const cc = a.data.reduce((acc, peaks) => acc + peaks.reduce((acc2, peak) => acc2 + (peak.ks_test_pvalue <= 0.05 && peak.$_len >= 40 && peak.value >= 0 ? 1 : 0), 0), 0);
	// 	return [line_head, cc].join("\t");
	// });

	// document.querySelector("textarea").value = aaa + "\n\n\n\n" + bbb;
}

/**
 * TSETA:Tr_ssDNA_parent
 * @param {boolean} [clear]
 */
 async function ssDNA_raw_limit_20211014(clear = false, row_height = 3, dataset_init = false) {
	
	await [
		"QM6a_SS_Rad51_1.sort.bam.scaled.bed",
		"QM6a_SS_Rad51_2.sort.bam.scaled.bed",
		"QM6a_SS_Rad51_3.sort.bam.scaled.bed",
		"QM6a_SS_Sae2_1.sort.bam.scaled.bed",
		"QM6a_SS_Sae2_2.sort.bam.scaled.bed",
		"QM6a_SS_Sae2_3.sort.bam.scaled.bed",
		"QM6a_ss_spo11_rad51_1.sort.bam.scaled.bed",
		"QM6a_ss_spo11_rad51_2.sort.bam.scaled.bed",
		"QM6a_ss_spo11_rad51_3.sort.bam.scaled.bed",
		"QM6a_ss_spo11_sae2_1.sort.bam.scaled.bed",
		"QM6a_ss_spo11_sae2_2.sort.bam.scaled.bed",
		"QM6a_ss_spo11_sae2_3.sort.bam.scaled.bed",
	].reduce(async function (promise, file_name) {
		await promise;

		const url = `20210071/results/bamscale/bam_scaled/${file_name}`;

		const displayName = file_name.replace("_ss_", " ").replace(".sort.bam.scaled.bed", "");

		const cfg = new module_Methyl_sampleData({
			ref: "QM6a",
			sample: displayName,
			name: displayName,
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: url,
			description: "",
			max_display_value: null,
			value_normalize: null,
			density_to_opacity: false,
			row_height: row_height,//row scale
			// row_height: 1.5,
			// seg_row_height: viewerState.seg_row_height
			// fontSize_scale: 0.8,
			fontSize_scale: 0.8,
			border_color: "#000000",
			rendering_condition: null,
			display_minus_value: false,

			fast: false,
			fast_binSize: 1,

			tags: ["sort.bam.scaled.bed"],
		});

		debugger

		const table_header = [
			"sChr",
			"start",
			"end",
			"id",
			"value",
		];
		const chr_idx = viewerState.nChr - 1;
		cfg.data = await g_methylRenderer._load_tsv(cfg, table_header, function onLoadRow(chrIdx, row, data) {
			return (
				chrIdx == chr_idx &&
				row.start >= (1029699 - 1000) &&
				row.end <= (1038762 + 1000)
			);
		});
		
		cfg.display_minus_value = false;
		cfg.max_display_value = 10;
		await ssDNA_peak_fold_CoordY(cfg, cfg.max_display_value ?? cfg.value_normalize, cfg.display_minus_value);

		await delayFrame();

		methyl_dataset_list.push(cfg);
	}, Promise.resolve());
}

async function ssDNA_20211216(nChr) {
	dataset.mode = "single";
	hide_all_marker();
	viewerState.$display_SNP_density = false;

	await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, true);

	el_input_chr.value = String(nChr);
	el_input_chr.oninput(null);
	await delayFrame();

	await promise_load_task;

	await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, false);

	// g_chrCoord.bp_start = 1029699;
	// g_chrCoord.bp_end = 1038762;

	g_chrCoord.bp_start = 974826;
	g_chrCoord.bp_end = 990211;

	methyl_dataset_list.forEach(a => a.hide = !a.hide);

	await drawFrame();
}

ssDNA_20211014.file_extname = "txt";// bed
ssDNA_20211014.lg = [
	2,
	// 3,
	// 4
];
ssDNA_20211014.show_pvalue = false;
/**
 * TSETA:Tr_ssDNA_parent
 * ssDNA_20211014.fork_feat
 * @param {boolean} [clear]
 */
async function ssDNA_20211014(clear = false, { row_height = 3, fork_feat = true, load_QM6a = true, load_CBS1 = true, load_dup3 = true }, dataset_init = false) {
	if (clear) {
		methyl_dataset_list.splice(0);// clear
	}

	if (dataset_init) {
		dataset.progeny_list.forEach(a => {
			RemoveGenomeHelper.rename(a, a.replace(/_f$/, "♀"));
			RemoveGenomeHelper.rename(a, a.replace(/_m$/, "♂"));
		});
		return;
	}

	// if (false && debug) {
	// 	methyl_dataset_list.filter((_, i) => i % 2 == 0).forEach(a => a.hide = true);
	// 	methyl_dataset_list.filter((_, i) => i % 2).forEach(a => a.max_display_value = 10)
	// 	await drawFrame();
	// }

	gff_type_priority["protein-coding gene"] = 1e+10;

	g_methylRenderer.sort_layout = true;

	// const LAYER_lg = 1e10;
	const LAYER_ref = 1e8;
	const LAYER_cmp = 1e6;
	const LAYER_cate = 1e4;
	const LAYER_dup = 1e2;

	await [
		// 2,
		// 3,
		// 4,
		...ssDNA_20211014.lg,
	].reduce(async (prev_promise, cutoff, lg_idx) => {
		await prev_promise;

		if (load_QM6a) {
			await load_group("QM6a", cutoff, lg_idx, [
				{ data_name: "Rad51_vs_Sae2", display_name: "<span>QM6a <i>rad51Δ</i> vs <i>sae2Δ</i></span>", },// main
				{ data_name: "spo11_rad51_vs_Sae2", display_name: "<span>QM6a <i>spo11Δrad51Δ</i> vs <i>sae2Δ</i></span>", },
				// "spo11_rad51_vs_spo11_sae2",// main
				// "spo11_sae2_vs_Sae2",
			], (cutoff, file_name) => {
				return `20210071/results/bamscale/diff_peaks/lg${cutoff}/Diff_peaks_${file_name}_KStest.${ssDNA_20211014.file_extname}`;
			});
		}

		if (load_CBS1 &&
			dataset.mode != "single" &&
			dataset.parental_list.includes("CBS1-1") && dataset.parental["CBS1-1"] != null
		) {
			await load_group("CBS1-1", cutoff, lg_idx, [
				{ data_name: "CBS1-1_rad51_vs_sae2", display_name: "<span>CBS1-1 <i>rad51Δ</i> vs <i>sae2Δ</i></span>", },// main
				{ data_name: "CBS1-1_spo11rad51_vs_sae2", display_name: "<span>CBS1-1 <i>spo11Δrad51Δ</i> vs <i>sae2Δ</i></span>", },
				// "CBS1-1_spo11rad51_vs_spo11sae2",// main
				// "CBS1-1_spo11sae2_vs_sae2",// "CBS1-1_sae2_vs_spo11sae2", //
			], (cutoff, file_name) => {
				return `./ssDNA_pipeline/Diff_peaks_${file_name}_KStest.${ssDNA_20211014.file_extname}`;
			});
		}

		await (load_dup3 ? [
			"QM6a",
			"CBS1-1",
		] : []).reduce(async (promise_ref, ref, ref_idx) => {
			await promise_ref;
			await [
				"rad51",
				"spo11rad51",
			].reduce(async (promise_type, cmp, cmp_idx) => {
				await promise_type;
				/**
				 * @type {module_Methyl_sampleData[][]}
				 */
				const dup_list = await [
					"1", "2", "3",
				].reduce(async (promise_rep, rep, rep_idx) => {
					const sub_list = await promise_rep;

					const file_name = `Diff_peaks_lg2_${ref}_${cmp}_${rep}_vs_sae2_KStest.txt`;
					// const file_name = `Diff_peaks_lg2_${ref}_${type}_${rep}_vs_sae2_KStest_p0.05.txt`;

					const url = `./ssDNA_pipeline/${file_name}`;

					const displayName = [
						`<div>${ref} ${{ rad51: "<i>rad51Δ</i>", spo11rad51: "<i>spo11Δrad51Δ</i>" }[cmp]} ${rep} vs <i>sae2Δ</i></div>`,
						"<div>scaled coverage (cutoff 2)</div>",
					].join("");

					sub_list.push(await load_by_preset(ref, displayName, displayName, url, [
						`${cmp}_vs_sae2`, rep,
					], false));

					await delayFrame();

					return sub_list;
				}, Promise.resolve([]));// 3dup

				// // skip if classify, cate_idx != 0
				// if (fork_feat == true && dup_list[0].length != 4) {
				// 	debugger
				// }

				// make order
				{
					const ref_layer = LAYER_ref * (ref_idx + 1);
					const cmp_layer = LAYER_cmp * (cmp_idx + 1);

					dup_list.forEach((sub_list, dup_idx) => {
						const dup_layer = LAYER_dup * (dup_idx + 0);

						sub_list.forEach((cate, cate_idx) => {
							const cate_layer = LAYER_cate * (cate_idx + 1);

							cate.order = ref_layer + cmp_layer + cate_layer + dup_layer + 1;
						});
					});
				}
				await delayFrame();
				
				const cate_idx = 0;

				const cfg_fat_peak = new module_Methyl_sampleData(dup_list[0][cate_idx]);
				cfg_fat_peak.name = [
					`<div>${ref} ${{ rad51: "<i>rad51Δ</i>", spo11rad51: "<i>spo11Δrad51Δ</i>" }[cmp]} vs <i>sae2Δ</i></div>`,
					"<div>scaled coverage (cutoff 2)</div>",
					"<div>≥2 samples</div>",
				].join("");
				cfg_fat_peak.tags = [
					`${cmp}_vs_sae2`
				];

				const chr_list = dataset.genome_info_list.find(a => a.name == dup_list[0][0].ref).chr_list;
				cfg_fat_peak.data = chr_list.map(a => []);// init
				
				await chr_list.reduce(async function (chr_promise, chr_info, chr_idx) {
					await chr_promise;

					const chr_len = chr_list[chr_idx].length;
					
					const alnmap = new Uint8Array(chr_len);
					const pvalue_trip_list = [
						new Float32Array(chr_len),
						new Float32Array(chr_len),
						new Float32Array(chr_len),
					];

					pvalue_trip_list.forEach(aa => aa.fill(NaN));

					dup_list.forEach((sub_list, dup_idx) => {
						// const dup_bit = 0b1 < dup_idx;
						sub_list[cate_idx].data[chr_idx].forEach(d => {
							// if (d.start == 2451809 && d.end == 2452941) {
							// 	console.log(sub_list[cate_idx].tags.join(" "), { dup_idx, cate_idx }, d);
							// }

							/**
							 * ks_test_pvalue ??
							 * ssDNA_20211014 |> load_by_preset |> rendering_condition
							 * @see {@link ssDNA_20211014}
							 * @see {@link ssDNA_20211014.load_by_preset.rendering_condition}
							 */
							if (d.value >= 0) {
								for (let pos = d.start; pos < d.end; ++pos) {
									const pos_idx = pos - 1;
									alnmap[pos_idx] += 1;
									pvalue_trip_list[dup_idx][pos_idx] = d.ks_test_pvalue;
								}
							}
						});
					});

					typedArray_to_rangeList(alnmap.map(v => v >= 2 ? v : 0), cfg_fat_peak.data[chr_idx]);
					cfg_fat_peak.data[chr_idx].forEach(peak => {
						// peak.value = 10;// 2 or 3
						peak.dup3 = 2;
						peak.ks_test_pvalue = 0.1;// fill
						// 20211223
						peak.pvalue_trip = [
							pvalue_trip_list[0][peak.start],
							pvalue_trip_list[1][peak.start],
							pvalue_trip_list[2][peak.start],
						];
					});
					
					await delayFrame();
				}, Promise.resolve());
				
				// clear no use data
				dup_list.forEach(sub_list => {
					// sub_list[cate_idx].data[chr_idx].forEach(d => {
					// 	d.dup3 = Math.max.apply(null, alnmap.slice(d.start, d.end));
					// });

					methyl_dataset_list.splice(methyl_dataset_list.findIndex(a => a == sub_list[cate_idx]), 1);
				});
				await delayFrame();

				// 20211220
				cfg_fat_peak.max_display_value = 3;
				await ssDNA_peak_fold_CoordY(cfg_fat_peak, cfg_fat_peak.max_display_value, false, 0);

				// delete cfg_fat_peak.hide;// cfg_fat_peak.hide = false;
				delete cfg_fat_peak.hide;

				const forked_list = await fork_by_feat(cfg_fat_peak, fork_feat);
				forked_list.forEach((forked_cfg, cate_idx) => {
					forked_cfg.order += (LAYER_cate * cate_idx) + 1;
				});
				// methyl_dataset_list.push(cfg_fat_peak);

				await delayFrame();
			}, Promise.resolve());// strain

			await delayFrame();
		}, Promise.resolve());// ref

		await delayFrame();
	}, Promise.resolve());

	/**
	 * @param {string} ref
	 * @param {number} cutoff
	 * @param {number} lg_idx
	 * @param {{ data_name: string; display_name: string; }[]} cmp_list
	 * @param {(cutoff: number; file_name: string) => string} getUrl
	 */
	async function load_group(ref, cutoff, lg_idx, cmp_list, getUrl) {
		const ref_layer = LAYER_ref * (dataset.genome_info_list.findIndex(a => a.name == ref) + 1);

		await cmp_list.reduce(async (prev_promise, { data_name: sampleName, display_name }, sample_idx) => {
			await prev_promise;

			const group_id = lg_idx * 1000 + sample_idx * 100;

			const file_name = `lg${cutoff}_${sampleName}`;

			const displayName = `
				<div>${display_name}</div>
				<div>scaled coverage (cutoff ${cutoff})</div>
				<div>average</div>
			`;
			const displayName_p005 = `${displayName}`;
			// const displayName_p005 = `${displayName} (p-value≦0.05)`;

			const url = getUrl(cutoff, file_name);
			const cate_list = await load_by_preset(ref, sampleName, displayName_p005, url, [sampleName.replace(/spo11_/g, "spo11").replace("CBS1-1_", "").toLowerCase()], fork_feat);

			const strain_layer = LAYER_cmp * (sample_idx + 1);
			for (let cate_idx = 0; cate_idx < cate_list.length; ++cate_idx) {
				const cate_layer = LAYER_cate * (cate_idx + 1);
				const dup_layer = LAYER_dup * 0;//
				cate_list[cate_idx].order = ref_layer + strain_layer + cate_layer + dup_layer;
			}

			// const cfg = await load_by_preset(ref, sampleName, displayName, `20210071/results/bamscale/diff_peaks/lg${cutoff}/diff_bins_scaled_${file_name}.bedgraph`, fork_feat);
			// cfg.border_color = ["#000", "#700", "#070", "#007"][sample_idx];

			await delayFrame();// auto GC
		}, Promise.resolve());
	}

	/**
	 * @param {string} ref
	 * @param {string} sampleName
	 * @param {string} displayName
	 * @param {string} url
	 * @param {string[]} tags
	 * @param {boolean} fork_feat
	 * @returns {Promise<module_Methyl_sampleData[]>}
	 */
	async function load_by_preset(ref, sampleName, displayName, url, tags, fork_feat)  {
		const cfg = new module_Methyl_sampleData({
			ref: ref,
			sample: sampleName,
			name: displayName,
			// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
			url: url,
			description: "",
			max_display_value: null,
			value_normalize: null,
			density_to_opacity: false,
			row_height: row_height,//row scale
			// row_height: 1.5,
			// seg_row_height: viewerState.seg_row_height
			// fontSize_scale: 0.8,
			fontSize_scale: 0.8,
			border_color: "#000000",
			// rendering_condition: [
			// 	new module_Methyl_sampleData_RenderingCondition({ color: "#0000FF", condition: d => d.ks_test_pvalue > 0.05 && d.ks_test_pvalue <= 0.1 }),
			// 	new module_Methyl_sampleData_RenderingCondition({ color: "#FF0000", condition: d => d.ks_test_pvalue <= 0.05 }),
			// ],
			display_minus_value: false,

			fast: false,
			fast_binSize: 1,

			tags: tags,
		});

		if (url.endsWith(".bed")) {
			cfg.name = `${cfg.name} peaks`;
			await g_methylRenderer.load_bed(cfg);
			cfg.data.forEach(aa => aa.forEach(v => v.value = 1));//simple peak
			methyl_dataset_list.push(cfg);
			return [cfg];
		}
		else if (url.endsWith(".bedgraph")) {
			cfg.name = `${cfg.name}`;
			await g_methylRenderer.load_bedgraph(cfg);
			await load_peak_data(cfg, fork_feat, async function (cfg) {
				cfg.display_minus_value = true;
				cfg.max_display_value = 10;
				await ssDNA_peak_fold_CoordY(cfg, cfg.max_display_value ?? cfg.value_normalize, true);
			}, async function (dst, display_minus_value) {
				dst.max_display_value = display_minus_value;
				await ssDNA_peak_fold_CoordY(dst, dst.max_display_value, true);
			});
			methyl_dataset_list.push(cfg);
			return [cfg];
		}
		else if (url.endsWith(".txt")) {
			const name = cfg.name;
			{
				// cfg.name = `${name} (scaled coverage)`;

				// 20211109
				const table_header = [
					"Peak",
					"sChr",
					"start",
					"end",
					"value",// Mean
					"Highest",
					"ks_test_pvalue", // ks_test_pvalue
				];
				cfg.data = await g_methylRenderer._load_tsv(cfg, table_header, function onLoadRow(chrIdx, row, data) {
					row.ks_test_pvalue = Number(data.ks_test_pvalue);
					// return data.value >= 0;
					return true;
				});
				
				/**
				 * ks_test_pvalue ??
				 * ssDNA_20211014 |> load_by_preset |> rendering_condition
				 * @see {@link ssDNA_20211014}
				 * @see {@link ssDNA_20211014.load_by_preset.rendering_condition}
				 */
				
				cfg.rendering_condition ??= [];
				cfg.rendering_condition.push(new module_Methyl_sampleData_RenderingCondition({
					color: "red",
					condition: (s) => !s.dup3 && s.value >= 0 && s.ks_test_pvalue <= 0.05,
					min_width: 1,
				}));
				cfg.rendering_condition.push(new module_Methyl_sampleData_RenderingCondition({
					color: "blue",
					condition: (s) => !s.dup3 && s.value >= 0 && s.ks_test_pvalue > 0.05 && s.ks_test_pvalue <= 0.1,
					min_width: 1,
				}));
				cfg.rendering_condition.push(new module_Methyl_sampleData_RenderingCondition({
					color: "black",
					condition: (s) => s.dup3 >= 2 && s.value >= 0 && s.ks_test_pvalue <= 0.1,
					min_width: 1,
				}));

				cfg.max_display_value = 10;
				const forked_list = await fork_by_feat(cfg, fork_feat);
				
				return forked_list;
			}
			if (ssDNA_20211014.show_pvalue) {
				const fork = new module_Methyl_sampleData(cfg);
				fork.name = `${name} 0.05 - p-value`;
				const table_header = [
					"Peak",
					"sChr",
					"start",
					"end",
					"Mean",// Mean
					"Highest",
					"value", // ks_test_pvalue
				];
				fork.data = await g_methylRenderer._load_tsv(fork, table_header);
				
				fork.data.forEach(rows => {
					rows.forEach(peak => {
						peak.value = 0.05 - peak.value;
					});
				});

				const [
					cfg_pcg_only,
					cfg_PCG_and_intergenic,
					cfg_intergenic_only,
					cfg_island,
				] = await load_peak_data(fork, fork_feat, async function (fork) {
					fork.display_minus_value = false;
					fork.max_display_value = 0.05;
					await ssDNA_peak_fold_CoordY(fork, fork.value_normalize, false, 3);
				}, async function (fork_dst, display_minus_value) {
					fork_dst.max_display_value = 0.05;
					await ssDNA_peak_fold_CoordY(fork_dst, 0.05, false, 3);
				});

				const forked_list = [
					fork,
					cfg_pcg_only,
					cfg_intergenic_only,
					cfg_PCG_and_intergenic,
					cfg_island,
				].filter(a => a != null);
				
				forked_list.forEach(a => methyl_dataset_list.push(a));

				return forked_list;
			}
		}
	}

	/**
	 * methyl_dataset_list.push(cfg, cfg_pcg_only, cfg_intergenic_only, cfg_PCG_and_intergenic, cfg_island)
	 * @param {module_Methyl_sampleData} cfg
	 * @param {boolean} fork_feat
	 * @returns {Promise<module_Methyl_sampleData[]>}
	 */
	async function fork_by_feat(cfg, fork_feat) {
		const [
			cfg_pcg_only,
			cfg_PCG_and_intergenic,
			cfg_intergenic_only,
			cfg_island,
		] = await load_peak_data(cfg, fork_feat, async function (cfg) {
			// cfg.display_minus_value = true;
			// cfg.max_display_value = 10;
			await ssDNA_peak_fold_CoordY(cfg, cfg.max_display_value ?? cfg.value_normalize, cfg.display_minus_value);
		}, async function (dst, display_minus_value) {
			dst.max_display_value = display_minus_value;
			await ssDNA_peak_fold_CoordY(dst, dst.max_display_value, dst.display_minus_value);
		});

		const forked_list = [
			cfg,
			cfg_pcg_only,
			cfg_PCG_and_intergenic,
			cfg_intergenic_only,
			cfg_island,
		].filter(a => a != null);

		forked_list.forEach(a => methyl_dataset_list.push(a));
		return forked_list;
	}

	/**
	 * @param {module_Methyl_sampleData} cfg
	 * @param {boolean} fork_feat
	 * @param {(cfg: module_Methyl_sampleData) => Promise<void>} fn_coordY
	 * @param {(cfg: module_Methyl_sampleData, display_minus_value: number) => Promise<void>} fn_fork_cooreY
	 */
	async function load_peak_data(cfg, fork_feat, fn_coordY, fn_fork_cooreY) {
		// const data_max_value = Math.max.apply(this, cfg.data.map(a => Math.max.apply(this, a.map(b => b.value))));
		// cfg.value_normalize = data_max_value;
		await fn_coordY(cfg);

		if (!fork_feat) {
			return [];// empty array
		}

		const aaa = dataset.genome_info_list[0].chr_list.map((_, chr_idx) => {
			// 20211105
			const aa = __chr_ssDNA_peak_overlap_gene_AT_island(chr_idx + 1, [cfg.ref], cfg, function cbfunc_init(point) {
				// delete point.gene_types; //empty
				// delete point.in_2_PCG;
				delete point.in_PCG;
				// point.not_intergenic_atisland = 0;
				
				// delete point.not_intergenic_atisland;
			}, function cbfunc_gene(point, parental_ref, gene_type, gene_idx, inside) {
				if (gene_type == "protein-coding gene") {
					point.in_PCG = inside ? 1 : -1;
					// if (point.gene_types == null) {
					// 	point.gene_types = {};
					// }
					// if (point.gene_types[parental_ref] == null) {
					// 	// point.gene_types[parental_ref] = [];
					// 	point.gene_types[parental_ref] = 0;
					// }
					// // point.gene_types[parental_ref].push(gene_type);
					// point.gene_types[parental_ref] += 1;
				}
			}, function cbfunc_island(point, parental_ref) {
				// point.not_intergenic_atisland += 1;
			});

			const all_pcg_only = cfg.data[chr_idx].filter(point => {
				return point.in_PCG > 0;
				// if (point.gene_types != null) {
				// 	if (Object.values(point.gene_types).every(a => a > 0)) {
				// 		point.in_2_PCG = true;
				// 		return true;
				// 	}
				// }
			});
			const all_intergenic_only = cfg.data[chr_idx].filter(point => {
				// return !point.in_2_PCG && (point.AT_island == null || point.AT_island.size == 0);
				return !point.in_PCG && (point.AT_island == null || point.AT_island.size == 0);
				// point.in_PCG == null
			});
			const all_PCG_and_intergenic = cfg.data[chr_idx].filter(point => {
				// return !point.in_2_PCG && (point.AT_island == null || point.AT_island.size == 0);
				return point.in_PCG < 0 && (point.AT_island == null || point.AT_island.size == 0);
			});
			// // const all_atisland = all_intergenic.filter(point => {
			// // 	return point.AT_island?.size >= 1;
			// // });

			return {
				all_pcg_only,
				all_intergenic_only,
				all_PCG_and_intergenic,
				all_AT_island: aa.AT_island,
			};
		});

		const forked_list = [
			await fork_row("protein-coding gene", cfg, aaa.map(a => a.all_pcg_only),           cfg.max_display_value, [...cfg.tags, "PCG"]),// cfg_pcg_only: 
			await fork_row("PCG+intergenic",      cfg, aaa.map(a => a.all_PCG_and_intergenic), cfg.max_display_value, [...cfg.tags, "PCG+intergenic"]),// cfg_PCG_and_intergenic: 
			await fork_row("intergenic",          cfg, aaa.map(a => a.all_intergenic_only),    cfg.max_display_value, [...cfg.tags, "intergenic"]),// cfg_intergenic_only: 
			await fork_row("AT-island",           cfg, aaa.map(a => a.all_AT_island),          cfg.max_display_value, [...cfg.tags, "ATisland"]),// cfg_island: 
		];

		cfg.is_merged = true;// 20221111
		cfg.hide = true;

		/**
		 * @param {string} append_name
		 * @param {module_Methyl_sampleData} src_cfg
		 * @param {Partial<module_Methyl_ratioData>[][]} data
		 * @param {number} display_minus_value
		 * @param {string[]} tags
		 * @returns {Promise<module_Methyl_sampleData>}
		 */
		async function fork_row(append_name, src_cfg, data, display_minus_value, tags) {
			const dst = new module_Methyl_sampleData(src_cfg);
			// dst.name = `${src_cfg.name} ${append_name}`;
			dst.name = [
				`<div>${src_cfg.name}</div>`,
				`<div>${append_name}<div>`, // style="text-align: center;"
			].join("\n");

			dst.data = [];
			dst.data = data;
			dst.tags = tags;
			dst.from_forked = true;

			await fn_fork_cooreY(dst, display_minus_value);

			dst.hide = false;

			return dst;
		}

		return forked_list;
	}
	// if (in_DCS) {
	// 	CDS
	// }
	// else {
	// 	intergenic
	// }

	// var cc = new module_Methyl_sampleData({
	// 	display_minus_value: true,
	// 	max_display_value: 10,
	// 	value_normalize: 123,
	// 	row_height: 3,//row scale
	// 	fontSize_scale: 0.8,
	// })
	// await ssDNA_peak_fold_CoordY(cc, 10, true);
	// var tx = ui_ctx.getTransform();
	// ui_ctx.setTransform(1, 0, 0, 1, Math.floor(tx.e), Math.floor(tx.f));
	// ui_ctx.drawImage(cc.canvas_value_desc, 0, 0)
}

/**
 * @example
 * await capture_all_ssDNA(false, 2, 2, true, false)
 * methyl_dataset_list.forEach(a => a.hide = a.name.indexOf("≥2 samples") >= 0)
 */
async function ssDNA_AgeI_all_chr() {
	await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
		await load_multi_alignment_results(chr_idx);
	}));

	const ref_idx = 0;
	const ref = dataset.parental_list[ref_idx];// QM6a

	// methyl_dataset_list.filter(a => a.tags[1]).map(a => [a.ref, a.tags[0], a.tags[1]] + "")
	// "QM6a,rad51_vs_sae2,PCG"
	
	// const targets = methyl_dataset_list.filter(a => a.hide && a.ref == ref && a.tags[1] == null);
	const targets = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") >= 0 && a.ref == ref && a.tags[1] == null);

	{
		/** @example sym_map["ATisland"].description == "ATisland" */
		const sym_map = {
			"PCG": Symbol("PCG"),
			"PCG+intergenic": Symbol("PCG+intergenic"),
			"intergenic": Symbol("intergenic"),
			"ATisland": Symbol("ATisland"),
		};
		methyl_dataset_list.filter(a => a.tags[1]).forEach(cfg => {
			cfg.data.forEach(peaks => {
				peaks.forEach(peak => {
					peak.where = sym_map[cfg.tags[1]];
				});
			});
		});
	}

	/**
	 * peak_trip.pvalue_avg = peak_avg.ks_test_pvalue
	 */
	ssDNA_peak_link_avg_and_2sample(ref);

	const results = dataset.results.map((seq_list, chr_idx) => {
		const seq = dataset.results[chr_idx][dataset.genome_info_list[ref_idx].chr_list[chr_idx].chr];
		const output_header = chr_idx == 0;
		return ssDNA_AgeI(chr_idx, seq, targets, output_header);
	}).flat(1);
	
	document.querySelector("textarea").value = results.join("\n");

	// return document.querySelector("textarea").value;
}

/**
 * @param {number} chr_idx
 * @param {string} seq
 * @param {module_Methyl_sampleData[]} targets
 * @param {boolean} output_header
 */
function ssDNA_AgeI(chr_idx, seq, targets, output_header) {
	const results = [];

	const sub_header = [
		"AgeI start ~ peak distance",
		"peak start",
		"peak end",
		"AgeI end ~ peak distance",
		"where",
		"p-value (average triplicate)",
		"p-value (#1)",
		"p-value (#2)",
		"p-value (#3)",
		"value",// need replace to tags[0]
	];

	[...seq.replace(/-+/g, "").matchAll(/ACCGGT[ACGT-]{5000,10000}ACCGGT/g)].map(m => {
	// [...seq.replace(/-+/g, "").matchAll(/ACCGGT[ACGT-]{3000,10000}ACCGGT/g)].map(m => {
		// const start = ref1_pos_uint32array[m.index];
		// const end = ref1_pos_uint32array[m.index + m[0].length];
		const start = m.index;
		const end = m.index + m[0].length;
		// const result = targets.map(a => {
		// 	return [
		// 		// "",//target
		// 		// "",//peak_idx
		// 		"",//start
		// 		"",//end
		// 		"",//value
		// 	];
		// });// dummy
		
		/**
		 * @type {string[][][]}
		 */
		const result = targets.map(a => []);

		// generator sub table
		targets.map((target, target_idx) => {
			// 20211222
			const peaks = target.data[chr_idx];

			// not work
			// const peaks = merge_nearest_range(target.data[chr_idx]);
			
			peaks.forEach((peak, peak_idx) => {
				/**
				 * ks_test_pvalue ??
				 * ssDNA_20211014 |> load_by_preset |> rendering_condition
				 * @see {@link ssDNA_20211014}
				 * @see {@link ssDNA_20211014.load_by_preset.rendering_condition}
				 */
				if (peak.value >= 0) {
					// paired-AgeI overlap ssDNA peak
					if (start <= peak.end && end >= peak.start) {
						if (peak.where == null) {
							console.log(chr_idx, target.name, peak);
						}
						result[target_idx].push([
							// target.tags.join(" "),
							// peak_idx.toString(),
							(peak.start - start).toString(),// AgeI start ~ peak distance
							peak.start.toString(),
							peak.end.toString(),
							(end - peak.end).toString(),// AgeI end ~ peak distance
							peak.where.description,// typeof where == 'symbol'
							peak.pvalue_avg ?? peak.ks_test_pvalue,// "p-value (average triplicate)",
							peak.pvalue_trip?.[0],// "p-value (#1)", // 20211223
							peak.pvalue_trip?.[1],// "p-value (#2)",
							peak.pvalue_trip?.[2],// "p-value (#3)",
							peak.value.toString(),
						]);
					}
				}
			});
		});

		const rows = [];
		const num_of_row = Math.max(...result.map(a => a.length));
		if (num_of_row) {
			for (let i = 0; i < num_of_row; ++i) {
				const cols = [];
				for (let j = 0; j < result.length; ++j) {
					const subcols = result?.[j]?.[i] ?? [];
					cols.splice(j * sub_header.length, sub_header.length, ...subcols);
				}
				rows.push([
					(chr_idx + 1).toString(),
					start.toString(),
					end.toString(),
					end - start,// AgeI ~ AgeI distance
					...cols,
				]);
			}

			results.push(rows.map(a => a.join("\t")).join("\n"));
		}
	});
	
	if (output_header) {
		results.unshift(
			[
				"",// ch
				"AgeI",// start
				"",// end
				"",// AgeI ~ AgeI distance
				targets.map(a => [a.tags.join(" "), ...sub_header.slice(1).map(_ => "")]),//target|start colspan=sub_header.length
			].flat(2).join("\t"),

			[
				"ch",
				"start",
				"end",
				"AgeI ~ AgeI distance",
				...targets.map(a => [...sub_header.slice(0, -1), a.tags[0]]),
			].flat(2).join("\t"),
		);
	}

	return results;

	// var aaa = [...seq_list[0].replace(/-+/g, "").matchAll(/ACCGGT[ACGT-]{3000,10000}ACCGGT/g)][0];
	// seq_list[0].replace(/-+/g, "").substr(aaa.index, aaa[0].length).slice(0, "ACCGGT".length);
	// seq_list[0].replace(/-+/g, "").substr(aaa.index, aaa[0].length).slice(-"ACCGGT".length);

	// const seq_AgeI = "ACCGGT";
	// const regexp_AgeI = new RegExp(seq_AgeI, "g");

	// 384
	// [...seq_list[0].matchAll(regexp_AgeI)];
	
	// // 386
	// [...seq_list[0].replace(/-+/g, "").matchAll(regexp_AgeI)].map(m => {
	// 	const p = ref1_pos_uint32array[m.index];
	// 	// return seq_list[0].substr(p, seq_AgeI.length) == seq_AgeI;
	// 	return p;
	// }).slice(0, -1).map((pos, idx, pos_list) => {
	// 	const next_pos = pos_list[idx + 1];
	// 	const dist = next_pos - pos;
	// 	return dist >= 3000 && dist <= 100000;
	// });
}

/**
 * link avg and 2sample if (avg ∩ 2sample)
 * 
 * avg sample (p-value <= 0.05)
 * ≧2sample (no p-value)
 * 
 * @see {@link ssDNA_AgeI_all_chr} fill data
 * @see {@link ssDNA_AgeI}
 */
function ssDNA_peak_link_avg_and_2sample(ref) {
	const targets_avg = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") < 0 && a.ref == ref && a.tags[1] == null);
	const targets_trip = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") >= 0 && a.ref == ref && a.tags[1] == null);

	targets_avg.forEach(target_avg => {
		const target_trip = targets_trip.find(v => v.tags[0] == target_avg.tags[0]);
		target_avg.data.forEach((peaks_avg, chr_idx) => {
			const peaks_trip = target_trip.data[chr_idx];
			peaks_avg.forEach(peak_avg => {
				const peak_trip = peaks_trip.find(vv => peak_avg.start <= vv.end && peak_avg.end >= vv.start);
				if (peak_trip != null) {
					// 20211224
					peak_avg.trip = peak_trip;
					peak_trip.avg = peak_avg;
				}
			});
		});
	});
}

/**
 * @param {number} min_pvalue
 * @param {number} min_peak
 */
function ssDNA_calc_number_of_peak_if_linked(min_pvalue, min_peak) {
	if (min_peak < 2 || min_peak > 3) {
		throw new Error("min_peak");
	}
	
	const table = methyl_dataset_list.filter(a => a.name.indexOf("≥2 samples") < 0 && a.tags[1] != null).map(a => {
		const line_head = [a.ref, ...a.tags].join(" ");
		const cc = a.data.reduce((acc, peaks) => acc + peaks.reduce((acc2, peak) => acc2 + (peak.ks_test_pvalue <= min_pvalue && peak.trip?.value >= 2 ? 1 : 0), 0), 0);
		return [line_head, cc].join("\t");
	});

	document.querySelector("textarea").value = table.join("\n");

	return table;
}

function test_20211224(min_pvalue, min_val) {
	// min_peak always is 2
	ssDNA_peak_venn(min_pvalue, min_val);
	
	const min_peak = 2;
	ssDNA_calc_number_of_peak_if_linked(min_pvalue, min_peak);
}

class RemoveGenomeHelper {
	/**
	 * @param {string} name
	 */
	static remove(name) {
		for (let sets of ["parental", "progeny"]) {
			if (dataset[sets][name] != null) {
				const genome_idx = dataset.genome_info_list.findIndex(a => a.name == name);
				dataset.genome_info_list.splice(genome_idx, 1);
				dataset.all_telomere?.splice?.(genome_idx, 1);

				delete dataset[sets][name];

				dataset[`${sets}_list`].splice(dataset[`${sets}_list`].indexOf(name), 1);

				break;
			}
		};
		return RemoveGenomeHelper;
	}
	/**
	 * @param {string} name
	 * @param {string} new_name
	 */
	static rename(name, new_name) {
		for (let sets of ["parental", "progeny"]) {
			if (dataset[sets][name] != null) {
				dataset.genome_info_list.find(a => a.name == name).name = new_name;

				const entries = Object.entries(dataset[sets]);
				entries.find(([k, v], i) => {
					if (k == name) {
						entries[i][0] = new_name;
					 }
				});
				dataset[sets] = Object.fromEntries(entries);

				dataset[`${sets}_list`].splice(dataset[`${sets}_list`].indexOf(name), 1, new_name);

				break;
			}
		};
		
		for (let prop of ["gc_content", "_gc_content"]) {
			if (prop in dataset) {
				if (name in dataset[prop]) {
					dataset[prop][new_name] = dataset[prop][name];
					delete dataset[prop][name];
				}
			}
		}
		
		return RemoveGenomeHelper;
	}
	/**
	 * @param {string} name
	 */
	static move_to_last(name) {
		for (let sets of ["parental", "progeny"]) {
			if (dataset[sets][name] != null) {
				const genome_idx = dataset.genome_info_list.findIndex(a => a.name == name);
				last(dataset.genome_info_list, genome_idx);
				last(dataset.all_telomere, genome_idx);

				const sub_idx = dataset[`${sets}_list`].indexOf(name);

				const entries = Object.entries(dataset[sets]);
				last(entries, sub_idx);
				dataset[sets] = Object.fromEntries(entries);

				last(dataset[`${sets}_list`], sub_idx);

				/**
				 * @param {any[]} arr
				 * @param {number} idx
				 */
				function last(arr, idx) {
					if (Array.isArray(arr)) {
						arr.push(...arr.splice(idx, 1));
					}
				}

				break;
			}
		};
		return RemoveGenomeHelper;
	}
	static rebuild() {
		return RemoveGenomeHelper;
	}
}

let g_methylRenderer = new module_MethylRenderer(methyl_dataset_list);

class GC_Content {
	constructor() {
		this.name = "";
		this.chr = "";
		this.start = 0;
		this.end = 0;
		this.gc = 0;
	}
}
/** @type {{[name:string]:{[chr:string]:GC_Content[]}}} */
let gc_content = {};
/** @type {{[name:string]:{[chr:string]:number}}} */
let gc_content_average = {};

const el_display_bp_pos = document.getElementById("el_display_bp_pos");
const el_display_ref1_bp_pos = document.getElementById("el_display_ref1_bp_pos");
const el_display_ref2_bp_pos = document.getElementById("el_display_ref2_bp_pos");
/** @type {HTMLSelectElement} */
// @ts-ignore
const el_gc_content_window_size = document.getElementById("el_gc_content_window_size");

//////////////////////////////////////////////////////////////////////////////

/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_start = document.getElementById("el_input_start");
/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_end = document.getElementById("el_input_end");

//////////////////////////////////////////////////////////////////////////////

/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_ref1_start = document.getElementById("el_input_ref1_start");
/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_ref1_end = document.getElementById("el_input_ref1_end");

//////////////////////////////////////////////////////////////////////////////

/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_ref2_start = document.getElementById("el_input_ref2_start");
/** @type {HTMLInputElement} */
// @ts-ignore
const el_input_ref2_end = document.getElementById("el_input_ref2_end");

//////////////////////////////////////////////////////////////////////////////

/** @type {HTMLCanvasElement} */
const main_canvas = (document.createElement("canvas"));
/** @type {HTMLCanvasElement} */
const ui_canvas = (document.getElementById("canvas"));

/** @type {CanvasRenderingContext2D} */
const main_ctx = main_canvas.getContext("2d");

const ui_ctx = ui_canvas.getContext("2d");
{
	// set small width for auto resize canvas by parentElement
	main_canvas.width = 1563;//ui_canvas.parentElement.clientWidth - 20;//??
	// main_canvas.style.width = main_canvas.width + "px";

	main_canvas.height = ui_canvas.parentElement.clientHeight;
	// main_canvas.style.height = main_canvas.height + "px";

	ui_canvas.width = main_canvas.width;
	ui_canvas.height = main_canvas.height;

	// @ts-ignore
	main_ctx.$fillText = main_ctx.fillText;
	// @ts-ignore
	main_ctx.$strokeText = main_ctx.strokeText;

	// viewerState.padding_right;
	// viewerState.refresh_padding_right;

	/**
	 * @param {string} text
	 * @param {number} x
	 * @param {number} y
	 * @param {number} [maxWidth]
	 */
	main_ctx.fillText = function fillText(text, x, y, maxWidth) {
		if (viewerState.auto_padding_right) {
			const w = main_ctx.measureText(text).width;
			viewerState.refresh_padding_right = Math.max(viewerState.refresh_padding_right, w);
		}

		// @ts-ignore
		return main_ctx.$fillText(text, x, y, maxWidth);
	};

	/**
	 * @param {string} text
	 * @param {number} x
	 * @param {number} y
	 * @param {number} [maxWidth]
	 */
	main_ctx.strokeText = function strokeText(text, x, y, maxWidth) {
		const w = main_ctx.measureText(text).width;
		viewerState.refresh_padding_right = Math.max(viewerState.refresh_padding_right, w);

		// @ts-ignore
		return main_ctx.$strokeText(text, x, y, maxWidth);
	};
}

// const canvas_regs = [
// 	document.createElement("canvas"),
// 	document.createElement("canvas"),
// ];
// const ctx_segs = canvas_regs.map(c => c.getContext("2d"));

let getViewLength = () => (g_chrCoord.bp_end - g_chrCoord.bp_start + 1);
let getBPperPixel = () => getViewLength() / main_canvas.width;
let getPixelPerBP = () => main_canvas.width / getViewLength();
let g_maxPixelPerBP = 1;

const merge_near_priority_small = Object.fromEntries(Object.entries(ColorID).map(a => [a[1], 1]));
const merge_near_priority = Object.fromEntries(Object.entries(ColorID).map(a => [a[1], 1]));
{
	merge_near_priority_small[ColorID.identical] = 0; // 2 // 20210616
	// merge_near_priority_small[ColorID.diff] = 0; // 5
	// merge_near_priority_small[ColorID.none] = 0; // 8
}

class PlotInfo {
	title = "";
	row_height = 0;
	row_separate = 0;
	rowspan = 1;
	p_fontSize = 100;

	/**
	 * @type {(info: PlotInfo) => Promise<void>}
	 */
	func = null;

	/**
	 * @param {Partial<PlotInfo>} info
	 */
	constructor(info) {
		Object.assign(this, info);
	}
}

class ViewerState {
	$display_SNP_density = dataset.mode != "tetrad" && dataset.parental_list.length > 1;
	/**
	 * @param {number} value
	 * @param {number} max
	 */
	$SNP_density_func = (value, max) => value / max;

	/**
	 * ViewerState.prototype.$SNP_density_compare_parental = dataset.mode != "tetrad";
	 */
	$SNP_density_compare_parental = dataset.mode != "tetrad";

	removed_color = null;//"#0000FF";

	// /** @type {boolean} */
	// _display_13 = true;
	/** @type {boolean} */
	_display_31 = true;
	// /** @type {boolean} */
	// _display_04 = true;
	/** @type {boolean} */
	_display_40 = true;

	/** @type {boolean} */
	_display_2n2 = true;
	/** @type {boolean} */
	_display_1n3 = true;
	/** @type {boolean} */
	_display_3n1 = true;
	/** @type {boolean} */
	_display_4n0 = true;

	/** @type {boolean} */
	_display_rip = true;
	/** @type {boolean} */
	_display_rip_2 = true;

	/** @type {boolean} */
	_display_illegitimate_mutation = true;
	_display_illegitimate_mutation_indel = true;

	// illegitimate_mutation_31 = true;
	// illegitimate_mutation_40 = true;
	sss_04 = true;
	sss_22 = true;
	sss_31 = true;
	sss_40 = true;

	/** @type {boolean} */
	_display_snp = true;
	/** @type {boolean} */
	_display_snv = true;

	/** @type {boolean} */
	_display_co = true;

	/** @type {boolean} */
	_display_nonCo = false;

	/** @type {boolean} */
	display_nonCo2 = false;

	/** @type {number} */
	_nChr = 0;

	_disable_max_length = true;

	_crossover_only = true;

	// _padding_left = 16 * 5;//see: this.init(ctx)
	_padding_right = 16 * 10;//see: this.init(ctx)

	/** @type {[number,number][]} */
	_range_markers = [];

	/** @type {number[]} */
	_position_markers = [];
	/** @type {number[]} */
	_position_ref1_markers = [];
	/** @type {number[]} */
	_position_ref2_markers = [];

	_rip_display_weight = 10000000;//100 is small

	_display_parent_rip = false;

	$animationFrameId = 0;

	_gc_content_clip_indel = true;

	_display_rdna_border = false;

	seg_row_height = 32;
	//
	_seg_row_separate = 5;//15//5;
	set seg_row_separate(value) {
		document.querySelectorAll(".gc-plot").forEach(el => {
			el.style.marginBottom = `${value}px`;
		});

		document.querySelectorAll(".markers").forEach(el => {
			el.style.marginBottom = `${value}px`;
		});

		// progeny
		document.querySelectorAll("#progney-group .row-header > div > *").forEach(el => {
			el.style.marginBottom = `${value}px`;
		});

		this._seg_row_separate = value;
	}
	get seg_row_separate() {
		return this._seg_row_separate;
	}

	fill_prev_color = true;

	_display_progeny = true;
	get display_progeny() {
		return this._display_progeny;
	}
	set display_progeny(val) {
		document.querySelector("#progney-group").style.display = val ? "block" : "none";

		this._display_progeny = val;
	}

	display_GC_content = true;
	display_telomere = true;

	display_meth_ratio = true;
	display_nBS_meth_ratio = false;
	display_BS_nBS_meth_ratio = false;
	display_RNA_coverage = false;
	display_meth_diff_wt = false;

	display_gff = false;
	display_gff_CDS = false;
	display_gff_exon = false;
	display_gff_mRNA = false;

	auto_scale_RNA = true;
	/**
	 * @type {number}
	 * @see {@link window.rna_reads_func}
	 */
	rna_reads_max_display_value = null;
	/**
	 * @type {number}
	 * @see {@link window.rna_reads_func}
	 */
	rna_reads_low_cutoff_display_value = null;
	/**
	 * @type {number}
	 * @see {@link window.rna_reads_func}
	 */
	rna_reads_high_cutoff_display_value = null;

	display_cluster = [
		true, true, true, true, true, true,
	];
	display_new_gff3 = true;

	display_repeat_segment = false;

	marker_width = 1;

	/** @type {Partial<RepeatSegment>} */
	repeat_segment_filter_min = {
		identity: 65,
		set Identity(a) {
			this.identity = a;
		},
		get Identity() {
			return this.identity;
		},

		alignment: 100,
		set Alignment(a) {
			this.alignment = a;
		},
		get Alignment() {
			return this.alignment;
		},

		s_len: 100,
		set sLen(a) {
			this.s_len = a;
		},
		get sLen() {
			return this.s_len;
		},

		q_len: 100,
		set qLen(a) {
			this.q_len = a;
		},
		get qLen() {
			return this.q_len;
		},
	};
	/** @type {Partial<RepeatSegment>} */
	repeat_segment_filter_max = {
		set Identity(a) {
			this.identity = a;
		},
		get Identity() {
			return this.identity;
		},

		set Alignment(a) {
			this.alignment = a;
		},
		get Alignment() {
			return this.alignment;
		},

		set sLen(a) {
			this.s_len = a;
		},
		get sLen() {
			return this.s_len;
		},

		set qLen(a) {
			this.q_len = a;
		},
		get qLen() {
			return this.q_len;
		},
	};

	mask_NCO_is_RIP = false;
	display_half_color = false;
	colored_RIP_in_gene = false;

	// D8_bg_RIP_marker = true && dataset.mode != "SNP";

	display_rDNA = true;

	fill_large_identical = false;

	pageX = 0;

	/** @deprecated */
	plot_height_limit = null;
	/** @deprecated */
	plot_height_hard_limit = null;

	/** @returns height */
	get fixedPlotHeight() {
		return this.plot_height_hard_limit ?? this.plot_height_limit ?? 0;
	}
	/** @param {number} height */
	set fixedPlotHeight(height) {
		this.plot_height_hard_limit = height;
	}

	get_plot_height() {
		let height = this.fixedPlotHeight;
		if (height > 0 && Number.isSafeInteger(height)) {
			return height;
		}

		let { seg_row_height, seg_row_separate } = this;

		height += (seg_row_height + seg_row_separate) * seq_list.length;//seq

		height += (seg_row_height * 2 + seg_row_separate) * 2;//GC%

		// height += (seg_row_height + seg_row_separate) * 1;//CO+ marker
		// height += (seg_row_height + seg_row_separate) * 1;//NCO1 marker
		// height += (seg_row_height + seg_row_separate) * 1;//NCO2 marker
		// height += (seg_row_height + seg_row_separate) * allMarker.list.length;//nco markers

		// height += (Math.trunc(seg_row_height * 1.5) + seg_row_separate) * (6 * 2);//methy ratio
		// height += (Math.trunc(seg_row_height * 1.5) + seg_row_separate) * (6 * 2);//nBS methy ratio
		// height += (Math.trunc(seg_row_height * 1.5) + seg_row_separate) * (6 * 2);//rna mapping
		// height += (seg_row_height * 2 + seg_row_separate) * ((6 - 1) * 2);//methy CG diff
		// height += (seg_row_height * 2 + seg_row_separate) * ((6 - 1) * 2);//methy CHG diff
		// height += (seg_row_height * 2 + seg_row_separate) * ((6 - 1) * 2);//methy CHH diff
		// height += (seg_row_height * 2 + seg_row_separate) * (6 * 2);//methy BS-nBs CG diff

		// height += (seg_row_height + seg_row_separate) * (6 * 2);//sample separate

		// height += seg_row_separate;// padding bottom

		if (methyl_dataset_list?.length) {
			height += (seg_row_height + seg_row_separate) * 1.5 * methyl_dataset_list.length;
		}

		return Math.trunc(Math.max(height, main_ctx.getTransform().f));
	}

	/**
	 * @param {string} plot_title
	 * @param {string} plot_filename
	 */
	setPlotTitle(plot_title = `${dataset.name} Ch${this.nChr}`, plot_filename = undefined) {
		let thead = document.getElementById("plot-title");

		if (!thead) {
			thead = document.createElement("thead");
			thead.id = "plot-title";

			const table = document.getElementById("plot-area");
			table.insertBefore(thead, table.childNodes[0]);
		}

		this.$plot_title = plot_title;

		thead.innerHTML = `<tr><td colspan="2" id="plot-title-text" contenteditable="true" style="text-align: center;">${plot_title}</td></tr>`;
		
		const ts = document.getElementById("plot-title-text");
		if (ts != null) {
			ts.dataset.download = plot_filename;
		}

		return plot_title;
	}

	get ref1() {
		return dataset.parental_list[0];
	}
	get ref2() {
		return dataset.parental_list[1];
	}
	get nChr() {
		return this._nChr;
	}
	set nChr(value) {
		this._nChr = value;
		this.setPlotTitle();
		if (this._nChr != value) {
			drawFrame();
		}
	}

	get disable_max_length() {
		return this._disable_max_length;
	}
	set disable_max_length(value) {
		this._disable_max_length = value;
		drawFrame();
	}

	resizeCanvas() {
		const height = viewerState.get_plot_height();

		if (height == main_canvas.height || height == ui_canvas.height) {
			return null;
		}

		main_canvas.height = height;
		// main_canvas.style.height = main_canvas.height + "px";

		// canvas_regs.forEach(c => {
		// 	c.width = viewerState.max_view_width;
		// 	c.height = viewerState.seg_row_height;
		// });

		ui_canvas.width = main_canvas.width;
		ui_canvas.height = main_canvas.height;

		return height;
	}

	_max_all_chr_length = null;// 1580000;

	/**
	 * max_all_chr_length
	 */
	get max_all_chr_length() {
		if (this._max_all_chr_length) {
			return this._max_all_chr_length;
		}
		else {
			const max_len = Math.max.apply(globalThis, dataset.genome_info_list.map(b => b.chr_list.map(a => a.length)).flat(1));
			const digi = 10 ** (Math.ceil(Math.log10(max_len)) - 1);
			this._max_all_chr_length = Math.ceil(max_len / digi) * digi;

			Object.defineProperty(this, "max_all_chr_length", {
				get: function () {
					return this._max_all_chr_length;
				},
				configurable: true,
				enumerable: true,
			});

			return this._max_all_chr_length;
		}
	}

	get max_view_width() {
		if (seq_list[0]) {
			if (this.max_all_chr_length == null ||
				this.max_all_chr_length <= 0 ||
				this.disable_max_length
			) {
				return main_canvas.width - this._padding_right;
			}
			else {
				// const bp_size = (main_canvas.width - this._padding_right) / this.max_all_chr_length;
				const bp_size = (main_canvas.width - 155) / this.max_all_chr_length;
				return bp_size * Math.max(...seq_list.map(a => a.length));

				// let max_view_width = main_canvas.width * this.chr_scale_by_genome;
				// return max_view_width - this._padding_right;
			}
		}
		else {
			return 0;
		}
	}

	get display_rdna_border() {
		return this._display_rdna_border;
	}
	set display_rdna_border(value) {
		if (this._display_rdna_border != value) {
			this._display_rdna_border = value;
			drawFrame();
		}
	}

	get display_2n2() {
		return this._display_2n2;
	}
	set display_2n2(value) {
		if (this._display_2n2 != value) {
			this._display_2n2 = value;
			// drawFrame();
		}
	}

	// get display_13() {
	// 	return this._display_13;
	// }
	// set display_13(value) {
	// 	if (this._display_13 != value) {
	// 		this._display_13 = value;
	// 		// drawFrame();
	// 	}
	// }
	get display_31() {
		return this._display_31;
	}
	set display_31(value) {
		if (this._display_31 != value) {
			this._display_31 = value;
			// drawFrame();
		}
	}

	get display_1n3() {
		return this._display_1n3;
	}
	set display_1n3(value) {
		if (this._display_1n3 != value) {
			this._display_1n3 = value;
			// drawFrame();
		}
	}

	get display_3n1() {
		return this._display_3n1;
	}
	set display_3n1(value) {
		if (this._display_3n1 != value) {
			this._display_3n1 = value;
			// drawFrame();
		}
	}

	// get display_04() {
	// 	return this._display_04;
	// }
	// set display_04(value) {
	// 	if (this._display_04 != value) {
	// 		this._display_04 = value;
	// 		// drawFrame();
	// 	}
	// }

	get display_40() {
		return this._display_40;
	}
	set display_40(value) {
		if (this._display_40 != value) {
			this._display_40 = value;
			// drawFrame();
		}
	}

	get display_rip() {
		return this._display_rip;
	}
	set display_rip(value) {
		if (this._display_rip != value) {
			this._display_rip = value;
			// drawFrame();
		}
	}

	get display_4n0() {
		return this._display_4n0;
	}
	set display_4n0(value) {
		if (this._display_4n0 != value) {
			this._display_4n0 = value;
			// drawFrame();
		}
	}

	get display_illegitimate_mutation() {
		return this._display_illegitimate_mutation;
	}
	set display_illegitimate_mutation(value) {
		if (this._display_illegitimate_mutation != value) {
			this._display_illegitimate_mutation = value;
			// drawFrame();
		}
	}

	get display_snv() {
		return this._display_snv;
	}
	set display_snv(value) {
		if (this._display_snv != value) {
			this._display_snv = value;
			drawFrame();
		}
	}

	get display_snp() {
		return this._display_snp;
	}
	set display_snp(value) {
		if (this._display_snp != value) {
			this._display_snp = value;
			drawFrame();
		}
	}

	get display_co() {
		return this._display_co;
	}
	set display_co(value) {
		if (this._display_co != value) {
			this._display_co = value;
			// drawFrame();
		}
	}

	get display_conco() {
		return this._display_conco;
	}
	set display_conco(value) {
		if (this._display_conco != value) {
			this._display_conco = value;
			// drawFrame();
		}
	}

	get display_nonCo() {
		return this._display_nonCo;
	}
	set display_nonCo(value) {
		if (this._display_nonCo != value) {
			this._display_nonCo = value;
			// drawFrame();
		}
	}

	get crossover_only() {
		return this._crossover_only;
	}
	set crossover_only(value) {
		if (this._crossover_only != value) {
			this._crossover_only = value;
			drawFrame();
		}
	}

	// get padding_left() {
	// 	return this._padding_left;
	// }
	// set padding_left(value) {
	// 	if (this._padding_left != value) {
	// 		this._padding_left = value;
	// 		drawFrame();
	// 	}
	// }

	row_padding_top = 1;
	row_padding_bottom = 1;

	auto_padding_right = false;
	refresh_padding_right = 0;

	get padding_right() {
		return this._padding_right;
	}
	set padding_right(value) {
		// debugger
		if (this._padding_right != value) {
			this._padding_right = value;
			drawFrame();
		}
	}

	get display_parent_rip() {
		return this._display_parent_rip;
	}
	set display_parent_rip(value) {
		this._display_parent_rip = value;
		drawFrame();
	}

	get rip_display_weight() {
		return this._rip_display_weight;
	}
	set rip_display_weight(value) {
		this._rip_display_weight = value;
		drawFrame();
	}

	get gc_content_clip_indel() {
		return this._gc_content_clip_indel;
	}
	set gc_content_clip_indel(value) {
		this._gc_content_clip_indel = value;
		drawFrame();
	}

	/**
	 * @param {number} start
	 * @param {number} end
	 */
	pushRangeMarker(start, end) {
		this._range_markers.push([start, end]);
		drawFrame();
	}
	popRangeMarker() {
		if (this._range_markers.length) {
			let range = this._range_markers.pop();
			drawFrame();
			return range;
		}
	}

	/**
	 * @param {number} position
	 */
	pushPositionMarker(position) {
		this._position_markers.push(position);

		if (pos_ref1_uint32array) {
			let ref1_pos = pos_ref1_uint32array[position - 1];
			if (ref1_pos != null) {
				this._position_ref1_markers.push(ref1_pos);
			}
		}

		if (pos_ref2_uint32array) {
			let ref2_pos = pos_ref2_uint32array[position - 1];
			if (ref2_pos != null) {
				this._position_ref2_markers.push(ref2_pos);
			}
		}

		drawFrame();
	}
	popPositionMarker() {
		if (this._position_markers.length > 0) {
			let pos = this._position_markers.pop();
			let ref1_pos = this._position_ref1_markers.pop();
			let ref2_pos = this._position_ref2_markers.pop();
			drawFrame();
			return {
				pos, ref1_pos, ref2_pos,
			}
		}
	}

	_global_font_size = "16pt";
	_global_digi_font_size = "13pt";
	_global_font_family = "Arial";//16pt Arial

	get global_font_style() {
		return `${this._global_font_size} ${this._global_font_family}`;
	}

	get global_digi_font_style() {
		return `${this._global_digi_font_size} ${this._global_font_family}`;
	}

	/** @type {string} */
	set global_font_size(value) {
		[...document.styleSheets].forEach(ss => {
			/**
			 * @type {CSSStyleRule}
			 */
			// @ts-ignore
			const sr = ([...ss.cssRules].find(a => a.selectorText == ".plot"));

			sr.style.fontSize = value;

			this._global_font_style = value;
		});

		drawFrame();
	}

	/**
	 * @param {CanvasRenderingContext2D} ctx
	 */
	init(ctx) {
		ctx.fillStyle = "black";
		ctx.textAlign = "left";
		ctx.font = this.global_font_style;//"16pt Arial";//15 + "px Arial";//14

		// this.padding_right
		this._padding_right = ctx.measureText("5mC/(5mC+C)").width + 16;
		this.refresh_padding_right = this._padding_right;

		// ctx.canvas.width += this.padding_right;
	}

	// /**
	//  * @see {ViewerState#disable_max_length}
	//  */
	// enable_chr_genome_scale = false;

	// /**
	//  * @see {ViewerState#max_length}
	//  */
	// getChrInGenomeScale() {
	// 	if (this.enable_chr_genome_scale) {
	// 		const max_all_chr_length = Math.max(...dataset.genome_info_list.map(aa => Math.max(...aa.chr_list.map(a => a.length))));
	// 		const max_chr_length = Math.max(...dataset.genome_info_list.map(aa => aa.chr_list[viewerState.nChr - 1].length));
	// 		const scale = max_chr_length / max_all_chr_length;
	// 		return scale;
	// 	}
	// 	else {
	// 		return 1;
	// 	}
	// }
	
	/** @type {PlotInfo[]} */
	plot_list = [];

	$promise_prevFrame = Promise.resolve();
}
let color_user_marker_stroke = "#000000DD";
let color_user_marker_bk = "#00FF0030";

const viewerState = new ViewerState();
viewerState.init(main_ctx);

// @ts-ignore
window.viewerState = viewerState;

/**
I_D
ColorID.none | ColorID.indel

snp indel
ColorID.diff | ColorID.indel
 */
const $color_set_print = {
	"dad_bk": "#00FFFF",
	"dad": "#00FFFF",
	//"dad": "#0000FF",

	"mom_bk": "#FF90CB",
	"mom": "#FF90CB",
	//"mom": "#FF0000",

	// "identical": "lime",//"#FFFFFF",//20200914
	// "identical": "#FFFFFF",//"lime",//20200914
	"_identical": "#FFFFFF",
	get identical() { return viewerState.fill_large_identical ? "lime" : this._identical; },
	"deletion": "#FFFFFF",

	// "rip": "#000000",//darkgray (more dark) // RIP-1
	// "rip_2": "#000000",//darkgray (more dark) // RIP-2
	"rip": "crimson",//20200915//"#000000",//darkgray (more dark) // GA:AAAA or CT:TTTT
	"rip_Q": "crimson",//20200915//"#000000",//darkgray (more dark) // RIP-1 // GA:AGAG or CT:TCTC
	"rip_C": "crimson",//20200915//"#000000",//darkgray (more dark) // RIP-1
	"rip_QC": "crimson",//20200915//"#000000",//darkgray (more dark) // RIP-1 QC
	"rip_2_Q": "crimson",//20200915//"#000000",//darkgray (more dark) // RIP-2 // G-:A-A- or C-:T-T-
	"rip_2_C": "crimson",//20200915//"#000000",//darkgray (more dark) // RIP-2

	"co": "#FF0000",

	// "31": "#0000FF",
	// "40": "#00FF00",
	// "2n2": "#FFA500",
	// "3n1": "#9400d3",//darkviolet
	// "4n0": "#A0A0A0",

	// //45, 90, 135, 180, 225, 270, 315
	// "31": "hsl(60, 90%, 40%)",
	// "40": "hsl(120, 90%, 40%)",
	"22": "hsl(120, 70%, 50%)",//"rgb(38, 217, 38)"
	"31": "hsl(240, 100%, 50%)",//"#0000FF"
	"40": "hsl(240, 100%, 50%)",//"#0000FF"

	// "1n3": "hsl(190, 90%, 40%)",
	// "2n2": "hsl(250, 90%, 40%)",
	// "3n1": "hsl(280, 90%, 50%)",
	// "4n0": "hsl(0, 0%, 60%)",
	"1n3": "hsla(32, 89%, 43%, 1)",
	"2n2": "hsla(32, 89%, 43%, 1)",
	"3n1": "hsla(32, 89%, 43%, 1)",
	"4n0": "hsla(32, 89%, 43%, 1)",

	"sss_22": "hsl(120, 70%, 50%)",
	//
	"sss_04": "hsl(240, 100%, 50%)",
	"sss_13": "hsl(240, 100%, 50%)",
	"sss_31": "hsl(240, 100%, 50%)",
	"sss_40": "hsl(240, 100%, 50%)",
	// "sss_22": "hsl(250, 40%, 40%)",
	// "sss_31": "hsl(280, 40%, 50%)",
	// "sss_40": "hsl(0, 10%, 60%)",

	"illegitimate_mutation": "hsla(32, 89%, 43%, 1)",
	"illegitimate_mutation_indel": "hsla(32, 89%, 43%, 1)",
	"illegitimate_mutation_deletion": "hsla(32, 89%, 43%, 1)",
	// "illegitimate_mutation_31": "hsl(265, 40%, 50%)",
	// "illegitimate_mutation_40": "hsl(0, 20%, 40%)",
	// "illegitimate_mutation_deletion": "hsl(222, 20%, 40%)",

	"diff": "#EEEEEE",

	//

	"rDNA": "#FF0000",

	"snp": "gray",
	"snv": "green",
};
if (typeof window.version_2 == "boolean" && window.version_2) {
	$color_set_print["rip"] = "#000";//20200915//"#000000",//darkgray (more dark) // GA:AAAA or CT:TTTT
	$color_set_print["rip_Q"] = "#000";//20200915//"#000000",//darkgray (more dark) // RIP-1 // GA:AGAG or CT:TCTC
	$color_set_print["rip_C"] = "#000";//20200915//"#000000",//darkgray (more dark) // RIP-1
	$color_set_print["rip_QC"] = "#000";//20200915//"#000000",//darkgray (more dark) // RIP-1 QC
	$color_set_print["rip_2_Q"] = "#000";//20200915//"#000000",//darkgray (more dark) // RIP-2 // G-:A-A- or C-:T-T-
	$color_set_print["rip_2_C"] = "#000";//20200915//"#000000",//darkgray (more dark) // RIP-2
}

const $color_set_view = Object.assign({}, $color_set_print, {
	"rip": "#00FF00",
});

let current_colorset = $color_set_print;

const args_colors = {
	get [ColorID.dad]() { return current_colorset["dad"]; },		// 0
	get [ColorID.mom]() { return current_colorset["mom"]; },		// 1
	get [ColorID.identical]() { return current_colorset["identical"]; },	// 2 id
	// get [ColorID.dad_rip]() { return current_colorset["rip"]; },		// 3 ref1 RIP
	// get [ColorID.mom_rip]() { return current_colorset["rip_2"]; },		// 4 ref2 RIP
	get [ColorID.dad_rip]() { return current_colorset["dad"]; },		// 3 ref1 RIP
	get [ColorID.mom_rip]() { return current_colorset["mom"]; },		// 4 ref2 RIP
	get [ColorID.diff]() { return current_colorset["diff"]; },		// 5 diff
	get [ColorID.none]() { return "#FFFFFF"; },	// 8 none
};

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
}

// viewerState.display_31 = false;
// viewerState.display_40 = false;
// viewerState.display_2n2 = false;
// viewerState.display_1n3 = false;
// viewerState.display_3n1 = false;
// viewerState.display_4n0 = false;
// viewerState.display_rip = false;
// viewerState.display_illegitimate_mutation = false;
// viewerState.display_snp = false;
// viewerState.display_snv = false;
// viewerState.display_co = false;

// viewerState.illegitimate_mutation_31 = true;
// viewerState.illegitimate_mutation_40 = true;
// viewerState.sss_22 = true;
// viewerState.sss_31 = true;
// viewerState.sss_40 = true;

// const analysis_options__gff_sIdx = 0;
const analysis_options = {
	get mode() { return dataset.mode; },
	get nChr() { return viewerState.nChr },

	get point_filter() { return this.true; },
	get telomere() { return dataset.all_telomere[0][viewerState.nChr]; },
	get rDNA_info() { return dataset.rDNA_info; },
	// get rDNA_info() {
	// 	return {
	// 		chr: dataset.rDNA.nChr,
	// 		region_start: dataset.rDNA.region[0],
	// 		region_end: dataset.rDNA.region[1],
	// 	};
	// },

	//get show_rDNA_snp() { return true; },
	get co_list() { return dataset.crossover_list ? dataset.crossover_list[viewerState.nChr - 1] : null; },
	get nco_list() { return dataset.non_crossover_list ? dataset.non_crossover_list[viewerState.nChr - 1] : null; },
	get fill_prev_color() { return viewerState.fill_prev_color; },

	/** @type {ChromosomeData[]} - gemome[*].chrInfo */
	get chrInfo_list() { return dataset.genome_info_list.map(gInfo => gInfo.chr_list[analysis_options.nChr - 1]); },

	// get show_rDNA_snp() { return true; },
	// // get fill_rDNA() { return false; },
	// // get show_rDNA_non_InDel() { return true; },

	// rDNA no CO/NCO, but show 1:3 and 4:0 marker
	get show_rDNA_snp() { return false; },
	get show_rDNA_non_InDel() { return true; },

	get fill_mode() { return "poly snp filter"; },

	rip_gff: false,
	rip_repeat_methyl_gff: true,

	// get repeat_segment() {
	// 	return repeat_segment;
	// },

	// gff_sIdx: analysis_options__gff_sIdx,
	get gff() {
		if (gff_data_map) {
			const results = dataset.parental_list.reduce((obj, gff_ref, gff_idx) => {
				const symbol = dataset.genome_info_list[gff_idx].chr_list[analysis_options.nChr - 1].symbol;
				if (gff_data_map[gff_ref] && gff_data_map[gff_ref][symbol]) {
					obj[gff_ref] = gff_data_map[gff_ref][symbol];
				}
				return obj;
			}, {});
			return results;
		}
		else {
			return null;
		}
	},
	gene_rip_color: {
		"gene": "lime",
		"mRNA": "orange",
		"CDS": "purble",
		"exon": "red",
	},
};
initAnalyser(analysis_options);

allMarker.list.forEach(m => {
	$color_set_print[m.property] = $color_set_print[m.property] ?? m.color;
	$color_set_view[m.property] = $color_set_print[m.property] ?? m.color;
});

///////////////////////////////////////////////////////////////////////////////

function setup_row_header() {
	if (!document.getElementById("css_20221017")) {
		const css = document.createElement("style");
		css.innerHTML = `
#tetrad-mode-rows > .markers:nth-last-child(1) {
	margin-bottom: 0px;
}
#snp-mode-rows > .markers:nth-last-child(1) {
	margin-bottom: 0px;
}
`;
		document.body.append(css);
	}

	// const progneyGroup = document.getElementById("progney-group");
	// const el_appendRow = document.getElementById("append-row");
	// el_appendRow.after(progneyGroup);
	// update_seqRow_header();
	if (analysis_options.mode == "tetrad") {
		document.querySelectorAll("[data-mode=snp-mode]").forEach(elem => elem.style.display = "none");

		document.querySelectorAll("[data-mode=tetrad-mode]").forEach(elem => elem.style.display = null);

		// tetrad-mode-1 => crossover
		document.querySelectorAll("[data-mode-type=tetrad-mode-1]").forEach(elem => elem.style.display = null); //default
		document.querySelectorAll(".flex[data-mode-type=tetrad-mode-1]").forEach(elem => elem.style.display = "flex");
		document.querySelectorAll("[data-mode-type=tetrad-mode-2]").forEach(elem => elem.style.display = "none");

		// // tetrad-mode-2 => ??
		// document.querySelectorAll("[data-mode-type=tetrad-mode-1]").forEach(elem => elem.style.display = "none");
		// document.querySelectorAll("[data-mode-type=tetrad-mode-2]").forEach(elem => elem.style.display = null);//default
		// document.querySelectorAll(".flex[data-mode-type=tetrad-mode-2]").forEach(elem => elem.style.display = "flex");


		["SNP", "InDel"].reverse().forEach(type => {
			const el = append_parentAPD(type);
			if (el) {
				document.getElementById("progney-group").before(el);
			}
		});
	}
	else if (dataset.mode == "SNP") {
		//document.querySelectorAll("[data-mode=snp-mode]").forEach(elem => elem.style.display = "block");
		document.querySelectorAll("[data-mode=tetrad-mode]").forEach(elem => elem.style.display = "none");

		const el_seq_rows = document.getElementById("snp-mode-rows");
		if (el_seq_rows) {
			el_seq_rows.innerHTML = ""; // clear

			// ref
			dataset.parental_list.forEach(parental_name => {
				const display_name = ref_to_display_name(parental_name);
				//
				{
					let elem_div = document.createElement("div");
					el_seq_rows.append(elem_div);
					elem_div.outerHTML = `<div class="markers snp-target parental-name" data-mode="SNP" style="vertical-align: middle;"><span contenteditable="true" spellcheck="false">${display_name}</span></div>`;
				}
				// GC
				{
					let elem_div = document.createElement("div");
					el_seq_rows.append(elem_div);
					elem_div.outerHTML = `<div contenteditable="true" spellcheck="false" class="gc-plot" data-mode="SNP" style=""><span>${display_name} GC %</span></div>`;
				}
			});

			["SNP", "InDel"].forEach(type => {
				const el = append_parentAPD(type);
				if (el) {
					el_seq_rows.append(el);
				}
			});
			
			// target
			dataset.progeny_list.forEach(progeny_name => {
				let elem_div = document.createElement("div");
				el_seq_rows.append(elem_div);

				const display_name = ref_to_display_name(progeny_name);

				elem_div.outerHTML = `<div class="markers progeny-name" data-mode="SNP"><span contenteditable="true" style="vertical-align: middle;" data-id="${progeny_name}">${display_name}</span></div>`;
			});
		}
	}
	else if (dataset.mode == "single") {
		if (!20210304) {
			document.querySelector(".row-header").style.width = "20em";
		}

		document.querySelectorAll(`[data-mode="123"]`).forEach(a => a.remove());

		document.querySelectorAll("[data-mode=snp-mode]").forEach(elem => elem.style.display = "none");
		document.querySelectorAll("[data-mode=tetrad-mode]").forEach(elem => elem.style.display = "none");
		allMarker.list = [];

		// const a = document.createElement("div");
		const b = document.createElement("div");
		// document.getElementById("append-row").before(a);
		document.getElementById("append-row").before(b);
// 		a.outerHTML = `
// <div data-mode="123" style="" data-placehold="123">
// 	<div class="markers snp-reference" style="vertical-align: middle;"><span contenteditable="true" spellcheck="false">　</span></div>
// </div>`;
		b.outerHTML = `<div contenteditable="true" spellcheck="false" class="gc-plot" data-mode="123" style=""><span>${dataset.ref} GC %</span></div>`;
	}

	if (dataset.gc_content == null) {
		console.warn("no dataset.gc_content");
		document.querySelectorAll('.gc-plot[data-mode="tetrad-mode"]').forEach(a => a.style.display = "none");
	}

	/**
	 * @param {"SNP"|"InDel"} type
	 */
	function append_parentAPD(type) {
		if (dataset.parental_list.length > 1) {
			return append_el(`parental-${type}`);
		}

		/** @param {string} el_id */
		function append_el(el_id) {
			if (document.getElementById(el_id) == null) {
				const seg_row_height = viewerState.seg_row_height;
				const seg_row_separate = viewerState.seg_row_separate;
				const div = appendLeftTitle(`${dataset.parental_list[0]} vs ${dataset.parental_list[1]} ${type}`, seg_row_height, seg_row_separate, 1, 1 * 100);// QM6a ≠ CBS1-1
				div.id = el_id;
				return div;
			}
		}
	}
}
const el_markers_table = document.getElementById("markers_table");
const el_display_buttons_group = document.getElementById("display_buttons_group");
function setup_layout() {
	setup_row_header();

	// resize canvas

	///////////

	el_display_buttons_group.innerHTML = "";//clear

	if (0) {
		window.TPM_diff = 10;

		const el = document.createElement("input");
		el.type = "number";
		el.min = "1";
		el.value = window.TPM_diff;
		el.oninput = function (evt) {
			window.TPM_diff = Number(el.value);
			drawFrame();
		};
		el_display_buttons_group.append(el);
	}

	appendDisplayGeneCheckbox("gene", null, function (evt, el) {
		viewerState.display_gff = el.checked;
	}).checked = viewerState.display_gff;

	appendDisplayGeneCheckbox("hide", "gene_el", function (evt, el) {
		viewerState.display_gff_mRNA = !el.checked;
		viewerState.display_gff_CDS = !el.checked;
		viewerState.display_gff_exon = !el.checked;
	}).checked = true;

	appendDisplayGeneCheckbox("mRNA", "gene_el", function (evt, el) {
		viewerState.display_gff_mRNA = el.checked;
		viewerState.display_gff_CDS = !el.checked;
		viewerState.display_gff_exon = !el.checked;
	}).checked = viewerState.display_gff_mRNA;

	appendDisplayGeneCheckbox("CDS", "gene_el", function (evt, el) {
		viewerState.display_gff_CDS = el.checked;
		viewerState.display_gff_mRNA = !el.checked;
		viewerState.display_gff_exon = !el.checked;
	}).checked = viewerState.display_gff_CDS;

	appendDisplayGeneCheckbox("exon", "gene_el", function (evt, el) {
		viewerState.display_gff_exon = el.checked;
		viewerState.display_gff_mRNA = !el.checked;
		viewerState.display_gff_CDS = !el.checked;
	}).checked = viewerState.display_gff_exon;

	/**
	 * @param {string} label_text
	 * @param {string} group
	 * @param {(evt: MouseEvent, el: HTMLInputElement) => void} onclick
	 */
	function appendDisplayGeneCheckbox(label_text, group, onclick) {
		const elem_input = document.createElement("input");

		if (!group) {
			elem_input.type = "checkbox";
		}
		else {
			elem_input.type = "radio";
			elem_input.name = group;
		}

		elem_input.onclick = function (evt) {
			onclick(evt, elem_input);
			drawFrame();
		};

		elem_input.title = label_text;

		// elem_input.append(document.createTextNode(label_text));

		el_display_buttons_group.append(elem_input);

		return elem_input;
	}

	appendDisplayAllButton("display all", function (evt) {
		show_all_marker();
	});
	appendDisplayAllButton("hide all", function (evt) {
		hide_all_marker();
	});
	appendDisplayAllButton("nw", async function (evt) {
		const save_as_png = true;
		const height_scale = 0.5;//0.2;
		const resize_height = 1366;

		// await capture_all_new_wnd(1, 4, save_as_png, height_scale, resize_height);
		// await capture_all_new_wnd(5, 8, save_as_png, height_scale, resize_height);
		// await capture_all_new_wnd(9, 12, save_as_png, height_scale, resize_height);
		// await capture_all_new_wnd(13, 16, save_as_png, height_scale, resize_height);

		// wnd1 = await capture_all_new_wnd(1, 4, true, 1);
		// wnd2 = await capture_all_new_wnd(5, 8, true, 1);
		// wnd3 = await capture_all_new_wnd(9, 12, true, 1);
		// wnd4 = await capture_all_new_wnd(13, 16, true, 1);

		if (viewerState.save_capture_all_new_wnd !== false) {
			await capture_all_new_wnd(1, 16, true, 1);
		}
		else {
			await capture_all_new_wnd(1, 16, false, 1);
		}
		viewerState.save_capture_all_new_wnd = false;
	});
	appendDisplayAllButton("seeCO", async function (evt) {
		if (window.see_co_chr != viewerState.nChr) {
			window.see_co_chr = viewerState.nChr;
			window.see_co = see_CO("CO", 1.5);
		}
		else {
			window.see_co = window.see_co ?? see_CO("CO", 1.5);
		}
		window.see_co.next();
	});
	appendDisplayAllButton("seeNCO", async function (evt) {
		if (window.see_co_chr != viewerState.nChr) {
			window.see_co_chr = viewerState.nChr;
			window.see_co = see_CO("NCO", 1.5);
		}
		else {
			window.see_co = window.see_co ?? see_CO("CO", 1.5);
		}
		window.see_co.next();
	});

	if (dataset.mode == "tetrad") {
		const el_q_rep = appendMarkersTable({
			name: `${dataset.parental_list[0]}_repeats`,//el_marker_QM6a_repeats
			property: "repeat_segment_QM6a",
			order: 0,
		});
		const el_c_rep = appendMarkersTable({
			name: "CBS1-1_repeats",
			property: `${dataset.parental_list[1]}_CBS1-1`,
			order: 0,
		});
		appendDisplayButton({
			name: "repeat_segment",
			property: "repeat_segment",
			order: 0,
			onchange: function (evt) {
				// @ts-ignore
				viewerState.display_repeat_segment = this.checked;
				// viewerState["display_" + marker.property] = this.checked;
				el_q_rep.style.display = this.checked ? "block" : "none";
				el_c_rep.style.display = this.checked ? "block" : "none";
			},
		});
		// viewerState.display_repeat_segment
	}

	if (analysis_options.mode == "tetrad") {
		const rDNA_marker = {
			name: "rDNA",
			property: "rDNA",
			order: 1,
		};
		appendDisplayButton(rDNA_marker);

		const co_marker = {
			name: "CO",
			property: "co",
			order: 2,
			title: "Crossover",
		};
		appendMarkersTable(co_marker);
		appendDisplayButton(co_marker);

		const conco_marker = {
			name: "CO(NCO)",
			property: "conco",
			order: 3,
			title: "Crossover",
		};
		appendMarkersTable(conco_marker);
		appendDisplayButton(conco_marker);

		const nco_marker = {
			name: "NCO1",
			property: "nonCo",
			order: 4,
		};
		appendMarkersTable(nco_marker);
		appendDisplayButton(nco_marker);

		const nco2_marker = {
			name: "NCO2",
			property: "nonCo2",
			order: 6,
		};
		appendMarkersTable(nco2_marker);
		appendDisplayButton(nco2_marker);
	}

	/**
	 * @param {string} label_text
	 * @param {(evt: MouseEvent) => void} onclick
	 */
	function appendDisplayAllButton(label_text, onclick) {
		let elem_input = document.createElement("button");
		// elem_input.classList.add("markers");
		elem_input.onclick = function (evt) {
			onclick(evt);
			drawFrame();
		};

		elem_input.append(document.createTextNode(label_text));

		el_display_buttons_group.append(elem_input);
	}

	// no display 2:2 SNV
	// allMarker.list.slice(1).forEach(marker => {
	allMarker.list.forEach(marker => {
		appendMarkersTable(marker);
		appendDisplayButton(marker);
	});

	{
		const el_l = document.createElement("label");
		el_l.classList.add("checkbox");
		const el_i = document.createElement("input");
		el_l.id = "CONFIG_LOAD_METH";
		el_l.append(el_i);
		el_l.append(" load meth");
		el_i.type = "checkbox";
		el_i.checked = CONFIG_LOAD_METH;
		//CONFIG_LOAD_METH = true;
		el_i.onchange = function (evt) {
			CONFIG_LOAD_METH = el_i.checked;
		};
		document.getElementById("append-row-button").append(el_l);
	}

	// (async function () {
	// 	if (dataset.name == "Tr_ssDNA_parent") {
	// 		if (!document.querySelector(".init_dialog")) {
	// 			await load_JQueryUI();

	// 			const $el = $(`<div title="init_dialog" id="init_dialog"></div>`);
	// 			$(document.body).append($el);
	// 			$el.dialog();

	// 			$el.append($(`
	// 				<div>
	// 					<input type="number" class="load_chr_ssDNA_20211014" value="7" min="1" max="${dataset.genome_info_list[0].chr_list.length}"></input>
	// 					<button class="btn_load_ssDNA_20211014">load ssDNA 20211014</button>
	// 				</div>
	// 			`));
	// 			{
	// 				let nChr = 0;
	// 				$el.children(".load_chr_ssDNA_20211014").click(async function () {
	// 					nChr = $("#el_input_chr").val();
	// 				});
	// 				$el.children(".btn_load_ssDNA_20211014").click(async function () {
	// 					if (nChr > 0 && nChr <= dataset.genome_info_list[0].chr_list.length) {
	// 						hide_all_marker();
	// 						document.querySelector("#el_input_chr").value = nChr;
	// 						document.querySelector("#el_input_chr").oninput(null);
	// 						await promise_load_task;
	// 						await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, false);
	// 						viewerState.resizeCanvas();
	// 						await drawFrame();
	// 					}
	// 				});
	// 			}
	// 		}
	// 	}
	// })();
}
function appendMarkersTable(marker) {
	let elem_div = document.createElement("div");
	elem_div.id = `el_marker_${marker.name}`;
	elem_div.title = marker.property;
	elem_div.classList.add("markers");
	elem_div.style.color = current_colorset[marker.property] || "black";
	elem_div.innerHTML = `<span contenteditable="true" spellcheck="false">${marker.name}</span>`;
	el_markers_table.append(elem_div);
	return elem_div;
}
function appendDisplayButton(marker) {
	// 200415 markers_table
	let elem_label = document.createElement("label");
	elem_label.classList.add("checkbox");

	let elem_input = document.createElement("input");
	elem_input.type = "checkbox";
	elem_input.id = `el_display_${marker.property}`;
	elem_input.checked = true;
	elem_input.onchange = function (evt) {
		const el_marker_ = document.getElementById(`el_marker_${marker.name}`);
		// if (evt) {
			if (el_marker_) {
				el_marker_.style.display = elem_input.checked ? "block" : "none";
			}
			if (marker.onchange) {
				marker.onchange.call(this, evt);
			}
			// @ts-ignore
			viewerState["display_" + marker.property] = elem_input.checked;
		// }

		drawFrame();
	};

	elem_label.append(elem_input);
	elem_label.append(document.createTextNode(` display ${marker.name}`));

	el_display_buttons_group.append(elem_label);

	if (elem_input.checked != viewerState["display_" + marker.property]) {
		elem_input.click();//fire onchange event
	}
	else {
		elem_input.onchange(null);
	}
}
//</setup_layout>

setup_layout();

///////////////////////////////////////////////////////////////////////////////

if (analysis_options.mode != "tetrad") {
	viewerState.fill_prev_color = false;
	viewerState._crossover_only = false;
}

/** @type {HTMLSelectElement} */
// @ts-ignore
let el_select_colorset = document.getElementById("el_select_colorset");

// function resizeCanvas_by_ChrScale() {
// 	main_canvas.dataset.width = main_canvas.width;
// 	main_canvas.width = (ui_canvas.parentElement.clientWidth * viewerState.getChrInGenomeScale()) - 20;
// }

function show_marker(marker) {
	viewerState["display_" + marker.property] = true;
	const input = document.getElementById(`el_display_${marker.property}`);
	if (viewerState["display_" + marker.property] != input.checked) {
		input.click();
	}
	else {
		input.onchange(null);
	}
}
function show_all_marker() {
	document.querySelectorAll(`input[id*="el_display_"`).forEach(input => {
		// @ts-ignore
		viewerState[input.id.replace(/^el_/, "")] = true;
		if (!input.checked) {
			input.click();
		}
		input.onchange(null);
	});
}
function hide_all_marker() {
	document.querySelectorAll(`input[id*="el_display_"`).forEach(input => {
		// @ts-ignore
		viewerState[input.id.replace(/^el_/, "")] = false;
		if (input.checked) {
			input.click();
		}
		input.onchange(null);
	});
}

function onChangeColorSet(evt, colorset_name, immediate_render = true) {
	if (colorset_name) {
		el_select_colorset.value = colorset_name;
	}
	if (el_select_colorset.value == "view") {
		current_colorset = $color_set_view;
	}
	else {
		current_colorset = $color_set_print;
	}
	Object.keys(current_colorset).forEach(name => {
		let el = document.getElementById("el_" + name);
		if (el) {
			el.style.background = current_colorset[name];
		}
	});
	document.body.querySelectorAll(".ref1").forEach(el => el.style.background = $color_set_view["dad_bk"]);
	document.body.querySelectorAll(".ref2").forEach(el => el.style.background = $color_set_view["mom_bk"]);

	if (immediate_render) {
		drawFrame();
	}
};

/** @type {HTMLInputElement} */
let el_input_chr = (document.getElementById("el_input_chr"));
{
	document.getElementById("el_set_ref_pos_from_ma").onclick = (evt) => set_ref_pos_from_ma(evt);

	/** @type {HTMLSelectElement} */
	el_select_colorset.onchange = (evt) => onChangeColorSet(evt);
	el_select_colorset.onchange(null);//init

	el_input_chr.max = dataset.genome_info_list[0].chr_list.length;

	// @ts-ignore
	el_input_chr.value = viewerState._nChr;
	// @ts-ignore
	el_input_chr.oninput = async (evt) => {
		viewerState.nChr = Number((evt?.target ?? el_input_chr).value);// ES next: optional chaining / nullish coalescing

		await loadData(true);

		// if (window.display_raw_reads_count) {// WARN: DON NOT USE IN FULL LENGTH
		// 	rna_window_size = 0;
		// 	window.rna_reads_func = (v, m) => v / 10000;
		// 	drawFrame();
		// }
	};

	let el_input_disable_max_length = document.getElementById("el_input_disable_max_length");
	// @ts-ignore
	el_input_disable_max_length.value = viewerState._disable_max_length;
	// @ts-ignore
	el_input_disable_max_length.oninput = (evt) => viewerState.disable_max_length = !evt.target.checked;
	viewerState.disable_max_length = viewerState._disable_max_length;

	let el_input_rip_display_weight = document.getElementById("el_input_rip_display_weight");
	if (el_input_rip_display_weight) {
		// @ts-ignore
		el_input_rip_display_weight.value = viewerState._rip_display_weight;
		// @ts-ignore
		el_input_rip_display_weight.oninput = (evt) => viewerState.rip_display_weight = Number(evt.target.value);
		viewerState.rip_display_weight = viewerState._rip_display_weight;
	}
}

/**
 * @type {HTMLSelectElement}
 * @see {@link CNChuang_spore_1}
 */
let el_select_appearance = (document.getElementById("el_select_appearance"));
if (el_select_appearance) {
	if (dataset.appearance) {
		el_select_appearance.value = dataset.appearance;
	}

	el_select_appearance.onchange = async function (evt) {
		// dataset.appearance = "diff";
		// dataset.appearance = "fill";
		dataset.appearance = el_select_appearance.value;

		// auto reload
		if (seq_list?.[0]?.length > 0 &&
			el_input_chr?.value > 0 &&
			viewerState.nChr > 0
		) {
			/**
			 * @see {@link TSETA_before_load}
			 */
			el_input_chr?.oninput?.(null);

			await promise_load_task;

			setTimeout(async () => {
				viewerState.resizeCanvas();
				await drawFrame();
			}, 1000);
		}
	}
}

async function loadFrag_to_rangeMarker() {
	try {
		const frag_list = await fetchData("frag_list.json", "json");

		// @ts-ignore
		window.frag_list = frag_list;
		// @ts-ignore
		window.$frag_list = frag_list;

		frag_list[analysis_options.nChr].forEach(coord => {
			if (coord.centromere) {
				const r1_start = Number(coord.start.r1);
				const r1_end = Number(coord.end.r1) + 1;
				const ref_start = ref1_pos_uint32array[r1_start];
				const ref_end = ref1_pos_uint32array[r1_end];
				viewerState.pushRangeMarker(ref_start, ref_end);
			}
		});
	}
	catch (ex) {
		console.error(ex);
	}
}

function set_ui_seq_pos(pageX) {
	g_chrCoord.mouse_bp_pos = screen_to_bp(pageX);

	el_display_bp_pos.innerText = "";// clear
	el_display_ref1_bp_pos.innerText = "";// clear
	el_display_ref2_bp_pos.innerText = "";// clear

	if (g_chrCoord.b_draw_selected_range) {
		const [lesser, greater] = [g_chrCoord.mouse_bp_pos_0, g_chrCoord.mouse_bp_pos].sort((a, b) => a - b);

		_set_ui_seq_pos(lesser);

		el_display_bp_pos.innerText += "..";
		el_display_ref1_bp_pos.innerText += "..";
		el_display_ref2_bp_pos.innerText += "..";
		_set_ui_seq_pos(greater);

		el_display_bp_pos.innerText += `(${greater - lesser + 1})`;
		el_display_ref1_bp_pos.innerText += `(${greater - lesser + 1})`;
		el_display_ref2_bp_pos.innerText += `(${greater - lesser + 1})`;
	}
	else {
		_set_ui_seq_pos(g_chrCoord.mouse_bp_pos);
	}

	function _set_ui_seq_pos(mbp) {
		el_display_bp_pos.innerText += mbp.toString();

		if (pos_ref1_uint32array) {
			const ref1_pos = pos_ref1_uint32array[mbp - 1];
			if (ref1_pos != null) {
				el_display_ref1_bp_pos.innerText += ref1_pos.toString();
			}
		}

		if (pos_ref2_uint32array && el_display_ref2_bp_pos) {
			const ref2_pos = pos_ref2_uint32array[mbp - 1];
			if (ref2_pos != null) {
				el_display_ref2_bp_pos.innerText += ref2_pos.toString();
			}
		}
	}
}

/**
 * @param {MouseEvent} evt
 */
window.onmousemove = function (evt) {
	if (
		(evt.target.id == "top_GUI") ||
		evt.path?.some?.(el => el?.id == "top_GUI")
	) {
		return;
	}
	let rect = ui_canvas.getBoundingClientRect();
	if (!(evt.pageY >= rect.top && evt.pageY <= rect.bottom)) {
		return;
	}

	viewerState.pageX = evt.pageX;

	set_ui_seq_pos(evt.pageX);
};

ui_canvas.onclick = function (evt) {
	let rect = ui_canvas.getBoundingClientRect();
	if (!(evt.pageY >= rect.top && evt.pageY <= rect.bottom)) {
		return;
	}

	set_ui_seq_pos(evt.pageX);
};

ui_canvas.onmousedown = function (evt) {
	set_ui_seq_pos(evt.pageX);

	g_chrCoord.mouse_bp_pos_0 = g_chrCoord.mouse_bp_pos;

	if (evt.which == 1) {
		if(evt.ctrlKey) {
			viewerState.pushPositionMarker(g_chrCoord.mouse_bp_pos);
			evt.preventDefault();
		}
		else if(evt.shiftKey) {
			viewerState.pushRangeMarker(g_chrCoord.bp_start, g_chrCoord.bp_end);
			evt.preventDefault();
		}
		else if (evt.altKey) {
			g_chrCoord.b_draw_selected_range = true;
			render_UI();
		}
	}
	else if (evt.which == 2) {
		if (evt.altKey) {
			evt.preventDefault();
			g_chrCoord.b_draw_selected_range = true;
			render_UI();
		}
		else if (evt.ctrlKey) {
			evt.preventDefault();

			const {
				el_table,
				el_thead,
				el_tbody,
			 } = (() => {
				const el_id = "all_chrCoord";
				let el = document.getElementById(el_id);
				if (el) {
				}
				else {
					el = document.createElement("table");
					el.id = el_id;
					el.style.zIndex = "10000";
					el.style.position = "fixed";
					el.style.background = "#fffe";
					el.style.borderRadius = "0 15px 15px 15px";
					el.style.border = "1px solid red";

					const el_he = document.createElement("thead");
					el.append(el_he);
					el.append(document.createElement("tbody"));

					el_he.innerHTML = [
						`<tr>`,
						`	<td colspan="3">`,
						`		<button id="btn_close_g_chrCoord">close</button>`,
						`	</td>`,
						`</tr>`,
						`<tr>`,
						`	<th>ref</th>`,
						`	<th>pos</th>`,
						`	<th>range</th>`,
						// `	<th>end</th>`,
						`</tr>`,
					].join("\n");

					document.body.append(el);

					document.getElementById("btn_close_g_chrCoord").onclick = function () {
						el.style.display = "none";
					};
				}
				el.style.display = "table";
				return {
					el_table: el,
					el_thead: el.querySelector("thead"),
					el_tbody: el.querySelector("tbody"),
				};
			})();

			el_table.style.left = `${evt.pageX}px`;
			el_table.style.top = `${evt.pageY}px`;

			el_tbody.innerHTML = "";// clear
			{
				const el_tr = document.createElement("tr");
				el_tr.innerHTML = [
					`<td title="TSETA">view</td>`,
					`<td>${el_display_bp_pos?.innerText}</td>`,
					`<td>${g_chrCoord.bp_start}-${g_chrCoord.bp_end}</td>`,
				].join("\n");
				el_tbody?.append?.(el_tr);
			}

			dataset.genomeNameList.forEach((refId, ref_idx) => {
				// const refId = k;
				const display_name = dataset.display_name[refId] ?? refId;
				const view_chr_pos = g_chrCoord.pos_to_ref(ref_idx, Number(el_display_bp_pos?.innerText));
				const view_chr_start = g_chrCoord.pos_to_ref(ref_idx, g_chrCoord.bp_start);
				const view_chr_end = g_chrCoord.pos_to_ref(ref_idx, g_chrCoord.bp_end);

				const el_tr = document.createElement("tr");
				el_tr.innerHTML = [
					`<td title="${display_name}">${refId}</td>`,
					`<td>${view_chr_pos}</td>`,
					`<td>${view_chr_start}-${view_chr_end}</td>`,
				].join("\n");
				el_tbody?.append?.(el_tr);
			});
		}
	}
	else if (evt.which == 3) {
		if(evt.ctrlKey) {
			viewerState.popPositionMarker();
			evt.preventDefault();
		}
		else if(evt.shiftKey) {
			viewerState.popRangeMarker();
			evt.preventDefault();
		}
		// else if (evt.altKey) {
		// 	g_chrCoord.b_draw_selected_range = true;
		// 	render_UI();
		// }
	}
};
ui_canvas.oncontextmenu = function (evt) {
	g_chrCoord.b_draw_selected_range = false;

	if (evt.which == 3) {
		if(evt.ctrlKey) {
			evt.preventDefault();
		}
		else if(evt.shiftKey) {
			evt.preventDefault();
		}
		else if (evt.altKey) {
			evt.preventDefault();
			render_UI();
		}
	}
};
ui_canvas.onmousemove = function (evt) {
	if (evt.which == 1) {
		set_ui_seq_pos(evt.pageX);
		if (evt.altKey) {
			g_chrCoord.b_draw_selected_range = true;
			render_UI();
		}
		else {
			let move_bp = g_chrCoord.mouse_bp_pos_0 - g_chrCoord.mouse_bp_pos;
			if (move_bp != 0) {
				set_view_range(g_chrCoord.bp_start + move_bp, g_chrCoord.bp_end + move_bp);
			}
		}
	}
	else if (evt.which == 2) {
		set_ui_seq_pos(evt.pageX);
		if (evt.altKey) {
			evt.preventDefault();
			g_chrCoord.b_draw_selected_range = true;
			render_UI();
		}
	}
	// else if (evt.which == 3) {
	// 	if (evt.altKey) {
	// 		g_chrCoord.b_draw_selected_range = true;
	// 		render_UI();
	// 	}
	// }
};
ui_canvas.onmouseup = function (evt) {
	if (evt.which == 1) {
		set_ui_seq_pos(evt.pageX);

		if (evt.altKey) {
			render_UI();

			if (g_chrCoord.b_draw_selected_range) {
				g_chrCoord.b_draw_selected_range = false;
				const [lesser, greater] = [g_chrCoord.mouse_bp_pos_0, g_chrCoord.mouse_bp_pos].sort((a, b) => a - b);
				set_view_range(lesser, greater);
			}
		}
		else {
			let move_bp = g_chrCoord.mouse_bp_pos_0 - g_chrCoord.mouse_bp_pos;
			if (move_bp != 0) {
				set_view_range(g_chrCoord.bp_start + move_bp, g_chrCoord.bp_end + move_bp);
			}
		}
	}
	else if (evt.which == 2) {
		if (evt.altKey) {
			g_chrCoord.b_draw_selected_range = false;
			render_UI();
		}
	}
	else if (evt.which == 3) {
		if (evt.altKey) {
			evt.preventDefault();
			set_view_range(0, seq_list[0].length);// show_all()
		}
	}
};

const scale_mod_ctrlKey = 1;
if (scale_mod_ctrlKey) {
	// @ts-ignore
	ui_canvas.onmousewheel = function (evt) {
		if (evt.buttons || evt.ctrlKey) {
			evt.preventDefault();
		}
	};

	document.getElementById("plot-area").onwheel = document.getElementById("plot-area").onmousewheel = function (evt) {
		if (evt.buttons || evt.ctrlKey) {
			evt.preventDefault();
		}
	};

	if (window.onwheel === null) {
		window.onwheel = function (evt) {
			if (evt.buttons || evt.ctrlKey) {
				//evt.preventDefault();//web_ui.js:1072 [Intervention] Unable to preventDefault inside passive event listener due to target being treated as passive. See https://www.chromestatus.com/features/6662647093133312
				onMouseWheel.call(this, evt);
			}
		};
	}
	else {
		window.onmousewheel = function (evt) {
			if (evt.buttons || evt.ctrlKey) {
				//evt.preventDefault();//web_ui.js:1072 [Intervention] Unable to preventDefault inside passive event listener due to target being treated as passive. See https://www.chromestatus.com/features/6662647093133312
				onMouseWheel.call(this, evt);
			}
		};
	}
}
else {
	// @ts-ignore
	ui_canvas.onmousewheel = function (evt) {
		evt.preventDefault();
	};

	if (window.onwheel === null) {
		window.onwheel = onMouseWheel;
	}
	else {
		window.onmousewheel = onMouseWheel;
	}
}
/**
 * @param {MouseWheelEvent} evt
 */
function onMouseWheel(evt) {
	let rect = ui_canvas.getBoundingClientRect();
	if (!(evt.screenY >= rect.top && evt.screenY <= rect.bottom)) {
		return;
	}

	set_ui_seq_pos(evt.pageX);

	let start, end;

	let deltaY = 0;

	// @ts-ignore
	deltaY = (evt.deltaY || evt.wheelDelta || evt.wheelDeltaY);
	deltaY = deltaY > 0 ? 1 : (deltaY < 0 ? -1 : 0);

	const zoom_scale = 0.25;
	let d_left = (g_chrCoord.mouse_bp_pos - g_chrCoord.bp_start);
	let d_right = (g_chrCoord.bp_end - g_chrCoord.mouse_bp_pos);
	if (deltaY > 0) {
		start = g_chrCoord.bp_start + Math.ceil(Math.max(1, d_left * zoom_scale) * (-deltaY));
		end = g_chrCoord.bp_end - Math.ceil(Math.max(1, d_right * zoom_scale) * (-deltaY));
	}
	else if (deltaY < 0) {
		start = g_chrCoord.bp_start + Math.ceil(d_left * zoom_scale * (-deltaY));
		end = g_chrCoord.bp_end - Math.ceil(d_right * zoom_scale * (-deltaY));
	}

	set_view_range(start, end);
};

/**
 * @param {number} start
 * @param {number} end
 */
function set_view_range(start, end) {
	if (start < end) {
		if (start < 0) {
			end += -start;
			start = 1;
			end = Math.min(end, seq_list[0].length);
		}
		if (end > seq_list[0].length) {
			start -= end - seq_list[0].length;
			end = seq_list[0].length;
			start = Math.max(1, start);
		}

		start = Math.max(1, start);
		end = Math.min(end, seq_list[0].length);

		if (start != g_chrCoord.bp_start || g_chrCoord.bp_end != end) {
			g_chrCoord.bp_start = start;
			g_chrCoord.bp_end = end;

			el_input_start.value = g_chrCoord.bp_start.toString();
			el_input_end.value = g_chrCoord.bp_end.toString();

			merge_gc_content();

			drawFrame();
		}
	}
}

/**
 * @param {number} [merge_size] - count of window, auto merge window if null
 * @param {boolean} [refresh] - refresh frame
 */
function merge_gc_content(merge_size, refresh = false) {
	if (dataset.gc_content == null) {
		return;
	}
	//seq_list[0].length
	if (merge_size == null) {
		merge_size = Math.ceil(((g_chrCoord.bp_end - g_chrCoord.bp_start + 1) / dataset.GC_Content_window) / viewerState.max_view_width);
	}

	if (merge_size < 1) {
		return;
	}

	merge_size = Math.min(merge_size, Math.trunc(max_gc_content_window / dataset.GC_Content_window));

	for (let refName of dataset.parental_list) {
		if (dataset._gc_content[refName][viewerState.nChr] == null) {
			continue;
		}

		let n_l = [];
		for (let i = 0; i < dataset._gc_content[refName][viewerState.nChr].length; i += merge_size) {
			let td = { ...dataset._gc_content[refName][viewerState.nChr][i] };//{ chr: 6, start: 1, end: 500, gc: 0 };
			let j = 1;
			for (; j < merge_size; ++j) {
				let t = dataset._gc_content[refName][viewerState.nChr][i + j];
				if (!t) {
					break;
				}
				td.end = t.end;
				td.gc += t.gc;
			}
			td.gc /= j;
			n_l.push(td);
		}
		dataset.gc_content[refName][viewerState.nChr] = n_l;
	}

	if (Number(el_gc_content_window_size.value) != merge_size) {
		el_gc_content_window_size.value = merge_size.toFixed(0);//new window size
	}

	if (refresh) {
		drawFrame();
	}
}

/**
 * @param {number} pos_x
 */
function screen_to_bp(pos_x) {
	const length = g_chrCoord.bp_end - g_chrCoord.bp_start;
	const scale = 1 / (length + 1);
	const bp_size = viewerState.max_view_width * scale;
	const rect = ui_canvas.getBoundingClientRect();

	if (pos_x >= rect.left) {
		let center_bp_pos;

		if (pos_x <= rect.right) {
			center_bp_pos = g_chrCoord.bp_start + Math.trunc((pos_x - rect.left) / bp_size);
			center_bp_pos = Math.max(0, Math.min(center_bp_pos, seq_list[0].length));
		}
		else {
			center_bp_pos = seq_list[0].length;
		}

		return center_bp_pos;
	}
	else {
		return 0;
	}
}

class GTF_ROW {
	constructor() {
		this.start = 0;
		this.end = 0;
		this.query = "";
		this.subject = "";

		/** @type {""|"transcript"|"exon"|"gene"|"mRNA"|"CDS"} */
		this.type = "";

		/** @type {-1|0|1} */
		this.strand = 0;
	}

	get min() {
		return Math.min(this.start, this.end);
	}
	get max() {
		return Math.max(this.start, this.end);
	}

	/**
	 * @param {GTF_ROW} other
	 */
	isOverlap(other) {
		return other.start <= this.end && other.end >= this.start;
	}

	/**
	 * @param {GTF_ROW} self
	 * @param {GTF_ROW} other
	 */
	static isOverlap(self, other) {
		return other.start <= self.end && other.end >= self.start;
	}
}

/**
 * @param {string} url
 * @param {string} [url_whiteList]
 */
async function loadGTF(url, url_whiteList) {
	/** @type {string} */
	const text = await fetchData(url, "text");

	/** @type {Set<string>} */
	const white_list = await (async () => {
		let text;
		try {
			text = await fetchData(url_whiteList, "text");
		}
		catch (ex) {
			throw new Error(`Not found: ${url_whiteList}`);
		}
		try {
			return url_whiteList ? new Set(text.split("\n")) : null;
		}
		catch (ex) {
			console.error("No load lncRNA list");
			return new Set();
		}
	})();

	/** @type {{ [gene_id:string]:GTF_ROW }} */
	const gene_map = {};

	const raw_table =  text.trim().split("\n").filter(a => !a.startsWith("#")).map(line => {
		const [seqname, source, feature, start, end, score, strand, frame, attribute] = line.split("\t");
		const query = attribute.match(/transcript_id "([^"]+)";/)[1];

		if (!white_list || white_list.has(query)) {
			const gene_id = attribute.match(/gene_id "([^"]+)";/)[1];

			const parent = gene_map[gene_id];

			const gtf = new GTF_ROW();
			gtf.query = query;//.replace("STRG.", "");// + "_" + seqname;
			gtf.subject = seqname;
			gtf.start = Number(start);
			gtf.end = Number(end);

			gtf.type = /** @type {"transcript"|"exon"} */(feature);
			gtf.strand = strand == "+" ? 1 : (strand == "-" ? -1 : (parent ? parent.strand : 0));

			if (feature == "transcript") {
				gene_map[gene_id] = gtf;
			}

			return gtf;
		}
	}).filter(a => a).sort(_func_sort_GTF);

	return raw_table;
}

/**
 * @param {GTF_ROW} a
 * @param {GTF_ROW} b
 */
function _func_sort_GTF(a, b) {
	const d = a.start - b.start;
	if (d) {
		return d;
	}
	else {
		return a.end - b.end;
	}
}

/**
 * @param {string} url
 * @param {string} url_whiteList
 */
async function loadPAF(url, url_whiteList) {
	/** @type {string} */
	const text = await fetchData(url, "text");

	/** @type {Set<string>} */
	const white_list = new Set((await fetchData(url_whiteList, "text")).split("\n"));

	return text.trim().split("\n").map(a => a.trim()).filter(a => a).map(line => {
		/**
0	string	Query sequence name
1	int	Query sequence length
2	int	Query start (0-based; BED-like; closed)
3	int	Query end (0-based; BED-like; open)
4	char	Relative strand: "+" or "-"
5	string	Target sequence name
6	int	Target sequence length
7	int	Target start on original strand (0-based)
8	int	Target end on original strand (0-based)
9	int	Number of residue matches
10	int	Alignment block length
11	int	Mapping quality (0-255; 255 for missing)
		 */

		let columns = line.split("\t");
		let obj = {
			query: columns[0],
			qstart: Number(columns[2]) + 1,
			qend: Number(columns[3]) + 1,
			strand: columns[4],
			subject: columns[5],
			sstart: Number(columns[7]) + 1,
			send: Number(columns[8]) + 1,
		};

		if (white_list.has(obj.query)) {
			const gtf = new GTF_ROW();

			gtf.query = obj.query.replace("STRG.", "");// + "_" + seqname;
			gtf.subject = obj.subject;
			gtf.start = Number(obj.sstart);
			gtf.end = Number(obj.send);

			gtf.type = "transcript";
			gtf.strand = obj.strand == "+" ? 1 : (obj.strand == "-" ? -1 : 0);

			return gtf;
		}
	}).filter(a => a);
}

async function load_featureCounts(url) {
	/** @type {string} */
	const text = await fetchData(url, "text");
	const rows = text.split("\n").slice(1).map(line => line.trim().split("\t")).filter(a => a && a.length);

	/** @type {{ [geneID: string]: number; }} */
	const tpm_data = {};
	let max_tpm_value = 0;

	rows.forEach(row => {
		const geneId = row[0];
		const sample_list = row.slice(1);
		// const tpm_value = Number(Math.min(...row.slice(1).map(a => Number(a))).toFixed(1)); // min
		const tpm_value = Number((sample_list.reduce((t , v) => t + Number(v), 0) / sample_list.length).toFixed(1)); // avg
		tpm_data[geneId] = tpm_value;
		max_tpm_value = Math.max(max_tpm_value, tpm_value);
	});

	return {
		tpm_data,
		max_tpm_value,
	};
}

/**
 * jqui
 */
async function init_jquery_ui() {
	await load_JQueryUI();

	const today = new Date();
	
	if (
		(dataset.name == "QM6a_CBS1-1_WT_dim2_rid1_dimrid" || dataset.name == "QM6a_CBS1-1_NP43_NP42_WT_dim2_rid1_dimrid" ) &&
		today.getTime() >= (new Date("2022/09/16").getTime())
	) {
		const $el = $(`<div if="el_load_by_preset" title="load by preset"></div>`);
		$(document.body).append($el);
		$el.dialog({
			width: "36em",
			closeOnEscape: true,
		});

		// @ts-ignore
		init_jquery_ui.$el = $el;

		// $el.append($(`<div>TDOD: delete RIP from 'ctrl/BS/oxBS rows' where color='magenta'</div>`));

		const $el_sw_p = $([
			`<div>`,
			`	<form>`,
			`		<label>`,
			`			<input type="radio" name="preset" value="D" checked />D`,
			`		</label>`,
			`		<label>`,
			`			<input type="radio" name="preset" value="M" />M`,
			`		</label>`,
			`	</form>`,
			`</div>`,
		].join("\n"));
		$el.append($el_sw_p);

		const max_length = 10_000_000;
		const el_reload_panel = $([
			`<div>`,
			`	<div>`,
			`		sub_title: <input type="text" class="input_sub_title" />`,
			`	</div>`,
			`	<div>range:`,
			`		<input type="number" value="398991" min="1" max="${max_length}" class="loading_bounding_lower" />`,
			`		~`,
			`		<input type="number" value="447710" min="1" max="${max_length}" class="loading_bounding_higher" />`,
			`		<span class="loading_bounding_length"></span>`,
			`		<button class="use_current_pos">get pos</button>`,
			`	</div>`,
			`	<div>`,
			`		extends: <input type="number" value="10000" min="1" max="${max_length}" step="1000" class="loading_bounding_extends" />`,
			`	</div>`,
			`	<div>`,
			`		<button class="unload_all_row">unload all</button>`,
			`		<button class="load_all_methyl">load 5mC</button>`,
			`		<button class="load_all_RNA">load RNA</button>`,
			`	</div>`,
			`	<hr />`,
			`	`,
			`</div>`,
		].join("\n"));
		$el.append(el_reload_panel);

		add_button("load DEG", async () => {
			const degs_list = [
				"DEGs_gene_Benjamini.txt",
				"DEGs_gene_Bonferroni.txt",
				"DEGs_gene_FDR.txt",
				"DEGs_gene_raw.txt",
			];
			const urls = degs_list.map(v => {
				const data_url = new URL(`data/20221111_DEG_seq/list/${v}`, location);
				data_url.port = 80;
				
				const api_url = new URL(location);
				api_url.port = Number(api_url.port) + 1;
				api_url.pathname = data_url.href;
				return api_url.href;
			});
		});

		add_button("load_jgi", async () => {
			await load_jgi_QM6a_gff();
		});

		/**
		 * @param {string} title
		 * @param {() => any} onclick
		 */
		function add_button(title, onclick) {
			const $btn = $(`<button>${title}</button>`);
			$el.append($btn);
			$btn.click(async (...args) => {
				const evt = args[0];
				evt.target.disabled = true;
				try {
					await onclick(...args);
				}
				finally {
					evt.target.disabled = false;
				}
			});
			return $btn;
		}
		
		$el.find(".use_current_pos").click(() => {
			g_methylRenderer.setLoadingBounding();
			el_input_ref1_start.value = g_methylRenderer._real_pos[0].start;
			el_input_ref1_end.value = g_methylRenderer._real_pos[0].end;
			el_input_ref2_start.value = g_methylRenderer._real_pos[1].start;
			el_input_ref2_end.value = g_methylRenderer._real_pos[1].end;

			$(".loading_bounding_lower").val(el_input_ref1_start.value);
			$(".loading_bounding_higher").val(el_input_ref1_end.value);

			$el.find(".loading_bounding_length").text(el_input_ref1_end.value - el_input_ref1_start.value + 1);
		});

		let loading_bounding_lower = $el.find(".loading_bounding_lower").val();
		let loading_bounding_higher = $el.find(".loading_bounding_higher").val();
		$el.find(".loading_bounding_length").text(loading_bounding_higher - loading_bounding_lower + 1);

		$el.find(".loading_bounding_lower").on("input", function (evt) {
			loading_bounding_lower = Number(evt.target.value);
			$el.find(".loading_bounding_length").text(loading_bounding_higher - loading_bounding_lower + 1);
		});
		$el.find(".loading_bounding_higher").on("input", function (evt) {
			loading_bounding_higher = Number(evt.target.value);
			$el.find(".loading_bounding_length").text(loading_bounding_higher - loading_bounding_lower + 1);
		});
		
		el_reload_panel.find(".unload_all_row").click(function () {
			methyl_dataset_list.length = 0;// clear
		});
		el_reload_panel.find(".load_all_methyl").click(async function () {
			el_reload_panel.find(".load_all_methyl").attr("disabled", true);
			try {
				// methyl_dataset_list.length = 0;// clear
				
				const title = $(".input_sub_title").val();
				const start = Number($(".loading_bounding_lower").val());
				const end = Number($(".loading_bounding_higher").val());
				const ext = Number($(".loading_bounding_extends").val());
				
				// g_methylRenderer.setLoadingBounding

				await _load_presset_20220916(title, start, end, ext, get_presetList());
				// await _load_presset_20220916(title, start, end, ext, ["CBS1-1", "Cdim", "Crid", "np59_Cdimrid"]);
			}
			finally {
				el_reload_panel.find(".load_all_methyl").attr("disabled", false);
			}
		});

		el_reload_panel.find(".load_all_RNA").click(async function () {
			el_reload_panel.find(".load_all_RNA").attr("disabled", true);
			try {
				// 20221103
				// await load_presset_20221025_QM6a_RNA("pair.btn_name", get_presetList(), ext);
				
				const title = $(".input_sub_title").val();
				const nChr = viewerState.nChr;
				const start = Number($(".loading_bounding_lower").val());
				const end = Number($(".loading_bounding_higher").val());
				const ext = Number($(".loading_bounding_extends").val());

				await _load_presset_20221025_QM6a_RNA(title, start, end, nChr , ext);
			}
			finally {
				el_reload_panel.find(".load_all_RNA").attr("disabled", false);
			}
		});

		$el.append($("<hr />"));

		const $el_panel = $([
			`<div>`,
			`	<details>`,
			`		<summary>`,
			`			load preset`,
			`		</summary>`,
			`		<!-- placehold -->`,
			`	</details>`,
			`</div>`,
		].join("\n"));
		$el.append($el_panel);
		await (async ($panel) => {
			function get_presetList() {
				/** @type {("WT" | "Qdim" | "Qrid" | "np128_Qdimrid")[]} */
				const D_preset = ["WT", "Qdim", "Qrid", "np128_Qdimrid"];
				/** @type {("CBS1-1" | "Cdim" | "Crid" | "np59_Cdimrid")[]} */
				const M_preset = ["CBS1-1", "Cdim", "Crid", "np59_Cdimrid"];
				
				/** @type {"D" | "M"} */
				const pn = $el_sw_p.find("input[name=preset]:checked").val();
				return {
					D: D_preset,
					M: M_preset,
				}[pn];
			}

			// TSETA-5mC/C-toT: The TSETA program source code used in Figure 3, Figure 4 and Figure 5.
			// TSETA-RNA-seq: The TSETA program source code used in Figure 6 and Figure S13-S15.
			const btn_pair = [
				{ btn_name: "Top2", loader: load_presset_20220916, },
				{ btn_name: "REC8", loader: load_presset_20220916, },
				{ btn_name: "MCD1", loader: load_presset_20220916, },
				{ btn_name: "gliP", loader: load_presset_20220916, },
				{ btn_name: null, loader: null, },

				{ btn_name: "ChV 398991-447710", loader: load_presset_20220916, },
				{ btn_name: "ChV usk1..cel61a",  loader: load_presset_20220916, },
				{ btn_name: "ChVI GTX",          loader: load_presset_20220916, },
				{ btn_name: "ChVII msh4",          loader: load_presset_20220916, },
				{ btn_name: null, loader: null, },
				
				{ btn_name: "ChV 398991-447710 RNA", loader: load_presset_20221025_QM6a_RNA, },
				{ btn_name: "ChV usk1..cel61a RNA",  loader: load_presset_20221025_QM6a_RNA, },
				{ btn_name: "ChVI GTX RNA",          loader: load_presset_20221025_QM6a_RNA, },
				{ btn_name: null, loader: null, },

				...(await all_preset_20221104()).map(gg => {
					return {
						btn_name: `${gg.title} BS`,
						loader: async function (_preset_name, ref_list, ext) {
							await loadTSETA_ChrData(gg.nChr);
							return await _load_presset_20220916(gg.title, gg.start, gg.end, ext, ref_list);
						},
					};
				}),
				{ btn_name: null, loader: null, },

				...(await all_preset_20221104()).map(gg => {
					return {
						btn_name: `${gg.title} RNA`,
						loader: async function (_preset_name, ref_list, ext) {
							await _load_presset_20221025_QM6a_RNA(gg.title, gg.start, gg.end, gg.nChr , ext);

							const ref_idx = 0;
							const gene_list = gff_data_map[dataset.genome_info_list[ref_idx].name][dataset.genome_info_list[ref_idx].chr_list[viewerState.nChr - 1].symbol].filter(a => a.type == "gene").sort((a, b) => a.start - b.end);

							const gene = gene_list.find(a => a.attributes.ID == gg.id);
							
							// // const idx = gene_list.sort((a, b) => a.start - b.end).findIndex(a => a.attributes.ID == gg.id);
							// // const ref_start = (gene_list[idx - 1] ?? gene_list[idx]).start + 500;
							// // const ref_end = (gene_list[idx + 1] ?? gene_list[idx]).end + 500;

							// // const ref_start = [...gene_list].reverse().find(a => a.end <= gene.start).start - 500;
							// // const ref_end = gene_list.find(a => a.start >= gene.end).end + 500;
							// const ref_start = [...gene_list].reverse().find(a => a.start <= gene.start).start - 500;
							// const ref_end = gene_list.find(a => a.end >= gene.end).end + 500;

							// const start = g_chrCoord.ref_to_pos(ref_idx, ref_start);
							// const end = g_chrCoord.ref_to_pos(ref_idx, ref_end);
							// g_chrCoord.bp_start = start;
							// g_chrCoord.bp_end = end;
							// await drawFrame();

							{
								const dd = gg.name ? `${gg.name} (${gg.id})` : gg.id;
		
								const title_text = [
									`QM6a`,

									`Ch${romanize(gg.nChr)}`,
									// `${ref_start}-${ref_end}`,

									dd,
									"RNA",
									// pair.btn_name.endsWith("RNA") ? "RNA" : null,
		
									// `±${Math.trunc(ext / 1000)}kb`
								].filter(a => a).join(" ");
		
								viewerState.setPlotTitle(title_text, title_text);

								document.getElementById("plot-title-text").dataset.suffix = [dd, "RNA"].join(" ");
							}
							// await captureScreen();// 20221108
						},
					};
				}),
				{ btn_name: null, loader: null, },

				// { btn_name: "ChVII msh4 RNA",        loader: load_presset_20221025_QM6a_RNA, },
			].forEach((pair, ppidx) => {
				if (pair.btn_name && pair.loader) {
					const $btn = $(`<button class="btn ui-btn">${ppidx}. load ${pair.btn_name}</button>`);
					$panel.append($btn);
					$btn.click(async evt => {
						evt.target.disabled = true;
						
						const ext = Number($(".loading_bounding_extends").val());
						
						await pair.loader(pair.btn_name, get_presetList(), ext);
						evt.target.disabled = false;

						// if (false) {
						// 	const qpb = dataset.genome_info_list.findIndex(v => v.name == "QM6a");
						// 	const np43 = dataset.genome_info_list.findIndex(v => v.name == "np43_QM6a");
						// 	const np42 = dataset.genome_info_list.findIndex(v => v.name == "np42_CBS1-1");
						// 	methyl_dataset_list.unshift(add_RIP_marker(qpb, np42, false));
						// 	methyl_dataset_list.unshift(add_RIP_marker(qpb, np43, false));
						// 	await drawFrame();
						// }
					});
					
					function ch() {
						/** @type {"D" | "M"} */
						const pn = $el_sw_p.find("input[name=preset]:checked").val();
						$btn.text(`${ppidx}. load ${pair.btn_name} >> ${pn}`);
					}
					$el_sw_p.find("form").change(ch);
					ch();// init
				}
				else {
					const $hr = $(`<hr></hr>`);
					$panel.append($hr);
				}
			});
		})($el_panel.find("details"));

		window.addEventListener("keydown", evt => {
			if (evt.key == "\\" || evt.code == "Backslash") {
				if ($el.dialog("isOpen")) {
					$el.dialog("close");
				}
				else {
					$el.dialog("open");
				}
			}
		});
	}
	else if (dataset.name == "Tr_ssDNA_parent") {
		// return example_ssDNA_peak_venn();

		const $el = $(`<div title="load data" id="beforeLoad_ssDNA"></div>`);
		$(document.body).append($el);
		$el.dialog();
		
		// let hide_forked_feature = true;
		$el.append($(`<label style="user-select: none;"><input type="checkbox" class="hide_forked_feature"></input> hide_forked_feature</label>`));
		$el.append($(`</ br>`));
		$el.append($(`<input class="beforeLoad_ssDNA-nChr" placehold="chr">load ssDNA</input>`));
		$el.append($(`<button class="btn_beforeLoad_ssDNA">load ssDNA</button>`));
		// $el.find(".hide_forked_feature").click(function (evt) {
		// 	hide_forked_feature = evt.target.checked;
		// 	console.log($el.find(".hide_forked_feature").prop("checked"));
		// });
		$el.children(".btn_beforeLoad_ssDNA").click(async function () {
			dataset.telomere = {}
			dataset.all_telomere = []
			dataset.rDNA = {}
			dataset.rDNA_info = {}

			// dataset.mode = "single";
			hide_all_marker();
			viewerState.$display_SNP_density = false;

			await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, true);// clear and rename

			const nChr = $(".beforeLoad_ssDNA-nChr").val();

			debugger;

			el_input_chr.value = String(nChr);
			el_input_chr.oninput(null);
			await delayFrame();

			await promise_load_task;

			// (rad51 - sae2) > 2
			await ssDNA_20211014(true, { row_height: 3, fork_feat: true }, false);// load chr data
			document.getElementById("status").innerText = "20221111: done ssDNA_20211014";
			
			// see raw peak
			if ($el.find(".hide_forked_feature").prop("checked")) {
				methyl_dataset_list.forEach(a => a.hide = a.tags.length != 1);
			}

			// low performance, hide all
			g_methylRenderer.display = false;

			{
				const last_row_idx = methyl_dataset_list.length - 1;

				// fold change >= 2
				await load_cmp_ss(false, false, false);
				methyl_dataset_list[last_row_idx + 1].ref = "";
				methyl_dataset_list[last_row_idx + 2].ref = "";
			}

			await ssDNA_raw_limit_20211014(true, 3, false);// load chr data
			
			if (nChr == 2) {
				// // raw coverage
				if (0) {
					await Promise.all([
						load_SS("QM6a"),
						// load_SS("CBS1-1"),
					]);
	
					Object.keys(rnaCoverage_float32a_map).forEach(refId => {
						Object.keys(rnaCoverage_float32a_map[refId]).forEach(sampleName => {
							const arr = rnaCoverage_float32a_map[refId][sampleName];
							const ttt = arr.reduce((aa, v) => aa + v, 0) / arr.length
							rnaCoverage_norm_map[refId][sampleName] = 100000000 / ttt;
						// Float32Array.from(arr).map(a => a * 100000000 / ttt).reduce((aa, v) => aa + v, 0)
						});
					});
					
					// viewerState.auto_scale_RNA = false;
					// viewerState.rna_reads_max_display_value = 150
					// await drawFrame()
				}

				/**
				 * @see {@link ssDNA_AgeI}
				 * DSB site
				 * ChII_QM6a:1029699..1038762
				 */
				g_chrCoord.bp_start = ref1_pos_uint32array[1029699];
				g_chrCoord.bp_end = ref1_pos_uint32array[1038762];

				// probe
				// ChII_QM6a:1036045..1036446
				viewerState.pushRangeMarker(ref1_pos_uint32array[1036045], ref1_pos_uint32array[1036446]);

				await drawFrame();
			}
		});
	}
}

/**
 * @see {@link TSETA_onload}
 */
async function TSETA_before_load() {
	await load_JQueryUI();

	if (CNChuang_spore_1) {
		if (dataset.appearance == "diff" || dataset.appearance == "difference") {
			analysis_options.reverse_fill_left_flank = false;

			// current_colorset
			$color_set_print.dad = "#0000FF";
			$color_set_print.mom = "#FF0000";
			$color_set_print._identical = "#EEEEEE";
			$color_set_print["diff"] = "#000000";
			// current_colorset
			$color_set_view.dad = "#0000FF";
			$color_set_view.mom = "#FF0000";
			$color_set_view._identical = "#EEEEEE";
			$color_set_view["diff"] = "#777777";
		}
		else if (true || dataset.appearance == "fill") {
			analysis_options.reverse_fill_left_flank = true;

			// current_colorset
			$color_set_print.dad = "#00FFFF";
			$color_set_print.mom = "#FF90CB";
			$color_set_print._identical = "#FFFFFF";
			$color_set_print["diff"] = "#EEEEEE";
			// current_colorset
			$color_set_view.dad = "#00FFFF";
			$color_set_view.mom = "#FF90CB";
			$color_set_view._identical = "#FFFFFF";
			$color_set_view["diff"] = "#EEEEEE";
		}
		return;
	}

	if (dataset.name == "Tr_ssDNA_parent") {
	}
	else if (dataset.name == "CNChuang_spore_1") {
		hide_all_marker();

		[
			allMarker.map[22].property,
			allMarker.map[31].property,
			allMarker.map["sss_22"].property,
			allMarker.map["sss_13"].property,
			allMarker.map["sss_31"].property,
		].forEach(property => {
			const id = `el_display_${property}`;
			const el = document.getElementById(id);
			if (el != null && !el.checked) {
				el.click();
			}
		});
	}
}

/**
 * @param {boolean} immediateDrawFrame
 */
async function loadData(immediateDrawFrame) {
	if (!promise_load_task) {
		const el_input_chr = document.getElementById("el_input_chr");
		el_input_chr.disabled = true;

		try {
			if (1) {
				promise_load_task = (async function () {
					viewerState.refresh_padding_right = 0;// reset

					await TSETA_before_load();

					if (viewerState.rename_html_seq) {
						rename_seq_row_header();
					}

					// resizeCanvas_by_ChrScale();

					try {
						if (typeof window?.html2canvas != "function") {
							const p = new Promise((resolve, reject) => {
								const script= document.createElement("script");
								script.async = true;
								script.onload = resolve;
								script.onerror = reject;
								script.src = "../src/html2canvas.js";
								document.body.append(script);
							});
							await p;
							if (typeof window?.html2canvas != "function") {
								debugger;
							}
						}
					}
					catch (ex) {
						console.log(ex);
					}

					try {
						await _load_data(immediateDrawFrame);
					}
					catch (ex) {
						console.error(ex);
					}

					await TSETA_onload();
				})();
				await promise_load_task;// wait all task
				promise_load_task = null;// all loaded

				// viewerState.display_co = true;
				// viewerState.display_nonCo = false;
				// viewerState.display_31 = true;
				// viewerState.display_40 = true;
				// viewerState.display_1n3 = true;
				// viewerState.display_2n2 = true;
				// viewerState.display_3n1 = true;
				// viewerState.display_4n0 = true;
				// viewerState.display_rip = true;
				// viewerState.display_illegitimate_mutation = true;
				viewerState.display_gff = false;//20200603
				// try {
				// 	document.getElementById(`el_display_${"WT(veg)"}`).checked = true;
				// 	document.getElementById(`el_display_${"WT(veg)"}`).onchange();
				// }
				// catch (ex) {
				// 	console.error(ex);
				// }
				// display_sample[0] = true;
				viewerState.display_meth_ratio = true;//20200603

				// drawFrame();//20200603
			}
			else {
				promise_load_task = (async function () {
					try {
						await _load_data(true);
					}
					catch (ex) {
						console.error(ex);
					}
				})();
				await promise_load_task;// wait all task
				promise_load_task = null;// all loaded

				const rip_indel = allMarker.map.rip.values.filter(a => seq_list[0][a.pos] == "-" || seq_list[1][a.pos] == "-");
				const rip_3indel = rip_indel.map(a => seq_list.map(ss => ss[a.pos])).filter(aaa => aaa.filter(a => a == "-").length == 3);
				if (rip_indel.length != rip_3indel.length) {
					alert("rip_indel.length != rip_3indel.length");
					console.error(rip_3indel);
				}
				allMarker.map[40].values.filter(a => seq_list.some(ss => ss[a.pos] == "-")).forEach(a => a.gc_in_out = seq_list[0][a.pos] == "-" ? 0xF0000000 : seq_list[1][a.pos] == "-" ? 0x0F000000 : undefined);

				// const snp_wnd = 3000;
				// snp_ratio = new Float32Array(Math.ceil(seq_list[0].length / snp_wnd));
				// for (let i = 0, w = 0; i < seq_list[0].length; i += snp_wnd, ++w) {
				// 	snp_ratio[w] = 0;
				// 	for (let j = 0; j < snp_wnd; ++j) {
				// 		let p = i + j;
				// 		let r1 = seq_list[0][p];
				// 		let r2 = seq_list[1][p];
				// 		if (r1 != r2) {
				// 			snp_ratio[w] += 1;
				// 		}
				// 	}
				// 	snp_ratio[w] = snp_ratio[w] / snp_wnd;
				// }

				// drawFrame();
			}
		}
		catch (ex) {
			console.error(ex);
		}
		finally {
			// viewerState._chr_scale_by_genome = viewerState._getChrScaleByGenome();

			promise_load_task = null;// all loaded
			el_input_chr.disabled = false;
			// drawFrame();//20200603
		}
	}
	else {
		await promise_load_task;
		promise_load_task = null;// all loaded
	}

	try {
		const check_len_results = analysis_options.chrInfo_list.map((v, i) => seq_list[i].replace(/-/g, "").length - v.length);
		if (check_len_results.some(a => a != 0)) {
			alert("check_len_results: " + check_len_results.join(":"));
		}
	}
	catch (ex) {
		console.error(ex);
	}

	await drawFrame();
	viewerState.resizeCanvas();
	await drawFrame();
}

/**
 * @version 20210616
 * @see {@link TSETA_before_load}
 */
async function TSETA_onload() {
	document.querySelector("#data-rows > div > div.gc-plot")?.remove?.();

	setup_row_header();
	rename_seq_row_header();
	
	if (dataset.parental[dataset.ref].endsWith("QM6a.genome.fa")) {
		try {
			// @ts-ignore
			console.info(timeElapsed(), "loading: gff", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
			document.getElementById("status").innerText = "loading: gff";

			await load_QM6a_CBS11_gff();
		}
		catch (ex) {
			console.error(ex);
		}

		try {
			if (dataset.parental_list.length > 1) {
				await loadGeneTable(gff_data_map);
			}
		}
		catch (ex) {
			console.error(ex);
		}

		// 20210608
		if (dataset.name == "Tr_ssDNA_parent") {
			if (!document.querySelector(".load_dialog")) {
				const $el = $(`<div title="load data" id="load_dialog"></div>`);
				$(document.body).append($el);
				$el.dialog();
				
				$el.append($(`<button class="btn_load_ssDNA">load ssDNA</button>`));
				$el.children(".btn_load_ssDNA").click(async function () {
					// $el.children(".btn_load_ssDNA").remove();
					// $el.attr("disabled", "disabled");
					$el.children(".btn_load_ssDNA").html("ssDNA loaded");
					await Promise.all([
						"rad51_m", "rad51_f", "sae2_m", "sae2_f",
						"spo11rad51_f",// spo11rad51_f
						"spo11sae2_f",
					].map(a => load_SS(a)));
				});

				drawTempData_afterMethyl = _drawTempData_afterMethyl;
				try {
					/**
					 * {@label ssDNA_coverage_display_config}
					 */
					if (0) {
						viewerState.rna_reads_max_display_value = 200;
					}
					else if (0) {
						rna_display_log_func = Math.log2;
						rna_display_log_base = 2;
						rna_display_exp = 10;
					}
					else {
						display_normRNA();
					}
					// viewerState.rna_reads_low_cutoff_display_value = 500;
					// viewerState.rna_reads_high_cutoff_display_value = 2000;
					// viewerState.rna_reads_max_display_value = viewerState.rna_reads_high_cutoff_display_value;

					// await load_cmp_ss(false, true, false);

					if (dataset.mode == "SNP" && dataset.parental_list.length == 2 && dataset.progeny_list == 8) {
						viewerState.display_progeny = true;
					}
				}
				catch (ex) {
					console.error("TSETA_onload", "load_cmp_ss", ex);
				}
			}
		}

		{
			// region_rect.length

			if (0) {
				globalThis.$_region_rect = [];
				globalThis.$_genome_info_list = [];
				globalThis.$_progeny_list = [];
				dataset.$_progeny = {};

				["spo11rad51_m", "spo11rad51_f", "spo11sae2_m", "spo11sae2_f"].forEach(ref => {
					const refIdx = get_genome_index(ref);

					globalThis.$_region_rect[refIdx] = dataset.genome_info_list[refIdx];
					globalThis.$_genome_info_list[refIdx] = region_rect[refIdx];
					globalThis.$_progeny_list[refIdx] = ref;
					dataset.$_progeny[ref] = dataset.progeny[ref];

					dataset.genome_info_list[refIdx] = null;
					region_rect[refIdx] = null;
					dataset.progeny_list[dataset.progeny_list.findIndex(a => a == ref)] = null;
					delete dataset.progeny[ref];
				});
				dataset.genome_info_list = dataset.genome_info_list.filter(a => a);
				region_rect = region_rect.filter(a => a);
				dataset.progeny_list = dataset.progeny_list.filter(a => a);
			}

			setup_row_header();
			hide_all_marker();
		}
	}
}

// function update_seqRow_header() {
// 	try {
// 		const el_snpRef = document.querySelector("#data-rows .snp-reference");
// 		el_snpRef.querySelector("span").innerHTML = dataset.parental_list[0];
// 	}
// 	catch (ex) {
// 		console.error(ex);
// 	}
// 	try {
// 		/** @type {HTMLDivElement} */
// 		const el_snpTarget = document.querySelector("#data-rows .markers.snp-target");
// 		const el_targets = el_snpTarget.parentElement;

// 		el_targets.innerHTML = ""; //clear

// 		[
// 			...dataset.parental_list,
// 			...dataset.progeny_list,
// 		].forEach(name => {
// 			/** @type {HTMLDivElement} */
// 			const template = (el_snpTarget.cloneNode(true));

// 			template.querySelector("span").innerHTML = name;

// 			el_targets.append(template);
// 		});
// 	}
// 	catch (ex) {
// 		console.error(ex);
// 	}
// }

/**
 * @param {boolean} immediateDrawFrame
 */
async function _load_data(immediateDrawFrame) {
	await delayFrame();
	// @ts-ignore
	console.info(timeElapsed(), "loading: clean", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: clean";
	{
		seq_list = [""];
		methy_diff_map = null;
		methy_BS_nBS_diff_map = null;
		a_methy_nBS_diff_map = null;
		methy_ratio_map = null;
		nBS_methy_ratio_map = null;

		/** @see {@link loadCoverageData} */
		// w_rnaCoverage_ui32a_map = {};
		// m_rnaCoverage_ui32a_map = {};
		// rnaCoverage_max_map = {};
		// rnaCoverage_norm_map = {};

		gff_data_map = null;
	}
	{
		align_start_index = null;
		align_end_index = null;

		seg_snp = [];
		// rip_list = [];

		region_rect = [];

		/** @type {Uint32Array} ref1 pos map to multialign */
		ref1_pos_uint32array = null;
		/** @type {Uint32Array} multialign pos map to ref1 */
		pos_ref1_uint32array = null;
		/** @type {Uint32Array} ref2 pos map to multialign */
		ref2_pos_uint32array = null;
		/** @type {Uint32Array} multialign pos map to ref2 */
		pos_ref2_uint32array = null;
	}
	console.info(timeElapsed(), "loading: clear", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: clear";
	await delayFrame();

	// @ts-ignore
	console.info(timeElapsed(), "download fa", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: download fa";

	// @ts-ignore
	console.info(timeElapsed(), "downloading fa", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	let loaded_fa = await (async function () {
		// preload all fa

		let fa_tasks = dataset.results.map(async function (obj, chrIdx) {
			if (analysis_options.nChr != 0 && analysis_options.nChr != (chrIdx + 1)) {
				return null;
			}
			//let obj = dataset.results[viewerState.nChr - 1];
			if (typeof obj == "string") {
				let raw_fa = await fetchData(obj, "text");
				let fa = await parseFasta(raw_fa);
				dataset.results[chrIdx] = fa;
			}
		});
		await Promise.all(fa_tasks);

		// 20201110
		if (dataset.mode == "single") {
			dataset.results.forEach((fa, chrIdx) => {
				if (analysis_options.nChr != (chrIdx + 1)) {
					return null;
				}
				const seqName = dataset.genome_info_list[0].chr_list[chrIdx].chr;
				dataset.results[chrIdx] = {
					[seqName]: fa[seqName].replace(/-/g, ""),
				};
			});
		}

		if (analysis_options.nChr == 0) {
			const keys = Object.keys(dataset.results[0]);

			keys.map(k => dataset.results.map(ss => ss[k]).join(""));

			const entries = dataset.genome_info_list.map((genome_info, genome_idx) => {
				return [
					genome_info.name,
					genome_info.chr_list.map(chr_info => dataset.results[chr_info.index - 1][chr_info.chr]).join(""),
				];
			});

			dataset.results[analysis_options.nChr - 1] = Object.fromEntries(entries);
		}
		else {
			return dataset.results[analysis_options.nChr - 1];
		}
	})();
	document.getElementById("status").innerText = "loading: downloaded fa";

	if (0 &&  dataset.mode != "tetrad") {
		const el_snpRef = document.querySelector(".markers.snp-reference");
		const el_snpTarget = document.querySelector("#data-rows .markers.snp-target")

		const fa = dataset.results[viewerState.nChr - 1];
		let ks = Object.keys(fa);//.map(a => a.replace("000000", ""));

		const mv_list = {
			"5": "ChI",
			"2": "ChII",
			"1": "ChIII",
			"7": "ChIV",
			"6": "ChV",
			"3": "ChVI",
			"8": "ChVII",
			"9": "MT",
		};

		ks = ks.map(a => {
			const nTig = a.match(/tig0000000(\d)/)[1];
			return a.replace(/tig0000000(\d)/, mv_list[nTig]);
		});

		el_snpRef.innerHTML = `<span>${ks.shift()}</span>`;
		el_snpTarget.innerHTML = `<span>${ks.shift()}</span>`;

		// for (let i = 0; i < ks.length - 1; ++i) {
		// 	el_snpTarget.after(el_snpTarget.cloneNode(true));
		// }
	}

	// @ts-ignore
	console.info(timeElapsed(), "init seq", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: init seq";

	if (dataset.gc_content) {
		let name_map = dataset.gc_content;

		let average = {};
		Object.keys(name_map).forEach(name => {
			average[name] = {};
			Object.keys(name_map[name]).forEach(chr =>{
				let first = name_map[name][chr][0];
				let w = Math.abs(first.end - first.start);
				average[name][chr] = name_map[name][chr].reduce((tt, curr) => {
					return tt + curr.gc / Math.abs(curr.end - curr.start);
				}, 0) / name_map[name][chr].length * w;
			});
		});
		// @ts-ignore
		gc_content = name_map;
		// @ts-ignore
		window.$gc_content = gc_content;
		// @ts-ignore
		gc_content_average = average;
		// @ts-ignore
		window.$gc_average = average;
	}

	//init
	g_chrCoord.bp_start = 1;
	g_chrCoord.bp_end = 20;
	el_input_start.value = g_chrCoord.bp_start.toString();
	el_input_end.value = g_chrCoord.bp_end.toString();

//////////////////////////////////////////////////////////////////////////////

	el_input_ref1_start.oninput = (evt) => {
		if (ref1_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			g_chrCoord.bp_start = ref1_pos_uint32array[newVal] | 0;
			g_chrCoord.bp_start = Math.max(1, g_chrCoord.bp_start);
			el_input_start.value = g_chrCoord.bp_start.toString();
			if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
				drawFrame();
			}
		}
	};
	el_input_ref1_end.oninput = (evt) => {
		if (ref1_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			g_chrCoord.bp_end = ref1_pos_uint32array[newVal] || seq_list[0].length;
			g_chrCoord.bp_end = Math.min(g_chrCoord.bp_end, seq_list[0].length);
			el_input_end.value = g_chrCoord.bp_end.toString();
			if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
				drawFrame();
			}
		}
	};

//////////////////////////////////////////////////////////////////////////////

	el_input_ref2_start.oninput = (evt) => {
		if (ref2_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			g_chrCoord.bp_start = ref2_pos_uint32array[newVal] | 0;
			g_chrCoord.bp_start = Math.max(1, g_chrCoord.bp_start);
			el_input_start.value = g_chrCoord.bp_start.toString();
			if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
				drawFrame();
			}
		}
	};
	el_input_ref2_end.oninput = (evt) => {
		if (ref2_pos_uint32array) {
			// @ts-ignore
			let newVal = Number(evt.target.value);
			g_chrCoord.bp_end = ref2_pos_uint32array[newVal] || seq_list[1].length;
			g_chrCoord.bp_end = Math.min(g_chrCoord.bp_end, seq_list[1].length);
			el_input_end.value = g_chrCoord.bp_end.toString();
			if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
				drawFrame();
			}
		}
	};

//////////////////////////////////////////////////////////////////////////////

	el_input_start.oninput = (evt) => {
		// @ts-ignore
		let newVal = Number(evt.target.value);
		g_chrCoord.bp_start = Math.max(1, newVal);
		el_input_start.value = g_chrCoord.bp_start.toString();
		if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
			drawFrame();
		}
	};
	el_input_end.oninput = (evt) => {
		// @ts-ignore
		let newVal = Number(evt.target.value);
		g_chrCoord.bp_end = Math.min(newVal, seq_list[0].length);
		el_input_end.value = g_chrCoord.bp_end.toString();
		if (g_chrCoord.bp_end > g_chrCoord.bp_start) {
			drawFrame();
		}
	};

	//seq_list clear
	if (reverse_parental) {
		let seq_id_list = Object.keys(loaded_fa);
		seq_id_list = [
			seq_id_list[1], seq_id_list[0],
			...seq_id_list.slice(2),
		];
		seq_id_list.forEach((id, i) => seq_list[i] = loaded_fa[id]);
	}
	else {
		if (analysis_options.nChr == 0) {
			dataset.genome_info_list.forEach((genome_info, genome_idx) => {
				const seq = dataset.results[analysis_options.nChr - 1][genome_info.name];
				if (seq) {
					seq_list[genome_idx] = seq;
				}
			});
		}
		else {
			// 20201110
			dataset.genome_info_list.map(info => info.chr_list[analysis_options.nChr - 1].chr).forEach((seq_id, i) => {
				const seq = dataset.results[analysis_options.nChr - 1][seq_id];
				if (seq) {
					seq_list[i] = seq;
				}
			});
			// let seq_id_list = dataset.genomeNameList;//Object.keys(loaded_fa);
			// seq_id_list.forEach((id, i) => seq_list[i] = loaded_fa[id]);
		}
	}
	
	g_chrCoord.setSeqList(seq_list);

	//////////////////////////////////////////////////////////////////////////////

	// @ts-ignore
	console.info(timeElapsed(), "loading: cmp", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: cmp";

	// if (analyser_mode == "tetrad") {
		if (false) {
			make_fa_image();
		}
		else if (seq_list.length > 1) {
			loadCmpData(seq_list, analysis_options);
		}
		else if (dataset.mode == "single") {
			// 20201110
			init_viewModel(seq_list);
			for (let i = 0; i < seq_list[0].length; ++i) {
				ref1_pos_uint32array[i] = i + 1;
				pos_ref1_uint32array[i] = i + 1;
			}

			const display_meth_ratio = viewerState.display_meth_ratio;
			const display_nBS_meth_ratio = viewerState.display_nBS_meth_ratio;
			const display_BS_nBS_meth_ratio = viewerState.display_BS_nBS_meth_ratio;
			const display_RNA_coverage = viewerState.display_RNA_coverage;
			const display_meth_diff_wt = viewerState.display_meth_diff_wt;
			const display_new_gff3 = viewerState.display_new_gff3;

			Object.keys(viewerState).map(k => k.match(/^_display|display/) && typeof viewerState[k] == "boolean" ? k : 0).filter(a => a).forEach(k => viewerState[k] = false);

			viewerState.display_meth_ratio = display_meth_ratio;
			viewerState.display_nBS_meth_ratio = display_nBS_meth_ratio;
			viewerState.display_BS_nBS_meth_ratio = display_BS_nBS_meth_ratio;
			viewerState.display_RNA_coverage = display_RNA_coverage;
			viewerState.display_meth_diff_wt = display_meth_diff_wt;
			viewerState.display_new_gff3 = display_new_gff3;
		}
	// }
	// else {
	// 	for (let i = 0; i < seq_list[0].length; i += 500) {
	// 		const column_d = Array(seq_list.length).fill(0);
	// 		for (let j = 0; j < 500; ++j) {
	// 			const idx = i + j;
	// 			const column = seq_list.map(seq => seq[idx]);
	// 			const ref_value = column[0];
	// 			let cmp = column.map(value => value == ref_value ? 0 : 1);
	// 			column_d.forEach((d, i) => column_d[i] += cmp[i]);
	// 		}
	// 	}
	// }

	make_subject_pos_to_ma(seq_list);
	make_subject_ma_to_pos(seq_list);

	// try {
	// 	// repeat_segment[0] = await LoadRepeatSegment(dataset.parental_list[0], analysis_options.nChr, 1, ref1_pos_uint32array);
	// 	const ref_pos_list = [
	// 		ref1_pos_uint32array,
	// 		ref2_pos_uint32array,
	// 	];
	// 	const tasks = dataset.parental_list.map(async (refId, ref_idx) => {
	// 		if (!repeat_segment[refId]) {
	// 			repeat_segment[refId] = [];
	// 		}
	// 		const rows = await LoadRepeatSegment(refId, analysis_options.nChr, 2, ref_pos_list[ref_idx]);
	// 		repeat_segment[refId].push(rows);
	// 	});
	// 	await Promise.all(tasks);
	// }
	// catch (ex) {
	// 	console.error(ex);
	// }


	//////////////////////////////////////////////////////////////////////////////

	// if (gff_data_map) {
	// 	const ref_pos_map_list = make_ref_pos_map_list();
	// 	dataset.parental_list.forEach((refName, si) => {
	// 		const gff_data = gff_data_map[refName];
	// 		const pos_map = ref_pos_map_list[si];
	// 		Object.values(gff_data).forEach(gene_list => {
	// 			gene_list.forEach(gene => {
	// 				// gene.$length = pos_map[gene.end - 1] - gene.$start + 1;
	// 				gene.$start = pos_map[gene.start - 1];
	// 				gene.$end = pos_map[gene.end - 1];
	// 				gene.$length = gene.$end - gene.$start + 1;
	// 			});
	// 		});
	// 	});
	// }

	await delayFrame();

	// @ts-ignore
	console.info(timeElapsed(), "loading: gv", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: gv";

	// // @ts-ignore
	// window.$seg = seg_snp;//no save
	// @ts-ignore
	//window.parental_cmp = parental_cmp_uint8array;
	// @ts-ignore
	window.allMarker = allMarker;
	// @ts-ignore
	window.ref1_pos_uint32array = ref1_pos_uint32array;
	// @ts-ignore
	window.pos_ref1_uint32array = pos_ref1_uint32array;
	// @ts-ignore
	window.ref2_pos_uint32array = ref2_pos_uint32array;
	// @ts-ignore
	window.pos_ref2_uint32array = pos_ref2_uint32array;

	//max length
	g_chrCoord.bp_end = Math.min(g_chrCoord.bp_end, seq_list[0].length);

	el_input_start.max = (seq_list[0].length - 1).toString();
	el_input_end.max = seq_list[0].length.toString();

	g_chrCoord.bp_start = 1;
	g_chrCoord.bp_end = seq_list[0].length;
	el_input_start.value = g_chrCoord.bp_start.toString();
	el_input_end.value = g_chrCoord.bp_end.toString();
	g_maxPixelPerBP = getPixelPerBP();

	//////////////////////////////////////////////////////////////////////////////

	// @ts-ignore
	console.info(timeElapsed(), "merge GC", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: merge GC";
	if (dataset._gc_content) {
		el_gc_content_window_size.innerHTML = "";//clear

		const max_merge_size = Math.ceil((seq_list[0].length / dataset.GC_Content_window) / viewerState.max_view_width);
		if (max_merge_size > 0 && Number.isFinite(max_merge_size)) {
			for (let i = 1; i <= max_merge_size; ++i) {
				const opt = document.createElement("option");
				const window_size = i * dataset.GC_Content_window;
				if (window_size <= max_gc_content_window) {
					opt.innerText = window_size.toFixed(0);// display window size
					opt.value = i.toFixed(0);// count of window
					el_gc_content_window_size.append(opt);
				}
			}

			el_gc_content_window_size.onchange = function (evt) {
				merge_gc_content(Number(el_gc_content_window_size.value), true);
			};

			merge_gc_content();
		}
	}
	await delayFrame();

	//////////////////////////////////////////////////////////////////////////////

	// allMarker = { list: [], };
	// console.error("clear", "allMarker.list");

	// return drawFrame();

	//////////////////////////////////////////////////////////////////////////////

	try {
		if (!CONFIG_LOAD_METH) {
			// @ts-ignore
			console.info(timeElapsed(), "loaded without meth", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
			document.getElementById("status").innerText = "loaded without meth";

			document.getElementById("CONFIG_LOAD_METH")?.remove();

			await delayFrame();
		}
		// else {
		// }

		// debugger;

		// viewerState.seg_row_separate = 15;

		await g_methylRenderer.load(methyl_dataset_list);

		// {
		// 	window.cel_f = new Float32Array(window.$data_map.get("BS-seq_20210225/QM6a+cel+Methyl.methratio.float32"))
		// 	window.glu_f = new Float32Array(window.$data_map.get("BS-seq_20210225/QM6a+glu+Methyl.methratio.float32"))
		// 	window.cel_sub_glu = cel_f.map((cv, ci) => {
		// 		let gv = window.glu_f[ci];
		// 		cv = Number.isNaN(cv) ? 0 : cv;
		// 		gv = Number.isNaN(gv) ? 0 : gv;
		// 		return Math.abs(Math.abs(cv) + Math.abs(gv))
		// 	})
		//
		// 	window._cel_sub_glu = _load_meth_ratio_float32("QM6a", window.cel_sub_glu, (v, i) => {
		// 		if (v) {
		// 			return v;
		// 		}
		// 	});
		// }

		g_methylRenderer.split_strand(methyl_dataset_list);

		if (!20210304) {
			g_methylRenderer.make_methyl_Chm_row();
		}
	}
	catch (ex) {
		console.error(ex);
	}

	// @ts-ignore
	console.info(timeElapsed(), "loading: RNA coverage", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "loading: RNA coverage";
	if (viewerState.display_RNA_coverage) {
		rnaCoverage_max_map = {};
		rnaCoverage_norm_map = {};
		w_rnaCoverage_ui32a_map = {};
		m_rnaCoverage_ui32a_map = {};

		const tasks = Object.keys(methy_refMap_sampleList).map(async function (refId) {
			// if (refId != "QM6a") {
			// 	console.warn("skip", refId, "RNA coverage");
			// 	return;
			// }
			rnaCoverage_max_map[refId] = {};//define prop // make ordered map
			rnaCoverage_norm_map[refId] = {};//define prop // make ordered map
			w_rnaCoverage_ui32a_map[refId] = {};//define prop // make ordered map
			m_rnaCoverage_ui32a_map[refId] = {};//define prop // make ordered map

			const _tasks = methy_refMap_sampleList[refId].map(/** @param {string} sampleName */async function (sampleName) {
				try {
					const fileName = `${refId}_${methy_sample_fileName_map[sampleName]}_coverage.uint32`;
					await loadCoverageData(refId, sampleName, fileName);
				}
				catch (ex) {
					console.error(ex);
				}
			});

			await Promise.all(_tasks);
		});
		await Promise.all(tasks);
		need_merge_rna_reads = false;

		// @ts-ignore
		window.rnaCoverage_max_map = rnaCoverage_max_map;
		// @ts-ignore
		window.rnaCoverage_norm_map = rnaCoverage_norm_map;
		// @ts-ignore
		window.rnaCoverage_float32a_map = w_rnaCoverage_ui32a_map;
		// @ts-ignore
		window.m_rnaCoverage_float32a_map = m_rnaCoverage_ui32a_map;

		console.log("loaded", "RNA coverage");
		await delayFrame();
	}
	else {
		rnaCoverage_max_map = {};
		rnaCoverage_norm_map = {};
		w_rnaCoverage_ui32a_map = {};
		m_rnaCoverage_ui32a_map = {};

		const tasks = Object.keys(methy_refMap_sampleList).map(async function (refId) {
			// if (refId != "QM6a") {
			// 	console.warn("skip", refId, "RNA coverage");
			// 	return;
			// }
			rnaCoverage_max_map[refId] = {};//define prop // make ordered map
			rnaCoverage_norm_map[refId] = {};//define prop // make ordered map
			w_rnaCoverage_ui32a_map[refId] = {};//define prop // make ordered map
			m_rnaCoverage_ui32a_map[refId] = {};//define prop // make ordered map
		});
		need_merge_rna_reads = false;

		// @ts-ignore
		window.rnaCoverage_max_map = rnaCoverage_max_map;
		// @ts-ignore
		window.rnaCoverage_norm_map = rnaCoverage_norm_map;
		// @ts-ignore
		window.rnaCoverage_float32a_map = w_rnaCoverage_ui32a_map;
		// @ts-ignore
		window.m_rnaCoverage_float32a_map = m_rnaCoverage_ui32a_map;

		console.log("skip", "RNA coverage");
	}
	// if (0) {
	// 	const rna_coverage_infoMap = {
	// 		"WT(veg)": "Q_wt_coverage.uint32",
	// 		"rid(veg)": "Q_rid_coverage.uint32",

	// 		"WT(D4)": "Q_d4_coverage.uint32",
	// 		"rid(D4)": "Q_d4r_coverage.uint32",

	// 		"WT(D8)": "Q_d8_coverage.uint32",
	// 		"rid(D8)": "Q_d8r_coverage.uint32",

	// 		// CBS1_coverage   QM6a1_coverage   ridD83_coverage  WTD81_coverage
	// 		// Crid1_coverage  ridD41_coverage  WTD41_coverage
	// 	};
	// 	rnaCoverage_ui32a_map = Object.assign({}, rna_coverage_infoMap);//copy key
	// 	m_rnaCoverage_ui32a_map = Object.assign({}, rna_coverage_infoMap);//copy key
	// 	//
	// 	const promise_rna_coverage = Object.keys(rna_coverage_infoMap).map(async function (coverageName) {
	// 		const coverageFileName = rna_coverage_infoMap[coverageName];

	// 		/** @type {ArrayBuffer} */
	// 		const buffer = await fetchData(coverageFileName, "arraybuffer");

	// 		/** @type {number[]} */
	// 		const chrLengthList = dataset.genome_info_list[0].chr_list.map(a => a.length);

	// 		const chrOffsetInByte = [];
	// 		chrOffsetInByte[0] = 0;
	// 		chrLengthList.reduce((pos, len, i) => {
	// 			let offset = pos + (len * 4);
	// 			chrOffsetInByte[i + 1] = offset;
	// 			return offset;
	// 		}, 0);

	// 		const chrIdx = analysis_options.nChr - 1;
	// 		const offsetStart = chrOffsetInByte[chrIdx];

	// 		const ui32a = new Uint32Array(buffer, offsetStart, chrLengthList[chrIdx]);
	// 		const m_ui32a = new Uint32Array(Math.ceil(ui32a.length / rna_window_size));

	// 		let max = 1;
	// 		for (let i = 0; i < ui32a.length; ++i) {
	// 			max = Math.max(max, ui32a[i]);
	// 		}
	// 		console.log(coverageName, "max value", max);

	// 		rnaCoverage_max_map[coverageName] = max;
	// 		rnaCoverage_norm_map[coverageName] = max;

	// 		rnaCoverage_ui32a_map[coverageName] = ui32a;
	// 		m_rnaCoverage_ui32a_map[coverageName] = m_ui32a;

	// 		document.getElementById("status").innerText = "loading: loaded:" + coverageFileName;
	// 		await delayFrame();
	// 	});
	// 	await Promise.all(promise_rna_coverage);

	// 	// @ts-ignore
	// 	window.rnaCoverage_float32a_map = rnaCoverage_ui32a_map;
	// 	// @ts-ignore
	// 	window.m_rnaCoverage_float32a_map = m_rnaCoverage_ui32a_map;

	// 	console.log("loaded", "RNA coverage");
	// 	await delayFrame();
	// }

	// {
	// 	//const [k1, k2] = Object.keys(rnaCoverage_float32a_map);
	// 	const length = rnaCoverage_float32a_map[k1].length;
	// 	rnaCoverage_cmp_u8a_map[coverageName]
	// 	function cov_cmp(name1, name2) {
	// 		let red = 0;
	// 		for (let i = 0; i < length; ++i) {
	// 			const v1 = rnaCoverage_float32a_map[name1][i];
	// 			const c2 = rnaCoverage_float32a_map[name2][i] * 2;
	// 			if (v1 > c2) {
	// 				red += 1;
	// 			}
	// 		}
	// 		console.log(`${name1} > ${name2}`, red);
	// 	}
	// }
	// Object.keys(rnaCoverage_float32a_map).forEach(name => {
	// 	const f32a = rnaCoverage_float32a_map[name];
	// 	f32a.forEach((v, i) => {
	// 		f32a[i] = v ** 0.6;
	// 	});
	// });

	// @ts-ignore
	console.info(timeElapsed(), "loaded all", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000_000_000);
	document.getElementById("status").innerText = "";

	await delayFrame();

	// {
	// 	calc_methy();
	// }
}

class Plot_Helper {
	/** @type {ImageData} */
	imgData = null;

	/** @type {HTMLCanvasElement} */
	wCanvas = null;

	/** @type {CanvasRenderingContext2D} */
	wCtx = null;

	_split_row = 0;

	constructor() {
		/** @type {HTMLCanvasElement} */
		this.wCanvas = (document.createElement("canvas"));
		this.wCtx = this.wCanvas.getContext("2d");
	}

	/**
	 * @param {number} num_of_row
	 * @param {number} num_of_col
	 * @param {(x: number, y: number) => number[]} fn_getColor
	 */
	make(num_of_row, num_of_col, fn_getColor) {
		const { wCanvas, wCtx } = this;
		const { width, height, split_row } = this._get_image_size(num_of_col, num_of_row);

		const imgData = wCtx.createImageData(width, height);

		this.imgData = imgData;

		this._split_row = split_row;
		wCanvas.width = width;
		wCanvas.height = height;

		// wCtx.imageSmoothingQuality = "high";
		wCtx.imageSmoothingEnabled = false;

		for (let x = 0; x < width; ++x) {
			for (let yy = 0; yy < split_row; ++yy) {
				const y0 = yy * num_of_row;
				for (let dy = 0; dy < num_of_row; ++dy) {
					const y = y0 + dy;
					imgData_setColor(x, y, fn_getColor(x + yy * split_row, dy));
				}
			}
		}

		wCtx.putImageData(this.imgData, 0, 0);

		/**
		 * @param {number} x
		 * @param {number} y
		 * @param {number[]} param2
		 */
		function imgData_setColor(x, y, [r, g, b, a]) {
			let i = (x + imgData.width * y) * 4;
			imgData.data[i + 0] = r;
			imgData.data[i + 1] = g;
			imgData.data[i + 2] = b;
			imgData.data[i + 3] = a;
		}
	}

	/**
	 * @param {number} num_of_col
	 * @param {number} num_of_row
	 * @returns {{ width: number, height: number, split_row: number }}
	 */
	_get_image_size(num_of_col, num_of_row) {
		let length = num_of_col * num_of_row;
		let split_row = Math.ceil(Math.ceil(length ** 0.5) / num_of_row);
		let height = split_row * num_of_row;
		let width = Math.ceil(length / height);
		return {
			width,
			height,
			split_row,//split height
		};
	}

	/**
	 *
	 * @param {CanvasRenderingContext2D} ctx
	 * @param {number} row_idx
	 * @param {number} view_width
	 * @param {number} row_height
	 * @param {number} row_separate
	 * @param {number} view_start
	 * @param {number} view_end
	 * @param {number} num_of_col
	 * @param {number} num_of_row
	 * @param {(x: number, y: number) => any} predicate get detail
	 */
	draw_row(ctx, row_start_idx, row_end_idx, view_width, row_height, row_separate, view_start, view_end, num_of_col, num_of_row, predicate) {
		ctx.save();
		try {
			const view_length = view_end - view_start + 1;
			// const scale = 1 / view_length;
			// const bp_size = view_width * scale;

			const width = this.wCanvas.width;
			const split_row = this._split_row;

			//ctx.scale(10, 32);
			ctx.scale(view_width / view_length, 1);
			ctx.translate(-view_start, 0);
			for (let yy = 0; yy < split_row; ++yy) {
				let start_pos = yy * width;
				const sx = 0;
				const sy = yy * num_of_row + row_start_idx;
				const sw = width;
				const sh = (row_end_idx - row_start_idx) + 1;
				const dx = start_pos;
				const dy = (row_height + row_separate) * row_start_idx;
				const dw = width;
				const dh = row_height;
				ctx.drawImage(this.wCanvas, sx, sy, sw, sh, dx, dy, dw, dh);
			}

			if (predicate) {
				const fw = view_width / (view_end - view_start + 1);
				if (Math.ceil(fw) > 9) {
					ctx.font = `${Math.min(fw, row_height)}px Arial`;

					ctx.textAlign = "center";
					ctx.textBaseline = "bottom";

					for (let i = 0; i < num_of_row; ++i) {
						const dy = (row_height + row_separate) * i;

						for (let pos = 0; pos < (view_end - view_start + 1); ++pos) {
							const t = predicate(view_start + pos, i);

							const x = pos * fw;

							ctx.beginPath();
							ctx.strokeStyle = "#000000";
							ctx.strokeRect(x, dy, fw, row_height);

							ctx.fillText(t, x + fw * 0.5, dy + row_height);
						}
					}
				}
			}
		}
		finally {
			ctx.restore();
		}
	}
}

/**
 * @param {string} refId
 * @param {string} fileName
 */
async function loadChrUint32(refId, fileName) {
	/** @type {number[]} */
	const chrLengthList = dataset.genome_info_list[get_genome_index(refId)].chr_list.map(a => a.length);

	/** @type {ArrayBuffer} */
	const buffer = await fetchData(fileName, "arraybuffer");

	const chrOffsetInByte = [];
	chrOffsetInByte[0] = 0;
	chrLengthList.reduce((pos, len, i) => {
		let offset = pos + (len * 4);
		chrOffsetInByte[i + 1] = offset;
		return offset;
	}, 0);

	if (buffer.byteLength != chrLengthList.reduce((t, v) => t + v, 0) * Uint32Array.BYTES_PER_ELEMENT) {
		debugger;
	}

	return chrOffsetInByte.map((offsetStart, chrIdx) => {
		return new Uint32Array(buffer, offsetStart, chrLengthList[chrIdx]);
	});
}

/**
 * loadRNAMapping
 * loadDNAMapping
 * @param {"QM6a"|"CBS1-1"} refId
 * @param {string} sampleName
 */
async function loadCoverageData(refId, sampleName, fileName) {
	try {
		rnaCoverage_max_map[refId] = rnaCoverage_max_map[refId] ?? {};
		rnaCoverage_norm_map[refId] = rnaCoverage_norm_map[refId] ?? {};
		w_rnaCoverage_ui32a_map[refId] = w_rnaCoverage_ui32a_map[refId] ?? {};
		m_rnaCoverage_ui32a_map[refId] = m_rnaCoverage_ui32a_map[refId] ?? {};

		rnaCoverage_max_map[refId][sampleName] = null; //define prop // make ordered map
		rnaCoverage_norm_map[refId][sampleName] = null; //define prop // make ordered map
		w_rnaCoverage_ui32a_map[refId][sampleName] = null; //define prop // make ordered map
		m_rnaCoverage_ui32a_map[refId][sampleName] = null; //define prop // make ordered map

		/** @type {number[]} */
		const chrLengthList = dataset.genome_info_list[get_genome_index(refId)].chr_list.map(a => a.length);
		chrLengthList[-1] = chrLengthList.reduce((aa, v) => aa + v, 0);

		if (fileName) {
			// @ts-ignore
			console.info(timeElapsed(), "loading: RNA coverage:" + sampleName, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);

			/** @type {ArrayBuffer} */
			const buffer = await fetchData(fileName, "arraybuffer");
			document.getElementById("status").innerText = "loading: downloaded RNA coverage:" + sampleName;

			const chrOffsetInByte = [];
			chrOffsetInByte[0] = 0;
			chrLengthList.reduce((pos, len, i) => {
				let offset = pos + (len * 4);
				chrOffsetInByte[i + 1] = offset;
				return offset;
			}, 0);

			if (buffer.byteLength != chrLengthList.reduce((t, v) => t + v, 0) * Uint32Array.BYTES_PER_ELEMENT) {
				debugger;
			}

			const chrIdx = analysis_options.nChr - 1;
			const offsetStart = analysis_options.nChr == 0 ? 0 : chrOffsetInByte[chrIdx];

			const ui32a = new Uint32Array(buffer, offsetStart, chrLengthList[chrIdx]);
			const m_ui32a = new Uint32Array(Math.ceil(ui32a.length / rna_window_size));

			let max = 1;
			for (let i = 0; i < ui32a.length; ++i) {
				max = Math.max(max, ui32a[i]);
			}

			rnaCoverage_max_map[refId][sampleName] = max;
			rnaCoverage_norm_map[refId][sampleName] = 1;

			w_rnaCoverage_ui32a_map[refId][sampleName] = ui32a;
			m_rnaCoverage_ui32a_map[refId][sampleName] = m_ui32a;
		}
		else {
			console.warn("loadCoverageData", "if (fileName) {", fileName);

			const chrIdx = analysis_options.nChr - 1;

			const ui32a = new Uint32Array(chrLengthList[chrIdx]);
			const m_ui32a = new Uint32Array(Math.ceil(ui32a.length / rna_window_size));

			w_rnaCoverage_ui32a_map[refId][sampleName] = ui32a;
			m_rnaCoverage_ui32a_map[refId][sampleName] = m_ui32a;
		}

		//if (need_merge_rna_reads) {
		_merge_rna_reads(refId, sampleName);
		//}
		// @ts-ignore
		console.info(timeElapsed(), "loading: loaded RNA coverage:" + sampleName, "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		document.getElementById("status").innerText = "loading: loaded RNA coverage:" + sampleName;
		await delayFrame();

		return true;
	}
	catch (ex) {
		console.error(ex);
	}

	return false;
}

async function load_jgi_QM6a_gff() {
	const refName = "jgi_Trire_Chr";

	await load_ref_GFF({
		[refName]: {
			pos_map: ref1_pos_uint32array,
			file: "data/Trire/Trire_Chr.ExternalModels.gff3",
		},
	});

	const chr_gff = gff_data_map[refName][`scaffold_${viewerState.nChr}`];

	/** @see {draw_meth_ratio_value_drawArrow} */
	chr_gff.forEach(gff_data => {
		if (gff_data.type == "gene") {
			gff_data.attributes ??= {};
			if (!gff_data.attributes.Name) {
				gff_data.attributes.Name = gff_data.attributes?.ID;
			}
		}
	});
	
	const gene_row = _show_gene_on_row(refName, chr_gff);
	gene_row.name += " jgi";
	methyl_dataset_list.splice(1, 0, gene_row);
}

function show_gene_on_row(refIdx) {
	const chr_idx = viewerState.nChr - 1;
	
	const genome_info = dataset.genome_info_list[refIdx];
	const refName = genome_info.name;
	const sChr = genome_info.chr_list[chr_idx].symbol;

	const chr_gff = gff_data_map[refName][sChr];

	return _show_gene_on_row(refName, chr_gff);
}
/**
 * @param {string} refName
 * @param {GFF_ROW[]} chr_gff
 */
function _show_gene_on_row(refName, chr_gff) {
	const chr_idx = viewerState.nChr - 1;

	const gene_list = chr_gff.filter(gg => {
		return gg.type == "gene" && gg.$start <= g_chrCoord.bp_end && gg.$end >= g_chrCoord.bp_start;
	});

	const cfg = new module_Methyl_sampleData({
		ref: refName,
		sample: `${refName} gene`,
		name: `${refName} gene`,
		url: null,
		region: false,
		mid_line: false,
		density_to_opacity: false,
		data: [],
		rendering_condition: [
			new module_Methyl_sampleData_RenderingCondition({ color: "#04F7", condition: v => v, min_width: 1, }),// passthrough NaN
		],
		getValueTransformer: function () { return function (v) { return v; }; },
		// getValueTransformer: null,

		row_height: 3,
		drawShape: draw_meth_ratio_value_drawArrow,
	});

	cfg.data[chr_idx] = [];

	gene_list.forEach(gg => {
		const dd = new module_Methyl_ratioData();
		Object.assign(dd, {
			start: gg.start,
			end: gg.end,
			strand: gg.strand,
			// value: value,
			value: gg,
		});
		cfg.data[chr_idx].push(dd);
	});

	// methyl_dataset_list.push(cfg);

	return cfg;
}

/**
 * 20221017
 * @param {number} x1
 * @param {number} x2
 * @param {number} min_width
 * @param {number} row_height
 * @param {GFF_ROW} gff_data value = strand * abs(value)
 * @param {boolean} display_minus_value no use
 */
function draw_meth_ratio_value_drawArrow(x1, x2, min_width, row_height, gff_data, display_minus_value) {
	return _draw_meth_ratio_value_drawArrow(x1, x2, min_width, row_height, gff_data, display_minus_value);
}

/**
 * 20221017
 * @param {number} x1
 * @param {number} x2
 * @param {number} min_width
 * @param {number} row_height
 * @param {GFF_ROW} gff_data value = strand * abs(value)
 * @param {boolean} display_minus_value no use
 */
function _draw_meth_ratio_value_drawArrow(x1, x2, min_width, row_height, gff_data, display_minus_value) {
	const ctx = main_ctx;
	const rhh = row_height * 0.5;
	const shape_height = rhh * 0.5;

	const cx = (x1 + x2) * 0.5;

	const width = Math.max(min_width, x2 - x1);
	const strand = gff_data.strand;
	const dname = gff_data.attributes.Name || "";// || gff_data.attributes.ID || "";

	const text_width = width * 1.1;// Math.min(width, (shape_height * 0.5) * dname.length);
	
	ctx.font = `${gff_data.fontStyle || ""} ${Math.trunc(shape_height * 0.9)}px arial`;
	ctx.textAlign = "center";
	
	const color = "black";
	const stroke_color = "transparent";

	if (strand > 0) {
		const y_up = rhh - rhh * 0.25;
		drawArrow(ctx, x1, y_up, strand, width, shape_height, color, stroke_color);

		// ctx.strokeRect(x1, y_up, width, shape_height);

		ctx.textBaseline = "bottom";
		ctx.fillStyle = "black";
		ctx.beginPath();
		ctx.fillText(dname, cx, y_up - shape_height * 0.5, text_width);
	}
	else if (strand < 0) {
		const y_down = rhh + rhh * 0.25;
		
		drawArrow(ctx, x1, y_down, strand, width, shape_height, color, stroke_color);

		ctx.textBaseline = "top";
		ctx.fillStyle = "black";
		ctx.beginPath();
		ctx.fillText(dname, cx, y_down + shape_height * 0.5, text_width);
	}
	else {
		const y_center = rhh;
		drawArrow(ctx, cx, y_center, strand, width, shape_height, color, stroke_color);
	}
}

/**
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} x1 - start
 * @param {number} y - arrow center
 * @param {-1|1} strand - arrow orientation
 * @param {number} width
 * @param {number} height
 * @param {string} color
 * @param {string} strokeStyle
 */
function drawArrow(ctx, x1, y, strand, width, height, color, strokeStyle) {
	// let arrow_height = Math.ceil(24 * g_row_height_scale);
	// let triangle_height = arrow_height / 2;
	
	const tr_h = (height * 0.5) * (3 ** 0.5);
	// if (width <= (tr_h * 0.5) || width < 5) {
	// 	return draw_triangle(ctx, x1, y, strand, width, height, color, strokeStyle);
	// }
	let x2 = x1 + width;
	let hw = width * 0.5;
	//
	let bw = width - tr_h;
	// if (bw <= 0) {
	// 	return draw_triangle(ctx, x1, y, strand, width, height, color, strokeStyle);
	// }
	if (bw <= tr_h) {
		bw = width * 0.5;
	}
	// let bw = width * 0.6;
	// let aw = width * 0.4;
	//
	// aw = Math.min(width * 0.4, hw * (3 ** 0.5), tr_h);
	// if (aw > 5) {
	// 	bw = width - aw;
	// }
	// else {
	// 	bw = width * 0.4;
	// }
	//
	let cx = (x2 + x1) * 0.5;
	ctx.beginPath();
	
	// ctx.moveTo(Math.round(cx + (-hw) * strand), Math.trunc(y - height / 4));
	// ctx.lineTo(Math.round(cx + (-hw + bw) * strand), Math.trunc(y - height / 4));
	// ctx.lineTo(Math.round(cx + (-hw + bw) * strand), Math.trunc(y - height / 2));
	// ctx.lineTo(Math.round(cx + hw * strand), Math.round(y));
	// ctx.lineTo(Math.round(cx + (-hw + bw) * strand), Math.ceil(y + height / 2));
	// ctx.lineTo(Math.round(cx + (-hw + bw) * strand), Math.ceil(y + height / 4));
	// ctx.lineTo(Math.round(cx + (-hw) * strand), Math.ceil(y + height / 4));
	// ctx.closePath();
	ctx.moveTo(cx + (-hw) * strand, y - height / 4);
	ctx.lineTo(cx + (-hw + bw) * strand, y - height / 4);
	ctx.lineTo(cx + (-hw + bw) * strand, y - height / 2);
	ctx.lineTo(cx + hw * strand, y);
	ctx.lineTo(cx + (-hw + bw) * strand, y + height / 2);
	ctx.lineTo(cx + (-hw + bw) * strand, y + height / 4);
	ctx.lineTo(cx + (-hw) * strand, y + height / 4);
	ctx.closePath();

	// notice stroke and fill real size
	
	// ctx.lineWidth = Math.max(1, g_lineWidth * 0.5);
	if (strokeStyle) {
		ctx.strokeStyle = strokeStyle;
	}
	else {
		ctx.strokeStyle = color;
	}
	ctx.stroke();
	// ctx.lineWidth = g_lineWidth;

	// if (g_displayGeneOutline) {
	// 	ctx.strokeStyle = "black";
	// 	ctx.stroke();
	// }
	
	ctx.fillStyle = color;
	ctx.fill();

	//ctx.globalAlpha = 1;
}

/**
 * @requires ref1_pos_uint32array
 * @requires ref2_pos_uint32array
 */
async function load_QM6a_CBS11_gff(version = "fun3") {
	try {
		console.info(timeElapsed(), "loading: gff", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);
		document.getElementById("status").innerText = "loading: gff";
		// if (analysis_options.mode != "SNP" && !gff) {
		{
			const pos_map_map = {
			};
			if (version == "fun3") {//20221117
				pos_map_map[dataset.parental_list[0]] = {
					pos_map: ref1_pos_uint32array,
					file: "/common_data/Trichoderma_reesei_QM6a.gff3",
					// // file: "/20200720_v3_QCt/20200720_v3_QCt/QM6a.gff3",
					// // // file: "Trichoderma_reesei_QM6a.gff3",
					seqid_to_chr: Object.fromEntries(dataset.genome_info_list[0].chr_list.map((info, chr_idx) => [`scaffold_${(chr_idx + 1)}`, info.symbol])),
					// file: "/annotate/new_QM6a.gff3",
				};
			}
			else {
				pos_map_map[dataset.parental_list[0]] = {
					pos_map: ref1_pos_uint32array,
					file: "/20200720_v3_QCt/20200720_v3_QCt/QM6a.gff3",
					// // // file: "Trichoderma_reesei_QM6a.gff3",
				};
			}
			// [dataset.parental_list[0] + "_new"]: {
			// 	pos_map: ref1_pos_uint32array,
			// 	// file: "QM6a.gff3",
			// 	file: `${LNCRNA_ROOT_DIR}/Trichoderma_reesei_QM6a.gff3`,
			// },
			if (dataset.parental_list[1]) {
				if (version == "fun3") {
					pos_map_map[dataset.parental_list[1]] = {
						pos_map: ref2_pos_uint32array,
						file: "/common_data/Trichoderma_reesei_CBS1.gff3",
						// // file: "/20200720_v3_QCt/20200720_v3_QCt/Trichoderma_reesei_CBS.gff3",
						seqid_to_chr: Object.fromEntries(dataset.genome_info_list[1].chr_list.map((info, chr_idx) => [`CBS1-1_${(chr_idx + 1)}`, info.symbol])),
						// file: "/annotate/new_CBS1.gff3",
					};
				}
				else {
					pos_map_map[dataset.parental_list[1]] = {
						pos_map: ref2_pos_uint32array,
						file: "/20200720_v3_QCt/20200720_v3_QCt/Trichoderma_reesei_CBS.gff3",
						seqid_to_chr: Object.fromEntries(dataset.genome_info_list[1].chr_list.map((info, chr_idx) => [`Ch${(chr_idx + 1)}_CBS1-1`, info.symbol])),
					};
				}
			}

			// @ts-ignore
			console.info(timeElapsed(), "loaded: gff", "totalJSHeapSize", performance?.memory?.totalJSHeapSize / 1000000000);

			await load_ref_GFF(pos_map_map);

			document.getElementById("status").innerText = "loaded: gff";
			await delayFrame();

			// try {
			// 	loadGeneTable(gff);
			// 	await delayFrame();
			// }
			// catch (ex) {
			// 	console.error(ex);
			// }
		}
	}
	catch (ex) {
		console.error(ex);
	}
}

/**
 * @param {{ [name: string]: { pos_map: (number[]|Uint32Array); file: string; seqid_to_chr?: { [seqid: string]: String }; }; }} pos_map_map
 */
async function load_ref_GFF(pos_map_map) {
	// if (typeof window.parseGFF != "function" || document.getElementById("/src/gff.js") == null || typeof GFF_ROW == "undefined") {
	if (document.getElementById("/src/gff.js") == null) {
		const sc = document.createElement("script");

		const task = new Promise((resolve, reject) => {
			sc.onload = resolve;
			sc.onerror = reject;
			sc.src = "/src/gff.js";
			sc.async = true;
			sc.id = "/src/gff.js";

			document.body.append(sc);
		});

		await task;
	}

	if (gff_data_map == null) {
		gff_data_map = {};
	}
	// if (gff_data_queryMap == null) {
	// 	gff_data_queryMap = {};
	// }

	const task_list = Object.keys(pos_map_map).map(async function (key) {
		// gff_data_queryMap[key] = {};
		try {
			const text = await fetchData(pos_map_map[key].file, "text");
			const gff_data = parseGFF(text, pos_map_map[key].seqid_to_chr);
			Object.entries(gff_data).forEach(([subjectName, gene_list]) => {
				
				const pos_map = pos_map_map[key].pos_map;
				if (pos_map) {
					gene_list.forEach(gene => {
						// gene.$length = pos_map[gene.end - 1] - gene.$start + 1;
						gene.$start = pos_map[gene.start - 1];
						gene.$end = pos_map[gene.end - 1];
						gene.$length = gene.$end - gene.$start + 1;
					});
				}

				// if (dataset.genome_info_list.find(a => a.name == key).chr_list[viewerState.nChr - 1].symbol == subjectName) {
				// 	gff_data_queryMap[key][subjectName] = Array.from(new Array(Math.trunc(pos_map.length / 100))).map(_ => []);
				// 	gene_list.forEach(gene => {
				// 		gff_data_queryMap[key][subjectName][Math.trunc(gene.$start / 100)].push(gene);
				// 	});
				// }
			});
			Object.keys(gff_data).forEach(subjectName => {
				const mm = new Map(gff_data[subjectName].map(a => [
					a.geneID, {
						geneID: a.geneID,
						type: "protein-coding gene",// ATG-STOP
						strand: a.strand,
						start: Infinity,
						end: 0,
						$start: Infinity,
						$end: 0,
					}
				]));
				gff_data[subjectName].forEach((feat, k) => {
					if (feat.type == "CDS") {
						const gg = mm.get(feat.geneID);

						gg.start = Math.min(gg.start, feat.start);
						gg.end = Math.max(gg.end, feat.end);

						gg.$start = Math.min(gg.$start, feat.$start);
						gg.$end = Math.max(gg.$end, feat.$end);
					}
				});
				[...mm.values()].filter(a => a.start < Infinity && a.end > 0).forEach(feat => {
					gff_data[subjectName].push(feat);
				});

				gff_data[subjectName] = gff_data[subjectName].sort(_func_sort_GTF);
			});
			gff_data_map[key] = gff_data;
			// return {
			// 	key,
			// 	gff: gff_data,
			// };
		}
		catch (ex) {
			console.error(ex);
		}
	});

	await Promise.all(task_list);

	// const gff_list = await Promise.all(task_list);
	// gff_data_map = gff_list.reduce((obj, pair_gff) => {
	// 	if (pair_gff) {
	// 		obj[pair_gff.key] = pair_gff.gff;
	// 	}
	// 	return obj;
	// }, {});
}

function calc_methy() {
	let t1 = new Date();
	console.log("methy start", t1);

	// const sp_size_list = [
	// 	500, 3_000, 5_000, 8_000, 10_000, 16_000, 20_000, 24_000, 30_000,
	// ];
	const methy_in_block_info_list = [];
	for (let sp_size = 500; sp_size <= 30_000; sp_size += 500) {
	//for (let sp_size of sp_size_list) {
		const block_list = [...seq_list[1].matchAll(new RegExp(`-{${ sp_size },}`, "g"))];
		const list = block_list.map(a => {
			const mList = methy_list[4].filter(b => {
				return (
					b.value < 0 &&	// has methy
					b.pos >= a.index && b.pos <= (a.index + a[0].length) // in range
				);
			});
			return {
				block_pos: a.index,
				block_size: a[0].length,
				mList: mList,
			};
		}).filter(c => c.mList.length > 0);

		list.sort((m1, m2) => m2.mList.length - m1.mList.length);

		const info = {
			sp_size: sp_size,
			block_list: block_list,
			methy_in_block_list: list,
		};
		methy_in_block_info_list.push(info);

		if (block_list.length == list.length) {
			break;
		}
	}
	let out_tab = "";
	out_tab += ["QM6a-specific size", "methy & block", "block", "%", "mmp", "mms", "mn"].join("\t") + "\r\n";
	methy_in_block_info_list.forEach(info => {
		const n_b = info.block_list.length;
		const n_mib = info.methy_in_block_list.length;
		const line = [
			info.sp_size,
			n_mib,
			n_b,
			n_mib / n_b,
			info.methy_in_block_list[0].block_pos,
			info.methy_in_block_list[0].block_size,
			info.methy_in_block_list[0].mList.length,
		].join("\t") + "\r\n";
		console.log(line);
		out_tab += line;
	});
	console.info(out_tab);

	let t2 = new Date();
	console.log("methy time", Math.trunc((t2.getTime() - t1.getTime()) / 1000));
}

document.getElementById("el_show_all").onclick = show_all;
async function show_all() {
	g_chrCoord.bp_start = 1;
	g_chrCoord.bp_end = seq_list[0].length;
	el_input_start.value = g_chrCoord.bp_start.toString();
	el_input_end.value = g_chrCoord.bp_end.toString();

	merge_gc_content();

	await delayFrame();
	await _drawFrame();
	viewerState.resizeCanvas();
	await _drawFrame();
}

async function _show_all() {
	g_chrCoord.bp_start = 1;
	g_chrCoord.bp_end = seq_list[0].length;
	el_input_start.value = g_chrCoord.bp_start.toString();
	el_input_end.value = g_chrCoord.bp_end.toString();

	merge_gc_content();

	await _drawFrame();
}

function render_UI() {
	ui_ctx.setTransform(1, 0, 0, 1, 0, 0);
	ui_ctx.clearRect(0, 0, main_canvas.width, main_canvas.height);
	// ui_ctx.clearRect(0, 0, ui_canvas.width, ui_canvas.height);
	ui_ctx.drawImage(main_canvas, 0, 0);

	ui_ctx.setTransform(1, 0, 0, 1, 0.5, 0.5);

	if (g_chrCoord.b_draw_selected_range) {
		const max_view_width = viewerState.max_view_width;

		const view_length = g_chrCoord.bp_end - g_chrCoord.bp_start + 1;
		const scale = 1 / view_length;
		const bp_size = max_view_width * scale;

		const [mouse_bp_pos_1, mouse_bp_pos_2] = [g_chrCoord.mouse_bp_pos, g_chrCoord.mouse_bp_pos_0].sort((a, b) => a - b);

		const x1 = ((1 - g_chrCoord.bp_start) + mouse_bp_pos_1 - 1) * bp_size;
		const x2 = ((1 - g_chrCoord.bp_start) + mouse_bp_pos_2) * bp_size;

		const y1 = 0;
		const y2 = main_ctx.getTransform().f + 100;

		ui_ctx.beginPath()
		ui_ctx.moveTo(x1, y1);
		ui_ctx.lineTo(x2, y1);
		ui_ctx.lineTo(x2, y2);
		ui_ctx.lineTo(x1, y2);
		ui_ctx.closePath();

		ui_ctx.fillStyle = "#00FF004F";
		ui_ctx.fill();
	}
}

/**
 * @param {number} time time in second
 */
async function _drawFrame(time) {
	try {
		main_ctx.imageSmoothingQuality = "high";
		main_ctx.imageSmoothingEnabled = false;

		if (promise_load_task) {
			return;
		}
		// if (analysis_options.mode == "tetrad") {
			await render_tetrad(analysis_options);
		// }
		// else {
		// 	render_SNP(analysis_options);
		// }

		render_UI();
	}
	catch (ex) {
		console.error(ex);
	}
	
	// 20211220
	// await viewerState.plot_list.reduce(async function (prev_promise, plot) {
	// 	await prev_promise;
	// 	appendLeftTitle(plot.title, plot.row_height, plot.row_separate, plot.rowspan, plot.p_fontSize);
	// 	await plot.func(plot);
	// 	await delayFrame();
	// }, Promise.resolve());
	viewerState.plot_list.forEach(function (plot) {
		appendLeftTitle(plot.title, plot.row_height, plot.row_separate, plot.rowspan, plot.p_fontSize);
		plot.func(plot);
	});

	let bUpdateCanvas;
	if (viewerState.auto_padding_right) {
		if (viewerState.padding_right != (viewerState.refresh_padding_right + 16)) {
			viewerState.padding_right = (viewerState.refresh_padding_right + 16);
			bUpdateCanvas = true;
		}
	}
	if (viewerState.resizeCanvas() || bUpdateCanvas) {
		await _drawFrame(time);
	}
}

async function drawFrame() {
	if (promise_load_task) {
		return false;
	}

	await viewerState.$promise_prevFrame;
	viewerState.$promise_prevFrame = null;

	if (!viewerState.$animationFrameId) {
		viewerState.$promise_prevFrame = new Promise((resolve, reject) => {
			const animationFrameId = requestAnimationFrame(async time => {
				viewerState.$animationFrameId = animationFrameId;
				await _drawFrame(time);
				resolve();
			});
		});

		await viewerState.$promise_prevFrame;
		viewerState.$promise_prevFrame = null;
		viewerState.$animationFrameId = 0;

		return true;
	}
	else {
		return false;
	}
}

async function delayFrame() {
	// await viewerState.$promise_prevFrame;
	// viewerState.$promise_prevFrame = null;

	// viewerState.$promise_prevFrame = new Promise(function (resolve, reject) {
	const task = new Promise(function (resolve, reject) {
		try {
			requestAnimationFrame(function () {
				resolve();
			});
		}
		catch (ex) {
			reject(ex);
		}
	});
	await task;

	// await viewerState.$promise_prevFrame;
	// viewerState.$promise_prevFrame = null;
}

/**
 * @param {AnalysisOptions} options
 */
async function render_tetrad(options) {
	if (!seq_list[0].length) {
		return;
	}

	const el_appendRow = document.getElementById("append-row");
	el_appendRow.innerHTML = "";//clear

	await _render_tetrad(options);

	let v_scale = getPixelPerBP() / g_maxPixelPerBP;

	v_scale = seq_list[0].length / (g_chrCoord.bp_end - g_chrCoord.bp_start + 1);

	let list = [
		"scale: " + (v_scale).toFixed(2).toString(),
	];

	document.getElementById("scale").innerText = list.join(",");
}

/**
 * @param {Event} evt
 */
function set_ref_pos_from_ma(evt) {
	const ref1_start = pos_ref1_uint32array[g_chrCoord.bp_start - 1];
	const ref1_end = pos_ref1_uint32array[g_chrCoord.bp_end - 1];
	el_input_ref1_start.value = ref1_start;
	el_input_ref1_end.value = ref1_end;

	const ref2_start = pos_ref2_uint32array[g_chrCoord.bp_start - 1];
	const ref2_end = pos_ref2_uint32array[g_chrCoord.bp_end - 1];
	el_input_ref2_start.value = ref2_start;
	el_input_ref2_end.value = ref2_end;
}

//ctx.canvas.width / (g_chrCoord.bp_end - g_chrCoord.bp_start + 1)
/**
 * @param {AnalysisOptions} options
 */
async function _render_tetrad(options) {
	const max_view_width = viewerState.max_view_width;
	if (max_view_width <= 0) {
		return;
	}

	const view_length = g_chrCoord.bp_end - g_chrCoord.bp_start + 1;
	const scale = 1 / view_length;
	const bp_size = max_view_width * scale;
	const bp_per_px = 1 / bp_size;

	const bp_per_px_max = seq_list[0].length / max_view_width;
	const bp_per_px_min1 = Math.max(1, Math.max(1, bp_per_px) / bp_per_px_max * viewerState.rip_display_weight);

	update_UI_viewLength(view_length);

	main_ctx.lineJoin = "round";

	main_ctx.setTransform(1, 0, 0, 1, 0, 0);
	//ctx.translate(0.5, 0.5);
	//ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
	main_ctx.fillStyle = "#FFFFFF"
	main_ctx.fillRect(0, 0, main_ctx.canvas.width, main_ctx.canvas.height);

	let allow_draw_small_object = bp_size / 4 * 3 >= 7;

	if (allow_draw_small_object) {
		onChangeColorSet(null, "view", false);
	}
	else {
		onChangeColorSet(null, "print", false);
	}

	let seg_row_height = viewerState.seg_row_height;//32
	let seg_row_separate = viewerState.seg_row_separate;//15//5
	let gc_max_height = seg_row_height * 2;

	let draw_rdna_arrow = false;
	if (options.rDNA_info && options.rDNA_info.chr == viewerState.nChr) {
		let min_len_rDNA = Math.min(...options.rDNA_info.data.map(d => {
			return Math.min(...d.repeats.map((a, b) => Math.abs(a[1] - a[0])));
		}));
		if ((min_len_rDNA * bp_size) >= 1) {
			draw_rdna_arrow = true;
		}
	}

	// main_ctx.save();
	// const ref_info_list = [
	// 	{
	// 		refName: viewerState.ref1,
	// 		ref_mapto_ma: ref1_pos_uint32array,
	// 	},
	// 	{
	// 		refName: viewerState.ref2,
	// 		ref_mapto_ma: ref2_pos_uint32array,
	// 	}
	// ];
	const ref_pos_map_list = make_ref_pos_map_list();
	/** @type {string[]} */
	const ref_name_list = dataset.genome_info_list.map(a => a.name);

	// fill line in pixel
	main_ctx.lineWidth = 1;
	main_ctx.translate(0.5, 0.5);

	if (options.mode == "tetrad") {
		// draw parental
		for (let seg_id = 0; seg_id < 2; ++seg_id) {//region_rect.length
			// const pos_x = 0;
			// const pos_y = seg_id * (seg_row_height + seg_row_separate);

			if (seg_id < 2) {
				drawSegRow(
					main_ctx, 0, 0,
					seg_id, seg_row_height, seg_row_separate, max_view_width,
					bp_size, bp_per_px_min1,
					allow_draw_small_object,
					draw_rdna_arrow ? options.rDNA_info : null
				);

				//main_ctx.drawImage(canvas, 0, 0);
				{
					main_ctx.save();
					drawSegRow_cover(
						main_ctx, 0, 0,
						seg_id, seg_row_height, seg_row_separate, max_view_width,
						bp_size, bp_per_px_min1,
						allow_draw_small_object,
						draw_rdna_arrow ? options.rDNA_info : null
					);
					main_ctx.restore();
				}
				main_ctx.translate(0, seg_row_height + seg_row_separate);

				// main_ctx.drawImage(canvas, 0, 0);
				// main_ctx.drawImage(canvas, 0, seg_row_height);

				if (viewerState.display_GC_content) {
					// parental
					main_ctx.save();
					main_ctx.translate(0, 1);
					const r = draw_GC_row(
						seg_id, ref_name_list[seg_id], ref_pos_map_list[seg_id],
						region_rect,
						bp_size,
						max_view_width, gc_max_height, seg_row_separate
					);
					main_ctx.restore();
					if (r) {
						main_ctx.translate(0, gc_max_height + seg_row_separate);
					}
					[...document.querySelectorAll(`#data-rows .gc-plot[data-mode="tetrad-mode"]`)].forEach(a => a.style.display = "block");
				}
				else {
					[...document.querySelectorAll(`#data-rows .gc-plot[data-mode="tetrad-mode"]`)].forEach(a => a.style.display = "none");
				}
			}
		}
	}
	else if (options.mode != "single") {
		// 20210616
		for (let seg_id = 0; seg_id < dataset.parental_list.length; ++seg_id) {
			drawProgney(seg_id, seg_id + 1, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options);

			main_ctx.save();
			main_ctx.translate(0, 1);
			const r = draw_GC_row(
				seg_id, ref_name_list[seg_id], ref_pos_map_list[seg_id],
				region_rect,
				bp_size,
				max_view_width, gc_max_height, seg_row_separate
			);
			main_ctx.restore();
			if (r) {
				main_ctx.translate(0, gc_max_height + seg_row_separate);
			}
		}
	}

	draw_all_APD(bp_size, max_view_width, seg_row_height, seg_row_separate);

	main_ctx.textAlign = "left";

	if (viewerState.draw_all_sample_CO) {
		document.getElementById("progney-group").style.display = "none";

		await drawSampleRows(max_view_width, seg_row_height, seg_row_separate, bp_size, view_length);

		drawMarkers(bp_size, max_view_width, seg_row_height, seg_row_separate);

		return;
	}
	else {
		// draw progney
		if (options.mode == "tetrad") {
			if (viewerState.display_progeny) {
				drawProgney(2, region_rect.length, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options);
				document.getElementById("progney-group").style.display = "block";
			}
			// if (viewerState.display_gff && gff) {
			// 	drawTetradMode_annotation(seg_row_height, seg_row_separate, bp_size, max_view_width);
			// }
		}
		else if (options.mode != "single") {
			drawProgney(dataset.parental_list.length, region_rect.length, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options);

			// for (let seg_id = 0; seg_id < dataset.progeny_list.length; ++seg_id) {
			// 	drawProgney(dataset.parental_list.length + seg_id, dataset.parental_list.length + seg_id + 1, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options);
			// 	const r = draw_GC_row(
			// 		seg_id, ref_name_list[seg_id], ref_pos_map_list[seg_id],
			// 		region_rect,
			// 		bp_size,
			// 		max_view_width, gc_max_height, seg_row_separate
			// 	);
			// 	main_ctx.restore();
			// 	if (r) {
			// 		main_ctx.translate(0, gc_max_height + seg_row_separate);
			// 	}
			// }

			if (options.mode != "SNP") {
				/** @see {@link update_seqRow_header} */
				document.getElementById("progney-group").style.display = "block";
			}
		}
		else if (options.mode == "single") {
			// drawProgney(0, region_rect.length, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options);
			drawSegRow_cover(main_ctx, 0, 0, 0, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow ? options.rDNA_info : null);

			//begin GC%
			const seg_id = 0;
			main_ctx.save();
			main_ctx.translate(0, 1);
			const r = draw_GC_row(
				seg_id, ref_name_list[seg_id], ref_pos_map_list[seg_id],
				region_rect,
				bp_size,
				max_view_width, gc_max_height, seg_row_separate
			);
			main_ctx.restore();
			if (r) {
				main_ctx.translate(0, gc_max_height + seg_row_separate);
			}
		}
	}

	await drawAfterProgney(bp_size, max_view_width, seg_row_height, seg_row_separate);

	if (viewerState.display_repeat_segment) {
		draw_repeat_segment(bp_size, seg_row_height, seg_row_separate);
	}

	// pre-defined marker

	drawMarkers(bp_size, max_view_width, seg_row_height, seg_row_separate);

	// methy
	// drawMethyDiff(max_view_width, seg_row_height, seg_row_separate);
	await drawSampleRows(max_view_width, seg_row_height, seg_row_separate, bp_size, view_length);
	
	drawScaleBar(bp_size, view_length, max_view_width);

	// end methy
	await delayFrame();

	// end all
}

function drawScaleBar(bp_size, view_length, max_view_width) {
	if (view_length <= 12000) {
		const ctx = main_ctx;
		ctx.font = "24px Arial";
		ctx.translate(0, 0);
		ctx.strokeStyle = "black";
		ctx.fillStyle = "black";
		for (let x1 = g_chrCoord.bp_start; x1 < g_chrCoord.bp_end; x1 += 1000) {
			const p1 = (x1 - g_chrCoord.bp_start) * bp_size;
			const p2 = (x1 - g_chrCoord.bp_start + 1000) * bp_size;
			if (p2 < max_view_width) {

				ctx.beginPath();
				ctx.moveTo(p1, 24);
				ctx.lineTo(p1, 0);
				ctx.lineTo(p2, 0);
				ctx.stroke();
				ctx.fillText(pos_ref1_uint32array[x1], p1, 24);
			}
		}
		ctx.translate(0, 30);
	}
}

function draw_all_APD(bp_size, max_view_width, seg_row_height, seg_row_separate) {
	if (dataset.parental_list.length > 1) {
		if (viewerState.$display_SNP_density) {
			drawAPDDensity(true, false, "red", bp_size, max_view_width, seg_row_height, seg_row_separate);
			drawAPDDensity(false, true, "blue", bp_size, max_view_width, seg_row_height, seg_row_separate);

			// appendLeftTitle(`${dataset.parental_list[0]} ≠ ${dataset.parental_list[1]}`, seg_row_height, seg_row_separate, 1, 1 * 100);
			document.getElementById("parental-SNP").style.display = "block";
			document.getElementById("parental-InDel").style.display = "block";
		}
		else {
			document.getElementById("parental-SNP").style.display = "none";
			document.getElementById("parental-InDel").style.display = "none";
		}
	}
}

/**
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} bp_size
 * @param {number} max_view_width
 */
function drawTetradMode_annotation(seg_row_height, seg_row_separate, bp_size, max_view_width) {
	const k_h = [
		(seg_row_height + seg_row_separate) + ((seg_row_height * 2) + seg_row_separate),
		(seg_row_height + seg_row_separate),
	];
	const y0 = (-k_h[0] * 2) + (-k_h[1] * 4); // + seg_row_height * 0.25;
	let y1 = y0;
	// Object.keys(gff).forEach(function (refId, refIdx) {
	dataset.genomeNameList.forEach(function (refId, refIdx) {

		// gff gene name
		_drawTetradMode_annotation(refId, refIdx, y1, seg_row_height, bp_size, max_view_width);

		y1 = y1 + k_h[refIdx <= 2 ? 0 : 1];
	});
}

/**
 *
 * @param {string} refId
 * @param {number} refIdx
 * @param {number} y1
 * @param {number} seg_row_height
 * @param {number} bp_size
 * @param {number} max_view_width
 * @returns
 */
function _drawTetradMode_annotation(refId, refIdx, y1, seg_row_height, bp_size, max_view_width) {
	if (refId == null) {
		return;
	}
	const sChr = analysis_options.chrInfo_list[refIdx].symbol; //chr.replace(/_unitig_.*consensus/, "");

	if (!gff_data_map[refId]) {
		return;
	}

	/** @type {GFF_ROW[]} */
	const gene_list = Object.entries(gff_data_map[refId]).find(a => sChr == a[0])?.[1]; //.filter(a => a.type == "gene");

	// const gene_list = Object.entries(gff[refId]).find(a => sChr.indexOf(a[0]) >= 0)[1];//.filter(a => a.type == "gene");
	if (!gene_list) {
		return;
	}

	main_ctx.textAlign = "left";
	main_ctx.textBaseline = "middle";
	main_ctx.font = Math.trunc(Math.max(8.5, seg_row_height * 0.5 - 2)) + "px Arial";

	gene_list.forEach(gene => {
		if (gene.type == "region") {
			return;
		}
		//const { $start, $length } = gene;
		const start = ((1 - g_chrCoord.bp_start) + gene.$start) * bp_size;
		const dist = gene.$length * bp_size;
		const end = start + dist;

		if (end >= 0 && start <= max_view_width) {
			let x1 = Math.max(0, start);
			let x2 = Math.min(end, max_view_width);
			let len = x2 - x1;

			const gene_text = [gene.attributes.ID, gene.name ?? gene.attributes.Name].join(", ");
			const left_overflow = x1 != start ? "..." : "";
			const right_overflow = x2 != end ? "..." : "";

			draw_GFF3_gene(gene, x1, y1, len, seg_row_height * 0.5, gene.attributes.ID, gene_text, left_overflow, right_overflow);
		}
	});
}

function draw_repeat_segment(bp_size, seg_row_height, seg_row_separate) {
	main_ctx.fillStyle = "#FF00007F";
	main_ctx.strokeStyle = "#0000007F";

	const c_min = viewerState.repeat_segment_filter_min;
	const c_max = viewerState.repeat_segment_filter_max;

	dataset.parental_list.forEach(refId => {
		repeat_segment[refId].forEach(rows => {
			let filtered_rows = rows;

			if (c_min) {
				if (viewerState.repeat_segment_filter) {
					filtered_rows = filtered_rows.filter(viewerState.repeat_segment_filter);
				}
				else {
					filtered_rows = filtered_rows.filter(row => {
						if (c_min.identity && row.identity < c_min.identity) {
							return;
						}
						if (c_min.alignment && row.alignment < c_min.alignment) {
							return;
						}
						if (c_min.q_len && row.q_len < c_min.q_len) {
							return;
						}
						if (c_min.s_len && row.s_len < c_min.s_len) {
							return;
						}
						// if (c_min.s_len && row.s_len < c_min.s_len &&
						// 	c_min.q_len && row.q_len < c_min.q_len
						// ) {
						// 	return;
						// }
						return true;
					});
				}
			}
			if (c_max) {
				filtered_rows = filtered_rows.filter(row => {
					if (c_max.identity && row.identity > c_max.identity) {
						return;
					}
					if (c_max.alignment && row.alignment > c_max.alignment) {
						return;
					}
					if (c_max.q_len && row.q_len > c_max.q_len) {
						return;
					}
					if (c_max.s_len && row.s_len > c_max.s_len) {
						return;
					}
					// if (c_max.s_len && row.s_len > c_max.s_len &&
					// 	c_max.q_len && row.q_len > c_max.q_len
					// ) {
					// 	return;
					// }
					return true;
				});
			}

			filtered_rows.forEach(segment => {
				if (segment.start <= g_chrCoord.bp_end && segment.end >= g_chrCoord.bp_start) {
					const p1 = ((1 - g_chrCoord.bp_start) + segment.start) * bp_size;
					const p2 = ((1 - g_chrCoord.bp_start) + segment.end) * bp_size;

					main_ctx.beginPath();

					main_ctx.rect(p1, 0, p2 - p1, seg_row_height);
					main_ctx.fill();
					main_ctx.stroke();
				}
			});
			main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
		});
	});
}

function update_UI_viewLength(view_length) {
	try {
		const el_view_length = document.getElementById("el_view_length");
		if (el_view_length) {
			const ref1_start = pos_ref1_uint32array[g_chrCoord.bp_start - 1];
			const ref1_end = pos_ref1_uint32array[g_chrCoord.bp_end - 1];

			const ref2_start = pos_ref2_uint32array[g_chrCoord.bp_start - 1];
			const ref2_end = pos_ref2_uint32array[g_chrCoord.bp_end - 1];

			const ref1_length = (ref1_end - ref1_start + 1) || analysis_options.chrInfo_list[0].length;
			const ref2_length = (ref2_end - ref2_start + 1) || analysis_options.chrInfo_list[1].length;

			el_view_length.innerHTML = [
				`<span title="view_length">${view_length}</span>`,
				`<span title="ref1_length" style="background: ${$color_set_view["dad_bk"]};">${ref1_length}</span>`,
				`<span title="ref2_length" style="background: ${$color_set_view["mom_bk"]};">${ref2_length}</span>`,
			].join("");
		}
	}
	catch (ex) {
		console.error(ex);
	}
}

/**
 * @param {number} specific_seq_len
 * @returns {[start: number, end: number][]|number[][]}
 */
function get_ref1_strain_specific(specific_seq_len) {
	const regexp = new RegExp(`-{${specific_seq_len},}`, "g");
	// sp: strain-specific
	const sp = [...seq_list[1].matchAll(regexp)].map(a => [a.index, a.index + a[0].length + 1]);
	return sp;
}
/**
 * @param {number} specific_seq_len
 * @returns {[start: number, end: number][]|number[][]}
 */
function get_ref2_strain_specific(specific_seq_len) {
	const regexp = new RegExp(`-{${specific_seq_len},}`, "g");
	const sp = [...seq_list[0].matchAll(regexp)].map(a => [a.index, a.index + a[0].length + 1]);
	return sp;
}


/**
 * @param {GTF_ROW} gene
 * @param {number} x1
 * @param {number} y1
 * @param {number} len
 * @param {number} height
 * @param {string} geneId
 * @param {string} gene_text
 * @param {string} left_overflow
 * @param {string} right_overflow
 */
function draw_GFF3_gene(gene, x1, y1, len, height, geneId, gene_text, left_overflow, right_overflow) {
	if (gene.strand > 0) {
		main_ctx.fillStyle = "#007F004F";
	}
	else if (gene.strand < 0) {
		main_ctx.fillStyle = "#00007F4F";
	}

	if (gene.type == "gene") {
		main_ctx.fillRect(x1, y1, len, height);
	}

	main_ctx.fillStyle = "#000000";
	if ((gene.type == "gene") &&
		(gene.attributes.ID.length * 8) <= len) {
		if (gene.strand > 0) {
			main_ctx.textAlign = "left";
			main_ctx.fillText(`${left_overflow} 🢂 ${gene_text} ${right_overflow}`, x1, y1 + height * 0.5, len);
		}
		else if (gene.strand < 0) {
			main_ctx.textAlign = "right";
			main_ctx.fillText(`${left_overflow} ${gene_text} 🢀 ${right_overflow}`, x1 + len, y1 + height * 0.5, len);
		}
	}
	if ((viewerState.display_gff_CDS && gene.type == "CDS") ||
		(viewerState.display_gff_exon && gene.type == "exon") ||
		(viewerState.display_gff_mRNA && gene.type == "mRNA")
	) {
		main_ctx.textAlign = "left";
		main_ctx.fillStyle = "#7F00004F";
		main_ctx.fillRect(x1, y1 + height, len, height);
		if (len > 5) {
			main_ctx.strokeStyle = "#0000007F";
			main_ctx.strokeRect(x1, y1 + height, len, height);
		}
	}
}

/**
 * @param {number} bp_size
 * @param {number} view_length
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 */
async function drawAllAnnotation(bp_size, view_length, max_view_width, seg_row_height, seg_row_separate) {
	const split_strand = true;

	main_ctx.save();
	try {
		const sChr = analysis_options.chrInfo_list[0].symbol; //chr.replace(/_unitig_.*consensus/, "");
		const new_sChr = chr_new_name[sChr];

		if (QM6a_transposome[sChr]) {
			/** @type {GTF_ROW[]} */
			const gene_list = QM6a_transposome[sChr];

			appendLeftTitle("QM6a transposon", seg_row_height, seg_row_separate, 1);
			const options = {
				strain_name: "",
				split_strand: false,
				border_color: "#000000"
			};
			drawAnnotationFrame_1(max_view_width, seg_row_height, gene_list, bp_size, options , draw_GTF_gene_s);

			main_ctx.translate(0, 1 * seg_row_height + seg_row_separate);
		}

		if (gff_data_map && gff_data_map["QM6a_new"]) {
			/** @type {GTF_ROW[]} */
			const new_gene_list = gff_data_map["QM6a_new"][new_sChr];

			/** @type {GTF_ROW[]} */
			const gene_list = gff_data_map["QM6a"][sChr];

			const display_gene_list = viewerState.display_new_gff3 ? new_gene_list : gene_list;
			draw_PCG(display_gene_list);
		}
		else {
			draw_PCG(gff_data_map[dataset.ref][sChr]);
		}
		function draw_PCG(display_gene_list) {
			//  (PCG)
			appendLeftTitle("protein-coding gene", seg_row_height, seg_row_separate, 1);
			const options = {
				strain_name: "",
				split_strand: split_strand,
				border_color: "#000000"
			};
			// drawAnnotationFrame_1(max_view_width, seg_row_height, display_gene_list, bp_size, options, draw_GFF3_gene);
			drawAnnotationFrame_1(max_view_width, seg_row_height, display_gene_list, bp_size, options, draw_GTF_gene_s);
			main_ctx.translate(0, 1 * seg_row_height * 1 + seg_row_separate);
		}

		if (gene_cluster && gene_cluster.length) {
			//  (PCG)
			appendLeftTitle("protein-coding gene", seg_row_height, seg_row_separate, gene_cluster.length);
			gene_cluster.forEach((subcluster, cluster_idx) => {
				const cluster_gene_list = [...subcluster.data.values()].filter(a => a);

				const options = {
					strain_name: subcluster.name,
					split_strand: split_strand,
					border_color: "#000000"
				};
				drawAnnotationFrame_1(max_view_width, seg_row_height, cluster_gene_list, bp_size, options, (gene, x1, y1, len, height, geneId, gene_text, left_overflow, right_overflow, name_1 = null, lncRNAData = null) => {
					draw_GTF_marker(gene, x1, y1, len, height, cluster_idx);
				});

				main_ctx.translate(0, 1 * seg_row_height + seg_row_separate);
			});
		}
	}
	finally {
		main_ctx.restore();
	}
}
/**
 * @param {string} inner_text
 * @param {number} row_height
 * @param {number} seg_row_separate
 * @param {number} [rowspan]
 * @param {number} [p_fontSize]
 */
function appendLeftTitle(inner_text, row_height, seg_row_separate, rowspan = 1, p_fontSize = 100) {
	const el_appendRow = document.getElementById("append-row");
	const el_div = document.createElement("div");
	el_appendRow.append(el_div);
// 	el_div.outerHTML = `
// <div title="methy ratio" style="
// 	${p_fontSize && p_fontSize != 100 ? "font-size: " + p_fontSize + "%;" : ""}
// 	line-height: ${row_height * rowspan}px;
// 	margin-bottom: ${(seg_row_separate) * rowspan}px;
// ">
// 	<span contenteditable="true" spellcheck="false" style="
// 		display: inline-block;
// 		vertical-align: middle;
// 		line-height: normal;
// 	">${inner_text}</span>
// </div>`;

el_div.title = "methy ratio"
el_div.style.fontSize = p_fontSize && p_fontSize != 100 ? p_fontSize + "%" : "";
el_div.style.lineHeight = `${row_height * rowspan}px`;
el_div.style.marginBottom = `${(seg_row_separate) * rowspan}px`;

el_div.innerHTML = `
<span contenteditable="true" spellcheck="false" style="
	display: inline-block;
	vertical-align: middle;
	line-height: normal;
	">${inner_text}</span>
`;

	return el_div;
}

/**
 * @param {GTF_ROW[]} gene_list
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {{ strain_name: string; split_strand: boolean; border_color: string; }} options
 * @param {draw_GTF_gene_1} func_draw_gene
 */
function drawAnnotationFrame_1(max_view_width, seg_row_height, gene_list, bp_size, options, func_draw_gene) {
	const strain_name = options?.strain_name ?? "";
	const split_strand = options?.split_strand ?? false;
	const border_color = options?.border_color ?? "#000000";

	main_ctx.strokeStyle = border_color;
	main_ctx.strokeRect(0, 0 - viewerState.row_padding_top, max_view_width, seg_row_height + viewerState.row_padding_bottom);

	drawAnnotationRow_1(gene_list, bp_size, max_view_width, seg_row_height, split_strand, func_draw_gene);

	main_ctx.textAlign = "left";
	main_ctx.textBaseline = "middle";
	main_ctx.fillStyle = "#000000";
	main_ctx.font = viewerState.global_font_style;
	main_ctx.fillText(strain_name, max_view_width + 4, seg_row_height * 0.5);
}

/**
 * @param {GTF_ROW[]} gene_list
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {boolean} split_strand
 * @param {draw_GTF_gene_1} func_draw_gene
 */
function drawAnnotationRow_1(gene_list, bp_size, max_view_width, seg_row_height, split_strand, func_draw_gene) {
	main_ctx.save();

	main_ctx.textAlign = "left";
	main_ctx.textBaseline = "middle";

	// gff gene name
	main_ctx.font = Math.trunc(Math.max(8.5, seg_row_height * 0.5 - 2)) + "px Arial";

	const list = gene_list.map(gene => {
		// if (gene.type == "transcript") {
			const $start = ref1_pos_uint32array[gene.start - 1] - 1;
			const $end = ref1_pos_uint32array[gene.end - 1] - 1;
			const $length = $end - $start + 1;

			const start = ((1 - g_chrCoord.bp_start) + $start) * bp_size;
			const length = $length * bp_size;
			const end = start + length;

			if (end >= 0 && start <= max_view_width) {
				return {
					gene, start, end, length,
				}
			}
		// }
	}).filter(a => a);

	list.forEach(({ gene, start, end, length }) => {
		let x1 = Math.max(0, start);
		let x2 = Math.min(end, max_view_width);
		let len = Math.max(1, x2 - x1);
		let y1 = 0;
		let height;

		if (split_strand) {
			if (gene.strand > 0) {
				y1 = 0;
			}
			else if (gene.strand < 0) {
				y1 = seg_row_height * 0.5;
			}
			else {
				y1 = seg_row_height * 0.33;
			}
			height = seg_row_height * 0.5;
		}
		else {
			height = seg_row_height;
		}

		const gene_text = gene?.query ?? gene?.attributes?.ID;
		const left_overflow = x1 != start ? "..." : "";
		const right_overflow = x2 != end ? "..." : "";

		func_draw_gene(gene, x1, y1, len, height, gene_text, gene_text, left_overflow, right_overflow);

		// if (
		// 	(viewerState.display_gff_CDS && gene.type == "transcript") ||
		// 	(viewerState.display_gff_exon && gene.type == "exon")
		// ) {
		// 	main_ctx.fillStyle = "#7F00004F";
		// 	main_ctx.fillRect(x1, y1 + seg_row_height * 0.5, len, seg_row_height * 0.5);
		// 	if (len > 5) {
		// 		main_ctx.strokeStyle = "#0000007F";
		// 		main_ctx.strokeRect(x1, y1 + seg_row_height * 0.5, len, seg_row_height * 0.5);
		// 	}
		// }
	});

	main_ctx.restore();
}

/**
 * @param {GTF_ROW} gene
 * @param {number} x1
 * @param {number} y1
 * @param {number} len
 * @param {number} height
 * @param {number} cluster_idx
 */
function draw_GTF_marker(gene, x1, y1, len, height, cluster_idx) {
	const alpha = cluster_idx / gene_cluster.length;
	// const ga = alpha ** 2.2;
	// main_ctx.fillStyle = `hsl(${Math.trunc(360 * alpha)}deg,100%,50%)`;
	// main_ctx.fillStyle = `hsl(${Math.trunc(30 + 360 * alpha)}deg,${50 + (1 - alpha) * 50}%,${40 + (1 - alpha) * 20}%)`;

	// main_ctx.fillStyle = `hsl(${Math.trunc(360 * (1 - alpha))}deg, ${80 + 20 * (1 - ga)}%, ${45 + 10 * (1 - ga)}%)`;

	// if (gene.overlap_lncRNA == "veg") {
	// 	main_ctx.strokeStyle = "#F00";
	// 	main_ctx.fillStyle = "#F008";
	// }
	// else if (gene.overlap_lncRNA == "cellulose") {
	// 	main_ctx.strokeStyle = "#080";
	// 	main_ctx.fillStyle = "#0808";
	// }
	// else if (gene.overlap_lncRNA == "glucose") {
	// 	main_ctx.strokeStyle = "#00F";
	// 	main_ctx.fillStyle = "#00F8";
	// }
	// else {
	// 	main_ctx.strokeStyle = "#5558";
	// 	main_ctx.fillStyle = "#5558";
	// }

	if (gene.overlap_lncRNA) {// 0xFFF => white
		// const r = (gene.overlap_lncRNA & 0xF00) >> 8;
		// const g = (gene.overlap_lncRNA & 0x0F0) >> 4;
		// const b = (gene.overlap_lncRNA & 0x00F) >> 0;
		// if (len < 10) {// blend_lncRNA_color
		// 	main_ctx.strokeStyle = "#" + [r, g, b].map(a => a.toString(16)).join("");
		// }
		// else {
		// 	main_ctx.strokeStyle = "#" + [r >> 1, g >> 1, b >> 1].map(a => Math.trunc(a).toString(16)).join("") + "8";
		// }
		// main_ctx.fillStyle = "#" + (gene.overlap_lncRNA | 0).toString(16).padStart(3, "0") + "8";

		main_ctx.strokeStyle = blend_lncRNA_color[gene.overlap_lncRNA].stroke;
		main_ctx.fillStyle = blend_lncRNA_color[gene.overlap_lncRNA].fill;
	}
	else {
		main_ctx.strokeStyle = "#5558";
		main_ctx.fillStyle = "#5558";
	}

	main_ctx.strokeRect(x1, y1, len, height);

	main_ctx.globalAlpha = 0.5;
	main_ctx.fillRect(x1, y1, len, height);

	main_ctx.globalAlpha = 1;
}

/**
 * @param {GTF_ROW} gene
 * @param {number} x1
 * @param {number} y1
 * @param {number} len
 * @param {number} height
 * @param {string} geneId
 * @param {string} gene_text
 * @param {string} left_overflow
 * @param {string} right_overflow
 * @param {number} [data_idx]
 * @param {string} [name_1] cmp RNA_1 RNA_2
 * @param {LncRNAData} [lncRNAData]
 */
function draw_GTF_gene_s(gene, x1, y1, len, height, geneId, gene_text, left_overflow, right_overflow, data_idx = null, name_1 = null, lncRNAData = null) {
	let tpm_text = null;
	let tpm_value = null;
	let has_diff = false;
	let simple = false;

	if (gene.type == "transcript") {
		// main_ctx.strokeStyle = "#00F";
		// main_ctx.fillStyle = "#00F8";
		// if (len < 10) {// blend_lncRNA_color
		// 	main_ctx.strokeStyle = "#" + [0xF00, 0x0F0, 0x00F][data_idx].toString(16).padStart(3, "0");
		// }
		// else {
		// 	main_ctx.strokeStyle = "#" + [0xF00 >> 1, 0x0F0 >> 1, 0x00F >> 1][data_idx].toString(16).padStart(3, "0") + "8";
		// }
		// main_ctx.fillStyle = "#" + [0xF00, 0x0F0, 0x00F][data_idx].toString(16).padStart(3, "0") + "8";

		main_ctx.fillStyle = blend_lncRNA_color[0b1 << data_idx]?.fill ?? "#66666666";
		main_ctx.strokeStyle = blend_lncRNA_color[0b1 << data_idx]?.stroke ?? "#66666677";

		if (lncRNAData) {
			/** @type {string} */
			const geneId = gene?.query;

			tpm_value = lncRNAData.tx_tpm_data[geneId];
			tpm_text = `<${tpm_value}>`;

			// let alpha = Math.min(100, lncRNAData.tx_tpm_data[geneId]) / 100;

			// main_ctx.strokeStyle = "#00F";
			// main_ctx.fillStyle = `#${Math.trunc((1 - alpha) * 15).toString(16)}${Math.trunc(alpha * 15).toString(16)}08`;
		}
	}
	else if (gene.type == "gene") {
		// main_ctx.strokeStyle = "#F00";
		// main_ctx.fillStyle = "#F008";
		main_ctx.strokeStyle = "#777";
		main_ctx.fillStyle = "#7778";
	}
	else if (gene.type.indexOf("transposon") >= 0) {
		main_ctx.strokeStyle = "#F008";
		main_ctx.fillStyle = "#0004";
	}
	else if (gene.type == "mRNA" && viewerState.display_gff_mRNA) {
		main_ctx.strokeStyle = "#F004";
		main_ctx.fillStyle = "#0004";
		simple = true;
	}
	else if (gene.type == "CDS" && viewerState.display_gff_CDS) {
		main_ctx.strokeStyle = "#0F04";
		main_ctx.fillStyle = "#0004";
		simple = true;
	}
	else if (gene.type == "exon" && viewerState.display_gff_exon) {
		main_ctx.strokeStyle = "#00F4";
		main_ctx.fillStyle = "#0004";
		simple = true;
	}
	else {
		return;
	}

	if (has_diff) {
		main_ctx.strokeRect(x1 - 2.5, y1 - 2.5, len + 5, height + 5);
	}

	main_ctx.strokeRect(x1, y1, len, height);
	main_ctx.fillRect(x1, y1, len, height);

	main_ctx.lineJoin = "round";
	main_ctx.lineWidth = 1;
	// main_ctx.shadowColor = "transparent";
	// main_ctx.shadowBlur = 0;

	// if (has_diff) {
	// 	main_ctx.strokeStyle = "black";
	// 	main_ctx.strokeRect(x1 - 5, y1 - 5, len + 10, height + 10);
	// }

	if (!simple && (gene_text.length * 8) <= len) {
		if (gene.strand > 0) {
			main_ctx.textAlign = "left";
			main_ctx.fillStyle = "#000000";
			main_ctx.fillText(`${left_overflow} 🢂 ${[gene_text, tpm_text].join(" ")} ${right_overflow}`, x1, y1 + height * 0.5, len);
		}
		else if (gene.strand < 0) {
			main_ctx.textAlign = "right";
			main_ctx.fillStyle = "#000000";
			main_ctx.fillText(`${left_overflow} ${[gene_text, tpm_text].join(" ")} 🢀 ${right_overflow}`, x1 + len, y1 + height * 0.5, len);
		}
		else {
			main_ctx.textAlign = "center";
			main_ctx.fillStyle = "#000000";
			main_ctx.fillText(`${left_overflow} * ${[gene_text, tpm_text].join(" ")} * ${right_overflow}`, x1 + len * 0.5, y1 + height * 0.5, len);
		}
	}
}

/**
 * @param {GTF_ROW} gene
 * @param {number} x1
 * @param {number} y1
 * @param {number} len
 * @param {number} height
 * @param {string} geneId
 * @param {string} gene_text
 * @param {string} left_overflow
 * @param {string} right_overflow
 */
function draw_GTF_gene_1(gene, x1, y1, len, height, geneId, gene_text, left_overflow, right_overflow) {
	if (gene.type == "transcript") {
		if (gene.strand > 0) {
			main_ctx.fillStyle = "#007F004F";
		}
		else if (gene.strand < 0) {
			main_ctx.fillStyle = "#00007F4F";
		}
		else {
			main_ctx.fillStyle = "#CC66664F";
		}

		main_ctx.fillRect(x1, y1, len, height);

		if ((gene_text.length * 8) <= len) {
			if (gene.strand > 0) {
				main_ctx.textAlign = "left";
				main_ctx.fillStyle = "#000000";
				main_ctx.fillText(`${left_overflow} 🢂 ${gene_text} ${right_overflow}`, x1, y1 + height * 0.5, len);
			}
			else if (gene.strand < 0) {
				main_ctx.textAlign = "right";
				main_ctx.fillStyle = "#000000";
				main_ctx.fillText(`${left_overflow} ${gene_text} 🢀 ${right_overflow}`, x1 + len, y1 + height * 0.5, len);
			}
			else {
				main_ctx.textAlign = "center";
				main_ctx.fillStyle = "#000000";
				main_ctx.fillText(`${left_overflow} * ${gene_text} * ${right_overflow}`, x1 + len * 0.5, y1 + height * 0.5, len);
			}
		}
	}
	// else if (gene.type.endsWith("transposon")) {
	// 	main_ctx.fillStyle = "#FF1111DD";
	// 	main_ctx.fillRect(x1, y1, len, height);
	// 	main_ctx.fillStyle = "#000000";
	// 	if ((gene_text.length * 8) <= len) {
	// 		main_ctx.fillText(`${left_overflow} * ${gene_text} ${right_overflow}`, x1, y1 + height * 0.5, len);
	// 	}
	// }
	else if (viewerState.display_gff_exon) {
		main_ctx.fillStyle = "#0003";
		main_ctx.fillRect(x1, y1, len, height);
		main_ctx.fillStyle = "#000000";
		// if ((gene_text.length * 8) <= len) {
		// 	main_ctx.fillText(`${left_overflow} * ${gene_text} ${right_overflow}`, x1, y1 + seg_row_height * 0.25, len);
		// }
	}
}

/**
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 */
function drawAPDDensity(show_SNP, show_InDel, color, bp_size, max_view_width, seg_row_height, seg_row_separate) {
	main_ctx.strokeStyle = "balck";
	main_ctx.strokeRect(0, 0, max_view_width, seg_row_height);

	const bp_len = g_chrCoord.bp_end - g_chrCoord.bp_start + 1;

	const slen = Math.trunc(bp_len / max_view_width);

	// viewerState.$SNP_density_isSNP;
	const spores = seq_list.slice(2);

	const b_compare_parental = viewerState.$SNP_density_compare_parental;

	if (bp_len <= max_view_width) {
		for (let pos = g_chrCoord.bp_start; pos < g_chrCoord.bp_end; ++pos) {
			if (b_compare_parental && seq_list[0][pos] == seq_list[1][pos]) {
				continue;
			}

			const is_indel = seq_list[0][pos] == "-" || seq_list[1][pos] == "-";
			if (is_indel) {
				if (show_InDel == false) {
					continue;
				}
			}
			else if (show_SNP == false) {
				continue;
			}

			// const is_snp = seq_list.some(sss => sss[pos] != seq_list[0][pos]);
			const is_snp_ref1 = spores.reduce((n_snp, sss, sss_idx, array) => sss[pos] != array[0][pos] ? n_snp + 1 : n_snp, 0);
			const is_snp_ref2 = spores.reduce((n_snp, sss, sss_idx, array) => sss[pos] != array[1][pos] ? n_snp + 1 : n_snp, 0);
			const is_snp_22 = is_snp_ref1 == 2 && is_snp_ref2 == 2;
			if (b_compare_parental || is_snp_22) {
				const x1 = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
				if (Math.trunc(x1 - 0.5) >= 0 && x1 <= Math.trunc(max_view_width + 0.5)) {
					main_ctx.fillStyle = viewerState?.$SNP_density_color ?? color;
					main_ctx.fillRect(x1, 0, viewerState?.density_w ?? bp_size, seg_row_height);
				}
			}
		}
	}
	else {
		for (let i = g_chrCoord.bp_start; i < g_chrCoord.bp_end; i += slen) {
			let sc = 0;
			for (let pos = i; pos < Math.min(i + slen, g_chrCoord.bp_end); ++pos) {
				if (b_compare_parental && seq_list[0][pos] == seq_list[1][pos]) {
					continue;
				}

				const is_indel = seq_list[0][pos] == "-" || seq_list[1][pos] == "-";
				if (is_indel) {
					if (show_InDel == false) {
						continue;
					}
				}
				else if (show_SNP == false) {
					continue;
				}

				// const is_snp = seq_list.some(sss => sss[pos] != seq_list[0][pos]);
				const is_snp_ref1 = spores.reduce((n_snp, sss, sss_idx, array) => sss[pos] != array[0][pos] ? n_snp + 1 : n_snp, 0);
				const is_snp_ref2 = spores.reduce((n_snp, sss, sss_idx, array) => sss[pos] != array[1][pos] ? n_snp + 1 : n_snp, 0);
				const is_snp_22 = is_snp_ref1 == 2 && is_snp_ref2 == 2;
				if (b_compare_parental || is_snp_22) {
					++sc;
				}
			}
			const x1 = ((1 - g_chrCoord.bp_start) + i) * bp_size;
			const x2 = ((1 - g_chrCoord.bp_start) + i + slen) * bp_size;
			if (x1 <= max_view_width && x2 >= 0) {
				let h = viewerState.$SNP_density_func(sc, slen) * seg_row_height;
				main_ctx.fillStyle = viewerState?.$SNP_density_color ?? color;
				main_ctx.fillRect(x1, seg_row_height - h, x2 - x1, h);
			}
		}
	}

	main_ctx.translate(0, seg_row_height);
	main_ctx.translate(0, seg_row_separate);
}

// await drawRNACoverage("QM6a", "ss_rad51", pos_ref1_uint32array, row_height, bp_size, view_length, max_view_width, function stroke_black(ctx) {
// 	ctx.globalAlpha = 1;
// 	ctx.strokeStyle = "#000000";
// });

// QM6a_ssDNA_Rad51.uint32

function ss_ss() {
	// loadCoverageData("QM6a", "rad51-sae2", null)

	for (let i = 0; i < w_rnaCoverage_ui32a_map["QM6a"]["rad51-sae2"].length; ++i) {
		const r1 = w_rnaCoverage_ui32a_map["QM6a"]["ss_rad51_1"][i];
		const r2 = w_rnaCoverage_ui32a_map["QM6a"]["ss_rad51_2"][i];
		const r3 = w_rnaCoverage_ui32a_map["QM6a"]["ss_rad51_3"][i];

		const s1 = w_rnaCoverage_ui32a_map["QM6a"]["ss_sae2_1"][i];
		const s2 = w_rnaCoverage_ui32a_map["QM6a"]["ss_sae2_2"][i];
		const s3 = w_rnaCoverage_ui32a_map["QM6a"]["ss_sae2_3"][i];

		const r = (r1 + r2 + r3) / 3;
		const s = (s1 + s2 + s3) / 3;

		w_rnaCoverage_ui32a_map["QM6a"]["rad51-sae2"][i] = Math.max(0, r - s);

		// Math.abs(Math.log2(r / s)) * 10;
	}

	_merge_rna_reads("QM6a", "rad51-sae2");

	drawFrame()
}

/**
 * @see {@link module_MethylRenderer._load_tsv}
 * @see {@link drawTempData_afterMethyl}
 */
 async function load_ss_peak_QM6a(clear) {
	if (clear) {
		methyl_dataset_list.splice(0);
	}

	const headers_matrix = [
		// ["sChr", "start", "end", "tagName", "value"], // *.bed

		// { valueDesc: "pileup height",   headers: ["sChr", "start", "end", "_3", "_4", "value", "_6",    "_7",    "_8"], },
		// { valueDesc: "pvalue",          headers: ["sChr", "start", "end", "_3", "_4", "_5",    "value", "_7",    "_8"], },
		{ valueDesc: "fold enrichment", headers: ["sChr", "start", "end", "_3", "_4", "_5",    "_6",    "value", "_8"], },
		// { valueDesc: "qvalue",          headers: ["sChr", "start", "end", "_3", "_4", "_5",    "_6",    "_7",    "value"], },

		// xls
		// 0 chromosome name
		// 1 start position of peak
		// 2 end position of peak
		// 3 length of peak region
		// 4 absolute peak summit position
		// 5 pileup height at peak summit
		// 6 -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
		// 7 fold enrichment for this peak summit against random Poisson distribution with local lambda,
		// 8 -log10(qvalue) at peak summit
	];

	const root = `ssDNA`;
	const root_qm6a = `${root}/20210324_ssDNA_mapto_QM6a_CBS1-1`;
	const root_parental = `${root}/20210623_ssDNA_mapto`;

	for (let pair of headers_matrix) {
		// await _load_ss_peak("QM6a", "rad51,sae2", `QM6a rad51_* vs sae2_* ${pair.valueDesc} (effective genome size=34Mb)`, pair.headers, `${root_qm6a}/macs_20210809/QM6a_rad51_sae2_peaks.xls`);
		// await _load_ss_peak("QM6a", "rad51,sae2", `QM6a rad51 vs sae2 ${pair.valueDesc} (effective genome size=34Mb)`, pair.headers, `${root_qm6a}/macs_20210809/QM6a_rad51_merge_sae2_merge_peaks.xls`);
		// await _load_ss_peak("QM6a", "rad51,sae2", `QM6a rad51 vs sae2 ${pair.valueDesc} (effective genome size=30Mb)`, pair.headers, `${root_qm6a}/macs_20210810/QM6a_rad51_merge_sae2_merge_peaks.xls`);

		// await _load_ss_peak("QM6a", "rad51", `QM6a rad51_merge (effective genome size=34Mb)`, pair.headers, `${root_qm6a}/macs_20210810/QM6a_rad51_merge_peaks.xls`);

		for (let i = 1; i <= 3; ++i) {
			await _load_ss_peak("QM6a", "rad51", `QM6a rad51_${i}`, pair.headers, `${root_qm6a}/macs_20210810/QM6a_rad51_${i}_peaks.xls`);
		}
		for (let i = 1; i <= 3; ++i) {
			await _load_ss_peak("QM6a", "sae2", `QM6a sae2_${i}`, pair.headers, `${root_qm6a}/macs_20210810/QM6a_sae2_${i}_peaks.xls`);
		}

		// // // macs3
		// // for (let i = 1; i <= 3; ++i) {
		// // 	await _load_ss_peak("QM6a", "rad51", `QM6a rad51_${i} macs3`, pair.headers, `${root_qm6a}/macs_20210811/QM6a_rad51_${i}_peaks.xls`);
		// // }
		// // for (let i = 1; i <= 3; ++i) {
		// // 	await _load_ss_peak("QM6a", "sae2", `QM6a sae2_${i} macs3`, pair.headers, `${root_qm6a}/macs_20210811/QM6a_sae2_${i}_peaks.xls`);
		// // }

		// // // macs3 -g 34e+6
		// // for (let i = 1; i <= 3; ++i) {
		// // 	await _load_ss_peak("QM6a", "rad51", `QM6a rad51_${i} (effective genome size=34Mb)`, pair.headers, `${root_qm6a}/macs_20210811_g34M/QM6a_rad51_${i}_peaks.xls`);
		// // }
		// // for (let i = 1; i <= 3; ++i) {
		// // 	await _load_ss_peak("QM6a", "sae2", `QM6a sae2_${i} (effective genome size=34Mb)`, pair.headers, `${root_qm6a}/macs_20210811_g34M/QM6a_sae2_${i}_peaks.xls`);
		// // }
	}

	// const bed_headers = ["sChr", "start", "end", "tagName", "value"];
	// for (let i = 1; i <= 3; ++i) {
	// 	for (let j = 1; j <= 3; ++j) {
	// 		await _load_ss_peak("QM6a", `rad51_${i},sae2_${j}`, `QM6a rad51_${i} vs sae2_${j} summits`, bed_headers, `${root_qm6a}/macs_20210809/QM6a_rad51_${i}_sae2_${j}_summits.bed`);
	// 	}
	// }

	await drawFrame();
	viewerState.resizeCanvas();
	await drawFrame();
}

/**
 * @param {string} seq
 */
function calc_gc_content(seq) {
	const _seq = seq.replace(/-/g, "");
	if (_seq.length <= 0) {
		debugger
	}

	const cnt = {
		A: 0,
		T: 0,
		C: 0,
		G: 0,
	};
	for (let a of _seq) {
		cnt[a.toUpperCase()] += 1;
	}

	return (cnt.C + cnt.G) / _seq.length;
}

/**
 * @see {@link genome_ssDNA_peak_overlap_gene_AT_island}
 * @param {number} parent_idx
 * @param {string} target_name
 */
function ssDNA_peak_GC_content(parent_idx, target_name = "rad51_m rad51_1 fold enrichment") {
	const chr_idx = viewerState.nChr - 1;
	const sample = methyl_dataset_list.find(v => v.name == target_name);

	const ref_pos_map_list = make_ref_pos_map_list();
	const ref_pos_map = ref_pos_map_list[dataset.genome_info_list.findIndex(a => a.name == sample.ref)];

	const parental_ref_pos = ref_pos_map_list[parent_idx];
	const parental_pos_ref = make_pos_ref_map_list()[parent_idx];

	sample.data[chr_idx].forEach(point => {
		const y1 = ref_pos_map[point.start];
		const y2 = ref_pos_map[point.end];
		const _y_mid = ref_pos_map[Math.round((point.start + point.end) / 2)];
		const parental_1_mid = parental_pos_ref[_y_mid];

		point.gc = calc_gc_content(seq_list[parent_idx].slice(y1 - 1, y2));
		point.gc_10 = calc_gc_content(seq_list[parent_idx].slice(parental_ref_pos[parental_1_mid - 10], parental_ref_pos[parental_1_mid + 10]));
		point.gc_20 = calc_gc_content(seq_list[parent_idx].slice(parental_ref_pos[parental_1_mid - 20], parental_ref_pos[parental_1_mid + 20]));
		point.gc_30 = calc_gc_content(seq_list[parent_idx].slice(parental_ref_pos[parental_1_mid - 30], parental_ref_pos[parental_1_mid + 30]));
		point.gc_40 = calc_gc_content(seq_list[parent_idx].slice(parental_ref_pos[parental_1_mid - 40], parental_ref_pos[parental_1_mid + 40]));
		point.gc_50 = calc_gc_content(seq_list[parent_idx].slice(parental_ref_pos[parental_1_mid - 50], parental_ref_pos[parental_1_mid + 50]));
	});
	const total = {
		gc: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc)).reduce((acc, v, i, arr) => acc + v.gc / arr.length, 0),
		gc_10: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc_10)).reduce((acc, v, i, arr) => acc + v.gc_10 / arr.length, 0),
		gc_20: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc_20)).reduce((acc, v, i, arr) => acc + v.gc_20 / arr.length, 0),
		gc_30: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc_30)).reduce((acc, v, i, arr) => acc + v.gc_30 / arr.length, 0),
		gc_40: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc_40)).reduce((acc, v, i, arr) => acc + v.gc_40 / arr.length, 0),
		gc_50: sample.data[chr_idx].filter(v => !Number.isNaN(v.gc_50)).reduce((acc, v, i, arr) => acc + v.gc_50 / arr.length, 0),
	};
	return total;
}

const gff_type_priority = {
	"gene": 1,
	"mRNA": 1e+1,
	"CDS": 1e+2,
	"exon": 1e+3,
	"five_prime_UTR": 1e+4,
	"three_prime_UTR": 1e+4,
};

function _ssDNA_peak_overlap_gene_AT_island(target_name = "rad51_m rad51_1 fold enrichment") {
	const sample = methyl_dataset_list.find(v => v.name == target_name);
	const aaa = __ssDNA_peak_overlap_gene_AT_island(["QM6a", "CBS1-1"], sample, null, null, null);
	return {
		exon: aaa.length,
		CDS: aaa.length,
		mRNA: aaa.length,
		gene: aaa.length,
		five_prime_UTR: aaa.length,
		three_prime_UTR: aaa.length,
		AT_island: aaa.length,
		other: aaa.length,
	}
}
/**
 * @param {string[]} target_ref_list
 * @param {module_Methyl_sampleData} sample
 * @param {(point: Partial<module_Methyl_ratioData>) => void} cbfunc_init
 * @param {(point: Partial<module_Methyl_ratioData>, parental_ref: string, gene_type: string, gene_idx: number) => void} cbfunc_gene
 * @param {(point: Partial<module_Methyl_ratioData>, parental_ref: string) => void} cbfunc_island
 */
function __ssDNA_peak_overlap_gene_AT_island(target_ref_list, sample, cbfunc_init, cbfunc_gene, cbfunc_island) {
	return __chr_ssDNA_peak_overlap_gene_AT_island(viewerState.nChr, target_ref_list, sample, cbfunc_init, cbfunc_gene, cbfunc_island);
}
/**
 * @param {number} nChr
 * @param {string[]} target_ref_list
 * @param {module_Methyl_sampleData} sample
 * @param {(point: Partial<module_Methyl_ratioData>) => void} cbfunc_init
 * @param {(point: Partial<module_Methyl_ratioData>, parental_ref: string, gene_type: string, gene_idx: number, inside: boolean) => void} cbfunc_gene
 * @param {(point: Partial<module_Methyl_ratioData>, parental_ref: string, inside: boolean) => void} cbfunc_island
 */
function __chr_ssDNA_peak_overlap_gene_AT_island(nChr, target_ref_list, sample, cbfunc_init, cbfunc_gene, cbfunc_island) {
	const chr_idx = nChr - 1;

	const ref_idx = dataset.genome_info_list.findIndex(a => a.name == sample.ref);
	const ref_pos_map_list = _make_chr_ref_pos_map_list(nChr);
	const ref_pos_map = ref_pos_map_list[ref_idx];

	sample.data[chr_idx].forEach(point => {
		delete point.gene_type;
		delete point.gene;
		delete point.AT_island;
		if (cbfunc_init) {
			cbfunc_init(point);
		}
	});

	target_ref_list.forEach(parental_ref => {
		const parental_idx = dataset.genome_info_list.findIndex(a => a.name == parental_ref);
		const parental_chr = dataset.genome_info_list[parental_idx].chr_list[chr_idx].symbol;
		const parental_to_pos = ref_pos_map_list[parental_idx];

		// gff_data_map[parental_ref][parental_chr].forEach((gene, gene_idx) => {
		
		// const gene_list = gff_data_queryMap[parental_ref][parental_chr][Math.trunc(y1 / 100)];
		const gene_list = gff_data_map[parental_ref][parental_chr];

		if (parental_idx == null || parental_chr == null || gene_list == null) {
			console.error({
				nChr,
				ref_idx,
				"sample.ref": sample.ref,
				
				parental_idx,
				parental_ref,
				chr_idx,
				parental_chr,
			});
		}

		sample.data[chr_idx].forEach(point => {
			const [
				y1, y2
			] = (function () {
				if (sample.ref == parental_ref) {
					return [point.start, point.end];
				}
				else {
					return [ref_pos_map[point.start], ref_pos_map[point.end]];
				}
			})();
			// const [y1, y2] = [ref_pos_map[point.start], ref_pos_map[point.end]];

			// gene_list.forEach((gene, gene_idx) => {
			// for (let gene_idx = 0; gene_idx < gene_list.length; ++gene_idx) {
			for (let gene_idx = 0; gene_idx < gene_list.length; ++gene_idx) {
				const gene = gene_list[gene_idx];

				const [
					x1, x2
				] = (function () {
					if (sample.ref == parental_ref) {
						return [gene.start, gene.end];
					}
					else {
						return [
							parental_to_pos[gene.start],
							parental_to_pos[gene.end],
						];
					}
				})();

				if (x1 <= y2 && x2 >= y1) {
					// point.gene ??= [];
					// point.gene.push({
					// 	gene_idx,
					// 	type: gene.type,
					// });

					if (point.gene_type == null ||
						gff_type_priority[gene.type] >= gff_type_priority[point.gene_type]
					) {
						point.gene_type = gene.type;
						point.gene = gene_idx;
						if (point.gene != gene_idx) {
							console.log(point.gene, gene_idx);
						}
					}

					if (cbfunc_gene != null) {
						const inside = x1 <= y1 && y2 <= x2;
						cbfunc_gene(point, parental_ref, gene.type, gene_idx, inside);
					}
				}
			}
		});
	});

	target_ref_list.forEach(parental_ref => {
		const at_list = get_tseta_AT_island(nChr, [parental_ref], ref_pos_map_list);

		at_list.forEach(island => {
			const [x1, x2] = (function () {
				if (sample.ref == parental_ref) {
					return [island._raw.start, island._raw.end];
				}
				else {
					return [island.start, island.end];
				}
			})();

			sample.data[chr_idx].filter(point => {
				const [y1, y2] = (function () {
					if (sample.ref == parental_ref) {
						return [point.start, point.end];
					}
					else {
						return [ref_pos_map[point.start], ref_pos_map[point.end]];
					}
				})();

				if (x1 <= y2 && x2 >= y1) {
					if (!point.AT_island) {
						point.AT_island = new Set();
					}
					point.AT_island.add(parental_ref);// island.ref == ref
					
					if (cbfunc_island != null) {
						const inside = x1 <= y1 && y2 <= x2;
						cbfunc_island(point, parental_ref, inside);
					}
				}
			});
		});
	});
	// sample.data[chr_idx].filter(a => a.AT_island).length;

	return {
		exon: sample.data[chr_idx].filter(a => a.gene_type == "exon"),
		CDS: sample.data[chr_idx].filter(a => a.gene_type == "CDS"),
		mRNA: sample.data[chr_idx].filter(a => a.gene_type == "mRNA"),
		gene: sample.data[chr_idx].filter(a => a.gene_type == "gene"),
		five_prime_UTR: sample.data[chr_idx].filter(a => a.gene_type == "five_prime_UTR"),
		three_prime_UTR: sample.data[chr_idx].filter(a => a.gene_type == "three_prime_UTR"),

		// AT_island: sample.data[chr_idx].filter(a => !a.gene_type && a.AT_island),
		AT_island: sample.data[chr_idx].filter(a => a.AT_island),// 20211223
		other: sample.data[chr_idx].filter(a => !a.gene_type && !a.AT_island),
	};
}
function ssDNA_peak_overlap_gene_AT_island() {
	return {
		rad51_1_m: _ssDNA_peak_overlap_gene_AT_island("rad51_m rad51_1 fold enrichment"),
		rad51_2_m: _ssDNA_peak_overlap_gene_AT_island("rad51_m rad51_2 fold enrichment"),
		rad51_3_m: _ssDNA_peak_overlap_gene_AT_island("rad51_m rad51_3 fold enrichment"),
		rad51_1_f: _ssDNA_peak_overlap_gene_AT_island("rad51_f rad51_1 fold enrichment"),
		rad51_2_f: _ssDNA_peak_overlap_gene_AT_island("rad51_f rad51_2 fold enrichment"),
		rad51_3_f: _ssDNA_peak_overlap_gene_AT_island("rad51_f rad51_3 fold enrichment"),

		sae2_1_m: _ssDNA_peak_overlap_gene_AT_island("sae2_m sae2_1 fold enrichment"),
		sae2_2_m: _ssDNA_peak_overlap_gene_AT_island("sae2_m sae2_2 fold enrichment"),
		sae2_3_m: _ssDNA_peak_overlap_gene_AT_island("sae2_m sae2_3 fold enrichment"),
		sae2_1_f: _ssDNA_peak_overlap_gene_AT_island("sae2_f sae2_1 fold enrichment"),
		sae2_2_f: _ssDNA_peak_overlap_gene_AT_island("sae2_f sae2_2 fold enrichment"),
		sae2_3_f: _ssDNA_peak_overlap_gene_AT_island("sae2_f sae2_3 fold enrichment"),

		spo11rad51_1_f: _ssDNA_peak_overlap_gene_AT_island("spo11rad51_f spo11rad51_1 fold enrichment"),
		spo11rad51_2_f: _ssDNA_peak_overlap_gene_AT_island("spo11rad51_f spo11rad51_2 fold enrichment"),
		spo11rad51_3_f: _ssDNA_peak_overlap_gene_AT_island("spo11rad51_f spo11rad51_3 fold enrichment"),

		spo11sae2_1_f: _ssDNA_peak_overlap_gene_AT_island("spo11sae2_f spo11sae2_1 fold enrichment"),
		spo11sae2_2_f: _ssDNA_peak_overlap_gene_AT_island("spo11sae2_f spo11sae2_2 fold enrichment"),
		spo11sae2_3_f: _ssDNA_peak_overlap_gene_AT_island("spo11sae2_f spo11sae2_3 fold enrichment"),
	};
}

/**
 * @see {@link ssDNA_peak_GC_content}
 */
async function genome_ssDNA_peak_overlap_gene_AT_island() {
	const chr_results = [];
	const genome_results = {};

	TSETA_before_load();

	for (let chr_idx = dataset.genome_info_list[0].chr_list.length - 1; chr_idx >= 0; --chr_idx) {
		console.log("genome_ssDNA_peak_overlap_gene_AT_island", chr_idx);

		seq_list = await load_multi_alignment_results(chr_idx);
		window.seq_list = seq_list;
		g_chrCoord.bp_start = 1;
		g_chrCoord.bp_end = seq_list[0].length - 1;

		pos_ref1_uint32array = multialign_to_chrPos_posMap(seq_list[0]);
		ref1_pos_uint32array = chrPos_to_multialign_posMap(seq_list[0]);
		pos_ref2_uint32array = multialign_to_chrPos_posMap(seq_list[1]);
		ref2_pos_uint32array = chrPos_to_multialign_posMap(seq_list[1]);
		make_subject_pos_to_ma(seq_list);
		make_subject_ma_to_pos(seq_list);

		await load_QM6a_CBS11_gff();

		await load_ss_peak(true, false);

		chr_results[chr_idx] = await ssDNA_peak_overlap_gene_AT_island();

		[
			"rad51_1_m",
			"rad51_2_m",
			"rad51_3_m",
			"rad51_1_f",
			"rad51_2_f",
			"rad51_3_f",
			"sae2_1_m",
			"sae2_2_m",
			"sae2_3_m",
			"sae2_1_f",
			"sae2_2_f",
			"sae2_3_f",
			"spo11rad51_1_f",
			"spo11rad51_2_f",
			"spo11rad51_3_f",
			"spo11sae2_1_f",
			"spo11sae2_2_f",
			"spo11sae2_3_f",
		].forEach(name => {
			chr_results[chr_idx][name];

			if (genome_results[name] == null) {
				genome_results[name] = {
					exon: 0,
					CDS: 0,
					mRNA: 0,
					gene: 0,
					five_prime_UTR: 0,
					three_prime_UTR: 0,
					AT_island: 0,
					other: 0,
				};
			}
			genome_results[name].exon = chr_results[chr_idx][name].exon
			genome_results[name].CDS = chr_results[chr_idx][name].CDS
			genome_results[name].mRNA = chr_results[chr_idx][name].mRNA
			genome_results[name].gene = chr_results[chr_idx][name].gene
			genome_results[name].five_prime_UTR = chr_results[chr_idx][name].five_prime_UTR
			genome_results[name].three_prime_UTR = chr_results[chr_idx][name].three_prime_UTR
			genome_results[name].AT_island = chr_results[chr_idx][name].AT_island
			genome_results[name].other = chr_results[chr_idx][name].other
		});

		await delayFrame();
	}

	genome_ssDNA_peak_overlap_gene_AT_island.chr_results = chr_results;
	genome_ssDNA_peak_overlap_gene_AT_island.genome_results = genome_results;

	(function (genome_results) {
		const header = [
			"gene",
			"mRNA",
			"CDS",
			"exon",
			"five_prime_UTR",
			"three_prime_UTR",
			"other",
			"AT_island",
		];
		const lines = Object.entries(genome_results).map(([name, sample]) =>{
			const cols = header.map(k => sample[k]);
			return [name, ...cols].join("\t");
		});

		return [
			[, ...header].join("\t"),
			...lines
		].join("\n");
	})(genome_ssDNA_peak_overlap_gene_AT_island.genome_results);

	return chr_results;
}
genome_ssDNA_peak_overlap_gene_AT_island.chr_results = [];
genome_ssDNA_peak_overlap_gene_AT_island.genome_results = {};

/**
 * output plot
 * @param {string} name
 * @param {string[]} src_list
 */
function _ssDNA_peak_overlap(name, src_list) {
	const merge_src = src_list.map(a => methyl_dataset_list.find(v => v.name == a));

	const sample = new module_Methyl_sampleData(methyl_dataset_list[0]);// clone

	sample.name = name;

	sample.row_height = 4;
	sample.value_normalize = methyl_dataset_list[0].value_normalize * 3;// 8 * 3 = 24
	ssDNA_peak_fold_CoordY(sample, sample.value_normalize);

	sample.rendering_condition = [
		new module_Methyl_sampleData_RenderingCondition({ color: "#F44", condition: v => v.col == 2 }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#4C4", condition: v => v.col == 1 }),
		new module_Methyl_sampleData_RenderingCondition({ color: "#47F", condition: v => v.col == 0 }),
	];

	// // 0xF_FFFFF_FFFFF_FFFFFn
	// /**
	//  * @param {bigint} n
	//  */
	// function decode_col(n) {
	// 	return n & 0xF_00000_00000_00000n >> 60n;
	// }
	// /**
	//  * @param {bigint} n1
	//  * @param {bigint} n2
	//  * @param {bigint} n3
	//  */
	// function encode_col(n1, n2, n3) {
	// 	return (n1 | n2 | n3) << 0xF_00000_00000_00000n;
	// }
	// const data_map = new BigUint64Array()
	// sample.data.forEach((list, chr_idx) => {
	// 	for (let i = 0; i < 3; ++i) {
	// 		methyl_dataset_list[i].data[chr_idx].forEach(point => {
	// 			point.start
	// 			point.end
	// 			point.value
	// 		});
	// 	}
	// });

	function do_chr(chr_idx) {
		const col_map = new Uint8Array(seq_list[0].length);
		// const idx_map = [
		// 	new Uint16Array(seq_list[0].length),
		// 	new Uint16Array(seq_list[0].length),
		// 	new Uint16Array(seq_list[0].length),
		// ];
		const acc_map = src_list.map(a => new Float32Array(seq_list[0].length));

		for (let i = 0; i < src_list.length; ++i) {
			const col_bit = 1 << i;
			merge_src[i].data[chr_idx].forEach((point, point_idx) => {
				const start = point.start - 1;
				const end = point.end;
				// idx_map[i].fill(point_idx, start, end);

				acc_map[i].fill(point.value, start, end);

				// for (let pos = start; pos < end; ++pos) {
				// 	acc_map[i][pos] += point.value;
				// }

				for (let jj = start; jj < end; ++jj) {
					col_map[jj] |= col_bit;
				}
			});
		}

		// acc_map[0];
		for (let i = 1; i < src_list.length; ++i) {
			for (let pos = 0; pos < acc_map[i].length; ++pos) {
				acc_map[i][pos] += acc_map[i - 1][pos];
			}
		}

		// debugger

		/** @type {module_Methyl_ratioData[][]} */
		const out_list = merge_src.map(a => []);

		let start_pos = 0;
		let prev_col = col_map[start_pos];
		col_map.forEach((current_col, current_pos) => {
			if (current_col != prev_col) {
				if (prev_col) {
					const col = prev_col;
					const end_pos = current_pos - 1;

					for (let pos = start_pos; pos < end_pos; ++pos) {
						for (let i = 0; i < src_list.length; ++i) {
							const col_bit = 1 << i;
							if (col & col_bit) {
								const point = new module_Methyl_ratioData();

								// const data_idx = idx_map[i][pos];
								// const orig_data = merge_src[i].data[chr_idx][data_idx];
								point.start = pos + 1; // array-index to pos
								point.end = point.start + 1; // array-index to pos
								point.value = acc_map[i][pos];
								point.col = i;

								out_list[i].push(point);
							}
						}
					}
				}
				prev_col = current_col;
				start_pos = current_pos;
			}
		});

		return out_list.flat(1).sort((a, b) => a.start - b.start);
	}

	sample.data = [];
	sample.data[viewerState.nChr - 1] = do_chr(viewerState.nChr - 1);

	// sample.data.forEach(do_chr);

	return sample;
}
function ssDNA_peak_overlap() {
	const rad51_m = _ssDNA_peak_overlap("rad51_m rad51_+ fold enrichment", [
		"rad51_m rad51_1 fold enrichment",
		"rad51_m rad51_2 fold enrichment",
		"rad51_m rad51_3 fold enrichment",
	]);
	const rad51_f = _ssDNA_peak_overlap("rad51_f rad51_+ fold enrichment", [
		"rad51_f rad51_1 fold enrichment",
		"rad51_f rad51_2 fold enrichment",
		"rad51_f rad51_3 fold enrichment",
	]);
	const sae2_m = _ssDNA_peak_overlap("sae2_m sae2_+ fold enrichment", [
		"sae2_m sae2_1 fold enrichment",
		"sae2_m sae2_2 fold enrichment",
		"sae2_m sae2_3 fold enrichment",
	]);
	const sae2_f = _ssDNA_peak_overlap("sae2_f sae2_+ fold enrichment", [
		"sae2_f sae2_1 fold enrichment",
		"sae2_f sae2_2 fold enrichment",
		"sae2_f sae2_3 fold enrichment",
	]);

	const spo11rad51_f = _ssDNA_peak_overlap("spo11rad51_f spo11rad51_+ fold enrichment", [
		"spo11rad51_f spo11rad51_1 fold enrichment",
		"spo11rad51_f spo11rad51_2 fold enrichment",
		"spo11rad51_f spo11rad51_3 fold enrichment",
	]);
	const spo11sae2_f = _ssDNA_peak_overlap("spo11sae2_f spo11sae2_+ fold enrichment", [
		"spo11sae2_f spo11sae2_1 fold enrichment",
		"spo11sae2_f spo11sae2_2 fold enrichment",
		"spo11sae2_f spo11sae2_3 fold enrichment",
	]);

	methyl_dataset_list.unshift(rad51_m, rad51_f, sae2_m, sae2_f, spo11rad51_f, spo11sae2_f);
}

/**
 * @see {@link module_MethylRenderer._load_tsv}
 * @see {@link drawTempData_afterMethyl}
 * @param {string} ref
 * @param {string} sampleName
 * @param {string} dataName
 * @param {string[]|["sChr", "start", "end", "tagName", "value"]} headers
 * @param {string} url
 */
async function _load_ss_peak(ref, sampleName, dataName, headers, url) {
	const sample = new module_Methyl_sampleData({
		ref: ref,
		sample: sampleName,
		name: dataName,
		// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
		url: url,
		description: "",
		max_display_value: null,
		value_normalize: null,
		density_to_opacity: false,
		// row_height: 3,
		row_height: 1.5,
		fontSize_scale: 0.8,
		border_color: "#000000",
		rendering_condition: [
			new module_Methyl_sampleData_RenderingCondition({ color: "#000000", condition: d => d.value < 2 }),
			new module_Methyl_sampleData_RenderingCondition({ color: "#E00000", condition: d => d.value <= 0 }),
			new module_Methyl_sampleData_RenderingCondition({ color: "#000FE0", condition: d => d.value >= 2 }),
		],
	});
	methyl_dataset_list.push(sample);

	// await g_methylRenderer.load();
	sample.data = await g_methylRenderer._load_tsv(sample, headers);

	const data_max_value = Math.max.apply(this, sample.data.map(a => Math.max.apply(this, a.map(b => b.value))));
	const coord_max_value = 8;

	// @ts-ignore
	sample.$data_max_value = data_max_value;

	sample.value_normalize = coord_max_value;

	// await sample.apply_valueDesc_HTML_to_canvas(`
	// <div style="display: inline-block; line-height: ${viewerState.seg_row_height * sample.row_height}px;">
	// 	<div title="max_value" style="vertical-align: top; line-height: normal;">${(max_value).toFixed(0)}</div>
	// 	<div title="mid_value" style="line-height: normal; display: inline-block; vertical-align: middle;">${(max_value * 0.5).toFixed(0)}</div>
	// 	<div title="min_value" style="vertical-align: bottom; line-height: normal;">${0}</div>
	// </div>
	// `);

	await ssDNA_peak_fold_CoordY(sample, coord_max_value);

	const ref_pos_map_list = make_ref_pos_map_list();
	const ref_to_pos_map = ref_pos_map_list[dataset.genome_info_list.findIndex(a => a.name == sample.ref)];

	const at_list = get_tseta_AT_island(viewerState.nChr, ["QM6a", "CBS1-1"]);
	sample.data[viewerState.nChr - 1].forEach(data => {
		const start = ref_to_pos_map[data.start];
		const end = ref_to_pos_map[data.end];
		if (at_list.find(island => island.start <= end && island.end >= start)) {
			data.value *= -1;
		}
	});

	// w_rnaCoverage_ui32a_map
	// m_rnaCoverage_ui32a_map
	// rnaCoverage_max_map
	// rnaCoverage_norm_map
	// rna_window_size
	// _merge_rna_reads(refId, sampleName);

	return sample;
}

/**
 * @param {module_Methyl_sampleData} sample
 * @param {number} coord_max_value
 * @param {boolean} [has_minus_value]
 * @param {number} [_fractionDigit]
 */
async function ssDNA_peak_fold_CoordY(sample, coord_max_value, has_minus_value = false, _fractionDigit) {
	const fractionDigits = _fractionDigit ?? Math.ceil(sample.value_normalize - Math.trunc(sample.value_normalize));
	let hi, mi, lo;

	if (!Number.isFinite(coord_max_value) || Number.isNaN(coord_max_value)) {
		debugger;
	}

	if (has_minus_value) {
		hi = coord_max_value * 1.0;
		mi = coord_max_value * 0.0;
		lo = coord_max_value * -1.0;
	}
	else {
		hi = coord_max_value * 1.0;
		mi = coord_max_value * 0.5;
		lo = coord_max_value * 0.0;
	}

	await sample.apply_valueDesc_HTML_to_canvas(`
	<table style="height: ${viewerState.seg_row_height * sample.row_height}px; font-size: ${sample.fontSize_scale}em; border-collapse: collapse; border-spacing: 0px;">
		<tr title="max_value" data-value="${hi}" style="${sample.fontSize_scale <= 1 ? "" : "vertical-align: top;"}">
			<td>${hi.toFixed(fractionDigits)}</td>
		</tr>
		<tr title="mid_value" data-value="${mi}">
			<td>${mi.toFixed(fractionDigits)}</td>
		</tr>
		<tr title="min_value" data-value="${lo}" style="${sample.fontSize_scale <= 1 ? "" : "vertical-align: bottom;"}">
			<td>${lo.toFixed(fractionDigits)}</td>
		</tr>
	</table>
	`);
}

async function capture_all_chr_peak() {
	for (let nChr = 1; nChr <= dataset.genome_info_list[0].chr_list.length; ++nChr) {
		await promise_load_task;
		await delayFrame();

		el_input_chr.value = String(nChr);
		el_input_chr.oninput(null);
		await delayFrame();

		await promise_load_task;
		await delayFrame();

		await load_ss_peak(true);

		await captureScreen();
	}
}

/**
 * @see {@link module_MethylRenderer._load_tsv}
 * @see {@link drawTempData_afterMethyl}
 * @param {boolean} clear
 * @param {boolean} immediate_render
 */
async function load_ss_peak(clear = false, immediate_render = true) {
	if (clear) {
		methyl_dataset_list.splice(0);
	}

	const headers_matrix = [
		// ["sChr", "start", "end", "tagName", "value"], // *.bed

		// { valueDesc: "pileup height",   headers: ["sChr", "start", "end", "_3", "_4", "value", "_6",    "_7",    "_8"], },
		// { valueDesc: "pvalue",          headers: ["sChr", "start", "end", "_3", "_4", "_5",    "value", "_7",    "_8"], },
		{ valueDesc: "fold enrichment", headers: ["sChr", "start", "end", "_3", "_4", "_5",    "_6",    "value", "_8"], },
		// { valueDesc: "qvalue",          headers: ["sChr", "start", "end", "_3", "_4", "_5",    "_6",    "_7",    "value"], },

		// xls
		// 0 chromosome name
		// 1 start position of peak
		// 2 end position of peak
		// 3 length of peak region
		// 4 absolute peak summit position
		// 5 pileup height at peak summit
		// 6 -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
		// 7 fold enrichment for this peak summit against random Poisson distribution with local lambda,
		// 8 -log10(qvalue) at peak summit
	];

	const root = `ssDNA`;
	const root_qm6a = `${root}/20210324_ssDNA_mapto_QM6a_CBS1-1`;
	const root_parental = `${root}/20210623_ssDNA_mapto`;

	const rad51_m_list = [
		{ genome: "rad51_m", name: "rad51_1" },
		{ genome: "rad51_m", name: "rad51_2" },
		{ genome: "rad51_m", name: "rad51_3" },
	];
	const rad51_f_list = [
		{ genome: "rad51_f", name: "rad51_1" },
		{ genome: "rad51_f", name: "rad51_2" },
		{ genome: "rad51_f", name: "rad51_3" },
	];
	const sae2_m_list = [
		{ genome: "sae2_m", name: "sae2_1" },
		{ genome: "sae2_m", name: "sae2_2" },
		{ genome: "sae2_m", name: "sae2_3" },
	];
	const sae2_f_list = [
		{ genome: "sae2_f", name: "sae2_1" },
		{ genome: "sae2_f", name: "sae2_2" },
		{ genome: "sae2_f", name: "sae2_3" },
	];
	const spo11rad51_f_list = [
		{ genome: "spo11rad51_f", name: "spo11rad51_1" },
		{ genome: "spo11rad51_f", name: "spo11rad51_2" },
		{ genome: "spo11rad51_f", name: "spo11rad51_3" },
	];
	const spo11sae2_f_list = [
		{ genome: "spo11sae2_f", name: "spo11sae2_1" },
		{ genome: "spo11sae2_f", name: "spo11sae2_2" },
		{ genome: "spo11sae2_f", name: "spo11sae2_3" },
	];

	for (let pair of headers_matrix) {
		each([
			{
				params: "genome size = QM6a",
				path: "macs_20210811 macs3",
				list: [
					...rad51_m_list,
					...sae2_m_list,
					...rad51_f_list,
					...sae2_f_list,
				],
			},
			{
				params: "genome size = 3.4e+7",
				path: "macs_20210816 macs3",
				list: [
					...spo11rad51_f_list,
					...spo11sae2_f_list,
				],
			}
		]);
		/**
		 * @param {{ path: string; list: { genome: string; name: string; }[]; params: string; }[]} pair_list
		 */
		async function each(pair_list) {
			for (let { path, list } of pair_list) {
				for (let { genome, name } of list) {
					await _load_ss_peak(genome, genome, `${genome} ${name} ${pair.valueDesc}`, pair.headers, `${root_parental}/${path}/${genome}_${name}_peaks.xls`);

					// fetchData("");
				}
			}
		}
	}

	{
		const root = `ssDNA`;
		const pair = { valueDesc: "fold enrichment", headers: ["sChr", "start", "end", "_3", "_4", "_5",    "_6",    "value", "_8"], };
		const path = "macs_20210819";
		const genome = "QM6a";
		const name = "spo11rad51_sae2";
		const root_qm6a = `${root}/20210324_ssDNA_mapto_QM6a_CBS1-1`;
		await _load_ss_peak("QM6a", "spo11rad51,sae2", `QM6a spo11rad51_* vs sae2_* ${pair.valueDesc}`, pair.headers, `${root_qm6a}/macs_20210809/QM6a_rad51_sae2_peaks.xls`);
	}

	if (immediate_render) {
		await drawFrame();
		viewerState.resizeCanvas();
		await drawFrame();
	}
}

/**
 * @see {@link load_SS}
 * @see {@link load_SS}
 * @see {@link drawTempData_afterMethyl}
 */
async function load_cmp_ss(load_ss = false, clear_all_meth = true, NP = false) {
	document.getElementById("status").innerText = "begin load_cmp_ss";
	await delayFrame();

	if (clear_all_meth) {
		while (methyl_dataset_list.pop());
	}

	// viewerState.seg_row_separate = 5;
	viewerState.display_progeny = false;
	if (load_ss) {
		await Promise.all([
			load_SS("QM6a"),
			// load_SS("CBS1-1"),
		]);

		await Promise.all([
			"rad51_m", "rad51_f", "sae2_m", "sae2_f",
			"spo11rad51_f",// spo11rad51_f
			"spo11sae2_f",
		].map(a => load_SS(a)));
	}

	async function tsv_to_tab(url) {
		const QM6a_text = await (await fetch(url)).text();
		const QM6a_rows = QM6a_text.trim().split("\n").map(a => a.trim().split("\t"));
		return QM6a_rows;
	}

	// sp: strain-specific
	const qm6a_ss = get_ref1_strain_specific(500);
	// const cbs11_ss = get_ref2_strain_specific();

	if (!20211206) {
		const parental_rad51_vs_sae2_20210713 = await tsv_to_tab("ssDNA/20210713_parental_rad51_vs_sae2.txt");
		const rad51_sae2_parental_fc_2 = cmp_ss_2("rad51 sae2", parental_rad51_vs_sae2_20210713, 2);
		const rad51_sae2_parental_fc_3 = cmp_ss_2("rad51 sae2", parental_rad51_vs_sae2_20210713, 3);
		const rad51_sae2_parental_fc_4 = cmp_ss_2("rad51 sae2", parental_rad51_vs_sae2_20210713, 4);
		// const rad51_sae2_parental_fc_4 = cmp_ss_2("rad51 sae2", parental_rad51_vs_sae2_20210713, 4);
		// const rad51_sae2_parental_fc_8 = cmp_ss_2("rad51 sae2", parental_rad51_vs_sae2_20210713, 8);

		const parental_rad51_vs_sae2_20210802_wnd100 = await tsv_to_tab("ssDNA/20210802_parental_rad51_vs_sae2_wnd100.txt");
		const parental_rad51_vs_sae2_wnd100_fc_2 = cmp_ss_2("rad51 sae2 w:100bp", parental_rad51_vs_sae2_20210802_wnd100, 2);
		const parental_rad51_vs_sae2_wnd100_fc_3 = cmp_ss_2("rad51 sae2 w:100bp", parental_rad51_vs_sae2_20210802_wnd100, 3);
		const parental_rad51_vs_sae2_wnd100_fc_4 = cmp_ss_2("rad51 sae2 w:100bp", parental_rad51_vs_sae2_20210802_wnd100, 4);
	}

	// const m_rad51_vs_f_rad51 = await tsv_to_tab("ssDNA/m_rad51_vs_f_rad51.txt");
	// cmp_ss_2("m-rad51 f-rad51", m_rad51_vs_f_rad51, 2);
	// cmp_ss_2("m-rad51 f-rad51", m_rad51_vs_f_rad51, 4);
	// cmp_ss_2("m-rad51 f-rad51", m_rad51_vs_f_rad51, 8);

	// const m_sae2_vs_f_sae2 = await tsv_to_tab("ssDNA/m_sae2_vs_f_sae2.txt");
	// cmp_ss_2("m-sae2 f-sae2", m_sae2_vs_f_sae2, 2);
	// cmp_ss_2("m-sae2 f-sae2", m_sae2_vs_f_sae2, 4);
	// cmp_ss_2("m-sae2 f-sae2", m_sae2_vs_f_sae2, 8);

	// const QM6a_rows = await tsv_to_tab("ssDNA/QM6a_SS.depth.txt");
	// const q2 = cmp_ss("QM6a", "", QM6a_rows, 2);
	// const q4 = cmp_ss("QM6a", "", QM6a_rows, 4);
	// const q8 = cmp_ss("QM6a", "", QM6a_rows, 8);

	const q2 = cmp_ss("QM6a", "r s", await tsv_to_tab("ssDNA/20210505_cmp_S_R.txt"), 2);
	if (!20211206) {
		extract_by_range("Q sp r s", 2, "QM6a", q2, qm6a_ss);
	}

	// const CBS11_text = await (await fetch("ssDNA/CBS1-1_SS.depth.txt")).text();
	// const CBS11_rows = CBS11_text.trim().split("\n").map(a => a.trim().split("\t"));
	// const c2 = cmp_ss("CBS1-1", CBS11_rows, 2);
	// const c4 = cmp_ss("CBS1-1", CBS11_rows, 4);
	// const c8 = cmp_ss("CBS1-1", CBS11_rows, 8);

	// const QM6a_sr_ss_rows = await tsv_to_tab("ssDNA/QM6a_sr_ss_depth.txt");
	// const qsssr2 = cmp_ss("QM6a", "sr ss", QM6a_sr_ss_rows, 2);
	// const qsssr4 = cmp_ss("QM6a", "sr ss", QM6a_sr_ss_rows, 4);
	// const qsssr8 = cmp_ss("QM6a", "sr ss", QM6a_sr_ss_rows, 8);

	const qsssr2 = cmp_ss("QM6a", "sr ss", await tsv_to_tab("ssDNA/20210505_cmp_SR_SS.txt"), 2);
	if (!20211206) {
		extract_by_range("Q sp sr ss", 2, "QM6a", qsssr2, qm6a_ss);
	}

	function has_number_value(a, b) {
		if (a != 0 && b != 0) {
			return a;
		}
		else {
			return 0;
		}
	}
	// const has_plus_value = (a, b) => {
	// 	if (a > 0 && b > 0) {
	// 		return a;
	// 	}
	// 	else {
	// 		return 0;
	// 	}
	// }

	// const qm6a_cbs11_cmp_2 = cmp(`QM6a ∩ CBS1-1 ss`, 2, "QM6a", q2, c2, has_value);
	// cmp(`QM6a ∩ CBS1-1 ss`, 4, "QM6a", q4, c4, has_value);
	// cmp(`QM6a ∩ CBS1-1 ss`, 8, "QM6a", q8, c8, has_value);

	// q2.hide = true;
	// q4.hide = true;
	// q8.hide = true;

	/**
	 * @param {string} name
	 * @param {number} FC_min
	 * @param {string} ref
	 * @param {module_Methyl_sampleData} aa
	 * @param {number[][]} range_list
	 */
	function extract_by_range(name, FC_min, ref, aa, range_list) {
		const chr_idx = viewerState.nChr - 1;

		/** @type {module_Methyl_sampleData} */
		const bb = new module_Methyl_sampleData();

		bb.data = [];
		bb.data[chr_idx] = range_list.map(a => {
			return {
				start: a[0],
				end: a[1],
				value: 1,
			};
		});

		// /** @type {Partial<module_Methyl_ratioData>[][]} */
		// let data = [];

		// data[chr_idx] = [];

		// aa.data.

		// const mmmmm = construct_ss_FC(name, ref, data, 2);

		// methyl_dataset_list.push(mmmmm);

		return cmp(name, FC_min, ref, aa, bb, has_number_value);
	}

	if (!20211206) {
		cmp(`QM6a ∩`, 2, "QM6a", q2, qsssr2, has_number_value);
		// cmp(`*QM6a ∩`, 2, "QM6a", _q2, _qsssr2, has_number_value);
	}

	if (NP) {
		// const { plus, minus } = await load_ss_l_t("QM6a", "test", "ssDNA/QM6a_ss_long_test_blastn.txt");
		const { plus: res_blastn_plus, minus: res_blastn_minus } = await load_ss_l_t("QM6a", "b res", "ssDNA/blast_res.csv");
		const { plus: res_minimap2_plus, minus: res_minimap2_minus } = await load_ss_l_t("QM6a", "m res", "ssDNA/minimap2_res.csv");

		// cmp(`QM6a ss ±`, 2, "QM6a", plus, minus, (a, b) => {
		// 	return Math.abs(a - b) >= 5 ? a : 0;
		// });

		cmp(`QM6a ilu ∩ np ss +`, 2, "QM6a", res_blastn_plus, q2, has_number_value);
		cmp(`QM6a ilu ∩ np ss -`, 2, "QM6a", res_blastn_minus, q2, has_number_value);

		// cmp(`QM6a *ilu ∩ np ss +`, 2, "QM6a", res_blastn_plus, _q2, has_number_value);
		// cmp(`QM6a *ilu ∩ np ss -`, 2, "QM6a", res_blastn_minus, _q2, has_number_value);
	}

	/**
	 * @param {module_Methyl_sampleData} aa
	 * @param {module_Methyl_sampleData} bb
	 * @param {(a: number, b: number) => number} cbfunc
	 */
	function cmp(name, FC_min, ref, aa, bb, cbfunc) {
		const chr_idx = viewerState.nChr - 1;

		const ccc = (function () {
			const af = new Float32Array(seq_list[0].length);
			const bf = new Float32Array(seq_list[0].length);

			const aaa = aa.data[chr_idx];
			const bbb = bb.data[chr_idx];

			aaa.forEach(v => {
				af.fill(v.value, v.start, v.end);
			});

			bbb.forEach(v => {
				bf.fill(v.value, v.start, v.end);
			});

			return af.map((v, i) => {
				return cbfunc(v, bf[i]);
			});
		})();

		/** @type {module_Methyl_ratioData[][]} */
		const data = [];

		data[chr_idx] = typedArray_to_rangeList(ccc, data[chr_idx]);

		const mmmmm = construct_ss_FC(name, ref, data, FC_min);

		methyl_dataset_list.push(mmmmm);

		return mmmmm;
	}

	await g_methylRenderer.__html_to_canvas(methyl_dataset_list, true);

	document.getElementById("status").innerText = "end load_cmp_ss";
}

/**
 * @template T
 * @param {ArrayLike<number>} typed_array
 * @param {T[]} output_array
 * @return {T[]}
 */
function typedArray_to_rangeList(typed_array, output_array = []) {
	const c_data = output_array;

	/** @type {Partial<module_Methyl_ratioData>} */
	let current = {
		start: 1,
		end: 0,
		value: 0,
	};
	typed_array.forEach((v, i) => {
		if (v != current.value) {
			if (current.value) {
				current.end = i - 1;
				c_data.push(current);
			}

			current = {
				start: i,
				// end: i,
				value: v,
			};
		}
	});

	return c_data;
}

/**
 * @deprecated not work
 * @template T
 * @param {Partial<{ start: number; end: number; value: number; }>[]} arr
 */
function _merge_nearest_range(arr) {
	const _arr = arr.slice(0);// copy

	for (let i = 1; i < _arr.length; ++i) {
		const d = _arr[i].start - _arr[i - 1].end;
		if (d == 1) {
			_arr[i - 1].end = _arr[i].end;

			_arr.splice(i, 1);// remove current
			i -= 1;
		}
	}

	return _arr;
	
	// a = typedArray_to_rangeList([0, 1, 1, 2, 2, 3, 3, 3, 0, 0, 0, 1, 1, 1, 0], [])
	// for (let i = 1; i < a.length; ++i) {
	// 	const d = a[i].start - a[i - 1].end;
	// 	if (d == 1) {
	// 		a[i - 1].end = a[i].end;
	// 		a.splice(i, 1);
	// 		i -= 1;
	// 	}
	// }
	// a
}

/**
 * @param {""|"QM6a"|"CBS1-1"} ref
 * @param {string} tag
 * @param {string} file_name
 */
async function load_ss_l_t(ref, tag, file_name) {
	const rows = await (async function (file_name) {
		const text = await (await fetch(file_name)).text();
		if (file_name.endsWith(".csv")) {
			return text.trim().split("\n").map(a => a.trim().split(",").map(v => v.trim()));
		}
		else {
			return text.trim().split("\n").map(a => a.trim().split("\t"));
		}
	})(file_name);

	const genome_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
	const genome_info = dataset.genome_info_list[genome_idx];

	const chr_idx_map = Object.fromEntries(genome_info.chr_list.map((v, i) => [v.symbol, i]));

	const pos_map = [
		ref1_pos_uint32array,
		ref2_pos_uint32array,
	][genome_idx];

	const min_gc_content = dataset?.min_gc_content?.[genome_idx] ?? 6;

	// gc_content["QM6a"][viewerState.nChr][0] = {chr: 3, start: 1, end: 500, gc: 9.6};
	const at_list = dataset.gc_content[ref][viewerState.nChr].filter(a => {
		return a.gc <= min_gc_content;
	});

	/** @type {Partial<module_Methyl_ratioData>[][]} */
	const plus_data = genome_info.chr_list.map(_ => []);
	/** @type {Partial<module_Methyl_ratioData>[][]} */
	const minus_data = genome_info.chr_list.map(_ => []);
	rows.forEach(cols => {
		const [sChr, s_p1, s_p2, s_plus, s_minus, s_diff] = (function (cols) {
			if (cols[1].charAt(0) == "[") {
				const [, v1, v2] = cols[1].match(/\[(\d+)-(\d+)\]/);
				return [cols[0], v1, v2, cols[2], cols[3], cols[4]];
			}
			return cols;
		})(cols);
		const chr_idx = chr_idx_map[sChr.match(/([^_]+_[^_]+)_*.*/)[1]];

		const p1 = Number(s_p1);
		const p2 = Number(s_p2);

		const plus = Number(s_plus);
		const minus = Number(s_minus);

		const diff = Math.abs(plus - minus) >= 2 && (Math.max(plus, minus) / Math.min(plus, minus) >= 2);

		let island;
		if (at_list.find(island => island.start <= p2 && island.end >= p1)) {
			island = -1;
		}
		else {
			island = 1;
		}

		if (diff) {
			if (plus > minus) {
				plus_data[chr_idx].push({
					$start: p1,
					$end: p2,
					start: pos_map[p1 - 1],
					end: pos_map[p2 - 1],
					value: plus * island,
					strand: 1,
				});
			}
			else if (minus > plus) {
				minus_data[chr_idx].push({
					$start: p1,
					$end: p2,
					start: pos_map[p1 - 1],
					end: pos_map[p2 - 1],
					value: minus * island,
					strand: -1,
				});
			}
		}
	});

	const results = {
		plus: construct_ss_FC(`${ref} np ss ${tag} +`, ref, plus_data, null),
		minus: construct_ss_FC(`${ref} np ss ${tag} -`, ref, minus_data, null),
	};
	methyl_dataset_list.push(results.plus);
	methyl_dataset_list.push(results.minus);
	results.plus.max_display_value = 1;
	results.minus.max_display_value = 1;

	/** @type {Partial<module_Methyl_ratioData>[][]} */
	const ss_data = [];
	results.plus.data.forEach((plus_list, chr_idx) => {
		ss_data[chr_idx] = [];

		plus_list.forEach(plus => {
			const p11 = plus.start;
			const p12 = plus.end;
			results.minus.data[chr_idx].forEach(minus => {
				const p21 = minus.start;
				const p22 = minus.end;
				if (p11 <= p22 && p12 >= p21) { // DS
				}
				else {// SS
					if ((p21 - p11) > 0 && (p21 - p11) <= 500) { // true SS
						/** @type {Partial<module_Methyl_ratioData>} */
						const ss = {
							start: p11,
							end: p22,
							value: ((plus.value < 0 || minus.value < 0) ? -1 : 1),
							strand: 1,
							// value: plus.value,
						};
						ss_data[chr_idx].push(ss);
					}
					else { // false SS
					}
				}
			});
		});
	});
	results.ss = construct_ss_FC(`${ref} np ss ${tag} → ←`, ref, ss_data, null);
	methyl_dataset_list.push(results.ss);
	results.ss.max_display_value = 1;

	return results;
}

/**
 * use TSETA coordinate
 * @param {string} tag
 * @param {string[][]} rows
 * @param {number} FC_min
 */
 function cmp_ss_2(tag, rows, FC_min = 0, max_display_value) {
	const data = load_ss_FC_2(rows, FC_min);

	const mmmmm = construct_ss_FC(`${tag} FC >= ${FC_min}`, "", data, FC_min);

	methyl_dataset_list.push(mmmmm);

	return mmmmm;
}

/**
 * @param {""|"QM6a"|"CBS1-1"} ref
 * @param {string} tag
 * @param {string[][]} rows
 * @param {number} FC_min
 */
function cmp_ss(ref, tag, rows, FC_min = 0, max_display_value) {
	// const text = await (await fetch("ssDNA/QM6a_SS.depth.txt")).text()
	// const rows = text.trim().split("\n").map(a => a.trim().split("\t"));

	const data = load_ss_FC(ref, rows, FC_min);

	const mmmmm = construct_ss_FC(`${ref} ss ${tag} FC >= ${FC_min}`, ref, data, FC_min);

	methyl_dataset_list.push(mmmmm);

	return mmmmm;

// methyl_dataset_list.pop()
// methyl_dataset_list.pop()
// methyl_dataset_list.pop()
// cmp_ss(rows, 2)
// cmp_ss(rows, 4)
// cmp_ss(rows, 8)
}

/**
 * use TSETA coordinate
 * @param {number} FC_min
 */
function load_ss_FC_2(rows, FC_min) {
	const pos_ref_map_list = [
		pos_ref1_uint32array,
		pos_ref2_uint32array,
	];

	// dataset.ref
	const min_gc_content = dataset?.min_gc_content?.[0] ?? 6;

	// gc_content["QM6a"][viewerState.nChr][0] = {chr: 3, start: 1, end: 500, gc: 9.6};
	const ref1_at_list = gc_content[dataset.parental_list[0]][viewerState.nChr].filter(a => {
		return a.gc <= min_gc_content;
	});
	const ref2_at_list = gc_content[dataset.parental_list[1]][viewerState.nChr].filter(a => {
		return a.gc <= min_gc_content;
	});

	/** @type {Partial<module_Methyl_ratioData>[][]} */
	const data = dataset.genome_info_list[0].chr_list.map(_ => []);
	rows.forEach(cols => {
		const [sChr, range, pValue, s_fc] = cols;
		if (Number(pValue) <= 0.05) {
			const [s_p1, s_p2] = range.split("-");
			const chr_idx = Number(sChr.match(/ch_(\d+)/)[1]) - 1;

			const fc = Number(s_fc);
			const p1 = Number(s_p1);
			const p2 = Number(s_p2);

			const p11 = pos_ref_map_list[0][p1];
			const p12 = pos_ref_map_list[0][p2];
			const p21 = pos_ref_map_list[1][p1];
			const p22 = pos_ref_map_list[1][p1];

			let island;
			if (
				ref1_at_list.find(island => island.start <= p12 && island.end >= p11) ||
				ref2_at_list.find(island => island.start <= p22 && island.end >= p21)
			) {
				island = -1;
			}
			else {
				island = 1;
			}

			if (fc >= FC_min) {
				data[chr_idx].push({
					$start: p1,
					$end: p2,
					start: p1,
					end: p2,
					value: fc * island,
					strand: 1,
				});
			}
		}
	});

	return data;
}

/**
 * @param {""|"QM6a"|"CBS1-1"} ref
 * @param {number} FC_min
 */
function load_ss_FC(ref, rows, FC_min) {
	const genome_idx = dataset.genome_info_list.findIndex(a => a.name == ref);
	const genome_info = dataset.genome_info_list[genome_idx];

	const pos_map = [
		ref1_pos_uint32array,
		ref2_pos_uint32array,
	][genome_idx];

	const chr_idx_map = Object.fromEntries(genome_info.chr_list.map((v, i) => [v.symbol, i]));

	const min_gc_content = dataset?.min_gc_content?.[genome_idx] ?? 6;

	// gc_content["QM6a"][viewerState.nChr][0] = {chr: 3, start: 1, end: 500, gc: 9.6};
	const at_list = dataset.gc_content[ref][viewerState.nChr].filter(a => {
		return a.gc <= min_gc_content;
	});

	/** @type {Partial<module_Methyl_ratioData>[][]} */
	const data = genome_info.chr_list.map(_ => []);
	rows.forEach(cols => {
		const [sChr, range, pValue, s_fc] = cols;

		const fc = Number(s_fc);

		if (Number.isFinite(fc) && !Number.isNaN(fc) &&
			Number(pValue) <= 0.05
		) {
			const [s_p1, s_p2] = range.split("-");
			const chr_idx = chr_idx_map[sChr.match(/([^_]+_[^_]+)_*.*/)[1]];

			const p1 = Number(s_p1);
			const p2 = Number(s_p2);

			let island;
			if (at_list.find(island => island.start <= p2 && island.end >= p1)) {
				island = -1;
			}
			else {
				island = 1;
			}

			if (fc >= FC_min) {
				data[chr_idx].push({
					$start: p1,
					$end: p2,
					start: pos_map[p1 - 1],
					end: pos_map[p2 - 1],
					value: fc * island,
					strand: 1,
				});
			}
		}
	});

	return data;
}

/**
 * @param {string} name
 * @param {string} ref
 * @param {Partial<module_Methyl_ratioData>[][]} data
 * @param {number} FC_min
 */
function construct_ss_FC(name, ref, data, FC_min) {
	let max;
	// const max = Math.max(...data.map(vv => vv.reduce((aa, v) => Math.max(aa, Math.abs(v.value)), 0)));
	Object.keys(data).forEach(si => {
		const i = Number(si);
		const _max = data[i].reduce((aa, v) => Math.max(aa, Math.abs(v.value)), 0);
		// for (let j = 0; j < data[i].length; ++j) {
		// 	const val = data[i][j].value / _max;
		// 	data[i][j].value = val;
		// }
		if ((i + 1) == viewerState.nChr) {
			max = _max;
			console.log(name, max);
		}
	});

	const html_value_desc = FC_min ? `
	<div class="methyl-html">
    <style>
      .methyl-plot {
        display: flex;
        width: 10em;
      }
      .methyl-plot > * {
        flex: 1;
      }
      .methyl-values {
        display: flex;
        flex-direction: column;
      }
      .methyl-values > * {
        flex: 1;
      }
      .methyl-desc {
        margin: auto;
      }
    </style>
    <div class="methyl-plot">
      <div class="methyl-values">
        <div class="methyl-top">
          ${("FC≧" + FC_min) || Math.ceil(max)}
        </div>
        <div class="methyl-mid">
          ${"" || Math.trunc(max / 2)}
        </div>
        <div class="methyl-bottom">
          ${("FC＜" + FC_min) || 0}
        </div>
      </div>
      <div class="methyl-desc">
        FC
      </div>
    </div>
  </div>` : null;

	const mmmmm = new module_Methyl_sampleData({
		ref: ref,
		// sample: "NP43",
		name: name,
		// split_strand: true,
		html_value_desc: html_value_desc,
		url: "",
		region: false,

		display_minus_value: true,
		max_display_value: FC_min,
		value_normalize: max,
		density_to_opacity: true,
		range_repeat_density: 30,
		data: data,
	});
	return mmmmm;
}

// await load_SS("QM6a")
// await load_SS("CBS1-1")

// Q = QM6a
// sp = QM6a strain-specific (length > 1k)
// r = rad51
// s = sae2
// sr = spo11 rad51
// ss = spo11 sae2
// QM6a_sr_ss.xlsx

// TSETA:Tr_ssDNA_parent >> data:20210623_ssDNA_mapto/

/**
 * @see {@link drawTempData_afterMethyl}
 * @param {"QM6a" | "CBS1-1"} refId
 */
async function load_SS(refId) {
	if (refId == "QM6a" || refId == "CBS1-1") {
		const t1 = loadCoverageData(refId, "rad51_1", `ssDNA/${refId}_SS_Rad51_1.uint32`);// ss_rad51_1
		const t2 = loadCoverageData(refId, "rad51_2", `ssDNA/${refId}_SS_Rad51_2.uint32`);// ss_rad51_2
		const t3 = loadCoverageData(refId, "rad51_3", `ssDNA/${refId}_SS_Rad51_3.uint32`);// ss_rad51_3
		const t4 = loadCoverageData(refId, "sae2_1", `ssDNA/${refId}_SS_Sae2_1.uint32`);// ss_sae2_1
		const t5 = loadCoverageData(refId, "sae2_2", `ssDNA/${refId}_SS_Sae2_2.uint32`);// ss_sae2_2
		const t6 = loadCoverageData(refId, "sae2_3", `ssDNA/${refId}_SS_Sae2_3.uint32`);// ss_sae2_3

		await Promise.all([t1, t2, t3, t4, t5, t6]);
	}

	if (refId == "QM6a") {
		const t4 = loadCoverageData(refId, "spo11rad51_1", `ssDNA/${refId}_ss_spo11_rad51_1.uint32`);// ss_spo11_rad51_1
		const t5 = loadCoverageData(refId, "spo11rad51_2", `ssDNA/${refId}_ss_spo11_rad51_2.uint32`);// ss_spo11_rad51_2
		const t6 = loadCoverageData(refId, "spo11rad51_3", `ssDNA/${refId}_ss_spo11_rad51_3.uint32`);// ss_spo11_rad51_3
		const t1 = loadCoverageData(refId, "spo11sae2_1", `ssDNA/${refId}_ss_spo11_sae2_1.uint32`);// ss_spo11_sae2_1
		const t2 = loadCoverageData(refId, "spo11sae2_2", `ssDNA/${refId}_ss_spo11_sae2_2.uint32`);// ss_spo11_sae2_2
		const t3 = loadCoverageData(refId, "spo11sae2_3", `ssDNA/${refId}_ss_spo11_sae2_3.uint32`);// ss_spo11_sae2_3

		await Promise.all([t1, t2, t3, t4, t5, t6]);
	}

	// 20210628
	if (refId.startsWith("rad51")) {
		await Promise.all([
			loadCoverageData(refId, "rad51_1", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_rad51_1.uint32`),
			loadCoverageData(refId, "rad51_2", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_rad51_2.uint32`),
			loadCoverageData(refId, "rad51_3", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_rad51_3.uint32`),
		]);
	}
	if (refId.startsWith("sae2")) {
		await Promise.all([
			loadCoverageData(refId, "sae2_1", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_sae2_1.uint32`),
			loadCoverageData(refId, "sae2_2", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_sae2_2.uint32`),
			loadCoverageData(refId, "sae2_3", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_sae2_3.uint32`),
		]);
	}
	if (refId.startsWith("spo11rad51")) {
		await Promise.all([
			loadCoverageData(refId, "spo11rad51_1", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_1.uint32`),
			loadCoverageData(refId, "spo11rad51_2", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_2.uint32`),
			loadCoverageData(refId, "spo11rad51_3", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_3.uint32`),
		]);
	}
	if (refId.startsWith("spo11sae2")) {
		await Promise.all([
			loadCoverageData(refId, "spo11sae2_1", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_spo11sae2_1.uint32`),
			loadCoverageData(refId, "spo11sae2_2", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_spo11sae2_2.uint32`),
			loadCoverageData(refId, "spo11sae2_3", `ssDNA/ssDNA/20210623_ssDNA_mapto/${refId}_spo11sae2_3.uint32`),
		]);
	}
	if (false) {// ["rad51_m", "rad51_f", "sae2_m", "sae2_f"].includes(refId)
		// "rad51_m", "rad51_f", "sae2_m", "sae2_f", "spo11rad51_m", "spo11rad51_f", "spo11sae2_m", "spo11sae2_f"
		await Promise.all([
			// loadCoverageData(refId, "rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_1.uint32`),
			// loadCoverageData(refId, "rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_2.uint32`),
			// loadCoverageData(refId, "rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_3.uint32`),

			// loadCoverageData(refId, "rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_1.uint32`),
			// loadCoverageData(refId, "rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_2.uint32`),
			// loadCoverageData(refId, "rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_3.uint32`),

			// loadCoverageData(refId, "sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_1.uint32`),
			// loadCoverageData(refId, "sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_2.uint32`),
			// loadCoverageData(refId, "sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_3.uint32`),

			// loadCoverageData(refId, "sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_1.uint32`),
			// loadCoverageData(refId, "sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_2.uint32`),
			// loadCoverageData(refId, "sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_3.uint32`),

			// loadCoverageData(refId, "spo11rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_1.uint32`),
			// loadCoverageData(refId, "spo11rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_2.uint32`),
			// loadCoverageData(refId, "spo11rad51_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_f_spo11rad51_3.uint32`),

			// loadCoverageData(refId, "spo11rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_1.uint32`),
			// loadCoverageData(refId, "spo11rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_2.uint32`),
			// loadCoverageData(refId, "spo11rad51_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11rad51_m_spo11rad51_3.uint32`),

			// loadCoverageData(refId, "spo11sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_1.uint32`),
			// loadCoverageData(refId, "spo11sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_2.uint32`),
			// loadCoverageData(refId, "spo11sae2_f", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_f_spo11sae2_3.uint32`),

			// loadCoverageData(refId, "spo11sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_1.uint32`),
			// loadCoverageData(refId, "spo11sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_2.uint32`),
			// loadCoverageData(refId, "spo11sae2_m", `ssDNA/ssDNA/20210623_ssDNA_mapto/spo11sae2_m_spo11sae2_3.uint32`),
		]);
	}

	return true;
}

async function drawTempData_beforeMethyl(seg_row_height, seg_row_separate, bp_size, view_length, max_view_width) {
}

/**
 * @param {number} nChr
 * @param {string[]} refs
 */
function get_tseta_AT_island(nChr, refs, ref_pos_map_list = make_ref_pos_map_list()) {
	const min_gc_content = dataset?.min_gc_content?.[0] ?? 6;
	const at_list = refs.map((ref, refIdx) => {
		return (dataset._gc_content ?? dataset.gc_content)[ref][nChr].filter(a => {
			return a.gc <= min_gc_content;
		}).map(a => {
			const pos_map = ref_pos_map_list[refIdx];
			return {
				gc: a.gc,
				start: pos_map[a.start],
				end: pos_map[a.end],
				ref: ref,
				_raw: a,
			};
		});
	}).flat(1).sort((a, b) => a.start - b.end);

	return at_list;
}

/**
 * @type {_drawTempData_afterMethyl}
 */
let drawTempData_afterMethyl = null;

/**
 * @see {@link load_SS}
 * @see {@link load_cmp_ss display config} {@link ssDNA_coverage_display_config}
 * @see {@link load_ss_peak}
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} bp_size
 * @param {number} view_length
 * @param {number} max_view_width
 */
async function _drawTempData_afterMethyl(seg_row_height, seg_row_separate, bp_size, view_length, max_view_width) {
	const row_height = seg_row_height * 3;

	const globalAlpha = 1;//Math.min(Math.max(Math.min(0.01, 1 / rna_window_size), max_view_width / (g_chrCoord.bp_end - g_chrCoord.bp_start + 1) * rna_window_size * 10), 1);

	const pos_ref_map_list = make_pos_ref_map_list();

	const ref_pos_map_list = make_ref_pos_map_list();

	const at_list = get_tseta_AT_island(viewerState.nChr, ["QM6a", "CBS1-1"], ref_pos_map_list);

	if (dataset.ref == "QM6a") {
		// await drawRow("QM6a", pos_ref1_uint32array, "rad51-sae2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		// await drawRow("QM6a", pos_ref1_uint32array, "m2_rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("QM6a", "rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_rad51_1
		await drawRow("QM6a", "rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_rad51_2
		await drawRow("QM6a", "rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_rad51_3
		await drawRow("QM6a", "sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_sae2_1
		await drawRow("QM6a", "sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_sae2_2
		await drawRow("QM6a", "sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_sae2_3

		await drawRow("QM6a", "spo11rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_rad51_1
		await drawRow("QM6a", "spo11rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_rad51_2
		await drawRow("QM6a", "spo11rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_rad51_3
		await drawRow("QM6a", "spo11sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_sae2_1
		await drawRow("QM6a", "spo11sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_sae2_2
		await drawRow("QM6a", "spo11sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);// ss_spo11_sae2_3
	}
	// if (dataset.mode == "tetrad") {
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// 	await drawRow("CBS1-1", pos_ref2_uint32array, "ss_sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	// }
	if (dataset.name == "Tr_ssDNA_parent") {
		// dataset.genome_info_list.map(a => a.name)
		// dataset.progeny_list

		await drawRow("rad51_m", "rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("rad51_m", "rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("rad51_m", "rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

		await drawRow("sae2_m", "sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("sae2_m", "sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("sae2_m", "sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

		await drawRow("rad51_f", "rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("rad51_f", "rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("rad51_f", "rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

		await drawRow("sae2_f", "sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("sae2_f", "sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("sae2_f", "sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

		//

		await drawRow("spo11rad51_f", "spo11rad51_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("spo11rad51_f", "spo11rad51_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("spo11rad51_f", "spo11rad51_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

		await drawRow("spo11sae2_f", "spo11sae2_1", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("spo11sae2_f", "spo11sae2_2", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
		await drawRow("spo11sae2_f", "spo11sae2_3", seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	}

	async function drawRow(refName, row_name, seg_row_height, seg_row_separate, bp_size, view_length, max_view_width) {
		const pos_map = pos_ref_map_list[get_genome_index(refName)];
		const innerText = `${refName} ${row_name}`;
		const el_div = appendLeftTitle(innerText, row_height, seg_row_separate, 1, 100);

		const bRet = await drawRNACoverage(refName, row_name, pos_map, row_height, bp_size, view_length, max_view_width, stroke_black);
		if (!bRet) {
			el_div.remove();
		}
		else {
			main_ctx.translate(0, row_height + seg_row_separate);
		}
	}

	/**
	 * @param {CanvasRenderingContext2D} ctx
	 * @param {number} value
	 * @param {number} max_value
	 * @param {number} ref_pos
	 * @param {number} percent
	 * @param {number} pos
	 * @param {number} array_cursor
	 */
	function stroke_black(ctx, value, max_value, ref_pos, percent, pos, array_cursor) {
		ctx.globalAlpha = globalAlpha;
		if (at_list.find(island => island.start <= pos && pos <= island.end)) {
			ctx.strokeStyle = "#FF0000";
		}
		else {
			ctx.strokeStyle = "#000000";//value > 300 ? "#FF0000" : "#000000";
		}
	}
}

/**
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 *
 * @param {number} bp_size
 * @param {number} view_length
 */
async function drawSampleRows(max_view_width, seg_row_height, seg_row_separate, bp_size, view_length) {
	const info_map = {
		"QM6a": {
			to_ma_posmap: pos_ref1_uint32array,
		},
		"CBS1-1": {
			to_ma_posmap: pos_ref2_uint32array,
		},
	};
	const row_width = max_view_width;

	const el_appendRow = document.getElementById("append-row");

	await drawTempData_beforeMethyl(seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);

	// veg QM6a/CBS1-1
	await draw_sample_group("veg", "QM6a", "#0123AB", $color_set_view.dad_bk);
	await draw_sample_group("veg", "CBS1-1", "#AB2301", $color_set_view.mom_bk);

	if (viewerState.D8_bg_RIP_marker) {
		await draw_sample_group("D4", "QM6a", "#0123AB", $color_set_view.dad_bk);
		[
			allMarker.map.rip_Q,
			allMarker.map.rip_2_Q,
		].forEach(markersData => {
			const row_height = seg_row_height * 1.5;

			_padding_top(false, 1);

			// border
			main_ctx.beginPath();
			main_ctx.strokeStyle = "black";
			main_ctx.strokeRect(0, 0, row_width, row_height);

			drawMarkerRow(markersData, bp_size, max_view_width, row_height, seg_row_separate, false);

			_padding_top(false, -1);
		});
		await draw_sample_group("D8", "QM6a", "#0123AB", $color_set_view.dad_bk);
	}
	else {
		// D1 ~ D8 QM6a
		await Object.keys(display_methy_refMap_sampleList).slice(1).reduce(async function (prev_task, groupName, groupIdx) {
			await prev_task;
			await draw_sample_group(groupName, "QM6a", "#0123AB", $color_set_view.dad_bk);
		}, Promise.resolve());
	}

	if (NP32_NP30_diff) {
		const diff_rows = NP32_NP30_diff;

		const row_height = Math.trunc(seg_row_height * 2);

		const el_div = document.createElement("div");
		el_appendRow.append(el_div);
		const innerText = `NP32 - NP30 CG`;
		el_div.outerHTML = `<div title="methy diff" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

		draw_methy_diff_row(row_width, row_height, diff_rows, "#FF0000", "#0000FF", "#000000");

		main_ctx.translate(0, 1 * (row_height + seg_row_separate));
	}

	if (viewerState.D8_bg_RIP_marker) {
		await draw_sample_group("D4", "CBS1-1", "#AB2301", $color_set_view.mom_bk);
		[
			allMarker.map.rip_C,
			allMarker.map.rip_2_C,
		].forEach(markersData => {
			const row_height = seg_row_height * 1.5;

			_padding_top(false, 1);

			// border
			main_ctx.beginPath();
			main_ctx.strokeStyle = "black";
			main_ctx.strokeRect(0, 0, row_width, row_height);

			drawMarkerRow(markersData, bp_size, max_view_width, row_height, seg_row_separate, false);

			_padding_top(false, -1);
		});
		await draw_sample_group("D8", "CBS1-1", "#AB2301", $color_set_view.mom_bk);
	}
	else {
		// D1 ~ D8 CBS1-1
		await Object.keys(display_methy_refMap_sampleList).slice(1).reduce(async function (prev_task, groupName, groupIdx) {
			await prev_task;
			await draw_sample_group(groupName, "CBS1-1", "#AB2301", $color_set_view.mom_bk);
		}, Promise.resolve());
	}

	await g_methylRenderer.draw_all_methyl_row(el_appendRow, row_width, seg_row_height, seg_row_separate);

	if (drawTempData_afterMethyl) {
		await drawTempData_afterMethyl(seg_row_height, seg_row_separate, bp_size, view_length, max_view_width);
	}

	need_merge_rna_reads = false;

	// Object.keys(methy_ratio_map).forEach(refId => {
	// });

	/**
	 * @param {"veg"|"D4"|"D8"} groupName
	 * @param {"QM6a"|"CBS1-1"} refId
	 * @param {string} color
	 * @param {string} border_color
	 */
	async function draw_sample_group(groupName, refId, color, border_color, padding_top = true) {
		if (!methy_ratio_map) {
			return;
		}

		const methyRatioMap = methy_ratio_map?.[refId];
		const nBS_methyRatioMap = nBS_methy_ratio_map?.[refId];
		const BS_nBS_match_methyRatioMap = BS_nBS_match_methy_map?.[refId];
		const methyDiffMap = methy_diff_map?.[refId];
		// const BS_nBS_methyDiffMap = methy_BS_nBS_diff_map[refId];

		//see line:295
		await display_methy_refMap_sampleList[groupName][refId].reduce(async function (prev_task, sample, sampleIndex) {
			await prev_task;

			const sampleName = sample.title;

			// const sampleTypeIndex = sampleIndex >>> 1;//(int)(sampleIndex / 2) // wt + rid
			const sampleTypeIndex = sampleIndex;//
			if (!sample.display) {
				return;
			}

			if (padding_top) {
				_padding_top();
			}

			if (viewerState.display_meth_ratio) {
				const ratio_rows = methyRatioMap[sampleName];
				if (ratio_rows) {
					const row_height = Math.trunc(seg_row_height * 1.5);

					const el_div = document.createElement("div");
					el_appendRow.append(el_div);
					// const innerText = `${sampleName} ${ref_id_map[refId]}`;
					const innerText = `${sampleName} ${name_mapto_html[methyl_ref_id_map[refId]]}`;
					el_div.outerHTML = `<div title="methy ratio" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

					alert("draw_methy_ratio_row(config, ratio_rows, ...");
					await draw_methy_ratio_row(ratio_rows, row_width, row_height, color, border_color, "5mC/(5mC+C)");

					main_ctx.translate(0, 1 * (row_height + seg_row_separate));
				}
			}

			if (viewerState.display_nBS_meth_ratio && nBS_methyRatioMap && nBS_methyRatioMap[sampleName]) {
				const ratio_rows = nBS_methyRatioMap[sampleName];
				if (ratio_rows) {
					const row_height = Math.trunc(seg_row_height * 1.5);

					const el_div = document.createElement("div");
					el_appendRow.append(el_div);
					const innerText = `${sampleName} ${name_mapto_html[methyl_ref_id_map[refId]]} (nBS)`;
					el_div.outerHTML = `<div title="methy ratio" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

					alert("draw_methy_ratio_row(config, ratio_rows, ...");
					await draw_methy_ratio_row(ratio_rows, row_width, row_height, color, border_color, "C/C");

					main_ctx.translate(0, 1 * (row_height + seg_row_separate));
				}
			}

			if (viewerState.display_BS_nBS_meth_ratio && BS_nBS_match_methyRatioMap && BS_nBS_match_methyRatioMap[sampleName]) {
				const ratio_rows = BS_nBS_match_methyRatioMap[sampleName];
				if (ratio_rows) {
					const row_height = Math.trunc(seg_row_height * 1.5);

					const el_div = document.createElement("div");
					el_appendRow.append(el_div);
					const innerText = `${sampleName} ${name_mapto_html[methyl_ref_id_map[refId]]} (BS ∩ nBS)`;
					el_div.outerHTML = `<div title="methy ratio" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

					alert("draw_methy_ratio_row(config, ratio_rows, ...");
					await draw_methy_ratio_row(ratio_rows, row_width, row_height, color, border_color, "5mC/(5mC+C)");

					main_ctx.translate(0, 1 * (row_height + seg_row_separate));
				}
			}

			if (viewerState.display_RNA_coverage && info_map[refId] && info_map[refId].to_ma_posmap) {
				const row_width = max_view_width;
				const row_height = Math.trunc(seg_row_height * 1.5);

				const el_div = document.createElement("div");
				el_appendRow.append(el_div);
				const innerText = `${sampleName} RNA coverage`;
				el_div.outerHTML = `<div style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

				if (need_merge_rna_reads) {
					_merge_rna_reads(refId, sampleName);
				}

				drawRNACoverageOuter(refId, sampleName, row_width, row_height, seg_row_separate, border_color);
				await drawRNACoverage(refId, sampleName, info_map[refId].to_ma_posmap, row_height, bp_size, view_length, row_width, null);
				main_ctx.translate(0, 1 * (row_height + seg_row_separate));
			}

			if (viewerState.display_meth_diff_wt) {
				const diff_mt_rows = methyDiffMap[sampleName];
				if (diff_mt_rows) {
					MethyContextTypeList.forEach(contextType => {
						const diff_rows = diff_mt_rows[contextType];

						const row_height = Math.trunc(seg_row_height * 2);

						const el_div = document.createElement("div");
						el_appendRow.append(el_div);
						const innerText = `${sampleName}-WT(veg) ${name_mapto_html[methyl_ref_id_map[refId]]} ${contextType}`;
						el_div.outerHTML = `<div title="methy diff" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;

						draw_methy_diff_row(row_width, row_height, diff_rows, "#FF0000", "#0000FF", "#000000");

						main_ctx.translate(0, 1 * (row_height + seg_row_separate));
					});
				}
			}

			// if (BS_nBS_methyDiffMap) {
			// 	const nBS_diff_rows = BS_nBS_methyDiffMap[sampleName];
			// 	if (nBS_diff_rows) {
			// 		const row_height = Math.trunc(seg_row_height * 2);
			// 		const el_div = document.createElement("div");
			// 		el_appendRow.append(el_div);
			// 		const innerText = `${sampleName}-${sampleName}(nBS) ${ref_id_map[refId]}`;
			// 		el_div.outerHTML = `<div title="nBS methy diff" style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1.5}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${innerText}</span></div>`;
			// 		draw_methy_diff_row(row_width, row_height, nBS_diff_rows, "#FF0000", "#0000FF", "#000000");
			// 		main_ctx.translate(0, 1 * (row_height + seg_row_separate));
			// 	}
			// }
		}, Promise.resolve());

		//
	}

	function _padding_top(append_html = true, scale = 1, row_height = seg_row_height * 0) {
		if (append_html) {
			const el_div = document.createElement("div");
			el_appendRow.append(el_div);
			const innerText = ``;
			el_div.outerHTML = `<div style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle; display: block; height: ${row_height + 1}px;">${innerText}</span></div>`;
			el_div.title = "sample separate";
		}

		main_ctx.translate(0, scale * 1 * (row_height + seg_row_separate));
	}
}

/**
 * @param {module_Methyl_sampleData} config
 * @param {Partial<module_Methyl_ratioData>[]} rows
 * @param {number} row_width
 * @param {number} row_height
 * @param {string} color
 * @param {string} border_color
 * @see {@link module_MethylRenderer#load_np_meth_tsv}
 * @see {@link module_MethylRenderer#load_meth_ratio_float32}
 */
async function draw_methy_ratio_row(config, rows, row_width, row_height, color, border_color) {
	const row_hheight = row_height * 0.5;

	main_ctx.save();

	const view_length = g_chrCoord.bp_end - g_chrCoord.bp_start + 1;
	const scale = 1 / view_length;
	const bp_size = row_width * scale;

	// border
	main_ctx.beginPath();
	main_ctx.rect(0, 0 - viewerState.row_padding_top, row_width, row_height + viewerState.row_padding_bottom);
	main_ctx.strokeStyle = config.border_color ?? border_color;
	main_ctx.stroke();
	
	main_ctx.clip();

	// ratio = 5mC/(5mC+C)
	if (config.value_desc) {//if (window.draw_methy_ratio_desc)
		main_ctx.fillStyle = "black";
		main_ctx.textAlign = "left";
		main_ctx.font = viewerState.global_font_style;//Math.trunc(typeof window.methylated == "number" ? window.methylated : 14) + "px Arial";//(row_height / 3))

		main_ctx.textBaseline = "middle";
		main_ctx.fillText(config.value_desc, row_width + 5, row_hheight);
	}
	else {
		if (config.canvas_value_desc && config.canvas_value_desc.height) {
			main_ctx.save();
			try {
				const tx = main_ctx.getTransform();
				main_ctx.setTransform(1, 0, 0, 1, Math.trunc(tx.e), Math.trunc(tx.f));
				main_ctx.drawImage(config.canvas_value_desc, row_width + 5, Math.trunc((row_height - config.canvas_value_desc.height) / 2));
			}
			finally {
				main_ctx.restore();
			}
		}
	}

	// middle line
	if (config.mid_line) {
		main_ctx.beginPath();
		main_ctx.moveTo(0, row_hheight);
		main_ctx.lineTo(row_width, row_hheight);
		main_ctx.strokeStyle = "#80808080";
		main_ctx.stroke();
	}

	// 20221003: red line
	if (config.beforeRenderList) {
		for (let ahl of config.beforeRenderList) {
			ahl(config, row_width, row_height);
		}
	}
	
	const value_transformer = config.getValueTransformer();

	// const pos_map = {
	// 	[dataset.genomeNameList[0]]: ref1_pos_uint32array,
	// 	[dataset.genomeNameList[1]]: ref2_pos_uint32array,
	// 	[dataset.genomeNameList[2]]: window.subject_1_pos_map,
	// 	[dataset.genomeNameList[3]]: window.subject_2_pos_map,
	// 	[dataset.genomeNameList[4]]: window.subject_3_pos_map,
	// 	[dataset.genomeNameList[5]]: window.subject_4_pos_map,
	// }[config.ref];

	// let cnt_inst = 0;

	if (config.region) {
		// main_ctx.beginPath();

		// const visibility = rows.reduce((prev, data) => {
		// 	const pos1 = data.start;//see: load_methyDiff
		// 	const pos2 = data.end;//see: load_methyDiff
		// 	if (pos1 > g_chrCoord.bp_start && pos1 < g_chrCoord.bp_end) {
		// 		return prev + 1;
		// 	}
		// 	else {
		// 		return prev;
		// 	}
		// }, 0);

		// // @ts-ignore
		// const $$$aaa = window.$$$aaa;
		// const display_region_border = (visibility / rows.length) < ($$$aaa ?? 0.03);// view_length <= 20_000;
		// if (display_region_border) {
		// 	main_ctx.strokeStyle = "#F00";
		// 	main_ctx.fillStyle = "#0004";
		// }
		// else {
		// 	main_ctx.strokeStyle = "#F008";
		// 	main_ctx.fillStyle = "#F004";
		// }

		// for (let i = 0; i < rows.length; ++i) {
		// 	const data = rows[i];
		// 	const pos1 = data.start;//pos_map[data.start - 1];//see: load_methyDiff
		// 	const pos2 = data.end;//pos_map[data.end - 1];//see: load_methyDiff

		// 	if (pos1 > g_chrCoord.bp_start && pos1 < g_chrCoord.bp_end) {
		// 		const x1 = ((1 - g_chrCoord.bp_start) + pos1) * bp_size;
		// 		if (x1 >= 0 && x1 <= row_width) {
		// 			const value = data.value;

		// 			const x2 = ((1 - g_chrCoord.bp_start) + pos2) * bp_size;

		// 			const y1 = row_height - value * row_height - 1;

		// 			let len = x2 - x1;
		// 			if (config.strong > 0) {
		// 				len = Math.max(x2 - x1, config.strong);
		// 			}

		// 			if (len > 1) {
		// 				main_ctx.beginPath();
		// 				main_ctx.rect(x1, y1, len, row_height - y1);

		// 				// if (display_region_border) {
		// 				// 	const gradient = main_ctx.createLinearGradient(0, 0, len, row_height);
		// 				// 	gradient.addColorStop(0, "#F00");
		// 				// 	gradient.addColorStop(1, "#FFF8");
		// 				// 	main_ctx.fillStyle = gradient;
		// 				// }
		// 				main_ctx.fill();
		// 				if (display_region_border) {
		// 					main_ctx.stroke();
		// 				}
		// 			}
		// 			else {
		// 				main_ctx.beginPath();

		// 				main_ctx.moveTo(x1, y1);
		// 				main_ctx.lineTo(x1, row_height);

		// 				main_ctx.moveTo(x2, y1);
		// 				main_ctx.lineTo(x2, row_height);

		// 				main_ctx.stroke();
		// 			}

		// 			// const data2 = rows[i];
		// 			// const pos3 = data2.start;//see: load_methyDiff
		// 			// const x3 = ((1 - g_chrCoord.bp_start) + pos3) * bp_size;
		// 			// const y3 = row_height - data2.value * row_height - 1;
		// 			// main_ctx.beginPath();
		// 			// main_ctx.moveTo(x2, y1);
		// 			// main_ctx.lineTo(x2, row_height);
		// 			// main_ctx.lineTo(x3, y3);
		// 			// main_ctx.fill();
		// 		}
		// 	}
		// }

	}
	else {
		// move rectangle to center
		// main_ctx.lineWidth = Math.min(Math.max(1, bp_size), 24);
		main_ctx.translate(-bp_size * 0.5, 0);

		if (config.density_to_opacity) {
			if (view_length >= 1_000_000) {
				main_ctx.globalAlpha = 0x33 / 0x255;
			}
			else if (view_length >= 500_000) {
				main_ctx.globalAlpha = 0x55 / 0x255;
			}
			else {
				main_ctx.globalAlpha = 1;
			}
		}
		else if (config.func_density_to_opacity) {
			let num_point = 0;
			// const ga = Math.min(Math.max(0.01, vp / view_length), 1);
			forEachValue(() => true, (pt, x1, x2) => {
				num_point += 1;
			});
			main_ctx.globalAlpha = config.func_density_to_opacity(view_length, num_point, row_width);
		}
		// else if (config.full_anti_aliasing) {
		// 	let vp = 0;
		// 	forEachValue(() => true, (pt, x1, x2) => {
		// 		vp += 1;
		// 	});
		// 	const ga = Math.min(Math.max(0.01, vp / view_length), 1);
		// 	// config.func_density_to_opacity = () => ga;
		// 	// config.density_to_opacity = false;
		// 	// config.full_anti_aliasing = false;
		// 	main_ctx.globalAlpha = ga;
		// }
		else {
			main_ctx.globalAlpha = 1;
		}

		if (config.rendering_condition?.length > 0) {
			for (let cond of config.rendering_condition) {
				await drawAllValue(cond.color, cond.condition, cond.min_width);
			}
		}
		else {
			await drawAllValue("#E00000", d => d.value > 0, 1);
			await drawAllValue("#0000E0", d => d.value < 0, 1);
			await drawAllValue("#000000", d => d.value == 0, 1);
		}

		main_ctx.globalAlpha = 1;//reset
	}
	
	// 20221003: red line
	if (config.afterRenderList) {
		for (let ahl of config.afterRenderList) {
			ahl(config, row_width, row_height);
		}
	}

	main_ctx.restore();

	/**
	 * @param {string} value_color
	 * @param {(point: module_Methyl_ratioData) => boolean} value_filter
	 * @param {number} min_width
	 */
	async function drawAllValue(value_color, value_filter, min_width) {
		// main_ctx.strokeStyle = value_color;
		main_ctx.fillStyle = value_color;

		if (config.fast && rows[0]) {
			const len = g_chrCoord.bp_end - g_chrCoord.bp_start + 1;
			const wz = config.fast_binSize ?? Math.round(len / row_width);

			if (config.fast_binSize > 1) {
				// if (rows[0].g_chrCoord.bp_start != g_chrCoord.bp_start || rows[0].g_chrCoord.bp_end != g_chrCoord.bp_end) {
				// 	rows[0].g_chrCoord.bp_start = g_chrCoord.bp_start;
				// 	rows[0].g_chrCoord.bp_end = g_chrCoord.bp_end;

					// const t1 = new Date().getTime();

					// let cnt = 0;
					let prev_i = 0;
					let prev_px = 0;
					let fin_value = 0;
					for (let i = 0; i < rows.length; ++i) {
						const data = rows[i];
						if (value_filter(data)) {
							data.skip = true;

							const px = Math.trunc(data.start / wz);

							fin_value = Math.max(Math.abs(data.value), Math.abs(fin_value));

							// cnt += 1;
							prev_i = i;

							if (px != prev_px) {
								rows[prev_i].skip = false;
								rows[prev_i].fin_value = Math.sign(rows[prev_i].value) * fin_value;
								// rows[prev_i].alpha = 1 / cnt;

								prev_px = px;
								// cnt = 0;
								fin_value = 0;
							}
						}
					}

					// const t2 = new Date().getTime();
					// console.log((t2 - t1) / 1000);
				// }
			}
			else if (config.fast_binSize === undefined) {
				for (let i = 0; i < rows.length; ++i) {
					const data = rows[i];
					delete data.skip;
					delete data.fin_value;
				}
			}
			else { // if (config.fast_binSize <= 0)
				for (let i = 0; i < rows.length; ++i) {
					const data = rows[i];
					data.skip = false;
					data.fin_value = data.value;
				}
			}
		}
		// else if (config.fast === false) {
		// }

		// 1_000_000 / Math.ceil(1_000_000 / 2345 * 1000)

		const start_time_stamp = new Date().getTime();

		main_ctx.beginPath();
		
		forEachValue(value_filter, function default_forEachValue(data, x1, x2) {
			const value = value_transformer(config.fast ? data.fin_value : data.value);
			if (config.drawShape) {
				config.drawShape(x1, x2, min_width, row_height, value, config.display_minus_value); //20221017
			}
			else {
				draw_meth_ratio_value(x1, x2, min_width, row_height, value, config.display_minus_value);
			}
		});

		// main_ctx.lineCap = "butt";
		// main_ctx.lineJoin = "bevel";
		main_ctx.fill();
		// main_ctx.lineCap = "round";
		// main_ctx.lineJoin = "round";

		config._elapsed_ms = (new Date().getTime()) - start_time_stamp;// render command

		if (config._elapsed_ms > 1000) {
			await delayFrame();
		}
		config._elapsed_ms = (new Date().getTime()) - start_time_stamp;// real time elapsed of render

		// max_display_value
	}

	/**
	 * @param {(point: module_Methyl_ratioData) => boolean} value_filter
	 * @param {(data: Partial<module_Methyl_ratioData>, x1: number, x2: number) => void} drawValue
	 */
	function forEachValue(value_filter, drawValue) {
		for (let i = 0; i < rows.length; ++i) {
			const data = rows[i];
			if (data.skip != true && value_filter(data)) {
				/**
				 * @see {module_MethylRenderer}
				 */
				let start_pos = config.ref_to_pos_map ? config.ref_to_pos_map[data.start] : data.start; //pos_map[data.start - 1];
				let end_pos = config.ref_to_pos_map ? config.ref_to_pos_map[data.end] : data.end; //pos_map[data.end - 1];
				const len = end_pos - start_pos;

				if (start_pos > end_pos) {
					[start_pos, end_pos] = [end_pos, start_pos];
				}

				// if (config.range_repeat_density && len > config.range_repeat_density) {
				// 	if (end_pos > start_pos) {
				// 		for (let pos = start_pos; pos <= end_pos; pos += config.range_repeat_density) {
				// 			draw_value_point(pos, data);
				// 		}
				// 	}
				// 	else {
				// 		for (let pos = end_pos; pos <= start_pos; pos += config.range_repeat_density) {
				// 			draw_value_point(pos, data);
				// 		}
				// 	}
				// }
				// else {
				// 	draw_value_point(start_pos, data);
				// }
				if (start_pos <= g_chrCoord.bp_end && end_pos >= g_chrCoord.bp_start) {
					const x1 = ((1 - g_chrCoord.bp_start) + start_pos - 0.5) * bp_size;
					const x2 = ((1 - g_chrCoord.bp_start) + end_pos + 1 - 0.5) * bp_size;
					if (x2 >= 0 && x1 <= row_width) {
						drawValue(data, x1, x2);
						// ++cnt_inst;
					}
				}
			}
		}
	}// forEachValue

	// console.log(config.name, cnt_inst);
}

/**
 * @param {number} x1
 * @param {number} x2
 * @param {number} min_width
 * @param {number} row_height
 * @param {number} value
 * @param {boolean} display_minus_value
 */
function draw_meth_ratio_value(x1, x2, min_width, row_height, value, display_minus_value) {
	// x1 = Math.trunc(x1);
	const w = Math.max(min_width, x2 - x1);

	if (display_minus_value) {
		const hh = row_height * 0.5;
		if (value > 0) {
			// 1 ~ Infinity
			const h = Math.max(1, value * hh);
			const top_y = hh - h;
			main_ctx.moveTo(x1, top_y);
			// main_ctx.lineTo(x1, hh);
			main_ctx.rect(x1, top_y, w, h);
		}
		else if (value < 0) {
			// -Infinity ~ -1
			const h = Math.min(-1, value * hh);
			const top_y = hh - h;
			main_ctx.moveTo(x1, top_y);
			// main_ctx.lineTo(x1, hh);
			main_ctx.rect(x1, top_y, w, h);
		}
		else {
			// // main_ctx.fillRect(x - hw, top_y - 0.5, w, 1);
			main_ctx.moveTo(x1, hh - 0.5);
			// main_ctx.lineTo(x1, hh + 0.5);
			main_ctx.rect(x1, hh - 0.5, w, 1);
		}
	}
	else {
		value = Math.abs(value);

		if (value > 0) {
			// 1 ~ Infinity
			const h = Math.max(1, value * row_height);
			const top_y = row_height - h;
			// main_ctx.moveTo(x1, top_y);
			// main_ctx.lineTo(x1, row_height);
			main_ctx.rect(x1, top_y, w, h);
		}
		else if (value < 0) {
			// -Infinity ~ -1
			const h = Math.min(-1, value * row_height);
			const top_y = row_height - h;
			// main_ctx.moveTo(x1, top_y);
			// main_ctx.lineTo(x1, row_height);
			main_ctx.rect(x1, top_y, w, h);
		}
		else {
			// // main_ctx.fillRect(x - hw, top_y - 0.5, w, 1);
			// main_ctx.moveTo(x1, row_height - 0.5);
			// main_ctx.lineTo(x1, row_height + 0.5);
			main_ctx.rect(x1, row_height - 0.5, w, 1);
		}
	}
}

/**
 * @param {number} row_width
 * @param {number} row_height
 * @param {{ pos: number, value: number }[]} data
 * @param {string} color_up
 * @param {string} color_down
 * @param {string} border_color
 */
function draw_methy_diff_row(row_width, row_height, data, color_up, color_down, border_color) {
	const row_hheight = row_height * 0.5;

	main_ctx.save();

	const scale = 1 / (g_chrCoord.bp_end - g_chrCoord.bp_start + 1);
	const bp_size = row_width * scale;

	main_ctx.beginPath();
	main_ctx.strokeStyle = border_color;
	main_ctx.strokeRect(0, 0 - viewerState.row_padding_top, row_width, row_height + viewerState.row_padding_bottom);

	if (window.draw_methy_diff_desc) {
		main_ctx.fillStyle = "black";
		main_ctx.textAlign = "left";
		main_ctx.font = viewerState.global_font_style;//"16pt Arial";//Math.trunc(typeof window.methylated == "number" ? window.methylated : 14) + "px Arial";//(row_height / 3))

		main_ctx.textBaseline = "middle";//"hanging";
		main_ctx.fillText("methylated", row_width + 5, 0 + row_height * 0.25);//100//methylated

		main_ctx.textBaseline = "middle";//"bottom";
		main_ctx.fillText("de-methylated", row_width + 5, row_height - row_height * 0.25);//50//de-methylated
	}

	// middle line
	main_ctx.beginPath();
	main_ctx.moveTo(0, row_hheight);
	main_ctx.lineTo(row_width, row_hheight);
	main_ctx.strokeStyle = "#80808080";
	main_ctx.stroke();

	main_ctx.lineWidth = Math.min(Math.max(1, bp_size), 24);
	main_ctx.translate(-bp_size * 0.5, 0);

	drawDiff(-1);
	main_ctx.strokeStyle = color_up;//darkblue, deepblue
	main_ctx.stroke();

	drawDiff(1);
	main_ctx.strokeStyle = color_down;//"#000000"
	main_ctx.stroke();

	/**
	 * @param {1|-1} flag
	 */
	function drawDiff(flag) {
		main_ctx.beginPath();
		for (let i = 0; i < data.length; ++i) {
			const pos = data[i].pos;//see: load_methyDiff
			if (pos >= g_chrCoord.bp_start && pos <= g_chrCoord.bp_end) {
				const x = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
				if (x >= 0 && x <= row_width) {
					const value = data[i].value;

					if (flag < 0 && value < 0) {
						main_ctx.moveTo(x, row_hheight + value * row_hheight);
						main_ctx.lineTo(x, row_hheight);
					}
					else if (flag > 0 && value > 0) {
						main_ctx.moveTo(x, row_hheight + value * row_hheight);
						main_ctx.lineTo(x, row_hheight);
					}
				}
			}
		}
	}

	// //ctx.beginPath();
	// //ctx.moveTo(0, plot_height);
	// for (let i = 0; i < data.length; ++i) {
	// 	const pos = data[i].pos;//see: load_methy
	// 	if (pos >= g_chrCoord.bp_start && pos <= g_chrCoord.bp_end) {
	// 		const x = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
	// 		const value = data[i].value;

	// 		main_ctx.beginPath();
	// 	//	ctx.lineTo(data[i].pos * row_width, plot_height * 0.5 + (0 + data[i].value) * plot_height * 0.5);
	// 		main_ctx.moveTo(x, row_hheight + value * row_hheight);
	// 		main_ctx.lineTo(x, row_hheight);

	// 		//ctx.arc(data[i].pos * row_width, plot_height * 0.5 + data[i].value * plot_height * 0.5, )
	// 	}
	// }
	// //ctx.moveTo(row_width, plot_height);

	// // main_ctx.fillStyle = "rgba(127, 255, 32, 0.5)";
	// // main_ctx.fill();

	// // main_ctx.strokeStyle = "rgba(0, 0, 0, 0.5)";
	// // main_ctx.stroke();

	main_ctx.restore();
}

/**
 * @deprecated
 *
 * @param {string} refId
 * @param {string} sampleName
 *
 * @param {number} row_width = max_view_width
 * @param {number} row_height
 * @param {number} seg_row_separate
 * @param {string} border_color = $color_set_view.dad_bk
 */
function drawRNACoverageOuter(refId, sampleName, row_width, row_height, seg_row_separate, border_color) {
	const row_hheight = row_height * 0.5;

	// border
	main_ctx.beginPath();
	main_ctx.strokeStyle = border_color;
	main_ctx.strokeRect(0, 0 - viewerState.row_padding_top, row_width, row_height + viewerState.row_padding_bottom);

	// middle line
	main_ctx.beginPath();
	main_ctx.moveTo(0, row_hheight);
	main_ctx.lineTo(row_width, row_hheight);
	main_ctx.strokeStyle = "#80808080";
	main_ctx.stroke();

	// if (0) {//left text
	// 	main_ctx.font = viewerState.global_font_style;//"12pt Arial";
	// 	main_ctx.textAlign = "left";
	// 	main_ctx.textBaseline = "top";

	// 	main_ctx.fillStyle = "rgba(111, 111, 111, 0.7)";
	// 	main_ctx.fillText(sampleName, 0, 0);

	// 	main_ctx.strokeStyle = "rgba(99, 99, 99, 0.5)";
	// 	main_ctx.strokeText(sampleName, 0, 0);
	// }
	// else {
	// 	const el_div = document.createElement("div");
	// 	const el_appendRow = document.getElementById("append-row");
	// 	el_appendRow.append(el_div);
	// 	el_div.outerHTML = `<div style="line-height: ${row_height}px; margin-bottom: ${seg_row_separate - 1}px;"><span contenteditable="true" spellcheck="false" style="vertical-align: middle;">${sampleName}</span></div>`
	// }
	{
		//const max_value = rnaCoverage_max_map[refId][sampleName];
		// const norm = rnaCoverage_norm_map[refId][sampleName];

		main_ctx.fillStyle = "black";
		main_ctx.textAlign = "left";
		main_ctx.font = viewerState.global_font_style;//"16pt Arial";//Math.trunc(row_height / 3 / 1.5) + "px Arial";

		if (typeof window.rna_reads_func == "function") {
			const norm = rnaCoverage_norm_map[refId][sampleName] ?? 1;
			const max_value = viewerState.rna_reads_max_display_value ?? rnaCoverage_max_map[refId][sampleName] * norm;
			const rna_reads_func = x => window.rna_reads_func(x, max_value) * max_value;

			main_ctx.textBaseline = "hanging";
			main_ctx.fillText((rna_reads_func(max_value)).toFixed(0), row_width + 5, 0);//100%

			main_ctx.textBaseline = "middle";
			main_ctx.fillText((rna_reads_func(max_value * 0.5)).toFixed(0), row_width + 5, row_hheight);//50%
		}
		else {//log
			main_ctx.textBaseline = "hanging";
			main_ctx.fillText((get_rna_plot_display_value(1)).toFixed(0), row_width + 5, 0);//100%

			main_ctx.textBaseline = "middle";
			main_ctx.fillText((get_rna_plot_display_value(0.5)).toFixed(0), row_width + 5, row_hheight);//50%
		}

		main_ctx.textBaseline = "ideographic";// alphabetic
		main_ctx.fillText((0).toFixed(0), row_width + 5, row_height);//0%
	}
}

/**
 * @param {number} max_value
 * @param {number} row_width = max_view_width
 * @param {number} row_height
 * @param {string} border_color = $color_set_view.dad_bk
 */
function draw_RNA_coverage_border(max_value, row_width, row_height, border_color) {
	main_ctx.save();

	const row_hheight = row_height * 0.5;

	// border
	main_ctx.beginPath();
	main_ctx.strokeStyle = border_color;
	main_ctx.strokeRect(0, 0, row_width, row_height);

	// middle line
	main_ctx.beginPath();
	main_ctx.moveTo(0, row_hheight);
	main_ctx.lineTo(row_width, row_hheight);
	main_ctx.strokeStyle = "#80808080";
	main_ctx.stroke();

	main_ctx.fillStyle = "black";
	main_ctx.textAlign = "left";
	main_ctx.font = viewerState.global_digi_font_style;//"16pt Arial";//Math.trunc(row_height / 3 / 1.5) + "px Arial";

	if (typeof window.rna_reads_func == "function") {
		const max_v = viewerState.rna_reads_max_display_value ?? max_value;

		const rna_reads_func = x => window.rna_reads_func(x, max_v) * max_v;

		main_ctx.textBaseline = "hanging";
		main_ctx.fillText((rna_reads_func(max_v)).toFixed(0), row_width + 5, 0);//100%

		main_ctx.textBaseline = "middle";
		main_ctx.fillText((rna_reads_func(max_v * 0.5)).toFixed(0), row_width + 5, row_hheight);//50%

		main_ctx.textBaseline = "ideographic";// alphabetic
		main_ctx.fillText((0).toFixed(0), row_width + 5, row_height);//0%
	}
	else {//log
		// const log_exp_font = viewerState.global_digi_font_style.replace(/\d+/, (parseInt(viewerState.global_digi_font_style) * 0.5).toFixed(0));
		// const log_base_string = rna_display_log_base.toString();

		for (let i = 0; i < rna_display_middleValue_list.length; ++i) {
			const midVal = rna_display_middleValue_list[i];
			// const exp = rna_display_exp * midVal;

			if (i == 0) {
				main_ctx.textBaseline = "hanging";
			}
			else if (i == rna_display_middleValue_list.length - 1) {// last
				main_ctx.textBaseline = "ideographic";// alphabetic
			}
			else {
				main_ctx.textBaseline = "middle";
			}

			main_ctx.fillText((get_rna_plot_display_value(midVal)).toFixed(0), row_width + 5, row_height - row_height * midVal);

			// main_ctx.font = viewerState.global_digi_font_style;
			// const { width, actualBoundingBoxAscent } = main_ctx.measureText(log_base_string);
			// main_ctx.fillText(log_base_string, row_width + 5, row_height - row_height * midVal);//50%

			// main_ctx.textBaseline = "ideographic";// alphabetic
			// main_ctx.font = log_exp_font;
			// main_ctx.fillText(exp.toString(), row_width + 5 + width, row_height - row_height * midVal);//50%
		}

		// main_ctx.textBaseline = "hanging";
		// main_ctx.fillText((get_rna_plot_display_value(1)).toFixed(0), row_width + 5, 0);//100%

		// main_ctx.textBaseline = "middle";
		// // main_ctx.fillText((get_rna_plot_display_value(0.5)).toFixed(0), row_width + 5, row_hheight);//50%

		// 0%
		// ...
	}

	main_ctx.restore();
}

/**
 * @returns {Promise<HTMLDivElement>}
 */
async function MathML2SVG(math) {
	if (typeof window.MathJax != "object") {
		alert("dynamic load MathJax");
		// await import("https://cdn.jsdelivr.net/npm/mathjax@3/es5/mml-svg.js");
	}

	// @ts-ignore
	const MathJax = window.MathJax;

// 	const math = `
// 	<math>
// 		<msup>
// 			<mn>2</mn>
// 			<mn>5</mn>
// 		</msup>
// 	</math>
// `;
	/**
	 * @type {HTMLDivElement}
	 */
	const el = await MathJax.mathml2svgPromise(math.trim(), {
		display: true,// important
	});

	return el;
}

function get_rna_plot_display_value(delta) {
	return rna_display_log_base ** (rna_display_exp * delta);
}

/**
 * @param {number} max_value
 */
function get_rna_reads_func(max_value) {
	const max_display_value = get_rna_plot_display_value(1);
	let rna_reads_func;

	if (typeof window.rna_reads_func == "function") {
		rna_reads_func = x => window.rna_reads_func(x, viewerState.rna_reads_max_display_value ?? max_value);
	}
	else {
		// rna_reads_func = x => Math.min(x, max_display_value) ** rna__exp;
		rna_reads_func = x => rna_display_log_func(Math.min(x, max_display_value)) / rna_display_exp;
		// rna_reads_func = x => rna_display_log_func(x) / rna_display_log_func(max_value);
	}

	return rna_reads_func;
}

/**
 * @param {string} refId
 * @param {string} sampleName
 */
function _merge_rna_reads(refId, sampleName) {
	const w_ui32a = w_rnaCoverage_ui32a_map[refId][sampleName];
	const m_ui32a = m_rnaCoverage_ui32a_map[refId][sampleName];
	if (rna_func_type == "mean") {
		for (let i = 0; i < m_ui32a.length; ++i) {
			for (let j = 0, p = i * rna_window_size; j < rna_window_size && p < w_ui32a.length; ++j, ++p) {
				m_ui32a[i] += w_ui32a[p];
			}
			m_ui32a[i] = m_ui32a[i] / rna_window_size;
		}
	}
	else {
		let def_val = { "max": 0, "min": Infinity, }[rna_func_type];
		let rna_func = Math[rna_func_type];
		for (let i = 0; i < m_ui32a.length; ++i) {
			m_ui32a[i] = def_val;
			for (let j = 0, p = i * rna_window_size; j < rna_window_size && p < w_ui32a.length; ++j, ++p) {
				if (w_ui32a[p] > 0) {
					m_ui32a[i] = rna_func(m_ui32a[i], w_ui32a[p]);
				}
			}
		}
	}
}

/**
 *
 * @param {"QM6a"|"CBS1-1"} refId
 * @param {string} sampleName
 * @param {number} start
 * @param {number} end
 */
function getRNAcoverage(refId, sampleName, start, end) {
	/** @type {Uint32Array} */
	const w_ui32a = w_rnaCoverage_ui32a_map[refId][sampleName];//m_rnaCoverage_ui32a_map[refId][sampleName];//

	return w_ui32a.slice(start, end);
}

/**
 * @param {string} refId
 * @param {string} sampleId
 * @param {string} geneId
 */
function getGeneRNAcoverage(refId = "QM6a", sampleId = "veg", geneId = "TRQ_010977") {
	const refIdx = get_genome_index(refId);

	const symbol = analysis_options.chrInfo_list[refIdx].symbol;

	const gene = gff_data_map[refId][symbol].find(a => a.type == "gene" && a.geneID == geneId);

	const sss = w_rnaCoverage_ui32a_map[refId][sampleId].slice(gene.start - 1, gene.end);

	// avg = sss.reduce((p, v) => p + v) / sss.length;

	return sss;
}

function enumGeneRNAavg(refId = "QM6a", geneId = "TRQ_010977") {
	return Object.keys(lncRNA_dataset).map(k => lncRNA_dataset[k].coverage ? k : null).filter(sampleId => sampleId).map(sampleId => {
		const sss = getGeneRNAcoverage(refId, sampleId, geneId);
		const avg = sss.reduce((p, v) => p + v) / sss.length;
		console.log(refId, sampleId, avg);
		return {
			sampleId: sampleId,
			avg: avg,
		};
	});
}

function normRNAbyGeneId(refId = "QM6a", geneId = "TRQ_010977") {
	const iter_pair = enumGeneRNAavg(refId, geneId);

	const max_avg = Math.max(...iter_pair.map(a => a.avg));

	iter_pair.forEach(({ sampleId, avg }) => {
		rnaCoverage_norm_map[refId][sampleId] = max_avg / avg;
		console.log(refId, sampleId, max_avg / avg);
	});
}

/**
 * @example
 * seq_list[0].forEach(pos => {
 *  ref_pos = pos_map[pos]);
 *  val = ui32a[ref_pos];
 * });
 *
 * @param {string} refId
 * @param {string} sampleName
 *
 * @param {Uint32Array} pos_map = pos_ref1_uint32array
 * @param {number} row_height
 * @param {number} bp_size
 * @param {number} view_length
 * @param {number} row_width
 * @param {(main_ctx: CanvasRenderingContext2D, value: number, max_value: number, ref_pos: number, percent: number, pos: number, array_cursor: number) => void} func_stroke_color
 * @returns {Promise<boolean>} return false if no data
 */
async function drawRNACoverage(refId, sampleName, pos_map, row_height, bp_size, view_length, row_width, func_stroke_color) {
	let scaleDown_mode;
	let scale_size;
	let delay_rendering_mode = false;

	if (rna_window_size > 0 &&
		view_length >= (seq_list[0].length / rna_window_size)
	) {
		scaleDown_mode = true;
		scale_size = rna_window_size;
	}
	else {
		scaleDown_mode = false;
		scale_size = 1;

		if (rna_window_size == 0 &&
			view_length >= (seq_list[0].length / rna_window_size)
		) {
			delay_rendering_mode = true;
			alert("rendering...");
		}
	}

	main_ctx.save();

	if (m_rnaCoverage_ui32a_map[refId] == null || w_rnaCoverage_ui32a_map[refId] == null) {
		return false;
	}

	/** @type {Uint32Array} */
	const w_ui32a = scaleDown_mode ? m_rnaCoverage_ui32a_map[refId][sampleName] : w_rnaCoverage_ui32a_map[refId][sampleName];
	if (!w_ui32a) {
		return false;
	}

	//const max_display_value = typeof window.max_display_value == "number" ? window.max_display_value : rna_display_log_base ** rna__exp;

	const raw_max_value = rnaCoverage_max_map[refId][sampleName];
	const norm = rnaCoverage_norm_map[refId][sampleName] ?? 1;

	let auto_scale_RNA = false;

	let globalAlpha;
	if (view_length >= 1000_000) {
		main_ctx.strokeStyle = "#000000";
		globalAlpha = 0x33 / 0xFF;
	}
	else if (view_length >= 500_000) {
		main_ctx.strokeStyle = "#000000";
		globalAlpha = 0x55 / 0xFF;
	}
	else if (view_length <= 50_000 && viewerState.auto_scale_RNA) {
		auto_scale_RNA = true;
		main_ctx.strokeStyle = "#008";
		globalAlpha = 1;
	}
	else {
		main_ctx.strokeStyle = "#000000";
		globalAlpha = 1;
	}

	let auto_max_value = 0;
	if (auto_scale_RNA) {
		for (let i = g_chrCoord.bp_start - 1; i < g_chrCoord.bp_end; i += scale_size) {
			const pos = Math.min(i, seq_list[0].length);
			if (seq_list[get_genome_index(refId)][pos] == "-") {
				continue;
			}

			const raw_bp_pos = pos_map[pos] - 1;
			const mappingPos = scaleDown_mode ? Math.trunc(pos_map[pos] / scale_size) : raw_bp_pos;

			const value = w_ui32a[mappingPos];

			auto_max_value = Math.max(auto_max_value, value);
		}
	}

	const max_value = auto_scale_RNA ? auto_max_value : raw_max_value;

	const low_cutoff_display_value = viewerState.rna_reads_low_cutoff_display_value ?? 0;
	const high_cutoff_display_value = viewerState.rna_reads_high_cutoff_display_value ?? max_value;

	// draw_RNA_coverage_border(max_value, row_width, row_height, "#04F"); // QM6a border
	draw_RNA_coverage_border(max_value, row_width, row_height, "#000");

	/** @type {(x:number) => number} */
	const rna_reads_func = get_rna_reads_func(max_value * norm);

	main_ctx.lineWidth = Math.min(Math.max(1, scaleDown_mode ? (bp_size * 100) : bp_size), 24);
	const half_lineWidth = main_ctx.lineWidth * 0.5;

	// # of RNA-seq reads mapped
	main_ctx.fillStyle = "#F00";

	let draw_value;
	if (main_ctx.lineWidth >= 10) {
		draw_value = true;

		main_ctx.font = Math.min(Math.trunc(main_ctx.lineWidth), 18) + "px Arial";
		main_ctx.textAlign = "right";
		main_ctx.textBaseline = "middle";
	}

	for (let i = g_chrCoord.bp_start - 1; i < g_chrCoord.bp_end; i += scale_size) {
		const pos = Math.min(i, seq_list[0].length);
		if (seq_list[get_genome_index(refId)][pos] == "-") {
			continue;
		}

		const raw_bp_pos = pos_map[pos] - 1;
		const mappingPos = scaleDown_mode ? Math.trunc(pos_map[pos] / scale_size) : raw_bp_pos;

		/** value: [0, max_display_value] */
		const value = w_ui32a[mappingPos] * norm;// / max_value;

		if (value > low_cutoff_display_value && value < high_cutoff_display_value) {
			const x = ((1 - g_chrCoord.bp_start) + pos) * bp_size + (bp_size * 0.5);
			if (x >= 0 && x <= row_width) {
				const dy = rna_reads_func(value);//value ** rna__exp // apply color in rna_reads_func

				if (dy > 0) {
					main_ctx.globalAlpha = globalAlpha;

					if (func_stroke_color) {
						func_stroke_color(main_ctx, value, max_value, raw_bp_pos, dy, pos, mappingPos);
					}

					main_ctx.beginPath();

					const top_y = row_height - Math.max(1, dy * row_height);
					main_ctx.moveTo(x, top_y);
					main_ctx.lineTo(x, row_height);
					main_ctx.stroke();

					if (viewerState.rna_reads_max_display_value && value > viewerState.rna_reads_max_display_value) {
						main_ctx.fillRect(x - half_lineWidth, top_y, main_ctx.lineWidth, row_height * 0.1);
					}

					main_ctx.globalAlpha = 1;

					if (draw_value && window.draw_RNA_value) {
						main_ctx.save();
						main_ctx.translate(Math.trunc(x) + 0.5, row_height * 0.95);
						main_ctx.rotate(Math.PI * 0.5);
						main_ctx.fillText((value).toFixed(0), 0, 0);
						main_ctx.restore();
					}
				}//dy > 0
			}
		}
	}

	// ctx.fillStyle = "rgba(32, 127, 255, 0.5)";
	// ctx.fill();

	main_ctx.restore();

	if (delay_rendering_mode) {
		await delayFrame();
	}

	return true;
}

/**
 * @returns {Uint32Array[]}
 * @see {@link multialign_to_chrPos_posMap}
 */
function make_pos_ref_map_list() {
	const pos_ref_map_list = [
		pos_ref1_uint32array,
		// pos_ref2_uint32array,
		// window.pos_subject_1_map,
		// window.pos_subject_2_map,
		// window.pos_subject_3_map,
		// window.pos_subject_4_map,
	];
	// is has ref2 ?
	if (dataset.parental_list.length > 1) {// pos_ref2_uint32array != null
		pos_ref_map_list.push(pos_ref2_uint32array);
	}
	if (seq_list.length > pos_ref_map_list.length) {
		seq_list.slice(pos_ref_map_list.length).forEach((_, idx) => {
			pos_ref_map_list.push(window[`pos_subject_${(idx + 1)}_map`]);
		});
	}
	return pos_ref_map_list;
}

/**
 * @returns {Uint32Array[]}
 */
function make_ref_pos_map_list() {
	const ref_pos_map_list = [
		ref1_pos_uint32array,
		// ref2_pos_uint32array,
		// window.subject_1_pos_map,
		// window.subject_2_pos_map,
		// window.subject_3_pos_map,
		// window.subject_4_pos_map,
	];
	if (ref2_pos_uint32array != null) {
		ref_pos_map_list.push(ref2_pos_uint32array);
	}
	if (seq_list.length > ref_pos_map_list.length) {
		seq_list.slice(ref_pos_map_list.length).forEach((_, idx) => {
			ref_pos_map_list.push(window[`subject_${(idx + 1)}_pos_map`]);
		});
	}
	return ref_pos_map_list;
}

/**
 * @param {number} nChr
 * @returns {Uint32Array[]}
 */
function _make_chr_pos_ref_map_list(nChr) {
	const list = Object.values(dataset.results[nChr - 1]);
	return list.map(seq => multialign_to_chrPos_posMap(seq));
}

/**
 * @param {number} nChr
 * @returns {Uint32Array[]}
 */
function _make_chr_ref_pos_map_list(nChr) {
	const list = Object.values(dataset.results[nChr - 1]);
	return list.map(seq => chrPos_to_multialign_posMap(seq));
}

/**
 * @param {number} seg_idx_start
 * @param {number} seg_idx_end
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} max_view_width
 * @param {number} bp_size
 * @param {number} bp_per_px_min1
 * @param {any} allow_draw_small_object
 * @param {any} draw_rdna_arrow
 * @param {any} options
 */
function drawProgney(seg_idx_start = 2, seg_idx_end = region_rect.length, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow, options) {
	const pos_map = make_pos_ref_map_list();

	for (let seg_id = seg_idx_start; seg_id < seg_idx_end; ++seg_id) {
		drawSegRow(main_ctx, 0, 0, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow ? options.rDNA_info : null);
		main_ctx.save();
		drawSegRow_cover(main_ctx, 0, 0, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, draw_rdna_arrow ? options.rDNA_info : null);
		main_ctx.restore();
		main_ctx.translate(0, seg_row_height + seg_row_separate);
	}
	if (options.rDNA_info && options.rDNA_info.chr == viewerState.nChr) {
		if (viewerState.display_rdna_border) {
			console.error("rdna_border rect=?");
			let y2 = seq_list.length * (seg_row_height + seg_row_separate);
			// let start = ((1 - g_chrCoord.bp_start) + (options.rDNA_info.region_start ?? options.rDNA_info.alignment_start)) * bp_size;
			// let end = ((1 - g_chrCoord.bp_start) + (options.rDNA_info.region_end ?? options.rDNA_info.alignment_end)) * bp_size;
			const [start, end] = [
				dataset.rDNA_info.region_start,
				dataset.rDNA_info.region_end,
			];
			let x1 = pos_t(Math.min(start, end));
			let x2 = pos_t(Math.max(start, end));
			main_ctx.strokeStyle = "#7F7F00";
			main_ctx.beginPath();
			main_ctx.rect(x1, 0, x2 - x1, y2);
			main_ctx.stroke();
		}
	}

	function pos_t(pos) {
		return ((1 - g_chrCoord.bp_start) + pos_map[seg_idx_start][pos - 1]) * bp_size;
	}
}

/**
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} pos_x
 * @param {number} pos_y
 * @param {number} seg_id
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} max_view_width
 * @param {number} bp_size
 * @param {number} bp_per_px_min1
 * @param {boolean} allow_draw_small_object
 * @param {*} rDNA_info
 */
function drawSegRow_s(ctx, pos_x, pos_y, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, rDNA_info) {
	const ss = region_rect[seg_id];

	ctx.globalAlpha = 1;
	ctx.filter = "none";

	for (let si = 0; si < ss.length; ++si) {
		const start = ((1 - g_chrCoord.bp_start) + ss[si].start) * bp_size;
		const end = ((1 - g_chrCoord.bp_start) + ss[si].end) * bp_size;
		const length = end - start;
		const color_id = ss[si].col;
		const x1 = Math.trunc(Math.max(0, start));
		const xlen = 1;// Math.min(end, max_view_width) - x1;

		if (color_id == null) {//no data
			continue;//transparent
		}

		if (!(start <= max_view_width)) {
			continue;
		}

		const colId = color_id & ColorID.mask;
		const col = args_colors[colId];
		if (col == "#FFFFFF") {
			continue;
		}

		ctx.beginPath();
		ctx.rect(x1, pos_y, xlen, seg_row_height);
		ctx.fillStyle = col;
		ctx.fill();
	}
}

/**
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} pos_x
 * @param {number} pos_y
 * @param {number} seg_id
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} max_view_width
 * @param {number} bp_size
 * @param {number} bp_per_px_min1
 * @param {boolean} allow_draw_small_object
 * @param {*} rDNA_info
 */
function drawSegRow(ctx, pos_x, pos_y, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, rDNA_info) {
	if (dataset.mode != "single" && !region_rect?.[0]?.length) {
		//viewerState.max_view_width = canvas.width;
		return _draw_seg_row_v8(ctx, seg_id);
	}

	if (viewerState.drawSegRow_func) {
		return viewerState.drawSegRow_func(ctx, pos_x, pos_y, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, rDNA_info);
	}

	let prev_22_colorId = [
		ColorID.dad,
		ColorID.mom
	];
	/**
	 * @param {number} pos_start g_chrCoord.bp_start
	 */
	function calc_prev_parentalSNP_colorId(pos_start) {
		if (prev_22_colorId.length != 2) {
			return;
		}
		for (let i = pos_start - 1; i >= 0; --i) {
			let ref1 = seq_list[0][i];
			let ref2 = seq_list[1][i];
			let a = seq_list[2][i];
			let b = seq_list[3][i];
			let c = seq_list[4][i];
			let d = seq_list[5][i];
			let ss = [a, b, c, d];

			let rc1 = ss.filter(s => s == ref1).length;
			let rc2 = ss.filter(s => s == ref2).length;

			// if (ref1 != ref2 && (rc1 + rc2) == 4) {//SNP 0:4 or 1:3 or 2:2 or 3:1 or 4:0
			if (ref1 != ref2 && rc1 == rc2) {//SNP 2:2 // 20210125
				prev_22_colorId[2] = (a == ref1 ? ColorID.dad : (a == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[3] = (b == ref1 ? ColorID.dad : (b == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[4] = (c == ref1 ? ColorID.dad : (c == ref2 ? ColorID.mom : ColorID.diff));
				prev_22_colorId[5] = (d == ref1 ? ColorID.dad : (d == ref2 ? ColorID.mom : ColorID.diff));
				return;
			}
		}
	}

	let ss = region_rect[seg_id];

	// init unknow state
	let is_first_identical_seg = true;

	let prev_col_id = null;
	for (let si = 0; si < ss.length; ++si) {
		let start = ((1 - g_chrCoord.bp_start) + ss[si].start) * bp_size;
		let end = ((1 - g_chrCoord.bp_start) + ss[si].end) * bp_size;
		let length = end - start;
		let color_id = ss[si].col;

		if (color_id == null) {//no data
			continue;//transparent
		}

		if (!(start <= max_view_width)) {
			continue;
		}

		function left_near_aling_snp() {
			let colId;
			for (let iii = si - 1; iii >= 0; --iii) {
				const ccc = ss[iii].col & ColorID.mask;
				if (ccc == ColorID.dad || ccc == ColorID.dad_rip) {
					colId = ColorID.dad;
					break;
				}
				else if (ccc == ColorID.mom || ccc == ColorID.mom_rip) {
					colId = ColorID.mom;
					break;
				}
			}
			return colId;
		}
		function next_near_aling_snp() {
			let colId;
			for (let iii = si + 1; iii < ss.length; ++iii) {
				const ccc = ss[iii].col & ColorID.mask;
				if (ccc == ColorID.dad || ccc == ColorID.dad_rip) {
					colId = ColorID.dad;
					break;
				}
				else if (ccc == ColorID.mom || ccc == ColorID.mom_rip) {
					colId = ColorID.mom;
					break;
				}
			}
			return colId;
		}

		if (allow_draw_small_object) {
			if (end >= 0 && length >= 1) {
				let x1 = Math.max(0, start);
				let xlen = Math.min(end, max_view_width) - x1;

				let colId = color_id & ColorID.mask;

				ctx.filter = "brightness(0.8) saturate(2)";

				// viewerState.fill_prev_snp22_after_ill &&
				if (is_first_identical_seg && colId == ColorID.identical) {
					// 20210126
					// calc_prev_parentalSNP_colorId(g_chrCoord.bp_start);
					// colId = prev_22_colorId[seg_id];

					colId = left_near_aling_snp() ?? colId;

					prev_col_id = colId;
					ctx.filter = "saturate(0.5) brightness(1.5)";

					is_first_identical_seg = false;// remove init state
				}
				else if ((color_id & ColorID.mask) == ColorID.identical) {
					if (prev_col_id != null) {
						colId = prev_col_id;
						ctx.filter = "saturate(0.5) brightness(1.5)";
					}
				}
				else  {
					prev_col_id = colId;
					if (viewerState.cast_indel_to_snp22) {
						let indel = color_id & ColorID.indel_mask;
						if (indel) {
							colId = prev_22_colorId[seg_id];
						}
						else {
							prev_col_id = colId;
						}
					}
				}

				let col = args_colors[colId];

				ctx.beginPath();
				ctx.rect(x1, pos_y, xlen, seg_row_height);
				ctx.fillStyle = col;
				ctx.fill();
				ctx.globalAlpha = 1;
				ctx.filter = "none";
				//++ct;
			}
		}
		else if (dataset.appearance == "diff") {// CNChuang_spore_1
			const x1 = Math.max(0, start);
			const i_x1 = Math.trunc(x1) + 0.5;
			const xlen = Math.min(end, max_view_width) - i_x1;
			const i_xlen = Math.max(1, xlen);
			const colId = color_id & ColorID.mask;
			const indel = color_id & ColorID.indel_bit;
			const identical = color_id & ColorID.identical;
			const col = args_colors[colId];
			if (identical) {
				draw_block(i_x1 + 0.25, pos_y, i_xlen - 0.5, seg_row_height, col);
			}
			else if (!indel) {
				draw_block(Math.trunc(i_x1), pos_y, i_xlen, seg_row_height, col);
			}
			else if (colId == ColorID.diff) {
				draw_block(Math.trunc(i_x1), pos_y, i_xlen, seg_row_height, col);
			}
			function draw_block(x, y, w, h, col) {
				ctx.beginPath();
				ctx.rect(x, y, w, h);


				ctx.fillStyle = col;
				// ctx.filter = "brightness(0.9) contrast(1.5) saturate(10) hue-rotate(60deg)";

				ctx.fill();
				// ctx.filter = "none";
			}
		}
		else if (length >= 1) {
			if (end >= 0) {
				if ((color_id & ColorID.indel_mask) == ColorID.diff) {
					debugger;
				}
				if (seg_id >= 2 && (color_id & ColorID.indel_mask)) {
					// nothing // large gap
					//prev_col_id = null;//reset to empty
					//skip

					// if ((ss[si].end - ss[si].start) <= 1) {
					if (1) {
						// const colId = next_near_aling_snp();
						let colId = color_id & ColorID.mask;

						// keep empty
						if (viewerState.cast_indel_to_snp22) {
							const col = args_colors[colId];
							ctx.fillStyle = col;
							ctx.beginPath();
							let x1 = Math.max(0, start);
							let xlen = Math.min(end + 1, max_view_width) - x1;
							ctx.rect(x1, pos_y, xlen, seg_row_height);
							ctx.fill();
						}

						// save state
						prev_col_id = colId;
					}

					// if (viewerState.cast_indel_to_snp22) {
					// 	const colId = prev_22_colorId[seg_id];
					// 	let col = args_colors[colId];
					// 	ctx.fillStyle = col;
					// 	ctx.beginPath();
					// 	let x1 = Math.max(0, start);
					// 	let xlen = Math.min(end + 1, max_view_width) - x1;
					// 	ctx.rect(x1, pos_y, xlen, seg_row_height);
					// 	ctx.fill();
					// }
				}
				else {
					let b_fill;

					let colId = color_id & ColorID.mask;
					if (is_first_identical_seg && colId == ColorID.identical) {
						// calc_prev_parentalSNP_colorId(ss[si].start);
						// colId = prev_22_colorId[seg_id];

						colId = left_near_aling_snp() ?? colId;

						let col = args_colors[colId];
						ctx.fillStyle = col;
						b_fill = true;
						is_first_identical_seg = false;
						prev_col_id = colId;
					}
					let col = args_colors[colId];

					let x1 = Math.max(0, start);
					let xlen = Math.min(end + 1, max_view_width) - x1;

					if ((color_id & ColorID.mask) == ColorID.identical) {
						if (prev_col_id != null &&
							prev_col_id != ColorID.none &&
							prev_col_id != ColorID.identical &&
							prev_col_id != ColorID.diff
						) {
							if (args_colors[prev_col_id] == current_colorset["dad"]) {
								ctx.fillStyle = current_colorset["dad_bk"];
								b_fill = true;
							}
							else if (args_colors[prev_col_id] == current_colorset["mom"]) {
								ctx.fillStyle = current_colorset["mom_bk"];
								b_fill = true;
							}
							if (seg_id == 0) {
								ctx.fillStyle = current_colorset["dad_bk"];
								b_fill = true;
							}
							else if (seg_id == 1) {
								ctx.fillStyle = current_colorset["mom_bk"];
								b_fill = true;
							}

							/**
							 * 20200803
							 * TODO: large identical region
							 */
							if (viewerState.fill_large_identical && seg_id >= 2 && b_fill) {
								// prev_col_id = null;//20200803
								ctx.fillStyle = "lime";
							}
						}
						else {
							colId = left_near_aling_snp() ?? colId;
							const col = args_colors[colId];
							ctx.fillStyle = col;
							b_fill = true;
						}
					}
					else {
						let indel = color_id & ColorID.indel_mask;
						if (!indel) {
							if (seg_id < 2) {
								ctx.fillStyle = col;
								b_fill = true;
							}
							else {
								if (viewerState.crossover_only) {
									if (col == current_colorset["dad"]) {
										ctx.fillStyle = current_colorset["dad_bk"];
										b_fill = true;
									}
									else if (col == current_colorset["mom"]) {
										ctx.fillStyle = current_colorset["mom_bk"];
										b_fill = true;
									}
									else {
										ctx.fillStyle = col;
										b_fill = true;
									}
								}
								else {
									ctx.fillStyle = col;
									b_fill = true;
								}
							}
						}
						if (colId == ColorID.dad || colId == ColorID.mom) {
							prev_col_id = color_id & ColorID.mask;
						}
					}

					if (b_fill) {
						ctx.beginPath();
						ctx.rect(x1, pos_y, xlen, seg_row_height);
						ctx.fill();
					}
					else {
						prev_col_id = null;
					}
				}
			}
		}
		else {
			let sj = si;
			let next_start = 0;
			let next_end = 0;
			let n_col = [
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0
			];
			let next_length;
			let merge_length;
			let has_rip = false;
			let merge_snp_cnt = 0;
			//
			//console.time("merge_seg");
			// let seg_seq = new Set();
			for (; sj < ss.length; ++sj) {
				// const _prev_next_start = next_start;
				const _prev_next_end = next_end;
				const _prev_merge_length = merge_length;

				next_start = ((1 - g_chrCoord.bp_start) + ss[sj].start) * bp_size;
				next_end = ((1 - g_chrCoord.bp_start) + ss[sj].end) * bp_size;

				// seg_seq.add(ss[sj].col);

				if (next_length >= (window.merge_near_size ?? 2)) {// || (merge_snp_cnt / next_length) < (15 / 1.5)
					break;
				}

				next_length = next_end - next_start;

				let ml = next_end - start;
				if (next_length >= 1) {// (window.seg_merge_max_length ?? 2)
					// next_start = _prev_next_start;
					next_end = _prev_next_end;
					merge_length = _prev_merge_length;
					break;
				}
				merge_length = ml;

				if (next_start > max_view_width) {
					next_end = _prev_next_end;
					merge_length = _prev_merge_length;
					break;
				}

				// n_col[ss[sj].col] = (n_col[ss[sj].col] || 0) + next_length * [1, 1, 1, bp_per_px_min1, bp_per_px_min1, 1][ss[sj].col & ColorID.mask];
				n_col[ss[sj].col] = (n_col[ss[sj].col] || 0) + next_length;

				// n_col[ss[sj].col] += next_length;

				has_rip = ss[sj].col == ColorID.dad_rip || ss[sj].col == ColorID.mom_rip;

				merge_snp_cnt += 1;

				if (merge_length >= 1) {//merge size in px
					break;
				}
			}

			if (next_end >= 0) {
				//skip no diff if rect_length <= 5
				if (merge_length <= 2) {
					Object.values(ColorID).forEach(k => {
						// n_col[k] *= (merge_near_priority_small[k] || 1);
						n_col[k] *= (merge_near_priority_small[k] ?? 1);
					});
					// n_col[ColorID.identical] = 0;
				}
				else {
					Object.keys(merge_near_priority).forEach(k => {
						n_col[k] *= merge_near_priority[k] ?? 1;
					});
				}

				let n_col_sorted = [...n_col].sort((a, b) => b - a);
				let max_col = n_col_sorted[0];

				if (max_col) {
					let color_id = n_col.findIndex(a => a == max_col);

					// // TODO: use last color color if count_of_SNP_1 == count_of_SNP_2
					// const scsc = n_col.map(a => a == max_col ? a : 0);
					// const _color_id = [...seg_seq.keys()].reverse().find(k => scsc[k]);
					// if (_color_id >= 0 && _color_id != color_id) {
					// 	color_id = color_id;
					// }
					// color_id = _color_id >= 0 ? _color_id : color_id;

					if (seg_id >= 2 && (color_id & ColorID.indel_mask)) {
						// nothing // large gap
						//prev_col_id = null;//reset to empty
						//skip

						// prev_col_id = null;//20200803

						if (1) {
							// const colId = next_near_aling_snp();
							let colId = color_id & ColorID.mask;
							const col = args_colors[colId];

							// keep empty
							if (viewerState.cast_indel_to_snp22) {
								ctx.fillStyle = col;
								ctx.beginPath();
								let x1 = Math.max(0, start);
								let xlen = Math.min(end + 1, max_view_width) - x1;
								ctx.rect(x1, pos_y, xlen, seg_row_height);
								ctx.fill();
							}

							// save state
							prev_col_id = colId;
						}
					}
					else {
						let b_fill;
						let colId = color_id & ColorID.mask;
						let indel = color_id & ColorID.indel_mask;
						if (is_first_identical_seg && colId == ColorID.identical) {
							if (20200803) {
								prev_col_id = ColorID.identical;
								ctx.fillStyle = current_colorset["identical"];
								b_fill = true;
							}
							else {
								// calc_prev_parentalSNP_colorId(ss[si].start);
								// colId = prev_22_colorId[seg_id];

								colId = left_near_aling_snp() ?? colId;

								is_first_identical_seg = false;
								prev_col_id = colId;
							}
						}
						let col = args_colors[colId];

						let x1 = Math.max(0, start);
						let xlen = Math.min(x1 + merge_length + 1, max_view_width) - x1;

						if ((color_id & ColorID.mask) == ColorID.identical) {
							if (prev_col_id != null) {
								if (args_colors[prev_col_id] == current_colorset["dad"]) {
									ctx.fillStyle = current_colorset["dad_bk"];
									b_fill = true;
								}
								else if (args_colors[prev_col_id] == current_colorset["mom"]) {
									ctx.fillStyle = current_colorset["mom_bk"];
									b_fill = true;
								}
								if (seg_id == 0) {
									ctx.fillStyle = current_colorset["dad_bk"];
									b_fill = true;
								}
								else if (seg_id == 1) {
									ctx.fillStyle = current_colorset["mom_bk"];
									b_fill = true;
								}
							}
						}
						else {
							if (viewerState.crossover_only) {
								if (col == current_colorset["dad"]) {
									ctx.fillStyle = current_colorset["dad_bk"];
									b_fill = true;
								}
								else if (col == current_colorset["mom"]) {
									ctx.fillStyle = current_colorset["mom_bk"];
									b_fill = true;
								}
								else {
									ctx.fillStyle = col;
									b_fill = true;
								}
							}
							else {
								ctx.fillStyle = col;
								b_fill = true;
							}
							if (colId == ColorID.dad || colId == ColorID.mom) {
								prev_col_id = color_id & ColorID.mask;
							}
						}
						if (b_fill) {
							// main_ctx.fillStyle = "white";
							// main_ctx.strokeStyle = "white";
							ctx.fillRect(x1, pos_y, xlen, seg_row_height);
							// if (colId != has_rip)  {
							// 	debugger;
							// }
						}
						else {
							prev_col_id = null;
						}
					}
				}
				else {
					//large gap
				}
			}
			si = sj - 1;
		}// draw data

		// if (start <= viewerState.pageX && viewerState.pageX <= end) {
		// 	ctx.fillStyle = "black";
		// 	ctx.fillText(si + "", start, pos_y);
		// }
	}
}

/**
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} pos_x
 * @param {number} pos_y
 * @param {number} seg_id
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {number} max_view_width
 * @param {number} bp_size
 * @param {number} bp_per_px_min1
 * @param {boolean} allow_draw_small_object
 * @param {*} rDNA_info
 */
function drawSegRow_cover(ctx, pos_x, pos_y, seg_id, seg_row_height, seg_row_separate, max_view_width, bp_size, bp_per_px_min1, allow_draw_small_object, rDNA_info) {
	const seg_row_hheight = Math.trunc(seg_row_height * 0.5);

	let font_size = Math.max(9, Math.min(bp_size / 4 * 3, 16));
	if (font_size == 9) {
		ctx.lineWidth = 0.5;
	}

	// bp text
	const font = Math.trunc(font_size) + "px Arial";
	ctx.font = font;

	if (allow_draw_small_object) {
		main_ctx.textAlign = "center";
		main_ctx.textBaseline = "middle";
		ctx.strokeStyle = "black";
		ctx.fillStyle = "black";

		for (let si = g_chrCoord.bp_start - 1; si < g_chrCoord.bp_end; ++si) {
			let x1 = ((1 - g_chrCoord.bp_start) + si) * bp_size;
			let x2 = ((1 - g_chrCoord.bp_start) + si + 1) * bp_size;

			ctx.beginPath();
			ctx.rect(x1, pos_y, bp_size, seg_row_height);
			ctx.stroke();

			//center text
			ctx.fillText(
				seq_list[seg_id][si],
				x2 - (bp_size * 0.5),
				pos_y + seg_row_hheight
				);
		}
	}

	// rDNA
	if (viewerState.display_rDNA && rDNA_info) {
		const cy = pos_y + seg_row_hheight;

		const ali_rep = rDNA_info.data[seg_id].alignment_repeats;
		const min_raw_start = Math.min(...ali_rep[0]);
		ali_rep.forEach((ar_data, rep_idx) => {
			let raw_start = Math.min(...ar_data);
			let raw_len = Math.abs(ar_data[1] - ar_data[0]);

			let strand = ar_data[1] < ar_data[0] ? -1 : 1;

			let x1 = ((1 - g_chrCoord.bp_start) + raw_start) * bp_size;
			let x2 = ((1 - g_chrCoord.bp_start) + raw_start + raw_len - 1) * bp_size;
			let len = x2 - x1;
			let hlen = len * 0.5;
			let ahlen = len * 0.3;
			let cx = (x1 + x2) * 0.5;

			if (!(x1 <= max_view_width && x2 >= 0)) {
				return;
			}

			const qy = seg_row_hheight * 0.5;

			ctx.save();
			{
				ctx.translate(cx, cy);
				ctx.scale(strand, 1);

				//arrow
				ctx.beginPath();
				ctx.moveTo(0 - hlen, 0 - qy);
				ctx.lineTo(0 + ahlen, 0 - qy);
				ctx.lineTo(0 + ahlen, 0 - seg_row_hheight);
				ctx.lineTo(0 + hlen, 0);
				ctx.lineTo(0 + ahlen, 0 + seg_row_hheight);
				ctx.lineTo(0 + ahlen, 0 + qy);
				ctx.lineTo(0 - hlen, 0 + qy);
				ctx.closePath();

				ctx.globalAlpha = Math.min(bp_size >= 1 ? (1.5 - Math.min(bp_size * 2, 32) / 32) : 1, 1);
				ctx.fillStyle = current_colorset.rDNA;
				ctx.fill();
				ctx.globalAlpha = 1;
			}
			ctx.restore();

			if (viewerState.display_rdna_border) {
				ctx.strokeStyle = "#7F7F00";
				ctx.beginPath();
				ctx.rect(x1, pos_y + seg_row_hheight * 0.5, len, seg_row_hheight);
				ctx.stroke();
			}
		});
	}

	if (viewerState.display_gff && gff_data_map) {
		// _drawTetradMode_annotation(Object.keys(gff)[seg_id], seg_id, 0, seg_row_height, bp_size, max_view_width);
		_drawTetradMode_annotation(dataset.genomeNameList[seg_id], seg_id, 0, seg_row_height, bp_size, max_view_width);
	}
}

/**
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 */
function drawMarkers(bp_size, max_view_width, seg_row_height, seg_row_separate) {
	drawUserInputMarker(main_ctx.getTransform().f, bp_size, seg_row_separate, true);

	drawStrainMarker(bp_size, seg_row_separate);

	drawCrossoverMarker(bp_size, seg_row_height, seg_row_separate, align_start_index, align_end_index);

	// no display 2:2 SNV
	// for (let marker_idx = 1; marker_idx < allMarker.list.length; ++marker_idx) {
	for (let marker_idx = 0; marker_idx < allMarker.list.length; ++marker_idx) {
		const markersData = allMarker.list[marker_idx];
		const marker_name = markersData.property;

		if (marker_name && viewerState[`display_${marker_name}`]) {
			// if (viewerState.D8_bg_RIP_marker) {
			// 	continue;
			// }
			// else {
				drawMarkerRow(markersData, bp_size, max_view_width, seg_row_height, seg_row_separate);
			// }
		}
	}
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_separate
 */
function drawStrainMarker(bp_size, seg_row_separate) {
	if (viewerState.nChr == 0) {
		drawStrainMarker_linked_chr(bp_size, seg_row_separate);
	}
	else {
		drawStrainMarker_single_chr(bp_size, seg_row_separate);
	}
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_separate
 */
function drawStrainMarker_linked_chr(bp_size, seg_row_separate) {
	/** @type {number[]} */
	const all = dataset.genome_info_list[0].chr_list.map(a => dataset.results[a.index - 1][a.chr].length);//.reduce((aa, v) => aa + v, 0);

	const position_markers = [];

	position_markers.push(0);
	all.reduce((aa, v) => {
		position_markers.push(aa);
		return aa + v;
	});

	_drawUserInputMarker(position_markers, [], main_ctx.getTransform().f, bp_size, seg_row_separate, false);
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_separate
 */
function drawStrainMarker_single_chr(bp_size, seg_row_separate) {
	const range_list = [];
	try {
		if (viewerState.display_telomere) {
			const ref_pos_map_list = make_ref_pos_map_list();

			let left = 0;
			let right = seq_list[0].length;

			const nChr = viewerState.nChr;

			if (dataset.all_telomere?.[0]?.[nChr]?.length) {
				dataset.all_telomere.forEach((telomere, genome_idx) => {
					if (dataset.mode == "single" && genome_idx) {
						return;
					}
					else if (ref_pos_map_list[genome_idx] != null) {
						extends_telomere(telomere[nChr], genome_idx);
					}
				});
			}
			// if (dataset.telomere?.[nChr]) {
			// 	extends_telomere(dataset.telomere[nChr], 0);
			// }
			range_list.push([1, left]);
			range_list.push([right, seq_list[0].length]);

			function extends_telomere(telomere, genome_idx) {
				if (!(telomere && telomere[0] && telomere[1])) {
					return false;
				}
				try {
					if (telomere[0][0] != telomere[0][1]) {
						const l_telo_end = ref_pos_map_list[genome_idx][telomere[0][1] - 1 + 1];
						// range_list.push([1, l_telo_end]);
						left = Math.max(left, l_telo_end);
					}
					if (telomere[1][0] != telomere[1][1]) {
						const r_telo_start = ref_pos_map_list[genome_idx][telomere[1][0] - 1];
						// range_list.push([r_telo_start, seq_list[genome_idx].length]);
						right = Math.min(right, r_telo_start);
					}
				}
				catch (ex) {
					console.error(ex);
				}
			}
		}

		if (dataset.rDNA && dataset.rDNA.region && (
			viewerState.nChr == dataset.rDNA.nChr ||
			viewerState.nChr == dataset.rDNA.chr
		)) {
			const ref1 = dataset.rDNA.region.map(p => ref1_pos_uint32array[p - 1]);
			const ref2 = dataset.rDNA.region_ref2.map(p => ref2_pos_uint32array[p - 1]);

			range_list.push([
				Math.min(ref1[0], ref2[0]),
				Math.max(ref1[1], ref2[1]),
			]);

			// range_list.push([dataset.rDNA_info.region_start, dataset.rDNA_info.region_end]);
		}
	}
	catch (ex) {
		console.error(ex);
	}
	_drawUserInputMarker([], range_list, main_ctx.getTransform().f, bp_size, seg_row_separate, false);
}

/**
 * @param {any} markersData
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 * @param {boolean} [padding_bottom] default: true
 */
function drawMarkerRow(markersData, bp_size, max_view_width, seg_row_height, seg_row_separate, padding_bottom = true) {
	let color = current_colorset[markersData.property];
	const markers = markersData.values;

	_draw_marker_row(markers, bp_size, max_view_width, seg_row_height, color, { ref_in_out: false });
	if (viewerState.display_half_color) {
		_draw_marker_row(markers, bp_size, max_view_width, seg_row_height, color, { ref_in_out: true });
	}

	if (padding_bottom) {
		main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}
}

/**
 *
 * @param {any} markers
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} seg_row_height
 * @param {string} color
 * @param {{ ref_in_out:boolean; }} draw_flags
 */
function _draw_marker_row(markers, bp_size, max_view_width, seg_row_height, color, draw_flags) {
	let marker_h_width = (viewerState.marker_width ?? 2) / 2;

	for (let i = 0; i < markers.length; ++i) {
		const marker = markers[i];
		if (marker.hide) {
			continue;
		}
		let { pos, value, gc_in_out } = marker;
		let cx = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
		//let cx = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
		let x1 = cx - marker_h_width;
		let _x2 = cx + marker_h_width;
		if (!(_x2 >= 0 && x1 <= max_view_width)) {
			continue;
		}
		let width = Math.min(x1 + Math.max(bp_size, _x2 - x1), max_view_width) - x1;

		if (draw_flags.ref_in_out) {
			if (gc_in_out & 0x000000F0) {
				main_ctx.fillStyle = "greenyellow";
				main_ctx.fillRect(x1, 0, width, seg_row_height);
			}
			if (gc_in_out & 0x0F000000) {
				main_ctx.fillStyle = "hotpink";
				main_ctx.fillRect(x1, 0, width, seg_row_height * 0.5);
			}
			else if (gc_in_out & 0xF0000000) {
				main_ctx.fillStyle = "cyan";
				main_ctx.fillRect(x1, 0, width, seg_row_height * 0.5);
			}
		}
		else {
			// if (marker.mut_ref == 1) {// ref1 QM6a
			// 	main_ctx.fillStyle = current_colorset.dad_bk;
			// 	main_ctx.fillRect(x1, 0, width, seg_row_height * 0.5);
			// 	main_ctx.fillStyle = marker.color ?? color;
			// 	main_ctx.fillRect(x1, seg_row_height * 0.5, width, seg_row_height * 0.5);
			// }
			// else if (marker.mut_ref == 2) {// ref2 CBS1-1
			// 	main_ctx.fillStyle = current_colorset.mom_bk;
			// 	main_ctx.fillRect(x1, 0, width, seg_row_height * 0.5);
			// 	main_ctx.fillStyle = marker.color ?? color;
			// 	main_ctx.fillRect(x1, seg_row_height * 0.5, width, seg_row_height * 0.5);
			// }
			// else {
			// main_ctx.fillStyle = marker.color ?? color;
			// main_ctx.fillRect(x1, 0, width, seg_row_height);
			// }
			if (viewerState.colored_RIP_in_gene) {
				main_ctx.fillStyle = marker.color || color;
			}
			else {
				main_ctx.fillStyle = color;
			}
			main_ctx.fillRect(x1, 0, width, seg_row_height);
			if (viewerState.display_half_color) {
				if (marker.mut_ref == 1) { // ref1 QM6a
					main_ctx.fillStyle = current_colorset.dad_bk;
					main_ctx.fillRect(x1, seg_row_height * 0.5, width, seg_row_height * 0.5);
				}
				else if (marker.mut_ref == 2) { // ref2 CBS1-1
					main_ctx.fillStyle = current_colorset.mom_bk;
					main_ctx.fillRect(x1, seg_row_height * 0.5, width, seg_row_height * 0.5);
				}
			}
		}
	}
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_height
 * @param {number} align_start_index
 * @param {number} align_end_index
 */
function drawCrossoverMarker(bp_size, seg_row_height, seg_row_separate, align_start_index, align_end_index) {
	// crossover CO
	if (viewerState.display_co && dataset.crossover_list) {
		const list = dataset.crossover_list[viewerState.nChr - 1].filter(co_data => co_data.type == "CO");
		_drawCrossover(list, bp_size, seg_row_height, align_start_index, align_end_index, "#EE0000F0");
		main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}

	// crossover CO(NCO)
	if (viewerState.display_conco && dataset.crossover_list) {
		const list = dataset.crossover_list[viewerState.nChr - 1].filter(co_data => co_data.type == "CO(NCO)");
		_drawCrossover(list, bp_size, seg_row_height, align_start_index, align_end_index, "#EE7700F0");
		main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}

	// non_crossover_1 NCO1
	if (viewerState.display_nonCo && dataset.non_crossover_list) {
		const list = dataset.non_crossover_list[viewerState.nChr - 1].filter(co_data => co_data.nco_type == "NCO1");
		_drawCrossover(list, bp_size, seg_row_height, align_start_index, align_end_index, "#0000EEF0");
		main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}

	// non_crossover_1 NCO2
	if (viewerState.display_nonCo2 && dataset.non_crossover_list) {
		const list = dataset.non_crossover_list[viewerState.nChr - 1].filter(co_data => co_data.nco_type == "NCO2");
		_drawCrossover(list, bp_size, seg_row_height, align_start_index, align_end_index, "#0000EEF0");
		main_ctx.translate(0, 1 * (seg_row_height + seg_row_separate));
	}
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_height
 * @param {number} align_start_index
 * @param {number} align_end_index
 * @param {string} co_color
 */
function _drawCrossover(list, bp_size, seg_row_height, align_start_index, align_end_index, co_color) {
	main_ctx.font = "15px Arial";

	list.forEach(_co_data => {
		/** @type {CO_pos} */
		const co_data = (_co_data);
		if (co_data.hide) {
			return;
		}
		if (co_data.chr == viewerState.nChr) {
			// let x1 = ((1 - g_chrCoord.bp_start) + ref1_pos_uint32array[co_data.snp_start_out - 1]) * bp_size;
			// let x2 = ((1 - g_chrCoord.bp_start) + ref1_pos_uint32array[co_data.snp_start_in - 1]) * bp_size;
			// let x3 = ((1 - g_chrCoord.bp_start) + ref1_pos_uint32array[co_data.snp_end_in - 1]) * bp_size;
			// let x4 = ((1 - g_chrCoord.bp_start) + ref1_pos_uint32array[co_data.snp_end_out - 1]) * bp_size - 1;//pos[2:2 marker] - 1bp

			const x1 = ((1 - g_chrCoord.bp_start) + co_data.$snp_start_out) * bp_size;
			const x2 = ((1 - g_chrCoord.bp_start) + co_data.$snp_start_in) * bp_size;
			const x3 = ((1 - g_chrCoord.bp_start) + co_data.$snp_end_in + 1) * bp_size;
			const x4 = ((1 - g_chrCoord.bp_start) + co_data.$snp_end_out + 1) * bp_size - 1;//pos[2:2 marker] - 1bp

			const CO_marker_width = window.CO_marker_width ?? 4;

			let start;
			let end;
			if (co_data.$snp_start_in == co_data.$snp_end_in) {
				if (bp_size >= CO_marker_width) {
					start = x2;
					end = x2 + bp_size;
				}
				else {
					start = x2 - CO_marker_width;
					end = x2 + CO_marker_width;
				}
			}
			else {
				if ((x4 - x2) >= CO_marker_width) {
					start = x2;
					end = x4;
				}
				else {
					const c = (x2 + x3) * 0.5
					start = c - CO_marker_width;
					end = c + CO_marker_width;
				}
			}

			const len = end - start;

			if (end >= 0 && start <= viewerState.max_view_width) {
				// let in_length = Math.max(CO_marker_width, x3 - x2);
				// let out_length = Math.max(CO_marker_width, x4 - x1);

				if (co_data.why_remove) {
					if (viewerState.removed_color) {
						// main_ctx.fillStyle = "#0000FF4F";
						// main_ctx.fillRect(x1, 0, out_length, seg_row_height);

						// main_ctx.fillStyle = viewerState.removed_color;
						// main_ctx.fillRect(x2, 0, inout_length, seg_row_height);
						main_ctx.strokeStyle = viewerState.removed_color;
						main_ctx.strokeRect(start, 0, len, seg_row_height);

						main_ctx.fillStyle = "#000000";
						main_ctx.fillText(co_data.why_remove, start, 0);

						main_ctx.textAlign = "end";
						main_ctx.fillText(co_data.before_ref1_snp.join(","), start, 15);
						main_ctx.fillText(co_data.before_ref2_snp.join(","), start, 30);
						main_ctx.textAlign = "start";
						main_ctx.fillText(co_data.after_ref1_snp.join(","), end, 15);
						main_ctx.fillText(co_data.after_ref2_snp.join(","), end, 30);
					}
				}
				else {
					// main_ctx.fillStyle = "#FF00004F";
					// main_ctx.fillRect(x1, 0, out_length, seg_row_height);

					// main_ctx.beginPath();
					// main_ctx.rect(start, 0, len, seg_row_height);

					const ext = Math.min(2, len * 0.01);

					main_ctx.fillStyle = co_color;
					main_ctx.fillRect(start + ext * 0.5, 0, len - ext, seg_row_height);

					main_ctx.strokeStyle = co_color;
					main_ctx.strokeRect(start + ext * 0.5, 0, len - ext, seg_row_height);

					if (viewerState.removed_color) {
						main_ctx.textAlign = "end";
						main_ctx.fillText(co_data.before_ref1_snp.join(","), start, 15);
						main_ctx.fillText(co_data.before_ref2_snp.join(","), start, 30);
						main_ctx.textAlign = "start";
						main_ctx.fillText(co_data.after_ref1_snp.join(","), end, 15);
						main_ctx.fillText(co_data.after_ref2_snp.join(","), end, 30);
					}
				}
			}
		}
		else {
			console.log(co_data);
		}
	});

	mask_no_seq_region(bp_size, seg_row_height, align_start_index, align_end_index);
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_height
 * @param {number} align_start_index
 * @param {number} align_end_index
 */
function mask_no_seq_region(bp_size, seg_row_height, align_start_index, align_end_index) {
	// dummy
}

window.CO_marker_width = 1;

/**
 * @param {number} bp_size
 * @param {number} seg_row_height
 * @param {number} align_start_index
 * @param {number} align_end_index
 */
function _old_drawNonCrossover(bp_size, seg_row_height, align_start_index, align_end_index) {
	// main_ctx.save();
	main_ctx.font = "15px Arial";

	dataset.non_crossover_list[viewerState.nChr - 1].forEach(_nco_data => {
		/** @type {CO_pos} */
		const nco_data = (_nco_data);
		if (window.nco_filter && !window.nco_filter(nco_data)) {
			return;
		}
		if (nco_data.chr == viewerState.nChr) {
			if (nco_data.is_rip && viewerState.mask_NCO_is_RIP) {
				return;// continue;
			}
			let x1 = ((1 - g_chrCoord.bp_start) + nco_data.$snp_start_out) * bp_size;
			let x2 = ((1 - g_chrCoord.bp_start) + nco_data.$snp_start_in) * bp_size;
			let x3 = ((1 - g_chrCoord.bp_start) + nco_data.$snp_end_in) * bp_size;
			let x4 = ((1 - g_chrCoord.bp_start) + nco_data.$snp_end_out) * bp_size;//pos[2:2 marker] - 1bp

			const CO_marker_width = window.CO_marker_width ?? 1;

			let start;
			let end;
			if (nco_data.$snp_start_in == nco_data.$snp_end_in) {
				if (bp_size >= CO_marker_width) {
					start = x2;
					end = x2 + bp_size;
				}
				else {
					start = x2 - 2;
					end = x2 + 2;
				}
			}
			else {
				if ((x4 - x2) >= CO_marker_width) {
					start = x2;
					end = x4;
				}
				else {
					const c = (x2 + x3) * 0.5
					start = c - 2;
					end = c + 2;
				}
			}

			const len = end - start;

			if (end >= 0 && start <= viewerState.max_view_width) {
				if (nco_data.why_remove) {
					if (viewerState.removed_color) {
						main_ctx.strokeStyle = viewerState.removed_color;
						main_ctx.strokeRect(start, 0, len, seg_row_height);
						main_ctx.fillStyle = "#000000";
						main_ctx.fillText(nco_data.why_remove, start, 0);
					}
				}
				else {
					main_ctx.fillStyle = "#FF0000";
					main_ctx.fillRect(start, 0, len, seg_row_height);
				}
			}
		}
		else {
			console.log(nco_data);
		}
	});

	mask_no_seq_region(bp_size, seg_row_height, align_start_index, align_end_index);

	// main_ctx.restore();
}

/**
 * merge by subject
 * @param {({ sstart: number; send: number; score: number; } | BlastnCoord)[]} list
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

/**
 * @param {BlastnCoord[]} list
 */
function merge_overlap_alignment_bb(list) {
	/** @type {{ start: number; end: number; align: number; }[]} */
	const region_list = [];

	const tolerance = 1;

	list.filter(aln => {
		const [start, end] = [aln.sstart, aln.send].sort((a, b) => a - b);

		let region = check_range(start, end, aln.align);

		while (region) {
			const idx = region_list.findIndex(a => a == region);
			region_list.splice(idx, 1);// remove

			region = check_range(region.start, region.end, region.align);
		}
	});

	/**
	 * @param {number} start
	 * @param {number} end
	 * @param {number} align_len
	 */
	 function check_range(start, end, align_len) {
		const region = region_list.find(region => {
			return (
				(start - tolerance) <= region.end &&
				(end + tolerance) >= region.start
			);
		});
		if (region) {
			region.start = Math.min(region.start, start);
			region.end = Math.max(region.end, end);
			region.align = region.align + align_len;
			return region;
		}
		else {
			region_list.push({
				start: start,
				end: end,
				align: align_len,
			});
		}
	}

	return region_list.sort((a, b) => b.align - a.align);
}

/**
 * @param {number} bp_size
 * @param {number} seg_row_height
 * @param {number} seg_row_separate
 */
async function drawAfterProgney(bp_size, max_view_width, seg_row_height, seg_row_separate) {
	for (let func of drawAfterProgney.list) {
		await func(bp_size, max_view_width, seg_row_height, seg_row_separate);
	}
}
drawAfterProgney.list = [
];

/**
 * @param {number} pos_y - top
 * @param {number} bp_size
 * @param {number} seg_row_separate
 * @param {boolean} draw_background
 */
function drawUserInputMarker(pos_y, bp_size, seg_row_separate, draw_background = true) {
	_drawUserInputMarker(viewerState._position_markers, viewerState._range_markers, pos_y, bp_size, seg_row_separate, draw_background = true);
}

_drawUserInputMarker.lineWidth = 3;

/**
 * @param {number[]} position_markers
 * @param {number[][]} range_markers
 * @param {number} pos_y - top
 * @param {number} bp_size
 * @param {number} seg_row_separate
 * @param {boolean} draw_background
 */
function _drawUserInputMarker(position_markers, range_markers, pos_y, bp_size, seg_row_separate, draw_background = true) {
	main_ctx.save()

	const lineWidth = _drawUserInputMarker?.lineWidth ?? 0.5;
	const h_lineWidth = lineWidth * 0.5;
	main_ctx.lineWidth = lineWidth;

	const top = pos_y;
	for (let i = 0; i < range_markers.length; ++i) {
		const { 0: start, 1: end } = range_markers[i];
		const x1 = Math.floor(((1 - g_chrCoord.bp_start) + (start - 1)) * bp_size + h_lineWidth);
		const x2 = Math.ceil(((1 - g_chrCoord.bp_start) + (end - 0)) * bp_size - h_lineWidth);// point in range
		if (x2 >= 0 && x1 <= viewerState.max_view_width) {
			const w = Math.max(1, x2 - x1);
			const y1 = Math.floor(-top + h_lineWidth);
			const h = top - h_lineWidth - 0.5 - seg_row_separate;
			main_ctx.beginPath();
			main_ctx.rect(x1, y1, w, h);
			main_ctx.strokeStyle = color_user_marker_stroke ?? "#000000DD";
			main_ctx.stroke();

			if (draw_background) {
				main_ctx.beginPath();
				main_ctx.fillStyle = color_user_marker_bk ?? "#00FF0080";
				main_ctx.fillRect(
					(x1 + 0.5),
					(-top + 0.5) + (top - seg_row_separate),
					Math.max(1, x2 - x1),
					main_ctx.canvas.height
				);
			}
		}
	}
	for (let i = 0; i < position_markers.length; ++i) {
		let pos = position_markers[i] - 1;
		let x1 = ((1 - g_chrCoord.bp_start) + pos) * bp_size;
		if (x1 >= 0 && x1 <= viewerState.max_view_width) {
			main_ctx.beginPath();
			main_ctx.rect(x1 + 0.5, -top + 0.5, Math.max(1, bp_size), top - seg_row_separate);
			main_ctx.strokeStyle = color_user_marker_stroke ?? "#000000DD";
			main_ctx.stroke();

			if (draw_background) {
				main_ctx.beginPath();
				main_ctx.fillStyle = color_user_marker_bk ?? "#00FF0080";
				main_ctx.fillRect(
					(x1 + 0.5),
					(-top + 0.5) + (top - seg_row_separate),
					Math.max(1, bp_size),
					main_ctx.canvas.height
				);
			}
		}
	}

	main_ctx.restore();
}

/**
 * draw gc content
 * @param {number} seg_id
 * @param {string} refName
 * @param {*} ref_mapto_ma
 * @param {*} region_rect
 * @param {number} bp_size
 * @param {number} max_view_width
 * @param {number} gc_max_height
 * @param {number} seg_row_separate
 */
function draw_GC_row(seg_id, refName, ref_mapto_ma, region_rect, bp_size, max_view_width, gc_max_height, seg_row_separate) {
	if (gc_content == null || gc_content[refName] == null) {
		return false;
	}
	const gc_list = gc_content[refName][viewerState.nChr];
	if (!gc_list || gc_list.length <= 0) {
		return false;
	}

	const avg_top = 0;
	const avg_mid = Math.trunc(gc_max_height - gc_max_height * 0.5);
	const avg_bottom = gc_max_height;

	const draw_grid_hr = true;
	if (draw_grid_hr) {
		const step = 10;
		const d = gc_max_height / step;

		main_ctx.beginPath();

		const first_line = 1;// no first line
		const last_line = step - 1;// no top line
		for (let i = first_line; i <= last_line; ++i) {
			const y = Math.trunc(avg_top + d * i);
			main_ctx.moveTo(0, y);
			main_ctx.lineTo(max_view_width, y);
		}
	
		main_ctx.strokeStyle = "#DDDDDD";
		//ctx.setLineDash([2, 1]);
		main_ctx.stroke();
		//ctx.setLineDash([0, 0]);//reset
	}
	{
		// mid line
		main_ctx.beginPath();
		main_ctx.moveTo(0, avg_mid);
		main_ctx.lineTo(max_view_width, avg_mid);

		main_ctx.strokeStyle = "darkgray";
		//ctx.setLineDash([2, 1]);
		main_ctx.stroke();
		//ctx.setLineDash([0, 0]);//reset
	}

	draw_percent(main_ctx, max_view_width, avg_top, avg_mid, avg_bottom);

	if (viewerState.gc_content_clip_indel && region_rect?.[0]?.length) {
		main_ctx.save();//begin clip
		//
		main_ctx.beginPath();
		{//remove indel
			main_ctx.fillStyle = "#FFFFFF";
			let ss = region_rect[seg_id];
			for (let si = 0; si < ss.length; ++si) {
				let color_id = ss[si].col;
				if (color_id & ColorID.indel_mask) {
				}
				else {
					let start = ((1 - g_chrCoord.bp_start) + ss[si].start) * bp_size;
					let end = ((1 - g_chrCoord.bp_start) + ss[si].end) * bp_size;
					if (end >= 0 && start <= max_view_width) {
						let row_y_bottom = 1 * (gc_max_height + seg_row_separate);
						let row_y_top = row_y_bottom - gc_max_height;
						main_ctx.rect(start, row_y_top, end - start, gc_max_height);
						//ctx.fill();
					}
				}
			}
		}//gc remove indel
		main_ctx.clip();
	}

	let $seq_len_no_match_chrLen = false;

	let last_x = max_view_width;
	main_ctx.beginPath();
	main_ctx.moveTo(0, gc_max_height);
	for (let i = 0; i < gc_list.length; ++i) {
		let data = gc_list[i];
		let start = ref_mapto_ma[data.start];
		let end = ref_mapto_ma[data.end];
		if (start == null || end == null) {
			break;
		}
		if (end == 0) {
			debugger
			$seq_len_no_match_chrLen = true;
			break;
		}
		let x1 = ((1 - g_chrCoord.bp_start) + start) * bp_size;
		let x2 = ((1 - g_chrCoord.bp_start) + end + 1) * bp_size;
		if (x2 >= 0 && x1 <= max_view_width) {
			x1 = Math.max(0, x1);
			x2 = Math.min(x2, max_view_width);
			let w = x2 - x1;

			//let cx = (end + start) >> 1;//=> Math.trunc((end + start) / 2)
			let gc = (data.gc || 0) * gc_max_height / 100;
			let y = Math.trunc(avg_bottom - gc);

			main_ctx.lineTo(x1, y);
			main_ctx.lineTo(x2, y);

			last_x = x2;
		}
	}
	if ($seq_len_no_match_chrLen) {
		document.body.style.background = "red";
		debugger;
		// $seq_len_no_match_chrLen = analysis_options.chrInfo_list.map((v, i) => v.length == seq_list[i].replace(/-/g, "").length)
	}
	main_ctx.lineTo(last_x, gc_max_height);

    main_ctx.shadowBlur = 1;
	main_ctx.shadowOffsetY = -1;// remove solid-gray
	main_ctx.shadowColor = "rgba(255, 255, 255, 1)";
	main_ctx.strokeStyle = "black";
	main_ctx.stroke();
	
    main_ctx.shadowBlur = 0;
	main_ctx.shadowOffsetY = 0;
	main_ctx.shadowColor = "rgba(0, 0, 0, 0)";
	main_ctx.fillStyle = viewerState.gc_gap_color || "#FFFFCC";//"#FFFF0033"
	main_ctx.fill();

	if (viewerState.gc_content_clip_indel && region_rect?.[0]?.length) {
		main_ctx.restore();//end clip => clear clip
	}

	main_ctx.textAlign = "left";//reset
	main_ctx.textBaseline = "alphabetic";//reset

	if (!viewerState.gc_content_clip_indel) {
		main_ctx.fillStyle = "#FFFFFF22";//cover white
		let ss = region_rect[seg_id];
		for (let si = 0; si < ss.length; ++si) {
			let start = ((1 - g_chrCoord.bp_start) + ss[si].start) * bp_size;
			let end = ((1 - g_chrCoord.bp_start) + ss[si].end) * bp_size;
			if (end >= 0 && start <= max_view_width) {
				let color_id = ss[si].col;
				if (color_id & ColorID.indel_mask) {
					let row_y_bottom = (seg_id + 1) * (gc_max_height + seg_row_separate);
					let row_y_top = row_y_bottom - gc_max_height;
					main_ctx.beginPath();
					main_ctx.rect(start, row_y_top, end - start, gc_max_height);
					main_ctx.fill();
				}
				else {
				}
			}
		}
	}

	//TDOD: GC% plot border rectangle
	main_ctx.strokeStyle = "black";
	main_ctx.beginPath();

	main_ctx.rect(0, 0 - viewerState.row_padding_top, max_view_width, gc_max_height + viewerState.row_padding_bottom);

	//ctx.moveTo(0, sid * height);
	//ctx.lineTo(max_view_width, sid * height);

	main_ctx.stroke();

	return true;
}

/**
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} max_view_width
 * @param {number} avg_top pos y
 * @param {number} avg_mid pos y
 * @param {number} avg_bottom pos y
 */
function draw_percent(ctx, max_view_width, avg_top, avg_mid, avg_bottom) {
	ctx.fillStyle = "black";
	ctx.textAlign = "left";
	ctx.font = viewerState.global_digi_font_style; //"16pt Arial";//Math.trunc(gc_max_height / 3 / 1.5) + "px Arial";
	{
		const p100 = `${(100).toFixed(0)}`;
		const p50 = `${(50).toFixed(0)}`;
		const p0 = `${(0).toFixed(0)}`;

		ctx.textBaseline = "hanging";
		ctx.fillText(p100, max_view_width + 5, avg_top); //100%
		const p100_width = ctx.measureText(p100).width;

		ctx.textBaseline = "middle";
		ctx.fillText(p50, max_view_width + 5, avg_mid); //50%

		ctx.textBaseline = "ideographic"; // alphabetic
		ctx.fillText(p0, max_view_width + 5, avg_bottom - 1); //0%

		ctx.textBaseline = "middle";
		ctx.fillText(" GC%", max_view_width + 5 + p100_width, avg_mid);
	}
}

/**
 * @param {string} file_name
 * @param {(row:any) => boolean} filter
 * @param {Uint32Array} pos_map
 * @param {number} max_len
 */
async function load_methyDiff(file_name, filter, pos_map, max_len) {
	const meth_diff_header = [
		"chr",
		"loci", "diff",
		//"pos", "value",
	];

	const tb = await loadTable(file_name, meth_diff_header);
	const rows = filter ? tb.filter(filter) : tb;

	let max = 0;
	let min = Infinity;

	const data = rows.map(row => {
		const pos = pos_map[Number(row.loci)];//ref1_pos_uint32array//Number(row.loci) / max_len;
		const value = Number(row.diff);
		if (value > max) {
			max = value;
		}
		else if (value < min) {
			min = value;
		}
		return {
			pos, value,
		};
	});

	return data;
}

/**
 * @param {string} url
 * @param {string[]} header
 * @param {number} [start_line]
 * @returns {Promise<{ [prop:string]: string }[]>}
 */
async function loadTable(url, header, start_line) {
	//const text = (await (await fetch(url)).text()).toString();
	//return table_to_object_list(tsv_parse(text), header);
	const text = await fetchData(url, "text");
	return parseTable(text, header, start_line);
}

/**
 * @param {string} text
 * @param {string[]} prop_name_list
 * @param {number} [start_line]
 * @returns {{ [prop:string]: string }[]}
 */
function parseTable(text, prop_name_list, start_line = 0) {
	return text.trim().split("\n").slice(start_line).map(line => {
		/** @type {{ [prop:string]: string }} */
		const obj = {};
		line = line.trim();
		if (line) {
			return line.split("\t").reduce((obj, col, j) => {
				obj[prop_name_list[j]] = col.trim();
				return obj;
			}, obj);
		}
		else {
			return null;
		}
	}).filter(a => a);// remove empty line / row
}

/**
 * @param {string} in_seq
 */
function parseFasta(in_seq) {
	let all = in_seq.split(">");

	let results = {};

	all.filter(a => a.length).forEach((sub_fa) => {
		let li = sub_fa.indexOf("\n");
		let out_name = sub_fa.slice(0, li).trim().replace(/:/g, "");
		let out_file_name = out_name.match(/^([^ ]+)/)[1];
		let out_seq = sub_fa.slice(li).trim();

		out_seq = out_seq.replace(/\n/g, "");
		out_seq = out_seq.replace(/\r/g, "");
		out_seq = out_seq.replace(/ /g, "");

		results[out_file_name] = out_seq.toUpperCase();
	});

	return results;
}

// /**
//  * @param {string} url
//  * @param {""|"arraybuffer"|"blob"|"json"|"text"} responseType
//  */
// async function fetchData(url, responseType) {
// 	window.$data_map = window.$data_map ?? new Map();

// 	let data = window.$data_map.get(url);
// 	if (data) {
// 		return data;
// 	}
// 	else {
// 		window.$data_map.set(url, null);// fill null
// 	}

// 	const rep = await fetch(url);

// 	if (typeof rep[responseType] == "function") {
// 		data = await rep[responseType]();
// 	}
// 	else {
// 		data = await rep.text();
// 	}

// 	window.$data_map.set(url, data);// fill data

// 	return data;
// }

async function draw_slice_region(clear = true) {
	if (clear) {
		while (viewerState.popRangeMarker());
	}

	const text = await fetchData(`tmp/multi_coord_ch${viewerState.nChr}.txt`, "text");
	const v = text.split("\n").map(line => line.split("\t").filter(a => a != "|"));
	for (let i = 4; i < 12; i += 2) {
		const ref = 0;
		const s = Number(v[i + 0][ref + 3]);
		const e = Number(v[i + 1][ref + 3]);
		const rs = ref1_pos_uint32array[s - 1];
		const re = ref1_pos_uint32array[e - 1];
		viewerState.pushRangeMarker(rs, re);
		// console.log("ref", rs, re, re - rs)
		console.log("000", s, e, e - s)
	}
}

/**
 * @param {string} url
 * @param {""|"arraybuffer"|"blob"|"document"|"json"|"text"} responseType
 */
async function fetchData(url, responseType) {
	const err_trace = new Error(url);

	window.$dataUrl_list = window.$dataUrl_list ?? [];
	window.$dataUrl_list.push(url);

	window.$data_map = window.$data_map ?? new Map();
	let data = window.$data_map.get(url);
	if (data) {
		return data;
	}
	else {
		window.$data_map.set(url, null);// fill null
	}

	const promise = new Promise(function (resolve, reject) {
		let xhr = new XMLHttpRequest();
		xhr.open("GET", url, true);
		xhr.setRequestHeader("X-Requested-With", "XMLHttpRequest");

		if (responseType) {
			xhr.responseType = responseType;
		}

		xhr.timeout = 10 * 60 * 1000;//20000;

		xhr.onload = function () {
			if (this.status == 404 || this.status == 500) {
				console.error(err_trace);

				if (!window.$dataUrl_list.includes(url)) {
					alert("file: " + url);
					debugger;
				}
				//resolve(null);

				window.$data_map.delete(url);

				reject(this.status + ": " + url);
			}
			else if (this.status == 200) {
				window.$data_map.set(url, this.response);// fill data
				resolve(this.response);
			}
		};

		xhr.ontimeout = function (e) {
			debugger;
			//resolve(null);
			reject("timeout: " + url);
		};

		xhr.onabort = function (e) {
			debugger;
			reject("abort: " + url);
		};

		xhr.send();
	});

	window.$data_map.set(url, promise);

	return await promise;
}

async function html2canvas_async(el, options) {
	if (typeof window.html2canvas != "function") {
		html2canvas = (await import(JSON.parse(JSON.stringify("/src/_html2canvas.esm.js")))).default
	}

	return await html2canvas(el, {
		scrollX: -window.scrollX,
		scrollY: -window.scrollY,
		...options,
	});

	// return new Promise((resolve, reject) => {
	// 	html2canvas(el, {
	// 		scrollX: -window.scrollX,
	// 		scrollY: -window.scrollY,
	// 		...options,
	// 	}).then(resolve, reject);
	// });
}

/**
 * @param {string} add_file_name
 * @param {string} _file_name
 * @returns {Promise<HTMLCanvasElement>}
 */
async function captureScreen(add_file_name = "", _file_name = null) {
	// await render_tetrad(analysis_options);

	/** @type {HTMLCanvasElement} */
	const canvas = await html2canvas_async(document.querySelector("#plot-area"));

	const promise_rendering = new Promise(async function (resolve, reject) {
		canvas.toBlob(function (blob) {
			if (window.$plot_blobUrl) {
				URL.revokeObjectURL(window.$plot_blobUrl);
				window.$plot_blobUrl = null;
			}
			const blobUrl = URL.createObjectURL(blob);
			window.$plot_blobUrl = blobUrl;
			try {
				const el_a = window.el_download ?? document.createElement("a");
				if (!window.el_download) {
					window.el_download = el_a;
					document.body.append(el_a);
				}
				el_a.href = blobUrl;
				// const display_list_name = display_sample.map((v, i) => v ? Object.values(methy_refMap_sampleList.Q)[i] : null).filter(a => a).join("_");
				const today = new Date();
				const file_name = `${[
					[
						today.getFullYear(),
						(today.getMonth() + 1).toString().padStart(2, "0"),
						(today.getDate()).toString().padStart(2, "0"),
					].join(""),
					dataset.mode != "single" ? (today.getTime() / 1000).toFixed(0) : "",
					// display_list_name,

					`Ch${romanize(viewerState.nChr)}`,
					`${g_chrCoord.pos_to_ref(0, g_chrCoord.bp_start)}-${g_chrCoord.pos_to_ref(0, g_chrCoord.bp_end)}`,

					document.getElementById("plot-title-text")?.dataset?.suffix,// 20221108
					
					// (function () {
					// 	if (dataset.mode != "single") {
					// 		const ts = document.getElementById("plot-title-text");
					// 		if (ts != null && ts.dataset.download) {
					// 			return ts.dataset.download;
					// 		}
					// 		else {
					// 			return dataset.name;
					// 		}
					// 	}
					// 	return "";
					// })(),

					// document.getElementById("plot-title-text")?.dataset?.download ? "" : `ch${viewerState.nChr}`,
					// window.marker_mod,
				].filter(a => a).join("_")}${add_file_name}.png`;
				el_a.download = _file_name ?? file_name;
				el_a.innerText = _file_name ?? file_name;
				el_a.click();

				resolve(canvas);
			}
			catch (ex) {
				console.error(ex);
				reject();
			}
			// finally {
			// 	URL.revokeObjectURL(blobUrl);
			// }
		}, "png");
	});

	return await promise_rendering;
}

/**
 * _getGeneById(gff["CBS1-1"]["Ch1_CBS1-1"], "000016")
 * @param {GFF_ROW[]} geneList - gff["CBS1-1"]["Ch1_CBS1-1"]
 * @param {string} id - gene ID (fuzzy search)
 */
function _getGeneById(geneList, id) {
	return geneList.find(gene => {
		const m = gene._attributes.match(RegExp(`ID=.*${id};`));
		return m ? m[0] : null;
	});
}

/**
 * _queryGeneAll(gff["CBS1-1"]["Ch1_CBS1-1"], "000016")
 * @param {GFF_ROW[]} geneList - gff["CBS1-1"]["Ch1_CBS1-1"]
 * @param {string|RegExp} regexp - this._attributes.match(regexp)
 */
function _queryGeneAll(geneList, regexp) {
	return geneList.map(gene => {
		const m = gene._attributes.match(RegExp(regexp));
		return (m && m[0]) ? gene : null;
	}).filter(a => a != null);
}

/**
 * @param {GFF_ROW[]} geneList - gff["CBS1-1"]["Ch1_CBS1-1"]
 * @param {Uint32Array} pos_map
 * @param {string} geneId
 * @param {boolean} [push_marker]
 * @param {number} [ext_left]
 * @param {number} [ext_right]
 */
function _gotoGeneById(geneList, pos_map, geneId, push_marker = false, ext_left = 0, ext_right = 0) {
	const gene = _getGeneById(geneList, geneId);

	g_chrCoord.bp_start = pos_map[gene.start - ext_left - 1];
	g_chrCoord.bp_end = pos_map[gene.end + ext_right - 1];

	if (push_marker) {
		viewerState.pushRangeMarker(gene.$start, gene.$start + gene.$length - 1);
	}
}

/**
 * @param {number} geneId
 * @param {"QM6a"|"CBS1-1"} refId
 * @param {Uint32Array} pos_map
 */
function gotoGeneById(geneId, refId, pos_map) {
	const refIdx = dataset.parental_list.indexOf(refId);

	// dataset.genome_info_list[refId].chr_list.forEach(chrInfo => {
		// const _sChr = analysis_options.chrInfo_list[refIdx].chr.replace(/_unitig_.+consensus/, "");
		// const q_sChr = _sChr.replace(dataset.parental_list[0] + "_", "");
		// const sChr = q_sChr.replace(dataset.parental_list[1] + "_", "");
		// // const sChr = chrInfo.chr.replace(/_unitig_.+consensus/, "");
		const sChr = analysis_options.chrInfo_list[refIdx].symbol;

		if (!(geneId <= 0)) {// (Number.isSafeInteger(geneId) && geneId > 0) || isNaN(geneId) // or Not a number
			_gotoGeneById(gff_data_map[refId][sChr], pos_map, geneId, true, window.ext_left ?? 1000, window.ext_right ?? 1000);
		}
	//});
}

[
	{
		el_id_search: "ref1_gene",
		el_id_select: "ref1_gene_select",
		refId: dataset.parental_list[0],
		get pos_map() { return ref1_pos_uint32array; },
	},
	{
		el_id_search: "ref2_gene",
		el_id_select: "ref2_gene_select",
		refId: dataset.parental_list[1],
		get pos_map() { return ref2_pos_uint32array; },
	},
].forEach(info => {
	const el_search = document.getElementById(info.el_id_search);
	el_search.oninput = function (evt) {
		search_ID_oninput(info, evt);
	};
});

function search_ID_oninput(info, evt) {
	const { el_id_select, refId } = info;

	const el_select = document.getElementById(el_id_select);

	// 20201113 lncRNA // const search_text = [...(evt?.target?.value ?? "")].join("*.*");// mer3 => m*.*e*.*r*.*3*.*
	const search_text = (evt?.target?.value ?? "");
	const refIdx = get_genome_index(refId);
	el_select.innerHTML = ""; //clear
	el_select.value = ""; //no select
	dataset.genome_info_list[refIdx].chr_list.forEach(chrInfo => {
		// const sChr = chrInfo.chr.replace(/_unitig_.+consensus/, "");
		const sChr = chrInfo.symbol;
		const geneList = gff_data_map[refId][sChr];
		const all = _queryGeneAll(geneList, search_text).filter(gene => ["gene", "mRNA"].includes(gene.type));
		if (all.length) {
			all.forEach(gene => {
				const geneId = gene.attributes.ID;

				const el_opt = document.createElement("option");
				// if (gene.$attributes.Name) {
				// 	debugger
				// }
				el_opt.innerText = [`Ch${chrInfo.index}`, geneId, gene.$attributes.Name].join(", "); // chrInfo.index <- chrInfo.nChr
				el_opt.value = JSON.stringify({
					nChr: chrInfo.index,
					geneId,
				});

				el_select.append(el_opt);
			});
			el_select.onchange = async function (_evt) {
				if (!el_select.value) {
					return;
				}
				try {
					const { nChr, geneId } = JSON.parse(el_select.value);

					if (nChr != viewerState.nChr) {
						viewerState.nChr = nChr;
						await loadData(true);

						el_input_chr.value = nChr; // no fire event
					}
					gotoGeneById(geneId, refId, info.pos_map);
				}
				catch (ex) {
					console.error(ex);
				}
			};
		}
	});
}

async function loadGeneTable(gff) {
	// /** @type {string[][]} */
	// const tab1 = (await fetchData("CBS1-2_QM6a.txt", "text")).trim().split("\n").slice(1).map(line => line.split("\t").map(col => col.trim()));
	// /** @type {string[][]} */
	// const tab2 = (await fetchData("CBS1-2_CBS1-1.txt", "text")).trim().split("\n").slice(1).map(line => line.split("\t").map(col => col.trim()));

	const header = [
		"id_1", "id_2", "name",
		"sChr_1", "nChr_1", "start_1", "end_1",
		"sChr_2", "nChr_2", "start_2", "end_2",
	];
	const tab1 = await loadTable("CBS1-2_QM6a.txt", header, 1);
	const tab2 = await loadTable("CBS1-2_CBS1-1.txt", header, 1);

	// [
	// 	{
	// 		refId: dataset.parental_list[0],
	// 		table: tab1,
	// 		pos_map: ref1_pos_uint32array,
	// 	},
	// 	{
	// 		refId: dataset.parental_list[1],
	// 		table: tab1,
	// 		pos_map: ref2_pos_uint32array,
	// 	},
	// ].forEach(({ refId, table, pos_map }) => {
	// 	table.forEach(row => {
	// 		const start_2 = row.start_2 = Number(row.start_2);
	// 		const end_2 = row.end_2 = Number(row.end_2);
	// 		row.$start = pos_map[start_2 - 1];
	// 		row.$length = pos_map[end_2 - 1] - row.$start + 1;
	// 	});
	// });

	const cbs2_map = new Map();

	// set gff gene name
	{
		const tab_map = {
			[dataset.parental_list[0]]: new Map(),
			[dataset.parental_list[1]]: new Map(),
		};
		// Object.keys(gff).forEach(function (refId, refIdx) {
		dataset.genomeNameList.forEach(function (refId, refIdx) {
			// const sChr = analysis_options.chrInfo_list[refIdx].chr.replace(/_unitig_.*consensus/, "");

			Object.values(gff[refId]).flat(1).forEach(gene => {
				tab_map[refId.replace("_new", "")].set(gene.attributes.ID, gene);//20210324
			});

			// /** @type {GFF_ROW[]} */
			// const gene_list = Object.values(gff[refId][sChr]).filter(a => a.type == "gene");
		});

		[
			{
				refId: dataset.parental_list[0],
				table: tab1,
				pos_map: ref1_pos_uint32array,
			},
			{
				refId: dataset.parental_list[1],
				table: tab2,
				pos_map: ref2_pos_uint32array,
			},
		].forEach(({ refId, table, pos_map }) => {
			table.forEach(row => {
				if (row.id_1 && row.id_1 != "-") {
					row.refId = refId;
					const list = cbs2_map.get(row.id_1);
					if (list) {
						list.push(row);
					}
					else {
						cbs2_map.set(row.id_1, [row]);
					}
				}
			});

			table.forEach(row => {
				if (row.name && row.name != "-") {
					if (row.id_2 && row.id_2 != "-") {
						const gene = tab_map[refId].get(row.id_2.replace(/-T\d$/, ""));
						if (gene) {
							gene.name = row.name;
						}
					}
				}
			});
		});
	}

	const el_search_gene = document.getElementById("search_gene");
	const el_search_gene_select = document.getElementById("search_gene_select");

	// const el_ref1_gene_select = document.getElementById("ref1_gene_select");
	// const el_ref2_gene_select = document.getElementById("ref2_gene_select");

	el_search_gene.disabled = true;
	el_search_gene_select.disabled = true;

	el_search_gene.oninput = function (evt) {
		const search_text = el_search_gene.value?.toLowerCase() ?? "";//[...(el_search_gene.value?.toLowerCase() ?? "")].join("*.*") + "*.*";// mer3 => m*.*e*.*r*.*3*.*
		if (!search_text) {
			return;
		}

		el_search_gene_select.innerHTML = "";// clear

		const search_results = new Map();
		[
			{
				refId: dataset.parental_list[0],
				table: tab1,
			},
			{
				refId: dataset.parental_list[1],
				table: tab2,
			},
		].forEach(({ table, refId }) => {
			table.filter(row => {
				return [
					row.id_1.replace(/-T\d$/, "").toLowerCase(),
					row.id_2.replace(/-T\d$/, "").toLowerCase(),
					row.name.toLowerCase(),
				].some(s => s.indexOf(search_text) >= 0);
			}).forEach(row => {
				const list = cbs2_map.get(row.id_1);//3-way search
				if (list) {
					list.forEach(rr => {
						// rr.refId = refId;// 20200522
						search_results.set(rr.id_2, rr);
					});
				}

				if (row.id_2 != "-") {
					row.refId = refId;
					search_results.set(row.id_2, row);
				}
			});
		});
		search_results.forEach(row => {
			if (row.sChr_2 == "-") {
				return;// no match gene
			}
			const nChr = Number(row.nChr_2);
			const el_opt = document.createElement("option");
			el_opt.innerHTML = [`${dataset.genome_info_list[dataset.parental_list.indexOf(row.refId)].chr_list[nChr - 1].symbol}:${row.id_2}`, row.name].join(", ");
			el_opt.value = JSON.stringify({
				refId: row.refId,
				geneId: row.id_2,
				geneName: row.name,
				nChr: nChr,
			});

			el_search_gene_select.append(el_opt);

			// el_ref1_gene_select.append(qm6a_gn)
			// el_search_gene_select.append(gn)
		});

		// list2 -> ref2_gene_select.append(qm6a_gn), el_search_gene_select.append(gn)
	};
	el_search_gene_select.onchange = async function (evt) {
		const { refId, geneId, geneName, nChr } = JSON.parse(el_search_gene_select.value);

		if (refId == dataset.parental_list[0]) {
			// el_ref1_gene_select.value = JSON.stringify({
			// 	nChr: nChr,
			// 	geneId: geneId.replace(/-T\d$/, ""),
			// });
			try {
				if (nChr != viewerState.nChr) {
					viewerState.nChr = nChr;
					await delayFrame();// await update UI
					await loadData(true);

					el_input_chr.value = nChr;// no fire event
				}
				gotoGeneById(geneId, dataset.parental_list[0], ref1_pos_uint32array);
			}
			catch (ex) {
				console.error(ex);
			}
		}
		else if (refId == dataset.parental_list[1]) {
			// el_ref2_gene_select.value = JSON.stringify({
			// 	nChr: nChr,
			// 	geneId: geneId.replace(/-T\d$/, ""),
			// });
			try {
				if (nChr != viewerState.nChr) {
					viewerState.nChr = nChr;
					await delayFrame();// await update UI
					await loadData(true);

					el_input_chr.value = nChr;// no fire event
				}
				gotoGeneById(geneId, dataset.parental_list[1], ref2_pos_uint32array);
			}
			catch (ex) {
				console.error(ex);
			}
		}
	};

	el_search_gene.disabled = false;
	el_search_gene_select.disabled = false;
}

// //008760
// document.getElementById("ref1_gene").oninput = function (evt) {
// 	const geneId = Number(evt.target.value);
// 	const refName = dataset.parental_list[0];
// 	const pos_map = ref1_pos_uint32array;
// 	gotoGeneById(geneId, refName, pos_map);
// };
// document.getElementById("ref2_gene").oninput = function (evt) {
// 	const geneId = Number(evt.target.value);
// 	const refName = dataset.parental_list[1];
// 	const pos_map = ref2_pos_uint32array;
// 	gotoGeneById(geneId, refName, pos_map);
// };

// {
// 	var db = window.openDatabase("mydata", "1.0", "資料庫描述", 20000);
// 	//window.openDatabase("資料庫名字", "版本","資料庫描述",資料庫大小);
// 	if (db)
// 		alert("新建資料庫成功！");
// 	db.transaction(function (tx) {
// 		tx.executeSql("CREATE TABLE test (id int UNIQUE, mytitle TEXT, timestamp REAL)");
// 	});
// 	db.transaction(function (tx) {
// 		tx.executeSql("INSERT INTO test (mytitle, timestamp) values(?, ?)", ["WEB Database", new Date().getTime()], null, null);
// 	});
// 	//db.transaction(function(tx) {
// 	//  tx.executeSql("DROP TABLE qqs");
// 	//})
// 	//db.transaction(function(tx) {
// 	//  tx.executeSql("update test set mytitle=? where mytitle = 'fsafdsaf'",['xp'],null,null);
// 	//});
// 	db.transaction(function (tx) {
// 		tx.executeSql("SELECT * FROM test", [], function (tx, result) {
// 			for (var i = 0; i < result.rows.length; i) {
// 				console.log(result.rows.item(i)['mytitle']);
// 			}
// 		}, function (...args) {
// 			console.error(args);
// 		});
// 	});
// }


// aaa.slice(0, 5000).split("\n").map(line => line.trim().split(/\s+/))
// line.startsWith("/") -> prevLine.join
// text surrounded by double quotes (")

function timeElapsed() {
	const now = new Date();
	const elapsed = now.getTime() - time_app_start.getTime();
	return (elapsed / 1000).toFixed(0) + "s";
}

function restore_error_co() {
	dataset.crossover_list = window.$co_list.map(list => [...list]);
}
function remove_error_co() {
	// backup
	window.$co_list = window.$co_list ?? dataset.crossover_list.map(list => [...list]);

	function is_equal_co(co, co2) {
		return co.every((v, i) => co2[i] == v);
	}
	function is_not_22(co) {
		return co.filter(a => a == co[0]).length != 2;
	}
	for (let chrIdx = 0; chrIdx < dataset.crossover_list.length; ++chrIdx) {
		dataset.crossover_list[chrIdx].forEach((co, i) => {
			co.$index = i;
			//co.$delete = false;
			delete co.$delete;
		});
		let filtered_list = dataset.crossover_list[chrIdx].filter((co, i) => {
			if (is_equal_co(co.before, co.after)) {
				co.$delete = true;
				console.log(co.$index, co);
				return false;
			}
			else {
				return true;
			}
		});

		for (let i = 0; i < filtered_list.length - 1; ++i) {
			let co = filtered_list[i];
			let co2 = filtered_list[i + 1];
			if (is_equal_co(co.before, co2.after)) {
				if (is_not_22(co.after) && is_not_22(co2.before)) {
					co.$delete = true;
					co2.$delete = true;

					console.log(co.$index, co);
					console.log(co2.$index, co2);
				}
			}
		}

		dataset.crossover_list[chrIdx] = filtered_list.filter(co => !co.$delete);

		// for (let i = 0; i < dataset.crossover_list[chrIdx].length - 1; ++i) {
		// 	let co = dataset.crossover_list[chrIdx][i];
		// 	let co2 = dataset.crossover_list[chrIdx][i + 1];

		// 	// if (is_equal_co(co.after, co2.before)) {
		// 	// 	if (is_not_22(co.after) && is_not_22(co2.before)) {
		// 	// 		console.log(i, co);
		// 	// 		console.log(i + 1, co2);
		// 	// 	}
		// 	// }
		// 	if (is_equal_co(co.before, co2.after)) {
		// 		if (is_not_22(co.after) && is_not_22(co2.before)) {
		// 			console.log(i, co);
		// 			console.log(i + 1, co2);
		// 		}
		// 	}
		// 	if (is_equal_co(co.before, co.after)) {
		// 		console.log(i, co);
		// 	}
		// 	if (is_equal_co(co2.before, co2.after)) {
		// 		console.log(i + 1, co2);
		// 	}

		// 	// let { before, after } = co;
		// 	// if (is_not_22(after)) {
		// 	// 	console.log("after", i, co);
		// 	// }
		// 	// if (is_equal_co(before, after)) {
		// 	// 	console.log("a==b", i, co);
		// 	// }
		// }
	}
}



// window.nco_filter = (nco) => nco.GCasso_marker > 1 && !(nco.snp_start_out >= dataset.rDNA_info.data[0].range[0] && nco.snp_start_out <= dataset.rDNA_info.data[0].range[1])

async function marker_nco() {
	dataset.non_crossover_list[0][0];

	for (let nChr = 1; nChr <= 7; ++nChr) {
		const raw_fa = await fetchData(`mafft_ch${nChr}.fa`, "text");
		const _seq_list = Object.values(await parseFasta(raw_fa));

		const ref1_pos = chrPos_to_multialign_posMap(_seq_list[0]);

		const nco_list = dataset.non_crossover_list[nChr - 1].slice(0);//clone

		const outlist = nco_list.filter((nco, idx, array) => {
			delete nco.rip_snv22_markers;;//20200629
			delete nco.maybe_rip_snv40_markers;;//20200629

			const g1s = ref1_pos[nco.snp_start_in];
			const g1e = ref1_pos[nco.snp_end_in];

			const seq = _seq_list.map(ss => ss.slice(g1s, g1e));

			array[idx].indel = null;
			seq.forEach((ss, si) => {
				if ([...ss].every(a => a == "-")) {
					if (!array[idx].indel) {
						array[idx].indel = (si + 1);
					}
					else {
						array[idx].indel = array[idx].indel + "," + (si + 1);
					}
				}
			});

			const a = seq.map(ss => [...ss].find(a => a != "-")).filter(a => a)[0];
			if (a) {
				if (seq.length >= 1) {
					if (!seq.every(ss => [...ss].every(v => v == "-" || v == a))) {
						array[idx].poly_mer = Math.max(...seq.map(ss => ss.length));
						return true;
					}
				}
			}
			return false;
		});

		const text = JSON.stringify(outlist);

		downloadTextFile(`nco_ch${nChr}_filtered.json`, text);
	}
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {Uint32Array}
 * @version 20220627 webView
 */
function chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 1;
	let posmap = new Uint32Array(ma_seq.length);
	// for (let index = 0; index < ma_seq.length; ++index) {
	// 	const element = ma_seq[index];
	// 	// if (element != "-") {
	// 	// 	posmap[nPos] = index + 1;
	// 	// 	++nPos;
	// 	// }
	// }
	// 20220627
	for (let pos = 0; pos < ma_seq.length; ++pos) {
		const element = ma_seq[pos];
		(function (ref_pos) {
			if (!posmap[ref_pos]) {
				posmap[ref_pos] = pos + 1;
			}
		})(Math.max(1, element != "-" ? nPos : (nPos - 1)));
		
		if (element != "-") {
			++nPos;
		}
		// if (nPos == 10) {
		// 	console.log(posmap.slice(0, 10));
		// 	break;
		// }
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @see {@link chrPos_to_multialign_posMap}
 * @param {string[]} ma_seq
 * @returns {Uint32Array[]}
 */
function ms_ma_to_pos(ma_seq) {
	const posmap = ma_seq.map(seq => new Uint32Array(seq.length));

	//

	return posmap;
}

/**
 * index: 0 ~ length - 1
 * posmap[index + 1] = ++nPos;
 * @param {string} ma_seq
 * @returns {Uint32Array}
 * @version 20220627 webView
 */
function multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 1;
	let posmap = new Uint32Array(ma_seq.length);

	posmap[0] = 0;

	for (let pos = 0; pos < ma_seq.length; ++pos) {
		const element = ma_seq[pos];

		// if (element != "-") {
		// 	++nPos;
		// }
		// // posmap[pos + 1] = nPos;
		// posmap[pos] = nPos - 1;// idx to idx

		(function (ref_pos) {
			posmap[pos] = ref_pos;
		})(Math.max(1, element != "-" ? nPos : (nPos - 1)));

		if (element != "-") {
			++nPos;
		}
	}
	
	// const a = posmap;
	// const b = pos_ref1_uint32array;
	// for (let i = 0; i < a.length; ++i) {
	// 	if (a[i] != b[i]) {
	// 		console.log(i, a[i], b[i]);
	// 		break;
	// 	}
	// }

	return posmap;
}

/**
 * @param {boolean} reload_current
 * @param {()=>Promise<void>} onload
 * @param {{ rep?: number; beforeCapture?: (rep: number, nChr: number)=>Promise<void>; }} options
 * @returns {Promise<HTMLCanvasElement|(HTMLCanvasElement[])>}
 */
async function captureScreenAll(reload_current = true, onload, start_chr, end_chr, download = true, options = null) {//height = 1920
	start_chr = start_chr != null ? start_chr : 1;
	end_chr = end_chr >= 1 && end_chr <= dataset.genome_info_list[0].chr_list.length ? end_chr : dataset.genome_info_list[0].chr_list.length;

	/**
	 * @type {HTMLCanvasElement|(HTMLCanvasElement[])}
	 */
	const save_canvas = [];

	for (let nChr = start_chr; nChr <= end_chr; ++nChr) {
		try {
			await promise_load_task;
			await delayFrame();
		}
		catch (ex) {
			console.error("captureScreenAll", ex);
		}

		if (reload_current == true || viewerState.nChr != nChr) {
			try {
				el_input_chr.value = String(nChr);
				el_input_chr.oninput(null);
				await delayFrame();
			}
			catch (ex) {
				console.error("captureScreenAll", ex);
			}

			try {
				await promise_load_task;
			}
			catch (ex) {
				console.error("captureScreenAll", ex);
			}

			if (typeof onload == "function") {
				await onload();
			}
		}

		// if (height) {
		// 	viewerState.get_plot_height = () => height;
		// 	// document.querySelector("#canvas").style.height = height + "px";
		// 	document.querySelector("#canvas").height = height;
		// }

		// load end
		await delayFrame();

		// resize canvas
		await show_all();
		viewerState.resizeCanvas();
		await delayFrame();

		await show_all();
		await delayFrame();

		if (download) {
			if (options && options.rep > 0 && options.beforeCapture != null) {
				const rep_canvas = [];
				for (let j = 0; j < options.rep; ++j) {
					await options.beforeCapture(j, nChr);
					await delayFrame();
					await drawFrame();
					const canvas = await captureScreen();
					rep_canvas.push(canvas);
				}
				save_canvas.push(rep_canvas);
			}
			else {
				const canvas = await captureScreen();
				save_canvas.push(canvas);
			}
		}
		console.log(nChr);
	}

	return save_canvas;
}

/**
 * @param {number} ref_idx
 * @returns {ChromosomeData[]}
 */
function get_chrList_order_by_length(ref_idx) {
	return [...dataset.genome_info_list[ref_idx].chr_list].sort((a, b) => b.length - a.length);//.map(a => a.index)
}

async function asdasdasdas() {
	const lines = [];

	for (let i = 0; i < 7; ++i) {
		viewerState.nChr = i + 1;
		const seq_list = Object.values(dataset.results[i]);
		pos_ref1_uint32array = multialign_to_chrPos_posMap(seq_list[0]);
		ref1_pos_uint32array = chrPos_to_multialign_posMap(seq_list[0]);

		await load_cmp_ss(false, true);

		const chr_line = methyl_dataset_list.filter(aaa => aaa.name.match(/→|ilu/)).map(aa => [
			aa.name,
			viewerState.nChr,
			aa.data[viewerState.nChr - 1].filter(v => v.value >= 0).length,
			aa.data[viewerState.nChr - 1].filter(v => v.value < 0).length
		].join("\t"));

		lines.push(...chr_line);
	}

	const textarea = document.querySelector("#top_GUI > table:nth-child(1) > tbody > tr > td:nth-child(3) > textarea");

	textarea.value = lines.join("\n");

	//methyl_dataset_list.filter(aaa => aaa.name.match(/→|ilu/)).map(aaa => aaa.data.map((aa, aa_idx) => aa.map(a => [aaa.name, aa_idx + 1, a.$start, a.$end, a.strand, a.value].join("\t")).join("\n")).join("\n")).join("\n")
}

async function ss_dna_tab() {
	const lines = [];

	for (let chr_idx = 0; chr_idx < dataset.genome_info_list[0].chr_list.length; ++chr_idx) {
		await load_multi_alignment_results(chr_idx);
		pos_ref1_uint32array = multialign_to_chrPos_posMap(seq_list[0]);
		ref1_pos_uint32array = chrPos_to_multialign_posMap(seq_list[0]);

		await load_cmp_ss(false, true);

		const _sp_s_r = methyl_dataset_list[1];
		const _sp_ss_sr = methyl_dataset_list[3];
		const list = [_sp_s_r, _sp_ss_sr];
		for (let li of list) {
			const chr_line = li.data[chr_idx].map(a => [
				li.name,
				chr_idx + 1,
				pos_ref1_uint32array[a.start - 1],
				pos_ref1_uint32array[a.end - 1],
				// a.strand,
				a.value,
			].join("\t"));

			lines.push(...chr_line);
		}
	}

	const textarea = document.querySelector("#top_GUI > table:nth-child(1) > tbody > tr > td:nth-child(3) > textarea");

	textarea.value = lines.join("\n");
}

async function load_multi_alignment_results(chr_idx) {
	// viewerState.nChr = chr_idx + 1;// auto requist frame
	// viewerState._nChr = chr_idx + 1;

	const seq_list = await (async function () {
		const obj = dataset.results[chr_idx];

		if (typeof obj == "string") {
			let raw_fa = await fetchData(obj, "text");
			let fa = await parseFasta(raw_fa);
			dataset.results[chr_idx] = fa;
		}

		return Object.values(dataset.results[chr_idx]);
	})();

	return seq_list;
}

async function cap_chr(start = 1, end = dataset.genome_info_list[0].chr_list.length, resize_height = 844) {
	hide_all_marker();

	for (let nChr = start; nChr <= end; ++nChr) {
		await promise_load_task;
		el_input_chr.value = String(nChr);
		el_input_chr.oninput(null);
		await promise_load_task;

		viewerState.fixedPlotHeight = resize_height;
		// viewerState.resizeCanvas(); await _drawFrame();

		await load_cmp_ss(false, true);

		await delayFrame();

		await show_all();

		await delayFrame();

		await captureScreen();
	}
}

// /** @type {HTMLCanvasElement[]} */
// capture_all_new_wnd.chr_canvas_list = [];

capture_all_new_wnd.scale_by_genome = true;

capture_all_new_wnd.save_per_chr = false;

// step_capture_all_new_wnd(1, true, 1)
async function step_capture_all_new_wnd(step = 2, save = false, height_scale = 1) {
	capture_all_new_wnd.scale_by_genome = false;

	capture_all_new_wnd.add_row_height = (
		(viewerState.seg_row_height + viewerState.seg_row_separate) * 4 - 2 + // progeny
		0
	);
	capture_all_new_wnd.markers = [
		"el_display_rDNA",
	];

	if (step == 1) {
		capture_all_new_wnd.save_per_chr = false;
	}

	viewerState.disable_max_length = true;

	for (let i = 1; i <= 16; i += step) {
		console.log(i, i - 1 + step);
		await capture_all_new_wnd(i, i - 1 + step, save, height_scale);

		await delayFrame();
	}
}

async function cap_v3() {
	hide_all_marker();

	capture_all_new_wnd.add_row_height = (
		(viewerState.seg_row_height + viewerState.seg_row_separate) * 4 - 2 + // progeny
		(viewerState.seg_row_height + viewerState.seg_row_separate) * 4 // markers
	);
	capture_all_new_wnd.markers = [
		"el_display_rDNA",
		"el_display_co", "el_display_conco", "el_display_nonCo", "el_display_nonCo2"
	];

	if (capture_all_new_wnd.markers) {
		capture_all_new_wnd.markers.forEach(id => {
			const el = document.getElementById(id);
			if (el != null && !el.checked) {
				el.click();
			}
		});
	}

	const resize_height = 535;

	viewerState.fixedPlotHeight = resize_height;
	// viewerState.resizeCanvas(); await _drawFrame();

	for (let nChr = 1; nChr <= 16; nChr += 1) {
		try {
			await promise_load_task;
			el_input_chr.value = String(nChr);
			el_input_chr.oninput(null);

			await promise_load_task;
		}
		catch (ex) {
			console.error(ex);
		}

		await capture_all_new_wnd.onload_chr?.();

		await _show_all();

		await delayFrame();

		await captureScreen();
	}
}

/** @type {HTMLCanvasElement[]} */
const capture_all_new_wnd_chr_canvas_list = [];

/**
 * @param {number} start_order
 * @param {number} end_order
 * @param {boolean} [save_as_png]
 * @param {number} [height_scale]
 * @param {number} [resize_height]
 */
async function capture_all_new_wnd(start_order = 1, end_order = dataset.genome_info_list[0].chr_list.length, save_as_png = false, height_scale = 0.2, resize_height = null) {
	viewerState.disable_max_length = !(capture_all_new_wnd.scale_by_genome ?? true);

	// hide_all_marker();
	// if (capture_all_new_wnd.markers) {
	// 	capture_all_new_wnd.markers.forEach(id => {
	// 		const el = document.getElementById(id);
	// 		if (el != null && !el.checked) {
	// 			el.click();
	// 		}
	// 	});
	// }

	const iframe = document.createElement("iframe");

	// const wnd_uid = "wnd_" + (new Date()).getTime().toString(36) + Math.trunc(Math.random() * Number.MAX_SAFE_INTEGER).toString(36);
	// iframe.id = wnd_uid;
	const iframe_class = "output_preview_wnd";
	iframe.classList.add(iframe_class);

	// const chr_list = dataset.genome_info_list[0].chr_list;

	const all_chr_list = get_chrList_order_by_length(0);
	const view_nChr_list = all_chr_list.map(a => a.index).slice(start_order - 1, end_order);

	// const max_all_chr_length = (function () {
	// 	return Math.max(...dataset.genome_info_list.map(aa => Math.max(...aa.chr_list.map(a => a.length))));
	// })();

	const max_all_chr_length = viewerState.max_all_chr_length;

	const {
		seg_row_height,
		seg_row_separate,
	} = viewerState;

	const src_x = 0;
	// const src_y = (
	// 	seg_row_height + seg_row_separate +   // parental
	// 	seg_row_height * 2 + seg_row_separate // gc
	// 	) * 2;
	// const src_height = Math.ceil(
	// 	(seg_row_height + seg_row_separate) * 4
	// 	// + (seg_row_height * 4 + seg_row_separate) * 2 // spo11-oligo DSB
	// ) - seg_row_separate + 1;//remove last seg_row_separate
	const src_y = 0;
	const src_height = viewerState.get_plot_height() ?? Math.ceil(
		// (
		// 	seg_row_height + seg_row_separate +   // parental
		// 	seg_row_height * 2 + seg_row_separate // gc
		// ) * 2 - 1 +
		// (capture_all_new_wnd.add_row_height ?? 0)
		((seg_row_height + seg_row_separate) + (seg_row_height * 2 + seg_row_separate)) +
		(seg_row_height + seg_row_separate) * 6
	);

	/** @type {HTMLCanvasElement[]} */
	let chr_canvas_list = capture_all_new_wnd_chr_canvas_list;

	//dataset.name.replace(/∕/g, "／")
	// const title = `${dataset.name}（${Object.keys(dataset.progeny).join(",")}）（ch${all_chr_list[start_order - 1].index}～${all_chr_list[end_order - 1].index}）`;
	const title = `${dataset.name}（${Object.keys(dataset.progeny).join(",")}）`;

	iframe.srcdoc = `<!DOCTYPE html>
<html lang="zh_TW">
	<head>
		<title>${title}</title>
	</head>
	<body>
		<canvas id="canvas" style="outline: 1px solid black; position: fixed; left: 0px; top: 0px; z-index: 1;"></canvas>
		<button id="btn_close" style="outline: 1px solid black; position: fixed; left: 0px; top: 0px; z-index: 10;">close</button>
		<button id="btn_close_any" style="outline: 1px solid black; position: fixed; left: 100px; top: 0px; z-index: 10;">close any</button>
		<a id="btn_download" style="outline: 1px solid black; position: fixed; left: 200px; top: 0px; z-index: 10;" download="${title}.png">download</a>
	</body>
</html>
	`;
	// iframe.srcdoc = (`<canvas id="canvas" width="${main_canvas.width}" height="${all_chr * src_height}" style="outline: 1px solid black;"></canvas>`);
	document.body.append(iframe);

	await new Promise((resolve, reject) => {
		iframe.onload = function() {
			resolve();
		};
	});

	iframe.contentWindow.document.getElementById("btn_close").onclick = () => iframe.remove();
	iframe.contentWindow.document.getElementById("btn_close_any").onclick = function () {
		window.document.querySelectorAll(`.${iframe_class}`).forEach(el => el.remove());
	};

	// load

	viewerState.auto_padding_right = false;
	//
	for (let nChr of view_nChr_list) {
		const chr_idx = nChr - 1;

		if (chr_canvas_list[nChr]) {
			continue;
		}

		// if (nChr != 4) {
		// 	continue;
		// }

		try {
			await promise_load_task;
			el_input_chr.value = String(nChr);
			el_input_chr.oninput(null);

			await promise_load_task;
			if (resize_height > 0) {
				viewerState.fixedPlotHeight = resize_height;
				// viewerState.resizeCanvas(); await _drawFrame();
			}
		}
		catch (ex) {
			console.error(ex);
		}

		await capture_all_new_wnd.onload_chr?.();

		await _show_all();
		// await delayFrame();

		// if (save_as_png) {
		// 	await captureScreen();
		// }

		/** @type {HTMLCanvasElement} */
		const canvas = chr_canvas_list[nChr] ?? (document.createElement("canvas"));
		const ctx = canvas.getContext("2d");

		chr_canvas_list[nChr] = canvas;

		canvas.width = viewerState.max_view_width + 1;
		canvas.height = src_height;
		// console.log({
		// 	scale: viewerState.getChrInGenomeScale(),
		// 	src_width,
		// });
		// debugger

		ctx.setTransform(1, 0, 0, 1, 0, 0);
		ctx.drawImage(main_canvas,
			src_x, src_y, canvas.width, canvas.height,
			0, 0, canvas.width, canvas.height
		);

		// const max_chr_length = Math.max(...dataset.genome_info_list.map(aa => aa.chr_list[chr_idx].length));
		// const scale = max_chr_length / max_all_chr_length;
		// const x = 0;
		// const y = chr_idx * src_height;
		// ctx.drawImage(main_canvas,
		// 	src_x, src_y, src_width, src_height,
		// 	x, y, src_width * scale, src_height
		// );

		if (save_as_png && capture_all_new_wnd.save_per_chr) {
			downloadCanvasAsPNG(canvas, viewerState.$plot_title);
		}

		console.log(viewerState.$plot_title);
	}

	// fin

	/** @type {HTMLCanvasElement} */
	const rows_head_canvas = await (async function () {
		const el = document.getElementById("data-rows");
		const c = await html2canvas_async(el);
		return c;
	})();

	const chr_plot_x = rows_head_canvas.width;// + 100;
	const chr_plot_separate = 20;
	const scale_bar_height = 70;
	const title_bar_height = 30;

	/** @type {HTMLCanvasElement} */
	const canvas = (iframe.contentWindow.document.getElementById("canvas"));
	const ctx = canvas.getContext("2d");

	ctx.imageSmoothingQuality = "high";
	ctx.imageSmoothingEnabled = false;

	// const height_scale = 0.5;
	const chr_height = src_height * height_scale + chr_plot_separate;

	const max_src_canvas_width = all_chr_list.reduce((max_width, chr_info) => {
		/** @type {HTMLCanvasElement} */
		const canvas = chr_canvas_list[chr_info.index];
		if (typeof canvas == "object") {
			return Math.max(canvas.width, max_width);
		}
		else {
			return max_width;
		}
	}, 0);
	const output_width = rows_head_canvas.width + max_src_canvas_width + chr_plot_x + 100;
	const output_height = title_bar_height + view_nChr_list.length * chr_height + scale_bar_height;

	canvas.width = output_width;
	canvas.height = output_height;

	ctx.setTransform(1, 0, 0, 1, 0, 0);
	ctx.fillStyle = "#FFFFFF";
	ctx.fillRect(0, 0, canvas.width, canvas.height);

	if (!viewerState.disable_max_length) {
		// scale bar
		ctx.save();
		{
			ctx.translate(0.5, 0.5);

			ctx.fillStyle = "black";
			ctx.font = "24px Arial";
			ctx.textBaseline = "bottom";
			ctx.textAlign = "center";
			ctx.fillText(title, canvas.width * 0.5, 30);
			ctx.translate(0, 30);

			const base_scale = 1000;
			const scale_bar_step_size = 10 ** Math.trunc(Math.log10(max_all_chr_length)) / 10;// 100_000
			for (let i = 0; i <= Math.floor(max_all_chr_length / scale_bar_step_size); ++i) {
				const p1 = (i + 0) * scale_bar_step_size;
				const p2 = (i + 1) * scale_bar_step_size;
				const x1 = Math.round(chr_plot_x + max_src_canvas_width * (p1 / max_all_chr_length));
				const x2 = Math.round(chr_plot_x + max_src_canvas_width * (p2 / max_all_chr_length));
				const y0 = output_height - scale_bar_height;
				const y1 = y0 - chr_plot_separate * 0.5;
				const y2 = y0 + 5;
				const y3 = y2 + 15;

				ctx.beginPath();

				ctx.moveTo(x1, y1);
				ctx.lineTo(x1, y3);

				ctx.moveTo(x1, y2);
				ctx.lineTo(x2, y2);

				ctx.moveTo(x2, y1);
				ctx.lineTo(x2, y3);

				ctx.strokeStyle = "#000000";
				ctx.stroke();

				ctx.font = "14px Arial";
				ctx.textAlign = "left";
				ctx.textBaseline = "bottom";
				ctx.fillStyle = "#000000";
				ctx.fillText(`${p1 / base_scale}kb`, Math.round(x1) + 5, y3);

				// console.log({ x1, x2, y1, y2 });
			}
		}
		ctx.restore();
	}

	ctx.translate(0, 30);//title

	const output_row_height = Math.round(src_height * height_scale);
	for (let [sub_plot_index, nChr] of Object.entries(view_nChr_list)) {
		const chr_idx = nChr - 1;

		/** @type {HTMLCanvasElement} */
		const canvas = chr_canvas_list[chr_idx + 1];
		if (typeof canvas == "object") {
			// const max_chr_length = Math.max(...dataset.genome_info_list.map(aa => aa.chr_list[chr_idx].length));
			// const scale = max_chr_length / max_all_chr_length;
			const y = (all_chr_list.findIndex(a => a.index == nChr) - (start_order - 1)) * chr_height;

			ctx.drawImage(rows_head_canvas,
				0, 0, rows_head_canvas.width, output_row_height,
				0, y, rows_head_canvas.width, output_row_height);

			{
				ctx.fillStyle = "#000000";
				ctx.textAlign = "right";
				ctx.textBaseline = "middle";

				// const fontSize = Math.round((src_height * height_scale * 0.5 - viewerState.seg_row_separate) * 0.75);
				// ctx.font = `${fontSize}px Arial`;

				ctx.font = `24px Arial`;
				ctx.fillText(`Ch${chr_idx + 1}`, 80, y + src_height * height_scale * 0.5);
			}

			// ctx.drawImage(canvas, chr_plot_x, y, Math.round(src_width * scale), output_row_height);
			if (canvas.width == 0 || output_row_height == 0) {
				debugger;
			}
			ctx.drawImage(canvas, chr_plot_x, y, canvas.width, output_row_height);
		}
	}

	iframe.style.position = "fixed";
	iframe.style.left = "0px";
	iframe.style.top = (start_order - 1) / view_nChr_list.length * canvas.height + "px";
	iframe.style.width = canvas.width + "px";
	iframe.style.height = canvas.height + "px";
	iframe.style.border = "none";
	iframe.style.outline = "auto";

	if (save_as_png) {
		console.log(start_order, end_order);
		downloadCanvasAsPNG(canvas, `${title}.png`, iframe.contentWindow.document);
	}

	const url = canvas.toDataURL("png");

	/** @type {HTMLAnchorElement} */
	const btn_download = (iframe.contentWindow.document.getElementById("btn_download"));
	btn_download.href = url;

	return iframe;
}

function downloadAllRawGenome() {
	const genome = [];
	Object.entries(dataset.results[0]).map((aa, si) => [aa[0].replace(dataset.genomeNameList[si], "").slice(1), aa]).forEach((aa, si) => {
		genome[si] = genome[si] ?? [];
		genome[si].push(aa);
	});
	genome.forEach((aaa, si) => {
		const t = aaa.map(aa => `>${aa[0]}\n${aa[1]}`).join("\n");
		downloadTextFile(`${dataset.genomeNameList[si]}.fa`, t);
	});
}

class DownloadHelper {
	_href = "";

	/**
	 * @param {string|ArrayBuffer} data
	 * @param {string} mime
	 */
	constructor(data, mime = "text/plain") {
		const blob = new Blob([data], {
			type: mime,
		});

		this._href = URL.createObjectURL(blob);

		return this;
	}

	/**
	 * @param {string|ArrayBuffer} data
	 * @param {string} mime
	 */
	buffer(data, mime = "application/octet-stream") {
		alert("not implement");
		throw new Error("not implement");
		return;

		const blob = new Blob([data], {
			type: mime,
		});

		this._href = URL.createObjectURL(blob);

		return this;
	}

	/**
	 * @param {string} filename
	 */
	download(filename) {
		const a = document.createElement("a");
		a.innerText = filename;
		a.download = filename;
		a.href = this._href;
		document.body.append(a);
		a.click();
		a.remove();

		URL.revokeObjectURL(this._href);
	}
}

/**
 * @param {string} download_name
 * @param {string} content
 */
function downloadTextFile(download_name, content) {
	const a = document.createElement("a");
	a.innerText = download_name;
	a.download = download_name;
	a.href = `data:text/plain;base64,${btoa(content)}`;
	document.body.append(a);
	a.click();
	a.remove();
}

/**
 * @param {HTMLCanvasElement} canvas
 * @param {string} file_name
 * @param {Document} [doc]
 */
async function downloadCanvasAsPNG(canvas, file_name, doc = canvas.ownerDocument || window.document) {
	const task = new Promise((resolve, reject) => {
		canvas.toBlob(function (blob) {
			let blobUrl;
			let el_a;
			try {
				blobUrl = URL.createObjectURL(blob);
				el_a = doc.createElement("a");
				doc.body.append(el_a);
				el_a.href = blobUrl;
				el_a.download = file_name;
				el_a.innerText = file_name;
				el_a.click();

				resolve();
			}
			catch (ex) {
				console.error(ex);
				reject();
			}
			finally {
				if (blobUrl) {
					URL.revokeObjectURL(blobUrl);
				}
				if (el_a) {
					el_a.remove();
				}
			}
		}, "png");
	});
	return await task;
}

// dataset.non_crossover_list[viewerState.nChr - 1].filter(nco => {
// 	const g1s = ref1_pos_uint32array[nco.snp_start_in];
// 	const g1e = ref1_pos_uint32array[nco.snp_end_in];
// 	const l2s = pos_ref2_uint32array[g1s];
// 	const l2e = pos_ref2_uint32array[g1e];
// 	const g2s = ref2_pos_uint32array[l2s];
// 	const g2e = ref2_pos_uint32array[l2e + 1];

// 	const ggs = g1s;//Math.min(g1s, g2s);
// 	const gge = g1e;//Math.max(g1e, g2e);

// 	const seq = seq_list.map(ss => ss.slice(ggs, gge));

// 	const a = seq.map(ss => [...ss].find(a => a != "-")).filter(a => a)[0];
// 	if (a) {
// 		console.log(ggs, seq.map((ss, i) => [i, ss.length, ss].join(",")).join(","));
// 		return !seq.every(ss => [...ss].every(v => v == "-" || v == a));
// 	}
// 	return false;

// 	const s1 = seq_list[0].slice(g1s, g1e);
// 	const s2 = seq_list[1].slice(g2s, g2e);
// 	console.log(g1s, s1, g2s, s2);
// 	return [...s1].every((v, i) => v == s1[i]) && [...s2].every((v, i) => v == s2[i]);
// });

// dataset.non_crossover_list[viewerState.nChr - 1].filter(nco => {
// 	const g1s = ref1_pos_uint32array[nco.snp_start_in];
// 	const g1e = ref1_pos_uint32array[nco.snp_end_in];
// 	const l2s = pos_ref2_uint32array[g1s];
// 	const l2e = pos_ref2_uint32array[g1e];
// 	const g2s = ref2_pos_uint32array[l2s];
// 	const g2e = ref2_pos_uint32array[l2e + 1] - 1;

// 	const ggs = g1s;//Math.min(g1s, g2s);
// 	const gge = g1e;//Math.max(g1e, g2e);

// 	const seq = seq_list.map(ss => ss.slice(ggs, gge));

// 	const a = seq.map(ss => [...ss].find(a => a != "-")).filter(a => a)[0];
// 	if (a) {
// 		console.log(ggs, seq.map((ss, i) => i + ss).join(","));
// 		return seq.every(ss => [...ss].every(v => v == "-" || v == a));
// 	}
// 	return false;

// 	const s1 = seq_list[0].slice(g1s, g1e);
// 	const s2 = seq_list[1].slice(g2s, g2e);
// 	console.log(g1s, s1, g2s, s2);
// 	return [...s1].every((v, i) => v == s1[i]) && [...s2].every((v, i) => v == s2[i]);
// });

{
	const el_clipboard = document.createElement("td");
	const el_input = document.createElement("textarea");
	el_input.id = "clipboard";
	const el_btn = document.createElement("button");
	el_clipboard.append(el_input);
	el_clipboard.append(el_btn);
	document.getElementById("status").after(el_clipboard);

	el_btn.innerText = "copy seq";
	el_btn.onclick = (evt) => copy_visable_seq(false);
}

/**
 * @param {boolan} [rev_comp]
 * @see {extract_seq}
 */
function copy_visable_seq(rev_comp = false) {
	const el_input = document.getElementById("clipboard");
	const fa = extract_seq(g_chrCoord.bp_start, g_chrCoord.bp_end, rev_comp);
	el_input.value = fa;
	el_input.select();
	document.execCommand("copy");
}

/**
 * 20220627
 * @param {number} seq_idx
 * @param {number} align_start
 * @param {number} align_end
 * @see {@link multialign_to_chrPos_posMap}
 */
function sub_seq_from_align(seq_idx, align_start, align_end) {
	return seq_list[seq_idx].slice(align_start - 1, align_end).replace(/-/g, "");
}

/**
 * 20220627
 * @param {number} seq_idx
 * @param {number} align_start
 * @param {number} align_end
 * @see {@link multialign_to_chrPos_posMap}
 */
function sub_seq_from_orig(seq_idx, start, end) {//404077, 404080
	return seq_list[seq_idx].replace(/-/g, "").slice(start - 1, end);
}

function test_sub_seq_all() {
	// 20220627
	const mamap = make_pos_ref_map_list();

	seq_list.map((_, i) => {
		const sa = sub_seq_from_align(i, g_chrCoord.bp_start, g_chrCoord.bp_end);
		const so = sub_seq_from_orig(i, mamap[g_chrCoord.bp_start] - 1, mamap[g_chrCoord.bp_end]);
		if (sa != so) {
			console.log(i, sa.slice(0, 10), so.slice(0, 10));
		}
	});
}

/**
 * @param {number} ma_start
 * @param {number} ma_end
 * @param {boolean} [rev_comp]
 * @param {number} [pos]
 * @see {copy_visable_seq}
 * @example extract_seq(g_chrCoord.bp_start, g_chrCoord.bp_end)
 */
function extract_seq(ma_start, ma_end, rev_comp = false, pos = Math.trunc((ma_start + ma_end) / 2), only_head = false) {
	// seq_list[0].slice(g_chrCoord.bp_start - 1, g_chrCoord.bp_end).replace(/-/g, "") == seq_list[0].replace(/-/g, "").slice(pos_ref1_uint32array[g_chrCoord.bp_start - 1] - 1, pos_ref1_uint32array[g_chrCoord.bp_end - 1])
	// seq_list[1].slice(g_chrCoord.bp_start - 1, g_chrCoord.bp_end).replace(/-/g, "") == seq_list[1].replace(/-/g, "").slice(   pos_subject_1_map[g_chrCoord.bp_start - 1],        pos_subject_1_map[g_chrCoord.bp_end - 1] + 1)
	const seq_names = dataset.genomeNameList;

	const ma_seq_start_idx = ma_start - 1;
	const ma_seq_end_idx = ma_end - 1;

	const range_info = [
		[dataset.parental_list[0], pos_ref1_uint32array[ma_seq_start_idx], pos_ref1_uint32array[ma_seq_end_idx]].join("-"),
		[dataset.parental_list[1], pos_ref2_uint32array[ma_seq_start_idx], pos_ref2_uint32array[ma_seq_end_idx]].join("-"),
		["viewer", ma_start, ma_end].join("-"),
	].join(" ");

	// 20220627
	const mamap = make_pos_ref_map_list();

	const fa = seq_list.map((ss, i) => {
		const pos_map = mamap[i];
		// const head = ">" + seq_names[i] + "_Ch" + viewerState.nChr + " " + range_info + "\n";
		const this_start_pos = pos_map[ma_seq_start_idx];// + (i >= dataset.parental_list.length ? 1 : 0);
		const this_end_pos = pos_map[ma_seq_end_idx];// + (i >= dataset.parental_list.length ? 1 : 0);
		const this_length = this_end_pos - this_start_pos + 1;

		const g_info = dataset.genome_info_list[i];
		const chr_info = g_info.chr_list[viewerState.nChr - 1];
		const gff_list = gff_data_map?.[g_info.name]?.[chr_info.raw_chr_name];

		const gff_in_region = gff_list ? gff_list.filter(a => {
			// const [this_start_pos, this_end_pos] = [pos_ref1_uint32array[g_chrCoord.bp_start - 1], pos_ref1_uint32array[g_chrCoord.bp_end - 1]];
			return a.type == "gene" && (a.start <= this_end_pos && a.end >= this_start_pos);
		}) : [];
		const gene_info_list = gff_in_region.map(gene => {
			const [
				rel_start,
				rel_end,
			] = [
				Number(gene.start) - this_start_pos + 1,
				Number(gene.end) - this_start_pos + 1,
			];
			if (rev_comp) {
				const strand = gene.strand > 0 ? "-" : (gene.strand < 0 ? "+" : "");
				return [
					gene.geneID,
					strand,
					Math.max(1, this_length - rel_end + 1),
					Math.min(this_length - rel_start + 1, this_length),
				].join(":");
			}
			else {
				const strand = gene.strand > 0 ? "+" : (gene.strand < 0 ? "-" : "");
				return [
					gene.geneID,
					strand,
					Math.max(1, rel_start),
					Math.min(rel_end, this_length),
				].join(":");
			}
		});
		if (rev_comp) {
			gene_info_list.reverse();
		}
		const gene_info = gene_info_list.join(" ");

		const head = `>${seq_names[i]}_Ch${viewerState.nChr}_s${this_start_pos}_e${this_end_pos} ${rev_comp ? "reverse complement " : ""}${range_info} ${gene_info}\n`;
		const body = ss.slice(ma_start - 1, ma_end).replace(/-/g, "");
		const _body = ss.replace(/-/g, "").slice(this_start_pos - 1, this_end_pos);
		if (body && body != _body) {// is has seq and equal two method
			console.log({
				name: seq_names[i],
				idx: i,
				"from align": [body, ma_start - 1, ma_end],
				"from orig": [_body, this_start_pos - 1, this_end_pos],
			});
			alert("error seq coord");
		}
		return head + (only_head ? "" : (rev_comp ? reverseComplement(body) : body));
	}).join("\n");

	return fa;
}

/**
 * @param {string} seq
 */
function reverseComplement(seq) {
	const complement_map = {};
	complement_map["A"] = "T";
	complement_map["T"] = "A";
	complement_map["C"] = "G";
	complement_map["G"] = "C";
	complement_map["N"] = "N";
	complement_map["a"] = "t";
	complement_map["t"] = "a";
	complement_map["c"] = "g";
	complement_map["g"] = "c";
	complement_map["n"] = "n";
	
	/** @type {Map<string, number} */
	const hasUnknow = new Map();
	
	const result = [...seq].reverse().map(a => {
		if (complement_map[a] == null || complement_map[a] == "N" || complement_map[a] == "n") {
			hasUnknow.set(a, (hasUnknow.get(a) || 0) + 1);
		}
		return complement_map[a];
	}).join("");

	if (result.length != seq.length || hasUnknow.size) {
		for (let ent of hasUnknow.entries()) {
			console.warn(`found ${ent[0]}: ${ent[1]}, ${(ent[1] * 100 / seq.length).toFixed(2)}%`);
		}
	}

	return result;
}

function find_rip_out_AT_island(view_range = 30, extract_range = 500) {
	const seq_names = dataset.genomeNameList;
	const ref1_name = seq_names[0];

	const fm = {
		"A": 0,
		"T": 0,
		"CBS1-1": 0,
		"G": 0,
	};
	[...seq_list[0]].forEach(a => ++fm[a]);
	// ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06)
	const block_min_gc = ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06) * 100;

	/** @type {"SNV RIP"|"RIP"} */
	const rip_snp_type = "RIP";
	const rip_markers_list = allMarker.list.filter(a => a.name.indexOf(rip_snp_type) >= 0);

	const map_pos_ref1 = pos_ref1_uint32array;
	const map_ref1_pos = ref1_pos_uint32array;

	const result = rip_markers_list.map(grp => {
		return grp.values.map(a => {
			const gc_wnd = dataset.gc_content[dataset.ref][viewerState.nChr].find(b => {
				return (
					b.start <= map_pos_ref1[a.pos] && map_pos_ref1[a.pos] <= b.end ||
					map_ref1_pos[b.start] <= a.pos && a.pos <= map_ref1_pos[b.end]
				);
			});
			return {
				pos: a.pos,
				// seq: seq_list.map(ss => ss.slice(a.pos - range, a.pos + range)).join("\n") + "\n" + " ".repeat(range) + "^",
				// p0s: seq_list.map(ss => ss[a.pos]).join(""),
				gc_wnd: gc_wnd,
				type: grp.type,
			};
		}).filter(a => a.gc > block_min_gc);//.forEach(a => console.log(a.pos, a.gc) + console.log(a.seq))
	});

	let fin_output_text = result.map(grp => {
		return grp.map(mark => {
			const pos = mark.pos;
			const ref1_pos = map_pos_ref1[mark.pos];
			const view = seq_list.map(ss => ss.slice(a.pos - view_range, a.pos + view_range)).join("\n") + "\n" + " ".repeat(view_range) + "^";
			const gc_value = mark.gc_wnd.gc;
			const gc_ref1_start = mark.gc_wnd.start;
			const gc_ref1_end = mark.gc_wnd.end;
			const ript_type = mark.type;

			const out_seq = extract_seq(pos - extract_range, pos + extract_range, false, ref1_pos);

			const info = {
				["Ch"]: viewerState.nChr,
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

			return [output_text, out_seq];
		}).join("\n" + "-".repeat(60) + "\n");
	});

	return fin_output_text;
}

// function find_rip_gff(nChr, seq_list, gff, analysis_results, gene_type) {
// 	const seq_names = dataset.genomeNameList;
// 	const ref1_name = seq_names[0];

// 	console.log(Object.keys(analysis_results).map(k => [k, typeof analysis_results[k]]));

// 	const fm = {
// 		"A": 0,
// 		"T": 0,
// 		"CBS1-1": 0,
// 		"G": 0,
// 	};
// 	[...seq_list[0]].forEach(a => ++fm[a]);
// 	// ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06)
// 	const block_min_gc = ((fm.C + fm.G) / (fm.A + fm.T + fm.C + fm.G) - 0.06) * 100;

// 	/** @type {"SNV RIP"|"RIP"} */
// 	const rip_snp_type = "RIP";
// 	const rip_markers_list = analysis_results.allMarker.list.filter(a => a.name.indexOf(rip_snp_type) >= 0);

// 	const map_pos_ref1 = analysis_results.pos_ref1_uint32array;
// 	const map_ref1_pos = analysis_results.ref1_pos_uint32array;

// 	console.log({
// 		// block_min_gc: block_min_gc,
// 		rip_markers_list: rip_markers_list.length,
// 	});

// 	const results = [];

// 	rip_markers_list.map(grp => {
// 		return grp.values.map(marker => {
// 			const gene = gff.find(gene => {
// 				// if (gene.type != "gene" || gene.type != "CDS") {
// 				if (gene.type != gene_type) {
// 					return null;
// 				}
// 				return (
// 					gene.start <= map_pos_ref1[marker.pos] && map_pos_ref1[marker.pos] <= gene.end ||
// 					map_ref1_pos[gene.start] <= marker.pos && marker.pos <= map_ref1_pos[gene.end]
// 				);
// 			});
// 			return {
// 				pos: marker.pos,
// 				marker: marker,
// 				gene: gene,
// 				rip_type: grp.name,
// 			};
// 		}).filter(a => a && a.gene).forEach(a => {
// 			results.push(a);
// 		});
// 	});
// 	results.sort((a, b) => a.pos - b.pos);

// 	return results;
// }


// Q
// allMarker.list[4].values.push(...allMarker.list[7].values.splice(0))

// C
// allMarker.list[3].values.push(...allMarker.list[6].values.splice(0))

// IM
// allMarker.list[17].values.push(...allMarker.list[18].values.splice(0), ...allMarker.list[19].values.splice(0))

if (typeof window.version_2 == "boolean" && window.version_2 ||
	dataset.ref.toUpperCase() != "QM6a".toUpperCase()
) {
	viewerState.D8_bg_RIP_marker = false;
}

function display_normRNA() {
	window.rna_reads_func = (x, max_value) => Math.min(x, max_value) / max_value;
	viewerState.rna_reads_max_display_value = null;
	normRNAbyGeneId();
}

/**
 * @param {string[]} seq_list
 */
function make_subject_pos_to_ma(seq_list) {
	seq_list.slice(dataset.parental_list.length).forEach((ma_seq, idx) => {
		// subject_1_pos_map
		const prop = `subject_${Number(idx) + 1}_pos_map`;
		// if (window[prop]) {
		// }
		// else {
			const pos_map = chrPos_to_multialign_posMap(ma_seq);
			window[prop] = pos_map;

			// if (dataset.genome_info_list[2 + idx].chr_list[viewerState.nChr - 1].length != window[prop][window[prop].length - 1]) {
			// 	console.log({
			// 		i: dataset.genome_info_list[2 + idx].chr_list[viewerState.nChr - 1].length,
			// 		a: window[prop][window[prop].length - 1],
			// 	});
			// }
		// }
	});
}

/**
 * @param {string[]} seq_list
 */
function make_subject_ma_to_pos(seq_list) {
	seq_list.slice(dataset.parental_list.length).forEach((ma_seq, idx) => {
		// pos_subject_1_map
		const prop = `pos_subject_${Number(idx) + 1}_map`;
		if (window[prop]) {
		}
		else {
			const pos_map = multialign_to_chrPos_posMap(ma_seq);
			window[prop] = pos_map;
		}
	});
}

/**
 * @abstract
 * @deprecated
 */
class AlignmentCoord {
	constructor() {
		this.start = 0;
		this.end = 0;
	}
}

/**
 * @abstract
 * @deprecated
 */
class TelomereCoord {
	constructor() {
		this.left = new AlignmentCoord();
		this.right = new AlignmentCoord();
	}
}

class BlastnCoord {
	qs = 0;
	qe = 0;
	ss = 0;
	se = 0;
	constructor() {
		this.query = "";
		this.subject = "";

		this.identity = 0;
		this.align = 0;
		this.mismatch = 0;
		this.gap = 0;

		this.qstart = 0;
		this.qend = 0;
		this.sstart = 0;
		this.send = 0;

		this.evalue = 0;
		this.score = 0;
	}

	get q_min() {
		return Math.min(this.qstart, this.qend);
	}
	get q_max() {
		return Math.max(this.qstart, this.qend);
	}

	get s_min() {
		return Math.min(this.sstart, this.send);
	}
	get s_max() {
		return Math.max(this.sstart, this.send);
	}

	get strand() {
		if (this.sstart > this.send || this.qstart > this.qend) {
			return -1;
		}
		else if (this.sstart < this.send && this.qstart < this.qend) {
			return 1;
		}
		else {
			//debugger;
			return 0;
		}
	}

	get slen() {
		return Math.abs(this.send - this.sstart);
	}

	get qlen() {
		return Math.abs(this.qend - this.qstart);
	}

	// /**
	//  * @param {number} start
	//  */
	// dist_sstart(start) {
	// 	return Math.abs(this.sstart - start);
	// }

	// /**
	//  * @param {number} start
	//  */
	// dist_send(start) {
	// 	return Math.abs(this.send - start);
	// }

	toArray() {
		return BlastnCoord.tableHeader.map(k => this[k]);
	}

	toString() {
		return JSON.stringify(this, null, "\t");
	}

	static get tableHeader() {
		return [
			"query", "subject",
			"identity", "align",
			"mismatch", "gap",
			"qstart", "qend", "sstart", "send",
			"evalue", "score"
		];
	}

	/**
	 * @param {Partial<BlastnCoord>} obj
	 * @returns {BlastnCoord}
	 */
	static fromObject(obj) {
		return Object.assign(new BlastnCoord(), obj);
	}

	/**
	 * @param {string} text
	 */
	static load_tab(text) {
		const tab = text.trim().split("\n").map(a => {
			const [
				query, subject,
				identity, align,
				mismatch, gap,
				qstart, qend, sstart, send,
				evalue, score
			] = a.trim().split("\t");
			return BlastnCoord.fromObject({
				query,
				subject,
				identity: Number(identity),
				align: Number(align),
				mismatch: Number(mismatch),
				gap: Number(gap),
				qstart: Number(qstart),
				qend: Number(qend),
				sstart: Number(sstart),
				send: Number(send),
				evalue: Number(evalue),
				score: Number(score),
			});
		});

		const qmax = Math.max.apply(this, tab.map(aln => {
			return Math.max(aln.qstart, aln.qend);
		}));
		const smax = Math.max.apply(this, tab.map(aln => {
			return Math.max(aln.sstart, aln.send);
		}));

		return {
			tab,
			qmax,
			smax,
		};
	}

	static async exec_blastn(query, subject_idx, { task = "megablast", evalue = "1e-5", options = "" }) {
		if (!(query && subject_idx == 0 || subject_idx == 1)) {
			throw new TypeError();
		}
		const subject = dataset.genome_info_list[subject_idx].chr_list[viewerState.nChr - 1].path;
		const outfmt = 6;
		const cmd = [
			`blastn`,
			task ? `-task ${task}` : "",
			query ? `-query "${query}"` : "",
			subject ? `-subject "${subject}"` : "",
			outfmt ? `-outfmt ${outfmt}` : "",
			evalue ? `-evalue ${evalue}` : "",
			options ? `${options}` : "",
		].join(" ");

		// const url = `//${location.hostname}:9527/blastn?cmd=${encodeURIComponent(cmd)}`;
		const url = `/blastn?cmd=${encodeURIComponent(cmd)}`;

		const text = await (await fetch(url, { "mode": "no-cors" })).text();
		push_range_marker_from_blastn(text, subject_idx);
	}
}

function push_range_marker_from_blastn(text, subject_idx) {
	const _tab = BlastnCoord.load_tab(text).tab;

	/** @type {Uint32Array} */
	const ref_pos = [ref1_pos_uint32array, ref2_pos_uint32array][subject_idx];

	// const tab = [_tab.sort((a, b) => Math.abs(b.send - b.sstart) - Math.abs(a.send - a.sstart))[0]];
	// const tab = [_tab.sort((a, b) => Math.min(a.send, a.sstart) - Math.min(b.send, b.sstart))[0]];
	const tab = _tab;

	tab.forEach(row => {
		let ss = ref_pos[row.sstart - 1];
		let se = ref_pos[row.send - 1];
		[ss, se] = [ss, se].sort((a, b) => a - b);
		viewerState.pushRangeMarker(ss, se);
	});
}
function download_all_final_table() {
	const group_ls = [
		// [
			"RAD51-WT+spo11-WT＃1", // #
			"RAD51-WT+spo11-WT＃2",
			"RAD51-WT+spo11-WT＃3", // #
		// ],
		// [
			"RAD51-WT+spo11-HA∕yf＃1",
			"RAD51-WT+spo11-HA∕yf＃2",
			"RAD51-WT+spo11-HA∕yf＃3",
		// ],
		// [
			"RAD51-3A+spo11-HA∕yf＃1", // #
			"RAD51-3A+spo11-HA∕yf＃2",
			"RAD51-3A+spo11-HA∕yf＃3",
		// ],
		// [
			"RAD51-3D+spo11-HA∕yf＃1",
			"RAD51-3D+spo11-HA∕yf＃2",
			"RAD51-3D+spo11-HA∕yf＃3",
		// ],
	];

	const div = document.createElement("div");
	document.body.append(div);
	for (let name of group_ls) {
		const url = `/${name}/${name}_final_table.txt`;
		const a = document.createElement("a");
		div.append(a);
		a.download = `${name}_final_table.txt`;
		a.href = url;
		a.innerHTML = a.download;
		a.target = "_blank";
		a.click();
	}
}

async function cap_all_CO(type = "CO", scale_len = 0.5, min_len = 70) {
	let start_in = 0;
	let end_in = 0;

	const list = type == "NCO" ? dataset.non_crossover_list : dataset.crossover_list;

	for (let ii = 0; ii < list[viewerState.nChr - 1].length; ++ii) {
		[g_chrCoord.bp_start, g_chrCoord.bp_end, start_in, end_in] = list[viewerState.nChr - 1].map(a => {
			const len = a.$snp_end_in - a.$snp_start_in + 1;
			const hlen = len >= min_len ? (len * scale_len) : (min_len * 0.5);// len * 0.5;
			const mm = (a.$snp_end_in + a.$snp_start_in) * 0.5;
			const ss = mm - hlen;
			const ee = mm + hlen;
			const range = [Math.floor(ss), Math.ceil(ee), a.$snp_start_in + 1, a.$snp_end_in + 1];
			return range;
		})[ii]
		while (viewerState.popRangeMarker());
		viewerState.pushRangeMarker(start_in, end_in == start_in ? start_in + 1 : end_in);
		drawFrame();
		await capturePlot(486, "_" + type + ii);
	}
	async function capturePlot(height = 486, name) {
			viewerState.get_plot_height = () => height;
			// document.querySelector("#canvas").style.height = height + "px";
			document.querySelector("#canvas").height = height;
			await captureScreen(name);
	}
}

function* see_CO(type = "CO", scale_len = 0.5, min_len = 70, out_scale = 5, clear = false) {
	const list = type == "NCO" ? dataset.non_crossover_list : dataset.crossover_list;

	// ii < list[viewerState.nChr - 1].length
	for (let ii = 0; ; ++ii) {
		ii = ii % list[viewerState.nChr - 1].length;
		if (clear) {
			while (viewerState.popRangeMarker());
		}
		yield _see_CO(list, ii, min_len, scale_len, true);
		yield _see_CO(list, ii, min_len * out_scale, scale_len * out_scale, false);
	}

	/**
	 * @param {*[]} list
	 * @param {number} ii
	 * @param {number} min_len
	 * @param {number} scale_len
	 * @param {boolean} rangeMarker
	 */
	function _see_CO(list, ii, min_len, scale_len, rangeMarker) {
		let start_in = 0;
		let end_in = 0;

		[g_chrCoord.bp_start, g_chrCoord.bp_end, start_in, end_in] = list[viewerState.nChr - 1].map(a => {
			const len = a.$snp_end_in - a.$snp_start_in + 1;
			const hlen = len >= min_len ? (len * scale_len) : (min_len * 0.5); // len * 0.5;
			const mm = (a.$snp_end_in + a.$snp_start_in) * 0.5;
			const ss = mm - hlen;
			const ee = mm + hlen;
			const range = [Math.floor(ss), Math.ceil(ee), a.$snp_start_in + 1, a.$snp_end_in + 1];
			return range;
		})[ii];

		el_input_start.value = String(g_chrCoord.bp_start);
		el_input_end.value = String(g_chrCoord.bp_end);

		if (rangeMarker) {
			viewerState.pushRangeMarker(start_in, end_in == start_in ? start_in + 1 : end_in);
		}
		drawFrame();
	}
}

/**
 * @see {@link drawTempData_afterMethyl}
 */
async function merge_ssDNA_parental_f_m() {
	/** @type {{ [name:string]: string }[]} */
	const loaded_fa = await (async function () {
		const fa_tasks = dataset.results.map(async function (obj, chrIdx) {
			if (typeof obj == "string") {
				let raw_fa = await fetchData(obj, "text");
				let fa = await parseFasta(raw_fa);
				dataset.results[chrIdx] = fa;
				return fa;
			}
			else {
				return dataset.results[chrIdx];
			}
		});
		return await Promise.all(fa_tasks);
	})();

	/** @type {Uint32Array[]} */
	let rad51_1_f;
	/** @type {Uint32Array[]} */
	let rad51_2_f;
	/** @type {Uint32Array[]} */
	let rad51_3_f;
	/** @type {Uint32Array[]} */
	let rad51_1_m;
	/** @type {Uint32Array[]} */
	let rad51_2_m;
	/** @type {Uint32Array[]} */
	let rad51_3_m;
	/** @type {Uint32Array[]} */
	let sae2_1_f;
	/** @type {Uint32Array[]} */
	let sae2_2_f;
	/** @type {Uint32Array[]} */
	let sae2_3_f;
	/** @type {Uint32Array[]} */
	let sae2_1_m;
	/** @type {Uint32Array[]} */
	let sae2_2_m;
	/** @type {Uint32Array[]} */
	let sae2_3_m;

	// rad51_1_f
	// rad51_2_f
	// rad51_3_f
	// rad51_1_m
	// rad51_2_m
	// rad51_3_m
	// sae2_1_f
	// sae2_2_f
	// sae2_3_f
	// sae2_1_m
	// sae2_2_m
	// sae2_3_m

	await Promise.all([
		loadChrUint32("rad51_f", `ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_1.uint32`).then(a => rad51_1_f = a),
		loadChrUint32("rad51_f", `ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_2.uint32`).then(a => rad51_2_f = a),
		loadChrUint32("rad51_f", `ssDNA/20210623_ssDNA_mapto/rad51_f_rad51_3.uint32`).then(a => rad51_3_f = a),
		loadChrUint32("rad51_m", `ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_1.uint32`).then(a => rad51_1_m = a),
		loadChrUint32("rad51_m", `ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_2.uint32`).then(a => rad51_2_m = a),
		loadChrUint32("rad51_m", `ssDNA/20210623_ssDNA_mapto/rad51_m_rad51_3.uint32`).then(a => rad51_3_m = a),
		loadChrUint32("sae2_f", `ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_1.uint32`).then(a => sae2_1_f = a),
		loadChrUint32("sae2_f", `ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_2.uint32`).then(a => sae2_2_f = a),
		loadChrUint32("sae2_f", `ssDNA/20210623_ssDNA_mapto/sae2_f_sae2_3.uint32`).then(a => sae2_3_f = a),
		loadChrUint32("sae2_m", `ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_1.uint32`).then(a => sae2_1_m = a),
		loadChrUint32("sae2_m", `ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_2.uint32`).then(a => sae2_2_m = a),
		loadChrUint32("sae2_m", `ssDNA/20210623_ssDNA_mapto/sae2_m_sae2_3.uint32`).then(a => sae2_3_m = a),
	]);

	const all_ma_len = loaded_fa.reduce((aa, chr_fa) => aa + Object.values(chr_fa)[0].length, 0);

	const c_rad51_1 = new Uint32Array(all_ma_len);
	const c_rad51_2 = new Uint32Array(all_ma_len);
	const c_rad51_3 = new Uint32Array(all_ma_len);
	const c_sae2_1 = new Uint32Array(all_ma_len);
	const c_sae2_2 = new Uint32Array(all_ma_len);
	const c_sae2_3 = new Uint32Array(all_ma_len);

	loaded_fa.reduce((file_pos, chr_fa, chr_idx) => {
		const all_ref_list = [
			...dataset.parental_list,
			...dataset.progeny_list,
		];
		// ["spo11rad51_m", "spo11rad51_f", "spo11sae2_m", "spo11sae2_f"];

		/** @type {{ [ref: string]: { seq: string, pos_map: Uint32Array; } }} */
		const chr_data = Object.fromEntries(all_ref_list.map((ref, refIdx) => {
			if (["rad51_m", "rad51_f", "sae2_m", "sae2_f"].includes(ref)) {
				const ma_seq = Object.values(chr_fa)[refIdx];
				return [
					ref,
					{
						seq: ma_seq,
						// pos_map: seq_list.map(ma_seq => multialign_to_chrPos_posMap(ma_seq)),
						pos_map: chrPos_to_multialign_posMap(ma_seq),
					},
				];
			}
		}).filter(a => a));

		const chr_ma_len = chr_data.rad51_m.seq.length;
		// const valid_map = new Uint32Array(chr_ma_len);
		// valid_map.map((_, pos) => fa_entries.every(a => a.seq[pos] != "-") ? 1 : 0);

		// const c_rad51_1 = new Uint32Array(chr_ma_len);
		// const c_rad51_2 = new Uint32Array(chr_ma_len);
		// const c_rad51_3 = new Uint32Array(chr_ma_len);
		// const c_sae2_1 = new Uint32Array(chr_ma_len);
		// const c_sae2_2 = new Uint32Array(chr_ma_len);
		// const c_sae2_3 = new Uint32Array(chr_ma_len);

		const seq_list = Object.values(chr_data).map(a => a.seq);
		for (let pos = 0; pos < chr_ma_len; ++pos) {
			if (seq_list.every(a => a[pos] != "-")) {
				const rad51_m_pos = chr_data.rad51_m.pos_map[pos];
				const rad51_f_pos = chr_data.rad51_f.pos_map[pos];
				const sae2_m_pos = chr_data.sae2_m.pos_map[pos];
				const sae2_f_pos = chr_data.sae2_f.pos_map[pos];

				c_rad51_1[file_pos + pos] = (rad51_1_m[chr_idx][rad51_m_pos] + rad51_1_f[chr_idx][rad51_f_pos]) * 0.5;
				c_rad51_2[file_pos + pos] = (rad51_2_m[chr_idx][rad51_m_pos] + rad51_2_f[chr_idx][rad51_f_pos]) * 0.5;
				c_rad51_3[file_pos + pos] = (rad51_3_m[chr_idx][rad51_m_pos] + rad51_3_f[chr_idx][rad51_f_pos]) * 0.5;
				c_sae2_1[file_pos + pos] = (sae2_1_m[chr_idx][sae2_m_pos] + sae2_1_f[chr_idx][sae2_f_pos]) * 0.5;
				c_sae2_2[file_pos + pos] = (sae2_2_m[chr_idx][sae2_m_pos] + sae2_2_f[chr_idx][sae2_f_pos]) * 0.5;
				c_sae2_3[file_pos + pos] = (sae2_3_m[chr_idx][sae2_m_pos] + sae2_3_f[chr_idx][sae2_f_pos]) * 0.5;
			}
		}

		return file_pos + chr_ma_len;
	}, 0);

	new DownloadHelper(c_rad51_1, "application/octet-stream").download("rad51_1");
	new DownloadHelper(c_rad51_2, "application/octet-stream").download("rad51_2");
	new DownloadHelper(c_rad51_3, "application/octet-stream").download("rad51_3");
	new DownloadHelper(c_sae2_1, "application/octet-stream").download("sae2_1");
	new DownloadHelper(c_sae2_2, "application/octet-stream").download("sae2_2");
	new DownloadHelper(c_sae2_3, "application/octet-stream").download("sae2_3");
}

/**
 * @param {number} n
 * @param {number} k
 * @returns {Generator<Number[], void>}
 */
 function* comb_C(n, k) {
	if (n < k || k <= 0) {
		return;
	}

	const comb = new Array(k).fill(0).map((v, i) => i);

	yield [...comb];//clone

	while (next_comb()) {
		yield [...comb];//clone
	}

	// do {
	// 	console.log(comb.map(v => v + 1));
	// }
	// while (next_comb());

	function next_comb() {
		let i = k - 1;
		const e = n - k;
		do {
			comb[i]++;
		}
		while (comb[i] > e + i && i--);

		if (comb[0] > e) {
			return false;
		}
		else {
			while (++i < k) {
				comb[i] = comb[i - 1] + 1;
			}
			return true;
		}
	}
}

/**
 * mappable genome size or effective genome size which is defined as the genome size
 * simple remove AT-rich island
 * @see {@link https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md#:~:text=on%20the%20chromosomes%2C-,the%20actual%20mappable%20genome%20size%20will%20be%20smaller%20than%20the%20original%20size%2C%20about%2090%25%20or%2070%25%20of%20the%20genome%20size,-.%20The%20default%20hs| MACS peak-caller}
 * @returns {number}
 */
function get_genome_size_without_AT_rich(specific_seq_len) {
	const genome_size = dataset.genome_info_list[0].chr_list.map(a => a.length).reduce((aa, v) => aa + v, 0);

	const arr = new Uint8Array(genome_size);
	arr.fill(1);

	const min_gc_content = dataset?.min_gc_content?.[0] ?? 6;
	const at_list = ["QM6a"].map((ref, refIdx) => {
		return Object.values(gc_content[ref]).flat(1).filter(a => {
			return a.gc <= min_gc_content;
		}).map(a => {
			arr.fill(0, a.start, a.end);
			return a.end - a.start + 1;
		});
	});//.flat(1);
	// const total_AT_rich_len = at_list.reduce((aa, v) => aa + v, 0);

	get_ref1_strain_specific(specific_seq_len).forEach(([start, end]) => {
		arr.fill(0, start, end);
	});

	// return genome_size - total_AT_rich_len;
	return arr.reduce((aa, v) => aa + v, 0);
}

/**
 * window size = ?
 *
 * @param {number} marker_pos
 * @param {number} start
 * @param {number} end
 * @param {number} kmer_size
 */
function kmer_analysis(marker_pos = (g_chrCoord.bp_start + g_chrCoord.bp_end) >> 1, start = g_chrCoord.bp_start, end = g_chrCoord.bp_end, kmer_size = 3) {
	/**
	 * @param {Map} a
	 * @param {Map} b
	 */
	var diff_map = (a, b) => {
		for (let [k, _] of a) {
			if (!b.has(k)) {
				return k;
			}
		}
		return true;
	};
	/**
	 * @param {Map} a
	 * @param {Map} b
	 */
	var compare_map = (a, b) => {
		if (a.size != b.size) return false;
		for (let [k, _] of a) {
			if (!b.has(k)) {
				return false;
			}
		}
		return true;
	};
	var mmm = seq_list.map(a => a.slice(start - 1, end)).map(ss => {
		ss = ss.replace(/-/g, "");
		// return new Set([...ss.slice(0, -(kmer_size - 1))].map((_, i) => ss.slice(i, i + kmer_size)));
		const m = new Map();
		[...ss.slice(0, -(kmer_size - 1))].map((_, i) => ss.slice(i, i + kmer_size)).forEach(mer => {
			m.set(mer, (m.get(mer) ?? 0) + 1);
		});
		return m;
	});

	var mmm_entries = mmm.map(m => [...m.entries()]).sort((a, b) => b[1] - a[1]).map(a => a[0]);
	if (!mmm_entries.every(a => mmm_entries[0][0] == a[0])) {
		console.error("not match kmer");
	}
	// if (MathEx.TTest([d, d1, d2], [m, m1, m2], 2, 3) = p_value) {
	// 	// is "2:2"
	// }

	const wnd_info = mmm_entries.sort((a, b) => b[1] - a[1])[0];
	// Math.max(...mmm.map(m => Math.max(...m.values())));

	wnd_info[1] += 2;
	const [mer, wnd_size] = wnd_info;

	console.log({
		mer,
		wnd_size
	});

	// const wnd_size = 100;

	mmm = seq_list.map(a => a.slice(marker_pos - wnd_size, marker_pos + wnd_size)).map(ss => {
		ss = ss.replace(/-/g, "");
		// return new Set([...ss.slice(0, -(kmer_size - 1))].map((_, i) => ss.slice(i, i + kmer_size)));
		const m = new Map();
		[...ss.slice(0, -(kmer_size - 1))].map((_, i) => ss.slice(i, i + kmer_size)).forEach(mer => {
			m.set(mer, (m.get(mer) ?? 0) + 1);
		});
		return m;
	});

	var p_rs = mmm.slice(0, 2).map(v => compare_map(v, mmm[0]) ? 0 : (compare_map(v, mmm[1]) ? 1 : NaN)).reduce((aa, v) => aa + v, 0)
	var rs = mmm.slice(2, 6).map(v => compare_map(v, mmm[0]) ? 0 : (compare_map(v, mmm[1]) ? 1 : NaN)).reduce((aa, v) => aa + v, 0)
	return rs == 2 || (p_rs == 0 && rs == 0);
}
// dataset.non_crossover_list[viewerState.nChr - 1].filter(a => !a.why_remove && kmer_analysis(a.$pos, a.$snp_start_out, a.$snp_end_out))
// [...allMarker.map["sss_13"].values, ...allMarker.map["sss_31"].values].filter(a => a.pos >= 15599 && a.pos <= 601595 && kmer_analysis(a.pos, a.pos - 50, a.pos + 50)).forEach(a => a.hide = true)

promise_load_task = null;// unlock layout


const rad51_seq = "ATGAGCGACGAGTATGACGAGGAGAACCAGGTGGCCGAGGAGGGCGGCATGACTGGCCCTGGAGCGCCGACGCCGCTCTCTGCTCTGGAGGGAGTTGCTGGATTGACAAAGCGCGACATCCAACTCGTCGTTGATGGCGGATTCAACACGGTCGAGTCGGTGGCTTACACCCCGCGCAGGGTGCTGGAGCAGATCAAGGGCATCTCAGAGCAGAAGGCGACCAAGATCTTGGCCGAGGCGTCAAAACTTGTGCCCATGGGGTTCACGACAGCCACTGAAATGCACCAGCGGCGAAGTGAGCTCATCTCCATTACCACTGGATCCAAGAACCTAGACACACTCCTGGCTGGAGGCATTGAAACGGGCTCCGTTACGGAGCTGTTTGGAGAGTTCAGGACAGGAAAGAGTCAGATCTGCCACACGCTGGCTGTGACGTGCCAGCTGCCTTTCGACATGGGCGGTGGTGAAGGCAAATGCCTGTACATTGACACCGAGGGTACCTTTCGACCCGTCCGACTGCTGGCCGTTGCCAATCGATTTGGGCTGTCTGGTGAAGAAGTCCTCGACAATGTCGCGTATGCGAGAGCGTACAACTCAGACCACCAGCTTCAGCTGCTGAACCAGGCAGCGGCCATGATGTGCGAGACAAGGTTTTCCCTGCTCATCGTCGACAGCGCTACTTCACTCTACCGAACGGACTTTACCGGCCGAGGTGAACTCTCGAATCGTCAGACACATTTGGCCAAGTTTATGAGGACACTGCAGCGGCTCGCAGACGAGTTCGGCATTGCCGTCGTCATCACCAACCAGGTTGTCGCGCAGGTCGATGGCGGACCGAGTGCCATGTTCAACCCTGATCCGAAGAAGCCCATTGGCGGCAATATTATTGCGCATGCCAGCACAACTCGAATCAGTTTGAAGAAGGGACGTGGAGAGACTCGAATCGCCAAGATTTACGATAGCCCTTGTCTGCCGGAGAGCGACACGCTGTTTGCCATTGGCGAGGACGGTATTGGTGACCCGGCGCCAAAGGACTTGGAGAAGGAGAAGGACTGA";

/** hygromycin */
const hph_seq = "ATGAAAAAGCCTGAACTCACCGCGACGTCTGTCGAGAAGTTTCTGATCGAAAAGTTCGACAGCGTCTCCGACCTGATGCAGCTCTCGGAGGGCGAAGAATCTCGTGCTTTCAGCTTCGATGTAGGAGGGCGTGGATATGTCCTGCGGGTAAATAGCTGCGCCGATGGTTTCTACAAAGATCGTTATGTTTATCGGCACTTTGCATCGGCCGCGCTCCCGATTCCGGAAGTGCTTGACATTGGGGAGTTCAGCGAGAGCCTGACCTATTGCATCTCCCGCCGTGCACAGGGTGTCACGTTGCAAGACCTGCCTGAAACCGAACTGCCCGCTGTTCTCCAGCCGGTCGCGGAGGCCATGGATGCGATCGCTGCGGCCGATCTTAGCCAGACGAGCGGGTTC";

async function load_JQueryUI() {
	if (window.jQuery != null && window.jQuery?.ui != null) {
		return;
	}
	await new Promise((resolve, reject) => {
		const jq = document.createElement("script");
		jq.onload = resolve;
		jq.onerror = reject;
		jq.src = "https://code.jquery.com/jquery-3.1.0.js";
		document.body.append(jq);
	});
	const tasks = Promise.all([
		new Promise((resolve, reject) => {
			const link = document.createElement("link");
			link.onload = resolve;
			link.onerror = reject;
			link.href = "https://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css";
			link.rel = "stylesheet";
			link.type = "text/css";
			document.body.append(link);
		}),
		new Promise((resolve, reject) => {
			const jui = document.createElement("script");
			jui.onload = resolve;
			jui.onerror = reject;
			jui.src = "https://code.jquery.com/ui/1.12.1/jquery-ui.js";
			document.body.append(jui);
		}),
	]);
	return await tasks;
}

async function download_ref1_indel_RIP() {
	const tasks = dataset.genome_info_list[0].chr_list.map(async function (_, chr_idx) {
		const nChr = chr_idx + 1;
		const list = JSON.parse(await (await fetch(`http://localhost:9099/20200720_v3_QCt/20200720_v3_QCt_ch${nChr}_ref1_2_rip.json`)).text())
		return list.map(a => [
			nChr,
			a.ref1_pos,
			a.seg[0],
		].join("\t")).join("\r\n")
	});
	const s = (await Promise.all(tasks)).join("\r\n");
	downloadTextFile("QM6a InDel RIP.txt", s);
}

async function download_ref1_strain_specific_A(specific_seq_len = 1000) {
	const ag = gen_ref1_strain_specific_A(specific_seq_len);
	let list = [];
	for await (let block of ag) {
		list = list.concat(block);
	}
	const s = list.join("\r\n");
	downloadTextFile("QM6a strain specific (min length 1kb).txt", s);
}
async function* gen_ref1_strain_specific_A(specific_seq_len = 1000) {
	for (let nChr = 1; nChr <= dataset.results.length; ++nChr) {
		if (dataset.results[nChr - 1] == null || typeof dataset.results[nChr - 1] == "string") {
			dataset.results[nChr - 1] = parseFasta(await fetchData(`mafft_ch${nChr}.fa`, "text"));
		}
		const fa = dataset.results[nChr - 1];

		/** @type {string[]} */
		const seq_list = Object.values(fa);

		const pos_ref1_uint32array = multialign_to_chrPos_posMap(seq_list[0]);
		// const ref1_pos_uint32array = chrPos_to_multialign_posMap(seq_list[0]);

		const ref1_sp_map = [...seq_list[0]].map((ref1, idx) => {
			if (ref1 != "-") {
				if (seq_list[1][idx] == "-") {
					return "-";
				}
			}
			return "0";
		});

		const regexp = new RegExp(`-{${specific_seq_len},}`, "g");
		// sp: strain-specific
		const sp = [
			...(ref1_sp_map.join("").matchAll(regexp))
		].map(a => {
			const start = a.index;
			const end = start + a[0].length + 1;
			// return {
			// 	chr: nChr,
			// 	start: pos_ref1_uint32array[start],
			// 	end: pos_ref1_uint32array[end],
			// };
			return [
				nChr,
				pos_ref1_uint32array[start],
				pos_ref1_uint32array[end],
			].join("\t");
		}).join("\r\n");

		yield sp;
	}
}

async function mark_methyl_point() {
	const repeat_rangeList = await get_repeat_rangeList();

	const nChr = viewerState.nChr;
	const chr_idx = nChr - 1;

	// const row = new module_Methyl_sampleData({
	// 	ref: dataset.ref,
	// 	sample: "QM6a",
	// 	name: "repeat",
	// 	// value_desc: "(M + H) / (M + H + C)",
	// 	// html_value_desc: make_methyl_value_desc("Cm + Chm", "Cm + Chm + C"),
	// 	url: null,
	// 	region: false,
	// 	mid_line: false,
	// 	density_to_opacity: false,
	// 	data: [],
	// 	rendering_condition: [
	// 		new module_Methyl_sampleData_RenderingCondition({ color: "#04F7", condition: v => v.col == 0b0010, min_width: 1, }),
	// 		new module_Methyl_sampleData_RenderingCondition({ color: "#f4F7", condition: v => v.col == 0b0011, min_width: 1, }),
	// 		new module_Methyl_sampleData_RenderingCondition({ color: "#F007", condition: v => v.col == 0b0001, min_width: 1, }),
	// 	],
	// });
	// row[chr_idx] = repeat_rangeList.map(b => {
	// 	if (b.value) {
	// 		return {
	// 			start: b.start,
	// 			end: b.end,
	// 			value: 1,
	// 			// strand: 1,
	// 			col: b.value,
	// 		};
	// 	}
	// }).filter(a => a)
	// methyl_dataset_list.push(row);

	/** @see {@link load_methly_preset_20210903} */
	const TAG_ZERO_READ = 0;

	methyl_dataset_list.forEach(md => {
		md.data[chr_idx].forEach(a => {
			if (a.t != TAG_ZERO_READ) {
				const b = repeat_rangeList.find(b => b.start <= a.start && b.end >= a.start);
				if (b) {
					a.col = b.value;
					return;
				}
			}
			delete a.col;
		});
		// md.rendering_condition.splice(0);
		// md.rendering_condition.splice(-3, 3,
		md.rendering_condition.push(
			// new module_Methyl_sampleData_RenderingCondition({ color: "#444", condition: v => v.col == 0b0000, min_width: 1, }),
			new module_Methyl_sampleData_RenderingCondition({ color: "#04F7", condition: v => v.col == 0b0010, min_width: 1, }),
			new module_Methyl_sampleData_RenderingCondition({ color: "#f4F7", condition: v => v.col == 0b0011, min_width: 1, }),
			new module_Methyl_sampleData_RenderingCondition({ color: "#F007", condition: v => v.col == 0b0001, min_width: 1, }),
		);
	});
}

async function captureRef1ChrRange(s, e) {
	g_chrCoord.bp_start = ref1_pos_uint32array[s];
	g_chrCoord.bp_end = ref1_pos_uint32array[e];

	viewerState.resizeCanvas()

	await drawFrame();

	await captureScreen("_" + [s, e].join("-"));
}
// 20211007
// await captureRef1ChrRange(398512, 415884)
// await captureRef1ChrRange(2564986, 2578317)

function asdas() {
	aaa_1 = new Uint32Array(seq_list[1].length)
	for (let i = 0, j = 0; i < seq_list[1].length; ++i) {
		if (seq_list[1][i] != "-") {
			// idx to idx
			aaa_1[i] = j++
		}
	}
	[aaa_1[404249], subject_1_pos_map[404080]]
	seq_list[1].replace(/-/g, "").slice(aaa_1[404248] - 1, aaa_1[404250])// AAG

	ss.replace(/-/g, "").slice(aaa_2[ma_seq_start_idx], aaa_2[ma_seq_end_idx] + 1)// AGG
}


function simple_find_RIP(ref_name = "QM6a", target_list = ["Qdim", "Qrid", "np128_Qdimrid"], save_tsv = true) {
	const bp_start_idx = g_chrCoord.bp_start - 1;
	const ref_idx = dataset.genome_info_list.findIndex(a => a.name == ref_name);
	const ref_seq = seq_list[ref_idx].slice(bp_start_idx, g_chrCoord.bp_end);

	const target_entries = seq_list.map((full_seq, i) => {
		const target_name = dataset.genome_info_list[i].name;
		if (target_list.includes(target_name)) {
			const seq = full_seq.slice(bp_start_idx, g_chrCoord.bp_end);
			return {
				idx: i,
				target_name,
				seq,
			};
		}
	}).filter(a => a != null);

	class RIP_Data {
		seq_idx = -1;
		
		pos = 0;
		ref_pos = 0;

		/** @type {"C->T"|"G->A"} */
		type = null;
	}

	/** @type {RIP_Data[]} */
	const rip_list = [];

	// 20220627
	const pos_map_list = make_pos_ref_map_list();

	for (let bp_idx = 0; bp_idx < ref_seq.length; ++bp_idx) {
		const ref_bp = ref_seq[bp_idx];
		
		const ma_pos = bp_start_idx + bp_idx + 1;
		const ref_pos = pos_map_list[ref_idx][ma_pos - 1];

		target_entries.forEach(target => {
			const pos_map = pos_map_list[target.idx];
			const target_pos = pos_map[ma_pos - 1];

			const bp = target.seq[bp_idx];
			if (ref_bp != "-" && ref_bp != bp) {
				if (ref_bp == "G" && bp == "A") {
					rip_list.push({
						seq_idx: target.idx,
						pos: target_pos,
						type: "G->A",
						ref_pos: ref_pos,
					});
				}
				else if (ref_bp == "C" && bp == "T") {
					rip_list.push({
						seq_idx: target.idx,
						pos: target_pos,
						type: "C->T",
						ref_pos: ref_pos,
					});
				}
			}
		});
	}

	if (save_tsv) {
		let tsv = "";

		const raw_chr_name = dataset.genome_info_list[ref_idx].chr_list[viewerState.nChr - 1].raw_chr_name;
		const ref_start = pos_map_list[ref_idx][g_chrCoord.bp_start - 1];
		const ref_end = pos_map_list[ref_idx][g_chrCoord.bp_end - 1];
		tsv += `${raw_chr_name}:${ref_start}..${ref_end}\n`;

		/** @type {RIP_Data[][]} */
		const rip_group = dataset.genome_info_list.map(a => []);
		rip_list.forEach(v => rip_group[v.seq_idx].push(v));

		rip_group.filter(a => a?.length > 0).forEach(vv => {
			tsv += [
				`${dataset.genome_info_list[ref_idx].name} pos`,
				`${dataset.genome_info_list[vv[0].seq_idx].name} pos`,
				"RIP type",
			].join("\t") + "\n";
	
			tsv += vv.map(v => {
				return [
					v.ref_pos,
					v.pos,
					v.type,
				].join("\t");
			}).join("\n") + "\n";
		});

		document.querySelector("textarea").value = tsv;
	}

	return rip_list;
}

class IPresetLoader {
	init() {
	}
	async load() {
	}
}

/**
 * 20221025
 * @alias load_presset_20220916
 * 
 * @param {"Top2"|"REC8"|"MCD1"|"gliP"} gene_name
 * @param {("WT"|"CBS1-1"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid")[]} ref_list
 * @param {number} ext 1000
 * @example
 * await load_presset_20220916("Top2", ["WT","Qdim","Qrid","np128_Qdimrid"])
 * await load_presset_20220916("REC8", ["WT","Qdim","Qrid","np128_Qdimrid"])
 * await load_presset_20220916("MCD1", ["WT","Qdim","Qrid","np128_Qdimrid"])
 * await load_presset_20220916("gliP", ["WT","Qdim","Qrid","np128_Qdimrid"])
 * 
 * await load_presset_20220916("MCD1", ["np128_Qdimrid"])
 * 
 * await load_presset_20220916("Top2", ["WT", "Cdim", "Crid", "np59_Cdimrid"])
 * await load_presset_20220916("REC8", ["WT", "Cdim", "Crid", "np59_Cdimrid"])
 * await load_presset_20220916("MCD1", ["WT", "Cdim", "Crid", "np59_Cdimrid"])
 * await load_presset_20220916("gliP", ["WT", "Cdim", "Crid", "np59_Cdimrid"])
 */
async function load_presset_20220916(gene_name, ref_list, ext = 1000) {
	class _PresetData {
		title = "";
		nChr = "";
		start = 0;
		end = 0;
	}
	/** @type {{ [gene_name: string]: _PresetData; }} */
	const _presets = {
		"Top2": { title: "Top2", nChr: "3", start: 2737345, end: 2744153, },
		"REC8": { title: "REC8", nChr: "4", start: 2318627, end: 2321939, },
		"MCD1": { title: "MCD1", nChr: "5", start: 2244735, end: 2247595, },
		"gliP": { title: "gliP", nChr: "6", start: 1366081, end: 1372952, },

		
		"ChV 398991-447710": { title: "ChV 398991-447710", nChr: 5, start: 398991, end: 447710, },
		"ChV usk1..cel61a": { title: "ChV usk1-SOR-BGC-axe1-cip1-cel61a", nChr: 5, start: 2173501, end: 2222743, },
		"ChVI GTX": { title: "ChVI GTX gene cluster", nChr: 6, start: 1343511, end: 1378812, },
		"ChVII msh4": { title: "ChVII msh4", nChr: 7, start: 3079958, end: 3083753, },
	};

	const preset = _presets[gene_name];
	if (!preset) {
		throw new Error("load_presset_20220916:", gene_name);
	}

	await loadTSETA_ChrData(preset.nChr);

	// const ext = 1000;
	return await _load_presset_20220916(preset.title, preset.start, preset.end, ext, ref_list);
}

async function loadTSETA_ChrData(nChr) {
	if (seq_list.length == 0 || !seq_list.every(ss => ss.length > 0)) {
		viewerState.nChr = nChr;
		el_input_chr.value = viewerState.nChr;
		await loadData(true);
		await drawFrame();// trigger GC
	}
}

/**
 * 20221025
 * @param {string} preset_name
 * @param {("WT"|"CBS1-1"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid")[]} ref_list no use
 */
async function load_presset_20221025_QM6a_RNA(preset_name, ref_list, ext = 10000) {
	const preset = {
		"ChV 398991-447710 RNA": {
			title: "ChV 398991-447710",
			start: 398991,
			end: 447710,
			nChr: 5,
		},
		"ChV usk1..cel61a RNA": {
			title: "ChV usk1-SOR-BGC-axe1-cip1-cel61a",
			start: 2173501,
			end: 2222743,
			nChr: 5,
		},
		"ChVI GTX RNA": {
			title: "ChVI GTX gene cluster",
			start: 1343511,
			end: 1378812,
			nChr: 6,
		},
		"ChVII msh4 RNA": {
			title: "ChVII msh4 RNA",
			start: 3079958,
			end: 3083753,
			nChr: 7,
		},
	}[preset_name];

	await _load_presset_20221025_QM6a_RNA(preset.title, preset.start, preset.end, preset.nChr , ext);
}

/**
 * @param {string} title
 * @param {number} start
 * @param {number} end
 * @param {number} nChr
 * @param {number} ext
 */
async function _load_presset_20221025_QM6a_RNA(title, start, end, nChr, ext) {
	const ref_idx = 0;

	await loadTSETA_ChrData(nChr);

	{
		current_colorset.mom = "#FF0000";
		current_colorset.dad = "#0000FF";

		viewerState.$display_SNP_density = false;

		_presset_20220916_setViewTarget(start, end, ext);

		_presset_20220916_load_gene(ref_idx);

		const title_text = [
			`QM6a`,
			title,
			`±${Math.trunc(ext / 1000)}kb`
		].filter(a => a).join(" ");
		viewerState.setPlotTitle(title_text, title_text);

		dataset.display_name = {
			"QM6a": "QM6a",
			"CBS1-1": "CBS1-1",

			"np43_QM6a": "QM6a WT (D-np43)",
			"np42_CBS1-1": "CBS1-1 WT (M-np42)",

			"np30": "QM6a dim2Δ (D-np30)",
			"np31": "CBS1-1 rid1Δ (M-np31)",
			"np32": "QM6a rid1Δ (D-np32)",
			"np128_Qdimrid": "dim2Δrid1Δ (D-np128)",
			"np59_Cdimrid": "dim2Δrid1Δ (M-np59)",

			"Qdim": "QM6a dim2Δ (D-np30)",
			"Crid": "CBS1-1 rid1Δ (M-np31)",
			"Qrid": "QM6a rid1Δ (D-np32)",
		};
		update_SNP_Mode_genome_display_name();
	}
	await delayFrame();

	await QM6a_RNA_depth_20221024();

	const genome_size = dataset.genome_info_list[0].chr_list.reduce((p, v) => p + v.length, 0);
	methyl_dataset_list.forEach(cc => cc.max_display_value = 1 / genome_size);

	await drawFrame();
}

/**
 * Figure 3, Figure 4 and Figure 5
 * @param {string} title
 * @param {number} start
 * @param {number} end
 * @param {number} ext
 * @param {("WT"|"CBS1-1"|"Qdim"|"Cdim"|"Qrid"|"Crid"|"np128_Qdimrid"|"np59_Cdimrid")[]} ref_list
 */
async function _load_presset_20220916(title, start, end, ext, ref_list) {
	if (seq_list.length == 0 || !seq_list.every(ss => ss.length > 0)) {
		throw new Error("No data");
	}

	current_colorset.mom = "#FF0000";
	current_colorset.dad = "#0000FF";

	_presset_20220916_setViewTarget(start, end, ext);
	
	if (ref_list == null || ref_list.length == 0) {
		ref_list = [
			"WT",
			"Qdim",
			"Qrid", "Crid",
			"np128_Qdimrid", "np59_Cdimrid",
		];
	}
	else {
		const ua = new Set(ref_list);
		if (ua.size != ref_list.length) {
			console.error("ua.size", "!=", "ref_list.length");
			ref_list = [...ua];
		}
	}
	
	// setup layout
	{
		document.querySelector("#data-rows > div").style.width = "18em";

		/**
		 * @see {load_methly_preset_20210924}
		 */
		dataset.display_name = {
			"QM6a": "QM6a",
			"CBS1-1": "CBS1-1",
			
			"np43_QM6a": "QM6a WT (D-np43)",
			"np42_CBS1-1": "CBS1-1 WT (M-np42)",

			"np30": "QM6a dim2Δ (D-np30)",//QdimΔ // np30
			"np31": "CBS1-1 rid1Δ (M-np31)",//CridΔ // np31
			"np32": "QM6a rid1Δ (D-np32)",//QridΔ // np32
			"np128_Qdimrid": "dim2Δrid1Δ (D-np128)", // np128_Qdimrid
			"np59_Cdimrid": "dim2Δrid1Δ (M-np59)", // np59_CdimΔridΔ

			"Qdim": "QM6a dim2Δ (D-np30)",
			"Crid": "CBS1-1 rid1Δ (M-np31)",
			"Qrid": "QM6a rid1Δ (D-np32)",
		};

		viewerState.$display_SNP_density = false;

		await delayFrame();
		await delayFrame();
	}

	const QM6a = dataset.genome_info_list.findIndex(a => a.name == "QM6a");
	const CBS1 = dataset.genome_info_list.findIndex(a => a.name == "CBS1-1");

	const np43 = dataset.genome_info_list.findIndex(a => a.name == "np43_QM6a");
	const np42 = dataset.genome_info_list.findIndex(a => a.name == "np42_CBS1-1");

	const Qdim = dataset.genome_info_list.findIndex(a => a.name == "Qdim");
	const Qrid = dataset.genome_info_list.findIndex(a => a.name == "Qrid");
	const Crid = dataset.genome_info_list.findIndex(a => a.name == "Crid");
	const np128_Qdimrid = dataset.genome_info_list.findIndex(a => a.name == "np128_Qdimrid");
	const np59_Cdimrid = dataset.genome_info_list.findIndex(a => a.name == "np59_Cdimrid");

	if (g_methylRenderer.isLimitLoad()) {//20221102
		add_RIP_marker(QM6a, CBS1);
		
		add_RIP_marker(QM6a, np43);
		add_RIP_marker(QM6a, np42);

		add_RIP_marker(QM6a, Qdim);
		// add_RIP_marker(QM6a, np59_Cdimrid);
		// add_RIP_marker(QM6a, Qdim);
		// add_RIP_marker(QM6a, np59_Cdimrid);
		add_RIP_marker(QM6a, Qrid);
		add_RIP_marker(QM6a, Crid);
		// add_RIP_marker(QM6a, Qrid);
		// add_RIP_marker(QM6a, Crid);
		add_RIP_marker(QM6a, np128_Qdimrid);
		add_RIP_marker(QM6a, np59_Cdimrid);
		// add_RIP_marker(QM6a, np128_Qdimrid);
		// add_RIP_marker(QM6a, np59_Cdimrid);
	}

	const add_display2 = false;// first BS-seq
	await ref_list.reduce(async (prev_promise, ref_type) => {
		await prev_promise;
		await delayFrame();
		
		// await load_methly_preset_20210924(ref_type, false);

		// QM6a sequencing by pacbio
		// QM6a-NP43 sequencing by nanopore

		// ×✖✕❌

		if (ref_type == "WT") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, CBS1);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("QM6a",   "QM6a",   "BS", "WT", ["veg"]                                                          );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a",   "QM6a",   "BS", "WT", ["veg"],                        "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "CBS1-1", "BS", "WT", ["veg"]                                                          );
			//
			await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "BS", "WT", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "BS", "WT", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "BS", "WT", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "BS", "WT", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
			//
			await __load_methratioFloat32_depthUint32_by_list("QM6a", "QM6a",   "ctrl", "WT", ["veg"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a", "QM6a",   "ctrl", "WT", ["veg"], "data/methyl_20221007/dsiplay-2"); // display-group-2
			//
			await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "ctrl", "WT", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "ctrl", "WT", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "ctrl", "WT", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a", "D × M",   "ctrl", "WT", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
		}
		else if (ref_type == "CBS1-1") {
			await __load_methratioFloat32_depthUint32_by_list("QM6a",   "QM6a",   "BS", "WT", ["veg"]                                                          );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("QM6a",   "QM6a",   "BS", "WT", ["veg"],                        "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "CBS1-1", "BS", "WT", ["veg"]                                                          );
			//
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "BS", "WT", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "BS", "WT", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "BS", "WT", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "BS", "WT", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
			//
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "CBS1-1", "ctrl", "WT", ["veg"                                                           ]);
			//
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "ctrl", "WT", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "ctrl", "WT", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "ctrl", "WT", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("CBS1-1", "D × M", "ctrl", "WT", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
		}
		else if (ref_type == "Qdim") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, Qdim);
				add_RIP_marker(QM6a, np59_Cdimrid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "Qdim",  "BS", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "Qdim",  "BS", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "D × M", "BS", "Pd", ["D2", "D4", "D5", "D6", "D8"]);
			//
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "Qdim",  "ctrl", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "Qdim",  "ctrl", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "D × M", "ctrl", "Pd", ["D2", "D4", "D5", "D6", "D8"]);
		}
		else if (ref_type == "Cdim") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, Qdim);
				add_RIP_marker(QM6a, np59_Cdimrid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "Qdim",  "BS", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "Qdim",  "BS", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "D × M", "BS", "Pd", ["D2", "D4", "D5", "D6", "D8"]);
			//
			await __load_methratioFloat32_depthUint32_by_list("Qdim",         "Qdim",  "ctrl", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "Qdim",  "ctrl", "Pd", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid", "D × M", "ctrl", "Pd", ["D2", "D4", "D5", "D6", "D8"]);
		}
		else if (ref_type == "Qrid") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, Qrid);
				add_RIP_marker(QM6a, Crid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "Qrid", "BS", "rid1", ["veg"]                                   );
			await __load_methratioFloat32_depthUint32_by_list("Crid", "Crid", "BS", "rid1", ["veg"]                                   );
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "BS", "rid1", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "BS", "rid1", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "BS", "rid1", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "BS", "rid1", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "Qrid", "ctrl", "rid1", ["veg"]                                   );
			await __load_methratioFloat32_depthUint32_by_list("Crid", "Crid", "ctrl", "rid1", ["veg"]                                   );
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "ctrl", "rid1", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "ctrl", "rid1", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "ctrl", "rid1", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Qrid", "D × M", "ctrl", "rid1", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
		}
		else if (ref_type == "Crid") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, Qrid);
				add_RIP_marker(QM6a, Crid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "Qrid", "BS", "rid1", ["veg"]                                   );
			await __load_methratioFloat32_depthUint32_by_list("Crid", "Crid", "BS", "rid1", ["veg"]                                   );
			//
			await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "BS", "rid1", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "BS", "rid1", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2//
			await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "BS", "rid1", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "BS", "rid1", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
			//
			await __load_methratioFloat32_depthUint32_by_list("Qrid", "Qrid", "ctrl", "rid1", ["veg"]                                   );
			await __load_methratioFloat32_depthUint32_by_list("Crid", "Crid", "ctrl", "rid1", ["veg"]                                   );
			//
			await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "ctrl", "rid1", ["D2", "D4"                  ]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "ctrl", "rid1", [      "D4"                  ], "data/methyl_20221007/dsiplay-2"); // display-group-2
			await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "ctrl", "rid1", [            "D5", "D6", "D8"]                                   );
			if (add_display2) await __load_methratioFloat32_depthUint32_by_list("Crid", "D × M", "ctrl", "rid1", [                        "D8"], "data/methyl_20221007/dsiplay-2"); // display-group-2
		}
		else if (ref_type == "np128_Qdimrid") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, np128_Qdimrid);
				add_RIP_marker(QM6a, np59_Cdimrid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "np128_Qdimrid", "BS", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "np59_Cdimrid",  "BS", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "D × M",         "BS", "ddr", ["D2", "D4", "D5", "D6", "D8"]);
			//
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "np128_Qdimrid", "ctrl", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "np59_Cdimrid",  "ctrl", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "D × M",         "ctrl", "ddr", ["D2", "D4", "D5", "D6", "D8"]);
		}
		else if (ref_type == "np59_Cdimrid") {
			if (!g_methylRenderer.isLimitLoad()) {//20221102
				add_RIP_marker(QM6a, np128_Qdimrid);
				add_RIP_marker(QM6a, np59_Cdimrid);
			}
			//
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "np128_Qdimrid", "BS", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "np59_Cdimrid",  "BS", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "D × M",         "BS", "ddr", ["D2", "D4", "D5", "D6", "D8"]);
			//
			await __load_methratioFloat32_depthUint32_by_list("np128_Qdimrid", "np128_Qdimrid", "ctrl", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "np59_Cdimrid",  "ctrl", "ddr", ["veg"]);
			await __load_methratioFloat32_depthUint32_by_list("np59_Cdimrid",  "D × M",         "ctrl", "ddr", ["D2", "D4", "D5", "D6", "D8"]);
		}
		
		await delayFrame();// trigger GC
		await drawFrame();
	}, Promise.resolve());
	await delayFrame();// trigger GC
	
	///////////////////////
	
	// methyl_dataset_list.forEach(cc => cc.hide = true);
	// cfg = methyl_dataset_list[1];
	methyl_dataset_list.forEach(cfg => {
		// cfg.hide = false;
		cfg.fast = false;
		cfg.fast_binSize = 10;
		cfg.density_to_opacity = false;
		cfg.func_density_to_opacity = () => 1;
		cfg.rendering_condition.forEach(v => v.min_width = 1);
		// cfg.value_normalize = 0.6; // display value: value / value_normalize
		cfg.max_display_value = 0.6;
		cfg.data[viewerState.nChr - 1].forEach(data => {
			delete data.skip;
			delete data.fin_value;
		});
	});
	// await drawFrame();// trigger GC

	const title_text = [
		`QM6a`,
		title,
		`±${Math.trunc(ext / 1000)}kb`
	].filter(a => a).join(" ");
	viewerState.setPlotTitle(title_text, title_text);

	update_SNP_Mode_genome_display_name();
	
	try {
		const avgline = new module_MethylRenderer_averageToLine(methyl_dataset_list.find(a => a.name == "QM6a dim2Δ veg +BS"));
		avgline.enable = true;
		avgline.updateAvgLine(methyl_dataset_list, g_chrCoord.bp_start - ext, g_chrCoord.bp_end + ext);

		methyl_dataset_list.filter(a => !a.url || a.name.indexOf("-BS") >= 0).forEach(a => a.afterRenderList.pop());
	}
	finally {
	}
	
	_presset_20220916_load_gene();
	await delayFrame();
	
	methyl_dataset_list.forEach(cc => {
		if (cc.name.indexOf("+BS") >= 0) {
			cc.rendering_condition.filter(rc => rc.desc == "TAG_RIP" || rc.desc == "TAG_C").forEach(rc => rc.color = "hsl(338deg 100% 50%)");// 20221102 color
		}
		else if (cc.name.indexOf("-BS") >= 0) {
			cc.rendering_condition.filter(rc => rc.desc == "TAG_RIP" || rc.desc == "TAG_C").forEach(rc => rc.color = "hsl(220deg 72% 50%)");// 20221102 color

			cc.name = cc.name + " (C-to-T)";
			cc.data[viewerState.nChr - 1].forEach(mm => {
				if (!("raw_value" in mm)) {
					mm.raw_value = mm.value;
				}
				mm.value = 1 - mm.raw_value;
			});
		}
	});

	await delayFrame();
	await drawFrame();// trigger GC

	methyl_ratio_cu_plot();
}

function _presset_20220916_load_gene(refIdx = 0) {
	// const refIdx = 0;
	const gene_row = show_gene_on_row(refIdx);

	const vn = {
		TrQ_007352: { name: "usk1",     fontStyle: "italic" },
		TrQ_007353: { name: "sor8",     fontStyle: "italic" },
		TrQ_007354: { name: "sor1",     fontStyle: "italic" },
		TrQ_007355: { name: "sor2",     fontStyle: "italic" },
		TrQ_007356: { name: "sor5",     fontStyle: "italic" },
		TrQ_007356: { name: "sor5",     fontStyle: "italic" },
		TrQ_007357: { name: "MSF/sor4", fontStyle: "italic" },
		TrQ_007358: { name: "sor4_2",   fontStyle: "italic" },
		TrQ_007359: { name: "sor3",     fontStyle: "italic" },
		TrQ_007360: { name: "sor7",     fontStyle: "italic" },
		TrQ_007361: { name: "ypr1",     fontStyle: "italic" },
		TrQ_007362: { name: "",         fontStyle: "italic" },
		TrQ_007363: { name: "axe1",     fontStyle: "italic" },
		TrQ_007364: { name: "",         fontStyle: "italic" },
		TrQ_007365: { name: "cel61a",   fontStyle: "italic" },

		TrQ_008316: { name: "T", },
		TrQ_008318: { name: "M", },
		TrQ_008323: { name: "G", },
		TrQ_008324: { name: "J", },
		TrQ_008325: { name: "K", },
		TrQ_008326: { name: "P", },
		TrQ_008327: { name: "I", },
		TrQ_008328: { name: "N", },

		TrQ_009909: { name: "msh4", },
	};
	gene_row.data[viewerState.nChr - 1].map(a => {
		const vv = vn[a.value.attributes.ID];
		a.value.attributes.Name = "";
		if (vv) {
			if (vv.fontStyle) {
				a.value.fontStyle = vv.fontStyle;
			}
			a.value.attributes.Name = vv.name;
		}
		// else {
		// 	a.value.attributes.Name = "";
		// }
	});

	methyl_dataset_list.unshift(gene_row);
}

function _presset_20220916_setViewTarget(start, end, ext1 = 1000, ext2 = 1000) {
	const e_start = start - ext1;
	const e_end = end + ext1;

	g_chrCoord.bp_start = ref1_pos_uint32array[e_start - 1] + 1;
	g_chrCoord.bp_end = ref1_pos_uint32array[e_end - 1] + 1;
	viewerState.pushRangeMarker(ref1_pos_uint32array[start - 1], ref1_pos_uint32array[end - 1])
	drawFrame();
	document.querySelector("#el_set_ref_pos_from_ma").click()

	color_user_marker_bk = "#DDFFCC7F";
	
	g_methylRenderer.setLoadingBounding(Math.max(1, g_chrCoord.bp_start - ext2), Math.min(g_chrCoord.bp_end + ext2, seq_list[0].length));
}

function methyl_ratio_cu_plot() {
	const t = methyl_dataset_list.filter(a => a.name.indexOf("+BS") >= 0).map(config => {
		const list = config.data[viewerState.nChr - 1];

		// =TRUNC((E19*100)/10*(30/10))
		const step = 50;
		const vec = Array(step);
		vec.fill(0);

		list.forEach(data => {
			// point
			const start_pos = config.ref_to_pos_map ? config.ref_to_pos_map[data.start] : data.start;//pos_map[data.start - 1];
			// const end_pos = config.ref_to_pos_map ? config.ref_to_pos_map[data.end] : data.end;//pos_map[data.end - 1];

			if (g_chrCoord.bp_start <= start_pos && start_pos <= g_chrCoord.bp_end) {
				const idx = Math.trunc((data.value ?? 0) * 100 / 10 * (step / 10));
				vec[Math.min(idx, vec.length - 1)] += 1;
			}
		});

		// vec = [
		// 	list.filter(a => a.value >=    0 && a.value < 0.05).length,
		// 	list.filter(a => a.value >= 0.05 && a.value < 0.10).length,
		// 	list.filter(a => a.value >= 0.15 && a.value < 0.20).length,
		// ];

		// a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
		// b = [a[0]];
		// a.slice(1).forEach((v, i) => b[i + 1] = b[i] + v)
		const cu_list = [vec[0]];
		vec.slice(1).forEach((v, i) => cu_list[i + 1] = cu_list[i] + v);

		const sum = vec.reduce((acc, v) => acc + v, 0);

		return [
			config.name,
			...cu_list.map(v => v / sum),
		];
	});
	document.querySelector("textarea").value = t.map(v => v.join("\t")).join("\n");
}

/**
 * 20221025
 */
async function QM6a_RNA_depth_20221024() {
	const dup_list = ["1", "2", "3"];

	// 20221102 color
	const color = {
		WT_veg: "hsl(330deg 100% 60%)",//332deg
		WT_D: "hsl(338deg 100% 50%)",
		rid1_veg: "hsl(200deg 85% 50%)",
		rid1_D: "hsl(220deg 72% 50%)",
	};
	// // find_color(x = "hsl(220deg") => methyl_dataset_list.find(c => c.rendering_condition.find(a => a.color.startsWith(x))).rendering_condition[0].color
	// const color = {
	// 	WT_veg: "hsl(33deg 100% 50%)",
	// 	WT_D: "hsl(30deg 100% 30%)",
	// 	rid1_veg: "hsl(200deg 0.85% 55%)",//hsl(200deg 0.85% 55%) // hsl(200deg 85% 50%)
	// 	rid1_D: "hsl(220deg 0.72% 30%)",//hsl(220deg 0.72% 30%) // hsl(220deg 72% 50%)
	// 	green: "#DDFFCC7F",
	// };
	// _drawUserInputMarker.lineWidth = 3;
	// color_user_marker_stroke

	class SampleDesc {
		/**
		 * @param {string} sample_id
		 * @param {string} color
		 */
		constructor(sample_id, color) {
			this.sample_id = sample_id;
			this.color = color;
		}
	}

	/** @type {SampleDesc[]} */
	const sample_list = [];

	// QM6a_WT_veg-1
	[
		"QM6a",
		"CBS1-1",
	].forEach(group_id => {
		dup_list.forEach(dup_id => {
			const sample_id = `${group_id}_veg-${dup_id}`;
			sample_list.push(new SampleDesc(sample_id, color.WT_veg));// pink
		});
	});
	
	//QM6a_WT_D1-1
	[
		"WT",
	].forEach(group_id => {
		["1", "2", "3", "4", "5", "6", "7", "8",].forEach(time_id => {
			dup_list.forEach(dup_id => {
				const sample_id = `${group_id}_D${time_id}-${dup_id}`;
				sample_list.push(new SampleDesc(sample_id, color.WT_D));// pink
			});
		});
	});

	sample_list.push(new SampleDesc("QM6a_rid1_veg-1", color.rid1_veg));
	sample_list.push(new SampleDesc("QM6a_rid1_veg-2", color.rid1_veg));
	sample_list.push(new SampleDesc("QM6a_rid1_veg-3", color.rid1_veg));
	sample_list.push(new SampleDesc("CBS1-1_rid1_veg-1", color.rid1_veg));
	sample_list.push(new SampleDesc("CBS1-1_rid1_veg-2", color.rid1_veg));
	sample_list.push(new SampleDesc("CBS1-1_rid1_veg-3", color.rid1_veg));

	[
		"rid",
	].forEach(group_id => {
		["1", "2", "3", "4", "5", "6", "7", "8",].forEach(time_id => {
			dup_list.forEach(dup_id => {
				const sample_id = `${group_id}_D${time_id}-${dup_id}`;
				sample_list.push(new SampleDesc(sample_id, color.rid1_D));
			});
		});
	});

	const ref_idx = dataset.genome_info_list.findIndex(v => v.name == "QM6a");
	const genome_size = dataset.genome_info_list[ref_idx].chr_list.reduce((aa, v) => aa + v.length, 0);

	const refId = dataset.genome_info_list[ref_idx].name;
	await sample_list.reduce(async (prev_promise, { sample_id, color, }) => {
		await prev_promise;

		const ui32_sample = new module_Methyl_sampleData({
			// ...sample,
			ref: refId,
			sample: `${refId} ${sample_id}`,
			// name: `# of ${sample_id} reads mapped to ${refId} gDNA`,//# of RNA reads mapped to QM6a gDNA
			name: `${sample_id} >> ${refId}`,
			// url: `/RNA_mapping/${refId}/${sample_id}.depth.float32`,
			url: `/RNA_mapping/QM6a/${sample_id}.depth.float32`,
			rendering_condition: [
				new module_Methyl_sampleData_RenderingCondition({ color: color, condition: v => v.value > 0, min_width: 0.3, }),
			],
			max_display_value: genome_size,
			// value_normalize: 1000,
		});
		const f32a = await g_methylRenderer.load_float32Array(ui32_sample);

		/** @type {module_Methyl_ratioData[]} */
		const data = [];
		ui32_sample.data[viewerState.nChr - 1] = data;

		f32a.forEach((value, idx) => {
			if (Number.isNaN(value)) {
				// val.t = TAG_UNMAP;
			}
			else if (value > 0) {
				const val = new module_Methyl_ratioData();
				const ref_pos = Number(idx) + 1;
				val.start = ref_pos;
				val.end = ref_pos;
				delete val.strand;
				val.value = value;
				data.push(val);
			}
		});
		
		methyl_dataset_list.push(ui32_sample);
		
		await delayFrame(); // release memory
	}, Promise.resolve());

	await drawFrame();

	if (0) {
		const mmx = array_max(methyl_dataset_list.map(cc => array_max(cc.data[viewerState.nChr - 1].map(v => v.value))));
		function each(cc) {
			// const dd = cc.data[viewerState.nChr - 1];
			// const tt = dd.reduce((acc, v) => acc + v.value, 0);
			// cc.max_display_value = tt / dd.length;
			cc.max_display_value = mmx;
			cc.rendering_condition[0].min_width = 0.1;
			cc.hide = false;
		}
		methyl_dataset_list.forEach(cc => {
			cc.hide = true;
			// cc.name = cc.name.replace("# of ", "").replace(" reads mapped to QM6a gDNA", "");
			each(cc)
		});
		function array_average(ls) {
			return ls.reduce((acc, v) => acc + Number(v), 0) / ls.length;
		}
		function array_max(ls) {
			let p = -Infinity;
			for (let v of ls) {
				if (!Number.isNaN(v) && v != null) {
					p = Math.max(v, p);
				}
			}
			return p;
		}
		await drawFrame();
	}
}

/**
 * @see {init_jquery_ui}
 * @see {load_presset_20220916}
 * @see {load_presset_20221025_QM6a_RNA}
 * @param {boolean} load_whole_genome
 */
async function all_preset_20221104(load_whole_genome) {
	if (load_whole_genome) {
		await load_QM6a_CBS11_gff();
		await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
			await load_multi_alignment_results(chr_idx);
		}));
	}

	const list = [
		// { title: "", id: "TrQ_000953", name: "",       nChr: 1, start: 3074212, end: 3075458, FoldChange: 0.288, logFC: -1.79, Pvalue: 9.27E-005, DGEtest: "Down", },
		// { title: "", id: "TrQ_000954", name: "",       nChr: 1, start: 3075498, end: 3076191, FoldChange: 0.280, logFC: -1.84, Pvalue: 7.14E-005, DGEtest: "Down", },
		// { title: "", id: "TrQ_001091", name: "Rsp",    nChr: 1, start: 3756146, end: 3757401, FoldChange: 0.268, logFC: -1.90, Pvalue: 0.0002550, DGEtest: "Down", },
		// { title: "", id: "TrQ_001092", name: "Rsp",    nChr: 1, start: 3757712, end: 3759289, FoldChange: 0.207, logFC: -2.27, Pvalue: 4.17E-005, DGEtest: "Down", },
		// { title: "", id: "TrQ_001775", name: "",       nChr: 1, start: 6056393, end: 6060896, FoldChange: 0.302, logFC: -1.73, Pvalue: 4.93E-005, DGEtest: "Down", },
		// { title: "", id: "TrQ_003328", name: "qde-2",  nChr: 2, start: 4416138, end: 4421178, FoldChange: 0.390, logFC: -1.36, Pvalue: 0.0016100, DGEtest: "Down", },
		// { title: "", id: "TrQ_005124", name: "rid-1",  nChr: 3, start: 4421064, end: 4422110, FoldChange: 0.239, logFC: -2.07, Pvalue: 0.0024500, DGEtest: "Down", },
		// { title: "", id: "TrQ_005125", name: "rid-1",  nChr: 3, start: 4422195, end: 4423115, FoldChange: 0.001, logFC: -10.6, Pvalue: 1.73E-189, DGEtest: "Down", },
		// { title: "", id: "TrQ_006066", name: "",       nChr: 4, start: 2318627, end: 2321939, FoldChange: 0.361, logFC: -1.47, Pvalue: 0.1640000, DGEtest: "Down", },
		// { title: "", id: "TrQ_006507", name: "stk-53", nChr: 4, start: 3782146, end: 3787637, FoldChange: 0.223, logFC: -2.17, Pvalue: 9.00E-005, DGEtest: "Down", },
		// { title: "", id: "TrQ_009638", name: "stk-21", nChr: 7, start: 2198617, end: 2200552, FoldChange: 0.142, logFC: -2.82, Pvalue: 5.28E-005, DGEtest: "Down", },
		{ title: "", id: "TrQ_009909", name: "msh4",   nChr: 7, start: 3079958, end: 3083753, FoldChange: 0.256, logFC: -1.96, Pvalue: 0.0011300, DGEtest: "Down", },
		{ title: "", id: "TrQ_009963", name: "ts",     nChr: 7, start: 3263000, end: 3263998, FoldChange: 0.487, logFC: -1.04, Pvalue: 0.0007720, DGEtest: "Down", },
		{ title: "", id: "TrQ_007361", name: "ypr1",   nChr: 5, start: 2206595, end: 2209633, FoldChange: 0.101, logFC: -3.30, Pvalue: 2.03E-05, DGEtest: "Down", },
		{ title: "", id: "TrQ_005193", name: "",       nChr: 3, start: 4637692, end: 4640101, FoldChange: 0.205, logFC: -2.29, Pvalue: 0.000772, DGEtest: "Down", },
		
	].sort((a, b) => a.FoldChange - b.FoldChange);

	list.forEach(gg => {
		const dd = gg.name ? `${gg.name} (${gg.id})` : gg.id;
		gg.title = `Ch${romanize(gg.nChr)} ${dd}`;
	});

	return list;
}

/**
 * @param {string} refName
 * @param {number} nChr
 * @param {null|{ feature?: string, b_out_fa?: boolean; auto_RC?: boolean; ext_func?: (start: number, end: number, is_rev: boolean) => { start: number; end: number; }; }} [options]
 * @example
function unit_test() {
	ggg = await getGeneSeqById("QM6a", 1, "TrQ_001091", false);
	gggg = await getGeneSeqById("QM6a", 1, "TrQ_001091", { ext_left: 500, ext_right: 100, });

	ext_left = gggg.seq.indexOf(ggg.seq);
	console.log("ext_left:", ext_left == 500);

	ext_right = gggg.seq.length - ext_left - ggg.seq.length;
	console.log("ext_left:", ext_right == 100);

	ggg.gene == gggg.gene

	gggg.seq.match(/.{500,500}(.*).{100,100}/)[1] == ggg.seq

	// await findGeneByGffAll("QM6a", 4, predicate = row => row.type == "CDS" && row.attributes?.ID?.startsWith("TrQ_006507"))
}

ls = await all_preset_20221104();

vv = await Promise.all(ls.map(a => getGeneSeqById("QM6a", a.nChr, a.id, {
	feature: "CDS",
    auto_RC: true,
    ext_func: (start, end, is_rev) => ({
        start: is_rev ? (end - 2) : (start - 0),
        end: is_rev ? (end + 0) : (start + 2),
    }),
})));
// vv = (await vv).flat(1);

vv = vv.map(arr => {
    if (arr[0].gene.strand < 0) {
        return [...arr].sort((a, b) => b.gene.end - a.gene.end)[0];
    }
    else if (arr[0].gene.strand > 0) {
        return arr[0];
    }
    else {
        throw new Error("arr[0].strand");
    }
});
document.querySelector("textarea").value = vv.map(a => a.fa_header + "\n" + a.seq).join("\n");
vv.map(a => a.gene.attributes.ID + ":" + a.gene.strand + ":" + a.seq + ":" + a.gene.start + "-" + a.gene.end)

 */
async function getGeneSeqById(refName = "QM6a", nChr = 1, gene_id = null, options = {}) {
	const feature = options?.feature ?? "gene";// mRNA

	const predicate = row => row.type == feature && row.attributes?.ID?.startsWith(gene_id);
	const gg = await findGeneByGffAll(refName, nChr, predicate);

	if (gg) {
		const chr_idx = nChr - 1;

		const ref_info = dataset.genome_info_list.find(a => a.name == refName);
		const chr_info = ref_info.chr_list[chr_idx];
		const chr_seq_name = chr_info.chr;
		
		if (typeof dataset.results[chr_idx] == "string") {
			await load_multi_alignment_results(chr_idx);
		}

		const chr_seq = dataset.results[chr_idx][chr_seq_name].replace(/-/g, "");

		const b_out_fa = options.b_out_fa;
		
		// const node_parent = await findGeneByGffAll(refName, nChr, row => row.type == "gene" && row.attributes?.ID?.startsWith(gene_id));
		// const node_parent_strand = node_parent[0].strand;
		// return (() => {
		// 	if (node_parent_strand > 0) {
		// 		return gg;
		// 	}
		// 	else {
		// 		return [...gg].sort((a, b) => b.start - a.start);
		// 	}
		// })().map(gene => {
		// 	const has_rc = options?.auto_RC && (gene.strand ?? node_parent_strand) < 0;
		return gg.map(gene => {
			const has_rc = options?.auto_RC && gene.strand < 0;
		
			// gff_data_map['QM6a']['ChVII_QM6a'].find(a => a.type == "gene" && a.attributes.Name)
			
			// const ext_ATG_fwd = options.ext_ATG_fwd ?? 0;
			// const ext_ATG_rev = options.ext_ATG_rev ?? 0;
			// const pos_start = gene.start;// has_rc ? (gene.start - ext_ATG_fwd) : (gene.start - ext_ATG_rev);
			// const pos_end = gene.end;// has_rc ? (gene.end + ext_ATG_rev) : (gene.end + ext_ATG_fwd);

			// const pos_start = has_rc ? (gene.start - 100) : (gene.start - 500);
			// const pos_end = has_rc ? (gene.start + 500) : (gene.start + 100);

			const {
				start: pos_start,
				end: pos_end,
			} = options.ext_func ? options.ext_func(gene.start, gene.end, has_rc) : gene;

			const seq = (s => has_rc ? reverseComplement(s) : s)(chr_seq.slice(pos_start - 1, pos_end));
			const strand = gene.strand > 0 ? "+" : (gene.strand < 0 ? "-" : "");
			const fa_header = ">" + [
				gene.attributes.ID,
				gene.attributes.Name,
				[
					refName,
					chr_info.symbol,
					`${gene.start}-${gene.end}`,
					strand
				].join(":"),
				has_rc ? "reverse complement" : null,
				pos_start,
				pos_end,
			].filter(v => v).join(" ");

			if (b_out_fa) {
				return [
					fa_header,
					seq,
				].join("\n");
			}
			else {
				return {
					gene,
					is_gene_seq: pos_start == gene.start && pos_end == gene.end ? true : false,
					pos_start,
					pos_end,
					has_reverse: has_rc,
					has_complement: has_rc,
					fa_header,
					seq,
				};
			}
		});// for ... of
	}
}

/**
 * @param {string} refName
 * @param {number} nChr
 * @param {predicate: (value: GFF_ROW, index: number, array: GFF_ROW[]) => unknown} predicate
 */
async function findGeneByGffAll(refName = "QM6a", nChr = 1, predicate = row => row.attributes.ID) {
	const chr_idx = nChr - 1;
	const sChr = dataset.genome_info_list.find(a => a.name == refName).chr_list[chr_idx].symbol;
	
	if (gff_data_map?.[refName]?.[sChr] == null) {
		await load_QM6a_CBS11_gff();
	}
	
	const ls = gff_data_map[refName][sChr];
	return ls.filter(predicate);
}

/**
 * @deprecated update_SNP_Mode_genome_display_name
 */
function rename_SNP_Mode_genome_display() {
	return update_SNP_Mode_genome_display_name();
}

function update_SNP_Mode_genome_display_name() {
	dataset.progeny_list.forEach(progeny_name => {
		const display_name = ref_to_display_name(progeny_name);
		document.querySelector(`span[data-id=${progeny_name}]`).innerHTML = display_name;
	});
}

// if (module.hot) {
	// module.hot.accept('./print.js', function() {
	// 	console.log('Accepting the updated printMe module!');
	// 	printMe();
	// });
// }

/**
 * 20221114
 * 20221116
 * @warn 20211225_avg_pvalue_005 fun3 ??
 * @see {20211225_avg_pvalue_005 venn plot}
 */
async function preset_ssDNA_venn_minimal(load_data = false) {
	{
		if (preset_ssDNA_venn_minimal.load_data) {
			load_data = false;
		}

		if (load_data) {
			// alert("20211225_avg_pvalue_005 fun3 ??");
			await load_QM6a_CBS11_gff();// fun3 ??
		}
		if (load_data) {
			await Promise.all(dataset.genome_info_list[0].chr_list.map(async (ch_info, chr_idx) => {
				await load_multi_alignment_results(chr_idx);
			}));
		}
	
		const load_dup3 = true;
		if (load_data) {
			ssDNA_20211014.file_extname = "txt"
			ssDNA_20211014.lg = [2];
			ssDNA_20211014.show_pvalue = false;
			await ssDNA_20211014(true, { row_height: 3, fork_feat: true, load_QM6a: true, load_CBS1: true, load_dup3: load_dup3, }, false);
			
			ssDNA_peak_link_avg_and_2sample("QM6a");
			ssDNA_peak_link_avg_and_2sample("CBS1-1");
		}
		
		const min_pvalue = 0.05;
		const min_val = 0;
		const bUseAvg = true;// 20211225_avg_pvalue_005 -> venn plot // paper data
		const cond = load_dup3 && !bUseAvg ? (peak => peak.trip?.value >= 2) : null;
		const min_cons_len = 40
		await ssDNA_peak_venn(min_pvalue, min_val, bUseAvg, cond);
	}

	if (load_data) {
		preset_ssDNA_venn_minimal.load_data = preset_ssDNA_venn_minimal.load_data || load_data;
	}
}

/**
 * @see {@link https://stackoverflow.com/a/9083076|stackoverflow}
 * @see {@link http://blog.stevenlevithan.com/archives/javascript-roman-numeral-converter}
 * @param {number} num
 * @returns {string}
 */
 function romanize(num) {
	if (isNaN(num))
		return String(num);
	var digits = String(+num).split(""),
		key = ["","C","CC","CCC","CD","D","DC","DCC","DCCC","CM",
			"","X","XX","XXX","XL","L","LX","LXX","LXXX","XC",
			"","I","II","III","IV","V","VI","VII","VIII","IX"],
		roman = "",
		i = 3;
	while (i--)
		roman = (key[+digits.pop() + (i * 10)] || "") + roman;
	return Array(+digits.join("") + 1).join("M") + roman;
}

/**
 * @see {@link http://blog.stevenlevithan.com/archives/javascript-roman-numeral-converter}
 * @param {string} str
 * @returns {number}
 */
function deromanize (str) {
	var	str = str.toUpperCase(),
		validator = /^M*(?:D?C{0,3}|C[MD])(?:L?X{0,3}|X[CL])(?:V?I{0,3}|I[XV])$/,
		token = /[MDLV]|C[MD]?|X[CL]?|I[XV]?/g,
		key = {M:1000,CM:900,D:500,CD:400,C:100,XC:90,L:50,XL:40,X:10,IX:9,V:5,IV:4,I:1},
		num = 0, m;
	if (!(str && validator.test(str)))
		return NaN;
	while (m = token.exec(str))
		num += key[m[0]];
	return num;
}

async function import_from_unpkg(package_name) {
	const package_json = (await import(`https://unpkg.com/${package_name}/package.json`, { assert: { type: "json" } })).default;
	const dependencies = Object.keys(package_json.dependencies);

	const importmap = {
		"imports": {
			// "loadash": "/libs/loadash/index.js",
		},
	};
	
	dependencies.forEach(k => {
		importmap.imports[k] = `https://unpkg.com/${k}`;
	});

	const el_importmap = document.createElement("script");
	el_importmap.type = "importmap";
	el_importmap.innerHTML = JSON.stringify(importmap, null, "\t");
	const promise = async_onload(el_importmap);

	console.log("append importmap");
	document.body.append(el_importmap);
	await promise;
	console.log("loaded importmap");
	
	return await import(`https://unpkg.com/${package_name}`);

	async function async_onload(el) {
		return new Promise((resolve, reject) => {
			el.onload = resolve;
			el.onerror = reject;
		});
	}
}

document.body.style.opacity = 0.5;
document.body.style.filter = "blue(1px)";
init_jquery_ui().then(_ => {
	document.body.style.opacity = 1;
	document.body.style.filter = "";// clear
})
