// @ts-check

const fs = require("fs");


/**
 * @param {any} value
 * @returns {boolean} !Number.isNaN(value) && Number.isFinite(value)
 */
function validateNumber(value) {
	return !Number.isNaN(value) && Number.isFinite(value);
}

/**
 * @param {any} value
 * @returns {boolean} !Number.isNaN(value) && Number.isFinite(value) && Number.isSafeInteger(value)
 */
function validateInteger(value) {
	return !Number.isNaN(value) && Number.isFinite(value) && Number.isSafeInteger(value);
}

/**
 * @see {@link ./validate_calmeth.js}
 * @param {string} log_file
 * @example
 * node ./validate_calmeth.js Qrid_rid1_D2_ctrl.run.log
 */
function validate_calmeth(log_file, print_log) {
	const str_log_text = fs.readFileSync(log_file).toString();
	return str_log_text.split("\n").map(line => {
		const [k, v] = line.split("\t").map(a => a.trim());
		switch (k) {
			case "Raw count of Met_C in CG:":
			case "Raw count of Non_Met_C in CG:":
			case "Raw count of Met_C in CHG:":
			case "Raw count of Non_Met_C in CHG:":
			case "Raw count of Met_C in CHH:":
			case "Raw count of Non_Met_C in CHH:":
				let numVal = parseInt(v, 10);
				if (print_log) {
					console.log(k, numVal);
				}
				return validateInteger(numVal);
			case "[CpG]":
			case "[mC]":
				let numArr = [...v.matchAll(/\d+/g)].map(a => parseInt(a[0], 10));
				// let [M, Mh, H, hU, U] = numArr;
				if (print_log) {
					console.log(k, numArr);
				}
				return numArr.every(a => validateInteger(a));
			case "mC/(C+T)":
			case "mCG/(CG+TG)":
			case "mCHG/(CHG+THG)":
			case "mCHH/(CHH+THH)":
				let va = v.match(/{(.*) \/ (.*)} = (.*)%/) ?? [];
				let v1 = parseInt(va[1], 10);
				let v2 = parseInt(va[2], 10);
				let v3 = parseFloat(va[3]);
				if (print_log) {
					console.log(k, v1, v2, v3);
				}
				return va && validateInteger(v1) && validateInteger(v2) && validateNumber(v3);
			// default:
			//   console.log(k);
			//   break;
		}
		return k ? k : true;// skip this line
	});
}

if (typeof module == "object" && typeof module.exports == "object") {
	module.exports.validate_calmeth = validate_calmeth;
}

// console.log(fs.realpathSync(process.argv[1]) == fs.realpathSync(__filename));

if (fs.realpathSync(process.argv[1]) == fs.realpathSync(__filename)) {
	const fail_list = process.argv.slice(2).filter(file_path => {
		const f = validate_calmeth(file_path, false);
		console.log(file_path, f.every(a => a));
		return !f;
	});

	console.log("fail:", fail_list.join(","));
}

const all_samples = [
    {
        "name": "QM6a WT veg +BS >> D",
        "url": "/BS_QM6a/QM6a_methyl.methratio.float32"
    },
    {
        "name": "CBS1-1 WT veg +BS >> M",
        "url": "/BS_20200411/CBS1-1/20200902_C_wt.methratio.float32"
    },
    {
        "name": "D x M WT D2 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_BS.methratio.float32"
    },
    {
        "name": "D x M WT D4 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D4_BS.methratio.float32"
    },
    {
        "name": "D x M WT D5 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D5_BS.methratio.float32"
    },
    {
        "name": "D x M WT D6 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D6_BS.methratio.float32"
    },
    {
        "name": "D x M WT D8 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D8_BS.methratio.float32"
    },
    {
        "name": "QM6a WT veg -BS >> D",
        "url": "/BS_QM6a/QM6a_control.methratio.float32"
    },
    {
        "name": "CBS1-1 WT veg -BS >> D",
        "url": "/nBS/QM6a_CBS_WT_ctrl.methratio.float32"
    },
    {
        "name": "D x M WT D2 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M WT D4 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M WT D5 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M WT D6 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M WT D8 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/WT_dim2_rid1_dimrid/QM6a/QM6a_WT_D8_ctrl.methratio.float32"
    },
    {
        "name": "QM6a dim2Δ veg +BS",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_veg_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D2 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D2_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D4 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D4_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D5 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D5_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D6 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D6_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D8 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D8_BS.methratio.float32"
    },
    {
        "name": "QM6a dim2Δ veg -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_veg_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D2 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D4 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D5 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D6 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D8 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qdim/Qdim_Pd_D8_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D2 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D2_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D4 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D4_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D5 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D5_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D6 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D6_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D8 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D8_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D2 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D4 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D5 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D6 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δ D8 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Cdimrid_Pd/np59_Cdimrid_Pd_D8_ctrl.methratio.float32"
    },
    {
        "name": "QM6a rid1 veg +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_BS.methratio.float32"
    },
    {
        "name": "CBS1-1 rid1 veg +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D2 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D2_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D4 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D4_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D5 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D5_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D6 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D6_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D8 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D8_BS.methratio.float32"
    },
    {
        "name": "QM6a rid1Δ veg -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_ctrl.methratio.float32"
    },
    {
        "name": "CBS1-1 rid1Δ veg -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D2 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D4 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D5 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D6 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D8 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_D8_ctrl.methratio.float32"
    },
    {
        "name": "QM6a rid1 veg +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_BS.methratio.float32"
    },
    {
        "name": "CBS1-1 rid1 veg +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D2 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D2_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D4 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D4_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D5 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D5_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D6 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D6_BS.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D8 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D8_BS.methratio.float32"
    },
    {
        "name": "QM6a rid1Δ veg -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Qrid_rid1_veg_ctrl.methratio.float32"
    },
    {
        "name": "CBS1-1 rid1Δ veg -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_veg_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D2 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D4 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D5 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D6 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M rid1Δ D8 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Batmeth2_Qrid_Crid/Crid_rid1_D8_ctrl.methratio.float32"
    },
    {
        "name": "dim2Δ rid1Δ veg +BS >> D",
        "url": "/data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_veg_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D2 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D2_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D4 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D4_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D5 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D5_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D6 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D6_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D8 BS +BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D8_BS.methratio.float32"
    },
    {
        "name": "dim2Δ rid1Δ veg -BS >> D",
        "url": "/data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_veg_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D2 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D4 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D5 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D6 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D8 -BS >> D",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np128_Qdimrid_ddr_D8_ctrl.methratio.float32"
    },
    {
        "name": "dim2Δ rid1Δ veg +BS >> M",
        "url": "/data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_veg_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D2 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D2_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D4 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D4_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D5 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D5_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D6 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D6_BS.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D8 BS +BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D8_BS.methratio.float32"
    },
    {
        "name": "dim2Δ rid1Δ veg -BS >> M",
        "url": "/data/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_veg_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D2 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D2_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D4 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D4_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D5 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D5_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D6 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D6_ctrl.methratio.float32"
    },
    {
        "name": "D x M dim2Δrid1Δ D8 -BS >> M",
        "url": "/BS_oxBS_ctrl/D2_D4_D5_D6_D8/Qdimrid_Cdimrid/np59_Cdimrid_ddr_D8_ctrl.methratio.float32"
    }
]

const child_process = require("child_process");

const lllss = child_process.execSync(`ls ../Batmeth2*/*.json`).toString().trim().split("\n");
const lllaa = lllss.filter(a => !a.endsWith(".info.json"));
console.log(lllaa)

