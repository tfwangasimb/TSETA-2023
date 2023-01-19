// @ts-check

const fs = require("fs");

let default_value = {
	blastn_bin: "blastn",
	mafft_bin: "mafft",
	version: null,
}

/**
 * @returns {default_value}
 */
function loadSetting() {
	try {
		let text = fs.readFileSync("./setting.json").toString();
		
		let setting = JSON.parse(text);

		setting.blastn_bin = setting.blastn_bin != null ? setting.blastn_bin : default_value.blastn_bin;
		setting.mafft_bin = setting.mafft_bin != null ? setting.mafft_bin : default_value.mafft_bin;

		setting.version = setting.version != null ? setting.version : null;

		return setting;
	}
	catch (ex) {
		return default_value;
	}
}

module.exports.loadSetting = loadSetting;

