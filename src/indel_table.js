// @ts-check

const fs = require("fs");

/**
 * @param {string} str
 * @param {string | RegExp} regexp
 */
function* string_matchAll(str, regexp) {
	const rx = new RegExp(regexp);
	let m = null;
	while ((m = rx.exec(str)) != null) {
		yield m;
	}
}

async function main() {
	const seq = "";
	const a = [...string_matchAll(seq, /(-+)/g)];
	
	a.map(aa => aa[1].length).filter(a => a == 1).length;
}

main();

