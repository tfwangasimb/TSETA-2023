// @ts-check

const readline = require("readline");

const {
	BlastnCoord,
	parseBlastnResults,
	isCollide, removeOverlap, groupByOverlap,
} = require("./BlastnCoord.js");


/**
 * @param {import("stream").Readable} stream
 * @returns {Promise<BlastnCoord[]>}
 */
async function parse_blastn_results(stream) {
	const rl = readline.createInterface({
		input: stream,
		crlfDelay: Infinity
	});
	// Note: we use the crlfDelay option to recognize all instances of CR LF
	// ('\r\n') in input.txt as a single line break.

	/** @type {BlastnCoord[]} */
	const list = [];

	for await (const line of rl) {
		const columns = line.split("\t");
		list.push(BlastnCoord.fromObject({
			query: columns[0],
			subject: columns[1],
			identity: Number(columns[2]),
			align: Number(columns[3]),
			mismatch: Number(columns[4]),
			gap: Number(columns[5]),
			qstart: Number(columns[6]),
			qend: Number(columns[7]),
			sstart: Number(columns[8]),
			send: Number(columns[9]),
			evalue: Number(columns[10]),
			score: Number(columns[11]),
		}));
	}

	// console.log("ddd");

	return list;
}

module.exports.parse_blastn_results = parse_blastn_results;
