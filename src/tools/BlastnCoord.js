// @ts-check


class BlastnCoord {
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
	
	/**
	 * @param {BlastnCoord[]} aln_list
	 */
	static to_tsv(aln_list) {
		const text = [
			BlastnCoord.tableHeader.join("\t"),
			aln_list.map(a => a.toArray().join("\t")).join("\n"),
		].join("\n");
		return text;
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
		// let coord = new BlastnCoord();
		// coord.qstart = obj.qstart;
		// coord.qend = obj.qend;
		// coord.sstart = obj.sstart;
		// coord.send = obj.send;
		// return coord;
		return Object.assign(new BlastnCoord(), obj);
	}
}

/**
 * @param {string} text
 * @returns {BlastnCoord[]}
 */
function parseBlastnResults(text) {
	let _table = text.trim().split("\n").filter(a => a && a.length).map(line => {
		let columns = line.split("\t");
		return BlastnCoord.fromObject({
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
		});
	});
	return _table;
}

/**
 * return duplicate group
 * @param {BlastnCoord[]} rows
 * @returns {BlastnCoord[][]}
 */
function groupByOverlap(rows) {
	/** @type {BlastnCoord[][]} */
	let groups = [];

	const tempRows = rows.slice(0);

	while (tempRows.length) {
		const left = tempRows.splice(0, 1)[0];
		if (!left) {
			throw new Error("tempRows.splice(0, 1)[0] -> undefined");
		}
		const overlap_group = [left];

		for (let j = 0; j < tempRows.length; ++j) {
			const right = tempRows[j];
			if (isCollide(left.sstart, left.send, right.sstart, right.send) ||
				isCollide(left.qstart, left.qend, right.qstart, right.qend)
			) {
				tempRows.slice(j, 0);
				overlap_group.push(right);
			}
		}

		groups.push(overlap_group.sort((a, b) => b.score - a.score));
	}

	return groups.sort((a, b) => b[0].score - a[0].score);
}

function isCollide(x11, x12, x21, x22) {
	return x11 <= x22 && x12 >= x21;
}

/**
 * @param {BlastnCoord[]} rows
 */
function removeOverlap(rows) {
	/** @type {Set<BlastnCoord>} */
	let removes = new Set();

	for (let i = 0; i < rows.length; ++i) {
		const left = rows[i];
		if (!removes.has(left)) {
			if (left.strand <= 0) {
				removes.add(left);
			}
			else {
				for (let j = i + 1; j < rows.length; ++j) {
					const right = rows[j];
					if (isCollide(left.sstart, left.send, right.sstart, right.send) ||
						isCollide(left.qstart, left.qend, right.qstart, right.qend)
					) {
						if (left.score < right.score) {
							removes.add(left);
						}
						else if (right.score < left.score) {
							removes.add(right);
						}
					}
				}
			}
		}
	}

	let filtered = rows.filter(left => !removes.has(left));
	
	// let prev_sstart = 0;
	// filtered = rows = rows.filter(left => {
	// 	let delta = left.sstart - prev_sstart;
	// 	if (delta >= 0 && delta <= init_max_delta) {
	// 		prev_sstart = left.sstart;
	// 		return true;
	// 	}
	// });

	return filtered;
}

module.exports.BlastnCoord = BlastnCoord;

module.exports.parseBlastnResults = parseBlastnResults;
module.exports.isCollide = isCollide;
module.exports.removeOverlap = removeOverlap;
module.exports.groupByOverlap = groupByOverlap;
