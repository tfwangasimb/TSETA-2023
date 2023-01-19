args = commandArgs(trailingOnly=TRUE)

###############################################################################

cmp_1 <- args[1]

sum_mean_1 <- args[2] # paste0(args[argc++], "_all.sort.bam.scaled.sum_mean.tsv") # "../bam_scaled/QM6a_SS_Rad51_all.sort.bam.scaled.sum_mean.tsv"
sum_mean_2 <- args[3] # paste0(args[argc++], "_all.sort.bam.scaled.sum_mean.tsv") # "../bam_scaled/QM6a_SS_Sae2_all.sort.bam.scaled.sum_mean.tsv"

min_lg <- as.numeric(args[4]) # 2
max_lg <- as.numeric(args[5]) # 4

min_peak_len <- as.numeric(args[6]) # 11 # >10 or >= 11
pValue <- as.numeric(args[7]) # 0.05

###############################################################################

dt1 <- read.table(sum_mean_1, sep="\t")
dt2 <- read.table(sum_mean_2, sep="\t")

for (l in min_lg:max_lg) {
	dat <- read.table(paste0("Diff_peaks_lg", l, "_", cmp_1, ".txt"),sep="\t", header=T)
	dat2 <- dat[ which((dat["End"]-dat["Start"])>=min_peak_len),]
	cmpdt <- dt1
	ctrdt <- dt2
	V1 <- cmpdt$V1
	V2 <- cmpdt$V2
	V3 <- cmpdt$V3
	V4 <- cmpdt$V7
	V5 <- ctrdt$V7
	sumdt <- data.frame(V1,V2,V3,V4,V5)
	dat2$ks_test_pvalue <- apply(dat2, 1, function(x) ks.test(as.integer(sumdt[ which(sumdt$V1 == x["Chr"] & sumdt$V2 >= as.integer(x["Start"]) & sumdt$V3 <= as.integer(x["End"])),]$V4),as.integer(sumdt[ which(sumdt$V1 == x["Chr"] & sumdt$V2 >= as.integer(x["Start"]) & sumdt$V3 <= as.integer(x["End"])),]$V5),exact=FALSE)$p.value)
	
	ks_file <- paste0("../", "Diff_peaks_lg", l, "_", cmp_1, "_KStest.txt")
	ks_file_p005 <- paste0("../", "Diff_peaks_lg", l, "_", cmp_1, "_KStest_p", pValue, ".txt")
	write.table(dat2, ks_file, quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(dat2[dat2$ks_test_pvalue <= pValue, ], ks_file_p005, quote = FALSE, sep = "\t", row.names = FALSE)
}
