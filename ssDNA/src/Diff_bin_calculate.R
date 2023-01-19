args = commandArgs(trailingOnly=TRUE)

sum_mean_1 <- args[1]
sum_mean_2 <- args[2]

Chr = strsplit(args[3], ",");
# Chr = c("ChI_QM6a", "ChII_QM6a", "ChIII_QM6a", "ChIV_QM6a", "ChV_QM6a", "ChVI_QM6a", "ChVII_QM6a")

min_lg <- as.numeric(args[4])
max_lg <- as.numeric(args[5])

output_name <- args[6]

# sae2
tbl <- read.table(sum_mean_2, sep="\t")

# Rad51 vs. Sae2
tbl1 <- read.table(sum_mean_1, sep="\t")
diff_mean <- round((tbl1$V8-tbl$V8),2)
chr <- tbl$V1
start <- tbl$V2
end <- tbl$V3
otbl <- data.frame(chr,start,end,diff_mean)
otbl <- otbl[order(factor(otbl$chr, levels = Chr)),]
write.table(otbl, paste0("diff_bins_scaled_", output_name, ".tsv"), sep="\t", col.names=F, row.names=F, quote = F)
for(lg in c(min_lg:max_lg)) {
	# Diff_peak_cutoff_lg.py -> clissify positive or negative peak
	write.table(otbl[abs(otbl$diff_mean) > lg,], paste0("diff_bins_scaled_lg", lg, "_", output_name, ".bedgraph"), sep="\t", col.names=F, row.names=F, quote = F)

	# 20211227 # ignore negative
	# write.table(otbl[otbl$diff_mean > lg,], paste0("diff_bins_scaled_lg", lg, "_", output_name, ".bedgraph"), sep="\t", col.names=F, row.names=F, quote = F)
	
	# write.table(otbl[abs(otbl$diff_mean) > 2,], "diff_bins_scaled_lg2_Rad51_vs_Sae2.bedgraph", sep="\t", col.names=F, row.names=F, quote = F)
	# write.table(otbl[abs(otbl$diff_mean) > 3,], "diff_bins_scaled_lg3_Rad51_vs_Sae2.bedgraph", sep="\t", col.names=F, row.names=F, quote = F)
	# write.table(otbl[abs(otbl$diff_mean) > 4,], "diff_bins_scaled_lg4_Rad51_vs_Sae2.bedgraph", sep="\t", col.names=F, row.names=F, quote = F)
	# #write.table(otbl[otbl$diff_mean > 2,][,c("chr","start","end")], "diff_bins_scaled_lg2_Rad51_vs_Sae2.bed", sep="\t", col.names=F, row.names=F, quote = F)
	# #write.table(otbl[otbl$diff_mean > 3,][,c("chr","start","end")], "diff_bins_scaled_lg3_Rad51_vs_Sae2.bed", sep="\t", col.names=F, row.names=F, quote = F)
	# #write.table(otbl[otbl$diff_mean > 4,][,c("chr","start","end")], "diff_bins_scaled_lg4_Rad51_vs_Sae2.bed", sep="\t", col.names=F, row.names=F, quote = F)
}
