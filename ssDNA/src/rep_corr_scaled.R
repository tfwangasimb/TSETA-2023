library(psych)
library(dplyr)

# Rad51
tbl <- read.table("CBS1-1_rad51_vs_sae2/CBS1-1_SS_Rad51_all.sort.bam.scaled.sum_mean.tsv", sep="\t")
colnames(tbl) <- c("Chr","Start","End","Rad51_1","Rad51_2","Rad51_3","Sum","Mean")
# pdf("rep_corr_scaled_CBS1-1_rad51.pdf",4,4)
png("rep_corr_scaled_CBS1-1_rad51.png")
set.seed(10)
pairs.panels(sample_n(tbl, 100000)[,4:6],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
dev.off()

# Sae2
tbl <- read.table("CBS1-1_sae2_vs_spo11sae2/CBS1-1_SS_Sae2_all.sort.bam.scaled.sum_mean.tsv", sep="\t")
colnames(tbl) <- c("Chr","Start","End","Sae2_1","Sae2_2","Sae2_3","Sum","Mean")
# pdf("rep_corr_scaled_CBS1-1_sae2.pdf",4,4)
png("rep_corr_scaled_CBS1-1_sae2.png")
set.seed(10)
pairs.panels(sample_n(tbl, 100000)[,4:6], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
dev.off()

# spo11_rad51
tbl <- read.table("CBS1-1_spo11rad51_vs_spo11sae2/CBS1-1_ss_spo11_rad51_all.sort.bam.scaled.sum_mean.tsv", sep="\t")
colnames(tbl) <- c("Chr","Start","End","spo11_rad51_1","spo11_rad51_2","spo11_rad51_3","Sum","Mean")
# pdf("rep_corr_scaled_CBS1-1_spo11rad51.pdf",4,4)
png("rep_corr_scaled_CBS1-1_spo11rad51.png")
set.seed(10)
pairs.panels(sample_n(tbl, 100000)[,4:6],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
dev.off()

tbl <- read.table("CBS1-1_spo11sae2_vs_sae2/CBS1-1_ss_spo11_sae2_all.sort.bam.scaled.sum_mean.tsv", sep="\t")
colnames(tbl) <- c("Chr","Start","End","spo11_sae2_1","spo11_sae2_2","spo11_sae2_3","Sum","Mean")
set.seed(10)
# pdf("rep_corr_scaled_CBS1-1_spo11sae2.pdf",4,4)
png("rep_corr_scaled_CBS1-1_spo11sae2.png")
pairs.panels(sample_n(tbl, 100000)[,4:6],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
dev.off()
