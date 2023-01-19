#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) <= 1) {
	print(args[1:length(args)]);
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

print(args[1])
print("----------")
print(args[2:5])
print("----------")
print(args[6:9])
print("----------")

# if no header
# use first column
# input_sets <- lapply(args[6:9], function (x) read.table(x, sep="\t", header=FALSE)$V1)

# has header
# use id column
input_sets <- lapply(args[6:9], function (x) read.table(x, sep="\t", header=TRUE)$id)

library(VennDiagram)
venn.diagram(
	x = input_sets,
	category.names = args[2:5],
	filename = paste0(args[1], '.png'),
	output = TRUE,
	
	# Output features
	imagetype = "png",
	# height = 1024,
	# width = 1024,
	# resolution = 300,
	# compression = "lzw",
	
	# Circles
	# lwd = 2,
	# lty = 'blank',
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),,
)

