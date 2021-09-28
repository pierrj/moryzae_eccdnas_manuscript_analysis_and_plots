library(qqman)

my_data <- read.delim('G3.normalized_binned.filtered.renamed.withbinnames.table')
snps_of_interest <- readLines('snps.tohighlight')

pdf('coverage_plot.pdf', width=15, height = 4)
manhattan(my_data, chr = "CHROMOSOME", bp = "BASE", snp = "SNP", p = "COUNT", highlight = snps_of_interest, logp=FALSE, ylab= "Split reads" , xlab= "Scaffold" , genomewideline=FALSE, suggestiveline = FALSE, main = NA, annotateTop = TRUE, cex.lab=1.5, cex.axis = 1.5, ylim = c(0,20000), col = c("black", "grey"))
dev.off()