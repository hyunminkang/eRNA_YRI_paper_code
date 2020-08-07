# pval_fcn = function(x, bk){
#     n <- length(x)
#     pvals <- numeric(length = n)
#     for(i in 1:n){
#         pvals[i] <-  1-(sum(bk<x[i]) / length(bk))
#     }
#     return(pvals)
# }

TRE_rpm <- read.table("bed/procap.rpm.lm.txt", header = F, sep = "\t", stringsAsFactors = F)
control_rpm <-read.table("bed/background_rpm.txt", header = F, sep = "\t", stringsAsFactors = F)

TRE_vector <- data.frame(rpm=TRE_rpm[, 3])
control_vector <- control_rpm[,3]

pvals <- 1-(apply(TRE_vector, 1, function(x) sum(control_vector<x))/length(control_vector))

fdrs <- p.adjust(pvals, method = "fdr")

stat_table <- data.frame(chr = TRE_rpm[,1], pos = TRE_rpm[,2], pval = pvals, fdr = fdrs, rpm = TRE_vector)

write.table(stat_table, file = "TRE_stats.txt", col.names = F, sep = "\t", quote = F)

filtered <- stat_table[which(stat_table$rpm>0.5),]

print(paste("The 0.5 rpm cut off corresponds to a p-value of", max(filtered$pval), "and an FDR of", max(filtered$fdr), sep = " "))

write.table(stat_table, file = "results/TRE_stats.txt", col.names = F, sep = "\t", quote = F)
write.table(filtered, file = "results/TRE_stats_filtered.txt", col.names = F, sep = "\t", quote = F)

pdf("results/FDR_plot.pdf")
plot(stat_table$fdr, stat_table$rpm)
dev.off()
