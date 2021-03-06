library(Rsubread)

# annotation file used is the merged GTF file created from cuffcompare
annotations <- file.choose()
read.table(annotations)

# input for each of the count operations for featurecounts is the sorted bam
day8_counts <- featureCounts(file.choose(), annot.ext = annotations,
                             isGTFAnnotationFile = TRUE,
                             GTF.featureType = "exon")
day8_counts$counts

day16_counts <- featureCounts(file.choose(), annot.ext = annotations,
                             isGTFAnnotationFile = TRUE)
day16_counts$counts

merged_counts <- merge(day8_counts$counts, day16_counts$counts, by=0, all=TRUE)
colnames(merged_counts) <- c("feature", "counts_day_8", "counts_day_16")
merged_counts$fold_change <- log2(merged_counts$counts_day_16/merged_counts$counts_day_8)
merged_counts <- merged_counts[order(-merged_counts$fold_change),]
merged_counts <- merged_counts[order(merged_counts$feature),]
head(merged_counts)

counts_cufflinks <- read.table(file.choose())
head(counts_cufflinks)

to_keep <- subset(counts_cufflinks, select = c(V1, V10))
colnames(to_keep) <- c("feature", "log_change_cufflinks")
to_keep <- to_keep[-1,]
head(to_keep)
to_keep <- to_keep[order(to_keep$feature),]
head(to_keep)

both_counts <- merge(merged_counts, to_keep, by.x = "feature", by.y = "feature",
                     all = TRUE)

# convert na values to 0 as measured in featurecounts
both_counts$fold_change[is.na(both_counts$fold_change)] <- 0
both_counts
length(both_counts$fold_change)
length(both_counts$log_change_cufflinks)
attach(both_counts)
plot(x = fold_change, y = log_change_cufflinks, pch = 19, col = 2, xlab = "log change, featureCounts",
     ylab = "log change, Cuffinks", main = "featureCounts vs. Cuffinks log change response")

plot(density(both_counts$fold_change), ylab = "density", xlab = "log2 fold change",
     main = "Distribution of log changes, Cufflinks vs. featureCounts", col = 1, lty = 2)
lines(density(as.numeric(both_counts$log_change_cufflinks)), col = 2, lty = 1)

legend(locator(1), legend = c("featureCounts", "Cufflinks"), col = c(2, 1), lty = c(1,2), cex = 0.8)


