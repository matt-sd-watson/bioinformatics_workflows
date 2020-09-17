library(Rsubread)

annotations <- file.choose()
annotations
day8_counts <- featureCounts(file.choose(), annot.ext = annotations,
                             isGTFAnnotationFile = TRUE)
day8_counts$counts

day16_counts <- featureCounts(file.choose(), annot.ext = annotations,
                             isGTFAnnotationFile = TRUE)
day16_counts$counts

merged_counts <- merge(day8_counts$counts, day16_counts$counts, by=0, all=TRUE)
colnames(merged_counts) <- c("feature", "counts_day_8", "counts_day_16")
merged_counts$fold_change <- log2(merged_counts$counts_day_8/merged_counts$counts_day_16)
merged_counts <- merged_counts[order(-merged_counts$fold_change),]
merged_counts
