# view the number of genes from the abundance quantification

cut -f9 cufflinks/day8/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l
cut -f9 cufflinks/day16/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l

# view the number of transcripts from the abundance quantification


cut -f9 cufflinks/day8/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l
cut -f9 cufflinks/day16/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l


# view the genes with a single transcript

cut –f9 cufflinks/day8/transcripts.gtf | cut –d ' ' –f2,4 | sort –u | cut –d ' ' –f1 | sort -u | uniq –c | grep –c " 1 "
cut –f9 cufflinks/day16/transcripts.gtf | cut –d ' ' –f2,4 | sort –u | cut –d ' ' –f1 | sort -u | uniq –c | grep –c " 1 "

