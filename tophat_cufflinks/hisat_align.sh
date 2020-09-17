# mapping using hisat2
# the outputs of the assay can be saved in the summary file
hisat2 -p 8 -x hisat_reference -U Day8.fastq -S align/day8.sam >& day8_summary.txt
hisat2 -p 8 -x hisat_reference -U Day16.fastq -S align/day16.sam >& day16_summary.txt

