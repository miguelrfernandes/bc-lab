kallisto index -i reference/transcripts.idx reference/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto quant -i reference/transcripts_new.idx -o output TCGA-DM-A288-01A_1.fastq TCGA-DM-A288-01A_2.fastq