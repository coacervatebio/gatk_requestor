samtools view ../alignments/HG00183.alt_bwamem_GRCh38DH.20150826.FIN.exome.cram \
chr1 -T ../resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-O cram -o ../alignments/view_hc_out_chr1.cram

samtools index ../alignments/view_hc_out_chr1.cram

gatk --java-options "-Xmx3g" HaplotypeCaller \
-R ../resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I ../alignments/view_hc_out_chr1.cram \
-O ../hc_out/view_hc_out_chr1.g.vcf.gz \
-ERC GVCF