gatk --java-options "-Xmx4g" HaplotypeCaller \
-R ../resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I ../alignments/HG00183.alt_bwamem_GRCh38DH.20150826.FIN.exome.cram \
-O ../hc_out/intervals_bm_out_chr1.g.vcf.gz \
-ERC GVCF \
-L chr1