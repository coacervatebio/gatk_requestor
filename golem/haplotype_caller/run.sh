#!/bin/bash

ALIGN_IN=$1
REG=$2
VCF_OUT=$3

echo "Unpacking Reference.."
tar -xzf /tmp/ref_hg38.tar.gz -C /golem/entrypoint

echo "Running HaplotypeCaller with:"
# echo "java -jar /home/gatk-local.jar HaplotypeCaller -I $ALIGN_IN -O $VCF_OUT -R /home/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -L $REG -ERC GVCF"
# java -jar /golem/entrypoint/gatk-local.jar HaplotypeCaller -I $ALIGN_IN -O $VCF_OUT -R /golem/entrypoint/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -L $REG -ERC GVCF