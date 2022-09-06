#!/bin/bash

ALIGN_IN=$1
REG=$2
VCF_OUT=$3

echo "Unpacking Reference.."
tar -xf /run/reference_HG38.tar.zst -C /golem/entrypoint

CMD="java -jar /run/gatk-local.jar HaplotypeCaller -I $ALIGN_IN -O $VCF_OUT -R /golem/entrypoint/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -L $REG -ERC GVCF"
echo "Running HaplotypeCaller with:"
echo $CMD
$CMD