#!/bin/bash

ALIGN_IN=$1
REG=$2
VCF_OUT=$3

echo "Running HaplotypeCaller.."
java -jar /home/gatk-local.jar HaplotypeCaller --version