FROM eclipse-temurin:8u332-b09-jre

COPY gatk_release/gatk-package-4.2.6.1-local.jar /home/gatk-local.jar
COPY ref_hg38.tar.gz /home/ref_hg38.tar.gz