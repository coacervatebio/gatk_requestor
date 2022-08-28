from pathlib import Path, PurePath
from common import ContainerTester, allowed_pattern
from runners import SnakemakeRunner
from checkers import VcfChecker

def test_genotype():

    data_path = PurePath("assets/genotype/")

    tmpdir = Path("/home/vagrant/tmp_guest")
    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path, tmpdir=tmpdir)
    tester.track_unexpected = False

    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()

debug = """

dr -v /home/vagrant/host_shared/coacervate/tests/requestor/assets/tmp_output/results:/mnt/results -v /home/vagrant/host_shared/coacervate/requestor/resources:/mnt/resources  broadinstitute/gatk /bin/bash
gatk --java-options '-Xmx4g' GenotypeGVCFs -R /mnt/resources/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V gendb:///mnt/results/combi_out/chr21_database -O /mnt/results/geno_out/combined.vcf.gz

/mnt/results/combi_out/chr21_database
/mnt/results/geno_out/combined_chr21.vcf.gz

"""
