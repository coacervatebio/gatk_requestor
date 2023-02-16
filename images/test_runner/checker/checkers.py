import logging
import subprocess as sp
import pysam as ps

LOGGER = logging.getLogger(__name__)

# Checkers
class SimpleChecker:
    """Byte for byte file comparison checker using `cmp`"""

    @staticmethod
    def compare_files(generated_file, expected_file):
        # Check that cmp returns no output (no difference between files)
        assert len(sp.check_output(["cmp", generated_file, expected_file])) == 0


class VcfChecker:
    """Compares VCFs based on records, excluding timestamped header"""

    @staticmethod
    def compare_files(generated_file, expected_file):
        gv = ps.VariantFile(generated_file, "r")
        ev = ps.VariantFile(expected_file, "r")
        g_recs = [str(r) for r in gv.fetch()]
        e_recs = [str(r) for r in ev.fetch()]
        assert g_recs == e_recs

# Utils
def allowed_pattern(tf):
    """
    Filter function to remove files from comparison matching certain patterns
    identifying files that change from run-to-run. E.g. failing comparison
    because the timestamp is different.

    :param tf: Target file being evaluated
    :returns False: File contains failure pattern
    :returns True: File does not contain failure pattern
    """
    expected_cmp_fail_patterns = ["__", "vcfheader.vcf", "gz.tbi"]

    s_tf = str(tf)
    for ecfp in expected_cmp_fail_patterns:
        if ecfp in s_tf:
            LOGGER.warn(f"Excluding {s_tf} for containing {ecfp}")
            return False
    return True
