import subprocess as sp
import pysam as ps


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
