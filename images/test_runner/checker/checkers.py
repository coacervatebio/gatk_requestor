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

class CramChecker:
    """
    Very crudely checks to ensure CRAM files do not *immediately* differ.
    Prevents needing to have the reference genome present in the test-checker
    image.
    """

    @staticmethod
    def compare_files(generated_file, expected_file):
        proc = sp.run(["cmp", generated_file, expected_file], capture_output=True)
        assert "line 1" not in str(proc.stdout)
