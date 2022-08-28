import subprocess as sp

class SimpleChecker:
    
    @staticmethod
    def compare_files(generated_file, expected_file):
        # Check that cmp returns no output (no difference between files)
        assert len(sp.check_output(["cmp", generated_file, expected_file])) == 0


class VcfChecker:
    def compare_files(generated_file, expected_file):
        ...