from pathlib import Path


def get_samples(dir):
    """Accepts directory containing input files and returns their filenames"""
    samples = set()
    for filepath in Path(dir).iterdir():
        while filepath.suffix in {".cram", ".crai", ".bam", ".bas", ".bai", ".sam"}:
            filepath = filepath.with_suffix("")
        samples.add(filepath.name)
    return samples
