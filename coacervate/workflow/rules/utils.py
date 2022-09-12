from pathlib import Path


def dir_to_samples(dir):
    samples = set()
    for filepath in Path(dir).iterdir():
        while filepath.suffix in {".cram", ".crai", ".bam", ".bas", ".bai", ".sam"}:
            filepath = filepath.with_suffix("")
        samples.add(filepath.name)
    return samples
