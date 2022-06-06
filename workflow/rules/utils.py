import os

def dir_to_samples(dir):
    samples = set()
    exts = ['.cram', '.crai', '.bam', '.bas', '.bai', '.sam']
    for filename in os.listdir(dir):
        for ext in exts:
            if ext in filename:
                sample = filename.replace(ext, '')
        samples.add(sample)
    return samples