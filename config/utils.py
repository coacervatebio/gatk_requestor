def file_to_sample(filename):
    exts = ['.cram', '.crai', '.bam', '.bas', '.bai', '.sam']
    for ext in exts:
        if ext in filename:
            filename = filename.replace(ext, '')
    return filename