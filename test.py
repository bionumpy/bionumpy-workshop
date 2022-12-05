import numpy as np
import bionumpy as bnp


def filter_reads(chunk):
    mask = np.mean(chunk.quality, axis=1) >= 20
    return chunk[mask]


f = bnp.open("ENCFF000RWH.fastq.gz")
out_file = bnp.open("ENCFF000RWH_filtered.fastq", "w")


for i, chunk in enumerate(f.read_chunks()):
    filtered_chunk = filter_reads(chunk)
    print(filtered_chunk)
    out_file.write(filtered_chunk)

    if i >= 5:
        break