import numpy as np
import bionumpy as bnp


@bnp.streamable()
def count_a(chunk):
    return np.sum(chunk.sequence == "A")


# first part
file = bnp.open("bionumpy-example-data/big.fq.gz")
chunk = file.read_chunk()
n = count_a(chunk)
print("Number of As:", n)


# second part
file = bnp.open("bionumpy-example-data/big.fq.gz")
chunks = file.read_chunks()
results = count_a(chunks)
print("Number of As: ", sum(results))