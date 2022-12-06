import bionumpy as bnp
from pyjaspar import jaspardb
from bionumpy.sequence.position_weight_matrix import PWM, get_motif_scores
import numpy as np
import plotly.express as plx

# Get the ctcf bed file
# ! wget https://www.encodeproject.org/files/ENCFF843VHC/@@download/ENCFF843VHC.bed.gz

peaks = bnp.open("ENCFF843VHC.bed.gz").read()
reference_genome = bnp.open_indexed("/home/knut/Data/hg38.fa")
contig_lenghts = reference_genome.get_contig_lengths()
# sorted_peaks = bnp.arithmetics.sort_all_intervals(
#     peaks, sort_order=list(contig_lenghts.keys()))


# We need to synch our peaks to the reference genome
# multistream = bnp.MultiStream(contig_lenghts,
#                               intervals=sorted_peaks,
#                               sequence=reference_genome)
# 

# Get the sequences of the peaks
#sequence_stream = bnp.sequence.get_sequences(multistream.sequence,
#                                              multistream.intervals)

# Merge all the sequences from the stream
# sequences = np.concatenate(list(sequence_stream))

sequences = reference_genome.get_interval_sequences(peaks)

# Get all human motifs :)
jdb_obj = jaspardb(release="JASPAR2020")
human_motifs = jdb_obj.fetch_motifs(collection="CORE", species=["9606"])

counts = []
lengths = []
# Calculate the motif scores for each sequence
for i, motif in enumerate(human_motifs):
    if i % 10 == 0:
        print(i)
    pwm = PWM.from_dict(motif.pwm)
    motif_scores = get_motif_scores(sequences, pwm)
    adjusted_scores = motif_scores-motif.length*np.log(0.25)
    counts.append((adjusted_scores.max(axis=-1) > np.log(0.9)).sum())
    lengths.append(motif.length)


# Check that lentghts are not an influencing factor
plx.scatter(lengths, counts)

names = [m.name for m in human_motifs]

# Sort by number of significant hit
args = np.argsort(counts)
sorted_names = [names[i] for i in args]
