import numpy as np
# Get k-mers from a sequence
def get_kmers(sequence, k):
    return [sequence[x:x+k].lower() for x in range(len(sequence) - k + 1)]

# Get the frequency of each k-mer in a sequence
def freq_kmers(seq, k):
    kmers = get_kmers(seq, k)
    freq = {}
    for kmer in kmers:
        if kmer not in freq: 
            freq[kmer] = 1
        else: 
            freq[kmer] += 1
    return freq


# change fasta file into a list of sequences
# fasta reader
def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines

# fasta parser


def parse_fasta(lines):
    seq_list = []
    seq_name = [line for line in lines if line.startswith('>')]
    seq = ''
    for line, index in zip(lines, range(len(lines))):
        if index == len(lines) - 1:
            seq += line.strip()
            seq_list.append(seq)
        if line.startswith('>'):
            seq_list.append(seq)
            seq = ''
            continue
        else:
            seq += line.strip()
    for i in seq_list:
        if i == '':
            seq_list.remove(i)
    return seq_list, seq_name


# Implement Spectrum Kernel
def spectrum_kernel(seq1, seq2, k):
    freq1 = freq_kmers(seq1, k)
    freq2 = freq_kmers(seq2, k)
    for key in freq1:
        if key not in freq2:
            freq2[key] = 0
    for key in freq2:
        if key not in freq1:
            freq1[key] = 0
    return np.dot([freq1[key] for key in freq1], [freq2[key] for key in freq1])

    ## XIao's work