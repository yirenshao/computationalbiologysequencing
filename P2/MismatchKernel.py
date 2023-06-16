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
    seq = ''
    seq_name = []
    for line, index in zip(lines, range(len(lines))):
        if index == len(lines) - 1:
            seq += line.strip()
            seq_list.append(seq)
        if line.startswith('>'):
            seq_list.append(seq)
            seq = ''
            name = line.strip()
            seq_name.append(name.split(" /")[0][1:])
            continue
        else:
            seq += line.strip()
    for i in seq_list:
        if i == '':
            seq_list.remove(i)
    return seq_list,seq_name

def diff(key1,key2):
    list1 = list(key1)
    list2 = list(key2)
    dif = 0
    for i in range(len(list1)):
        if list1[i] != list2[i]:
            dif += 1
    return dif

# Implement Mismatch Kernel
def mismatch_kernel(seq1, seq2, k):
	freq1 = freq_kmers(seq1, k)
	freq2 = freq_kmers(seq2, k)
	for key in freq1:
		if key not in freq2:
			freq2[key] = 0
	for key in freq2:
		if key not in freq1:
			freq1[key] = 0
   
	freq1 = dict( sorted(freq1.items(), key=lambda x: x[0]) )
	freq2 = dict( sorted(freq2.items(), key=lambda x: x[0]) )

	newfreq1 = freq1.copy()
	newfreq2 = freq2.copy()

	for key1 in freq1:
		for key2 in freq1:
			dif = diff(key1,key2)
			if dif == 1 and freq1[key2] == 0:
				newfreq1[key1] += 1
			if dif == 1 and freq1[key2] > 0:
				newfreq1[key1] += freq1[key2]
				
	for key1 in freq2:
		for key2 in freq2:
			dif = diff(key1,key2)
			if dif == 1 and freq2[key2] == 0:
				newfreq2[key1] += 1	
			if dif == 1 and freq2[key2] > 0:
				newfreq2[key1] += freq2[key2]

	return np.dot(list(newfreq1.values()), list(newfreq2.values()))

## Yiren's work