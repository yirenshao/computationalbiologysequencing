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
	#freq1 = {k: v for k, v in sorted(freq1.items(), key=lambda item: item[1], reverse=True)}
	#freq2 = {k: v for k, v in sorted(freq2.items(), key=lambda item: item[1], reverse=True)}
	return np.dot([freq1[key] for key in freq1], [freq2[key] for key in freq1])

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

#KNN
def dist(seq1, seq2, k):
    return np.sqrt(spectrum_kernel(seq1, seq1, k)-2*spectrum_kernel(seq1, seq2, k)+spectrum_kernel(seq2, seq2, k))

def KNN(train, test, knn_num, kmer_size):
    seq_train = [] 
    seq_name_train = []
    for path in train:
        seq, seq_name = parse_fasta(read_fasta(path))
        seq_train += seq
        seq_name_train += seq_name
    seq_test, seq_name_test = parse_fasta(read_fasta(test))
    train_class_list = ["exons" if "exon" in seq_name_train[i] else "introns" for i in range(len(seq_name_train))]
    test_true_class_list = ["exons" if "exon" in seq_name_test[i] else "introns" for i in range(len(seq_name_test))]
    test_pred_class_list = []
    for i in range(len(seq_test)):
        dist_list = []
        for j in range(len(seq_train)):
            dist_list.append(dist(seq_test[i], seq_train[j], kmer_size))
        dist_list = np.array(dist_list)
        indexes = np.argsort(dist_list)
        exons_num = 0
        introns_num = 0
        for index in indexes[:knn_num]:
            if train_class_list[index] == "exons":
                exons_num += 1
            else:
                introns_num += 1
        if exons_num > introns_num:
            test_pred_class_list.append("exons")
        else:
            test_pred_class_list.append("introns")
    correct = 0
    for i, j in zip( test_true_class_list, test_pred_class_list):
        if i == j:
            correct += 1
    
    accuracy = correct/len(test_true_class_list)

    return accuracy


train = ["train-exons100.fasta", "train-exons100.fasta", "train-exons30.fasta", "train-introns10.fasta", "train-introns30.fasta", "train-introns100.fasta"]
test = "test.fasta"

knn_neighbors = [1,3,5,7]
kmer_sizes = [2,4,6,8]
knn_results_list = []
for knn_num in knn_neighbors:
    for kmer_size in kmer_sizes:
        result = {}
        accuracy = KNN(train, test, knn_num, kmer_size)
        result['knn neighbors num'] = knn_num
        result['kmer size'] = kmer_size
        result['accuracy'] = accuracy
        knn_results_list.append(result)
        
        
d = {'knn neighbors num':[],'kmer size':[],'accracy',[]}
for i in range(len(knn_results_list)):
    d['knn neighbors num'].append(knn_results_list[i]['knn neighbors num'])
    d['kmer size'].append(knn_results_list[i]['kmer size'])
    d['accuracy'].
    append(knn_results_list[i]['accuracy'])
df = pd.DataFrame(data=d)

## Xiao's work