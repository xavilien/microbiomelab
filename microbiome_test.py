# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import alignment
import random
import matplotlib.pyplot as plt


def Load16SFastA(path, fraction = 1.0, database_size=200, query_size=50):
    # from a file, read in all sequences and store them in a dictionary
    # sequences can be randomly ignored for testing by adjusting the fraction'''

    random.seed(11)
    
    infile = open(path, 'r')
    database = {}
    queries = {}
    
    for line in infile:
        if ">" in line:
            my_id = line[1:-1]
            if random.random() < fraction:
                if len(database) < database_size:
                    database[my_id] = ""
                elif len(queries) < query_size:
                    queries[my_id] = ""
        else:
            if my_id in database:
                database[my_id] += line[:-1]
            elif my_id in queries:
                queries[my_id] += line[:-1]
    
    return database, queries


def ConvertLibaryToKmerSets(library, K=2):
    new_lib = {}
    for k in library.keys():
        new_lib[k] = set()
        seq = library[k]
        for i in range(len(seq) - K + 1):
            new_lib[k].add(seq[i:i+K])
        
    return new_lib

def JaccardIndex(s1, s2):
    numerator = float(len(s1.intersection(s2)))
    denominator = float(len(s1.union(s2)))
    return numerator/denominator

def KmerMatch(sequence_kmer_set, library_kmer_set):
    best_score = 0.0
    best_match = None
    
    #add your code here to find the best kmer match
    for k in library_kmer_set:
        score = JaccardIndex(sequence_kmer_set, library_kmer_set[k])
        if score > best_score:
            best_score = score
            best_match = k

    return best_score, best_match


def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None

    """
    sequence : query sequence
    library : list of 16S sequences
    
    returns the 16S sequence with the highest alignment score with sequence
    """

    for k in library: 
        (score, _, _) = alignment.local_align(sequence, library[k])
        if best_match is None or score > best_score:
            best_score = score 
            best_match = k
    
    return best_score, best_match   # stuff only to run when not called via 'import' here


def plot_agreement_curve_1(database, queries):
    num_thresholds = 10
    thresholds = [2 * i + 1 for i in range(num_thresholds)]
    overall_scores = []
    for K in thresholds:
        print(f"Testing agreement for kmers of length {K}")
        database_kmers = ConvertLibaryToKmerSets(database, K)
        queries_kmers = ConvertLibaryToKmerSets(queries, K)
        scores = []
        for k, sequence in queries.items():
            _, alignment_best_match = AlignmentMatch(sequence, database)
            _, kmer_best_match = KmerMatch(queries_kmers[k], database_kmers)
            if kmer_best_match == alignment_best_match:
                scores.append(1)
            else:
                scores.append(0)
        
        overall_scores.append(sum(scores) / len(scores))

    plt.plot(thresholds, overall_scores, marker='o')
    plt.xlabel('K (k-mer size)')
    plt.ylabel('Agreement (fraction matched)')
    plt.title('K-mer vs Alignment Agreement')
    plt.ylim(0, 1)
    plt.show()

    out_path = "kmer_vs_alignment_agreement.png"
    plt.gcf().savefig(out_path, bbox_inches='tight', dpi=300)

def calculateMinimizers(window_size, k, sequence):
    """
    returns: set of minimizers for each sequence
    """
    num_windows = int(len(sequence)/(window_size-k))
    minimizers = set()
    for i in range(num_windows):
        window_start_index = window_size * i
        minimizer = None
        for seq_index in range(len(sequence) - k):
            current_kmer = sequence[window_start_index + seq_index : window_start_index + seq_index + k]
            if minimizer is None or current_kmer < minimizer: 
                minimizer = current_kmer
        minimizers.add(minimizer)

    return minimizers

def calculate_database_minimizers(window_size, k, database):
    """
    returns: list of tuples where each tuple is a database sequence and its list of minimizers
    """
    database_minimizers = {}

    for seq in database:
        seq_minimizers = calculateMinimizers(window_size, k, database[seq])
        database_minimizers[seq] = (database[seq], seq_minimizers)

    return database_minimizers


def MinimizerMatch(sequence, minimizer, library_minimizers):
    best_score = -10000000000
    best_match = None

    """
    sequence : query sequence
    library : list of 16S sequences

    run alignment free sequencing on each query sequence against the database sequences (only run if they share at least one minimizer)
    
    returns the 16S sequence with the highest alignment score with sequence
    """

    for k in library_minimizers:
        if len(minimizer.intersection(library_minimizers[k][1])) >= 1:
            (score, _, _) = alignment.local_align(sequence, library_minimizers[k][0])
            if best_match is None or score > best_score:
                best_score = score 
                best_match = k
    
    return best_score, best_match   # stuff only to run when not called via 'import' here


def plot_agreement_curve_2(database, queries):
    best_k = 15
    num_window_sizes = 10
    window_sizes = [15 * i + 20 for i in range(num_window_sizes)]
    overall_scores = []
    for window_size in window_sizes:
        print(f"Testing agreement for window size of length {window_size}")
        database_minimizers = calculate_database_minimizers(window_size, best_k, database)
        queries_minimizers = calculate_database_minimizers(window_size, best_k, queries)
        database_kmers = ConvertLibaryToKmerSets(database, best_k)
        queries_kmers = ConvertLibaryToKmerSets(queries, best_k)
        scores = []
        for k, sequence in queries.items():
            _, minimizer_best_match = MinimizerMatch(sequence, queries_minimizers[k][1], database_minimizers)
            _, kmer_best_match = KmerMatch(queries_kmers[k], database_kmers)
            if kmer_best_match == minimizer_best_match:
                scores.append(1)
            else:
                scores.append(0)
        
        overall_scores.append(sum(scores) / len(scores))

    plt.plot(window_sizes, overall_scores, marker='o')
    plt.xlabel('Window size')
    plt.ylabel('Agreement (fraction matched)')
    plt.title('K-mer vs Alignment Agreement')
    plt.ylim(0, 1)
    plt.show()

    out_path = "kmer_vs_alignment_agreement_2.png"
    plt.gcf().savefig(out_path, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
   fn = "bacterial_16s_genes.fa"
   database, queries = Load16SFastA(fn, fraction=0.5, database_size=10, query_size=5)
   
   print ("Loaded %d 16s database sequences." % len(database))
   print ("Loaded %d 16s query sequences." % len(queries))

#    plot_agreement_curve_1(database, queries)

   plot_agreement_curve_2(database, queries)
   

   
   
