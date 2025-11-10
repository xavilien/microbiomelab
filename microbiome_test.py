# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import alignment
import random
import matplotlib.pyplot as plt
random.seed(11)


def Load16SFastA(path, fraction = 1.0, database_size=200, query_size=50):
    # from a file, read in all sequences and store them in a dictionary
    # sequences can be randomly ignored for testing by adjusting the fraction'''
    
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


# Agreement curve as a function of k-mer size
def plot_agreement_curve_1(database, queries):
    num_thresholds = 10
    thresholds = [2 * i + 1 for i in range(num_thresholds)]
    overall_scores = []

    best_K = None
    best_score = None

    agreement_results = {}
    for k, sequence in queries.items():
        _, alignment_best_match = AlignmentMatch(sequence, database)
        agreement_results[k] = alignment_best_match

    for K in thresholds:
        print(f"Testing agreement for kmers of length {K}")
        database_kmers = ConvertLibaryToKmerSets(database, K)
        queries_kmers = ConvertLibaryToKmerSets(queries, K)
        scores = []
        for k, sequence in queries.items():
            alignment_best_match = agreement_results[k]
            _, kmer_best_match = KmerMatch(queries_kmers[k], database_kmers)
            if kmer_best_match == alignment_best_match:
                scores.append(1)
            else:
                scores.append(0)
        
        score = sum(scores) / len(scores)
        overall_scores.append(score)
        if best_score is None or score > best_score:
            best_K = K
            best_score = score

    plt.figure()
    plt.plot(thresholds, overall_scores, marker='o')
    plt.xlabel('K-mer length')
    plt.ylabel('Agreement')
    plt.title('K-mer vs Alignment Agreement')
    plt.ylim(0, 1)
    out_path = "task4.png"
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()

    return best_K

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


def MinimizerMatch(sequence, minimizer, library_minimizers, query_kmers, database_kmers, minimum_minimizer_overlap=1):
    """
    sequence : query sequence
    library : list of 16S sequences

    run alignment free sequencing on each query sequence against the database sequences (only run if they share at least one minimizer)
    
    returns the 16S sequence with the highest alignment score with sequence
    """

    for k in library_minimizers:
        if len(minimizer.intersection(library_minimizers[k][1])) >= minimum_minimizer_overlap:
            return KmerMatch(query_kmers, database_kmers)


def plot_agreement_curve_2(database, queries, BEST_K=15):
    num_window_sizes = 10
    window_sizes = [15 * i + 20 for i in range(num_window_sizes)]
    overall_scores = []

    # Calculate local align ground truth
    agreement_results = {}
    for k, sequence in queries.items():
        _, alignment_best_match = AlignmentMatch(sequence, database)
        agreement_results[k] = alignment_best_match

    # Calculate kmers for queries and database based on best_K from previous part
    database_kmers = ConvertLibaryToKmerSets(database, BEST_K)
    queries_kmers = ConvertLibaryToKmerSets(queries, BEST_K)

    for window_size in window_sizes:
        print(f"Testing agreement for window size of length {window_size}")
        database_minimizers = calculate_database_minimizers(window_size, BEST_K, database)
        queries_minimizers = calculate_database_minimizers(window_size, BEST_K, queries)
        scores = []
        for k, sequence in queries.items():
            _, minimizer_best_match = MinimizerMatch(sequence, queries_minimizers[k][1], database_minimizers, queries_kmers[k], database_kmers)
            alignment_best_match = agreement_results[k]
            if minimizer_best_match == alignment_best_match:
                scores.append(1)
            else:
                scores.append(0)
        
        overall_scores.append(sum(scores) / len(scores))

    plt.figure()
    plt.plot(window_sizes, overall_scores, marker='o')
    plt.xlabel('Window size')
    plt.ylabel('Agreement')
    plt.title('Window size vs Alignment Agreement')
    plt.ylim(0, 1)
    out_path = "task6.png"
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()


def get_incorrect_base(base):
    bases = ["A", "T", "G", "C"]
    alternative_bases = [i for i in bases if i != base]
    return random.choice(alternative_bases)


def illumina_mutation(queries):
    mutated_queries = {}

    for k, seq in queries.items():
        mutated_seq = []
        for i in range(min(250, len(seq))):
            if random.random() > 0.99:
                mutated_seq.append(get_incorrect_base(seq[i]))
            else:
                mutated_seq.append(seq[i])

        mutated_queries[k] = "".join(mutated_seq)
        
    return mutated_queries


def nanopore_mutation(queries):
    mutated_queries = {}
    for k, seq in queries.items():
        mutated_seq = []
        for i in range(len(seq)):
            if random.random() > 0.9:
                mutated_seq.append(get_incorrect_base(seq[i]))
            else:
                mutated_seq.append(seq[i])
        
        mutated_queries[k] = "".join(mutated_seq)

    return mutated_queries


if __name__ == "__main__":
    fn = "bacterial_16s_genes.fa"
    database, queries = Load16SFastA(fn, fraction=0.5, database_size=10, query_size=5)

    print ("Loaded %d 16s database sequences." % len(database))
    print ("Loaded %d 16s query sequences." % len(queries))

    best_K = plot_agreement_curve_1(database, queries)
    print(f"Best K-mer size is {best_K}")
    plot_agreement_curve_2(database, queries, BEST_K = best_K)

    illumina_queries = illumina_mutation(queries)
    nanopore_queries = nanopore_mutation(queries)
   

   
   
