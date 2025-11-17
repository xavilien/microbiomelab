# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import alignment
import random
import matplotlib.pyplot as plt
import time
import json
import os
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


def check_benchmark(name):
    with open("benchmark.json", 'r') as f:
        benchmark = json.load(f)
    return name in benchmark


def save_benchmark(name, duration, max_agreement):
    benchmark = {}
    with open("benchmark.json", 'r') as f:
        benchmark = json.load(f)

    if name not in benchmark:
        benchmark[name] = (duration, max_agreement)
    with open("benchmark.json", 'w') as f:
        json.dump(benchmark, f)


def get_local_alignment_results(database, queries, query_type):
    file_path = f"saved_results/local_alignment_{query_type}.json"
    if os.path.isfile(file_path):
        with open(file_path) as f:
            agreement_results = json.load(f)
        return agreement_results

    start = time.time()
    agreement_results = {}
    for k, sequence in queries.items():
        _, alignment_best_match = AlignmentMatch(sequence, database)
        agreement_results[k] = alignment_best_match
    end = time.time()

    duration = end - start
    save_benchmark("local_alignment", duration, 1)
    
    with open(file_path, 'w') as f:
        json.dump(agreement_results, f)
    
    return agreement_results


# Agreement curve as a function of k-mer size
def run_alignment_free(database, queries, query_type):
    num_thresholds = 10
    thresholds = [2 * i + 1 for i in range(num_thresholds)]

    file_path = f"saved_results/alignment_free_{query_type}.json"
    if not os.path.isfile(file_path):
        overall_scores = []

        best_K = None
        best_score = None

        local_alignment_results = get_local_alignment_results(database, queries, query_type)

        for K in thresholds:
            # print(f"Testing agreement for kmers of length {K}")
            database_kmers = ConvertLibaryToKmerSets(database, K)
            queries_kmers = ConvertLibaryToKmerSets(queries, K)
            scores = []
            for k in queries:
                alignment_best_match = local_alignment_results[k]
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
        
        with open(file_path, 'w') as f:
            results = {
                "overall_scores": overall_scores,
                "best_K": best_K
            }

            json.dump(results, f)

    else:
        with open(file_path) as f:
            results = json.load(f)
            overall_scores = results["overall_scores"]
            best_K = results["best_K"]

    plt.figure()
    plt.plot(thresholds, overall_scores, marker='o')
    plt.xlabel('K-mer length')
    plt.ylabel('Agreement')
    plt.title('K-mer vs Alignment Agreement')
    plt.ylim(0, 1)
    plt.savefig(f"graphs/agreement_against_kmer_size_{query_type}.png", bbox_inches='tight', dpi=300)
    plt.close()

    return best_K


def benchmark_alignment_free(database, queries, query_type, best_K):
    if check_benchmark("alignment_free"):
        return
    
    local_alignment_results = get_local_alignment_results(database, queries, query_type)

    start = time.time()
    database_kmers = ConvertLibaryToKmerSets(database, best_K)
    queries_kmers = ConvertLibaryToKmerSets(queries, best_K)
    scores = []
    for k in queries:
        alignment_best_match = local_alignment_results[k]
        _, kmer_best_match = KmerMatch(queries_kmers[k], database_kmers)
        if kmer_best_match == alignment_best_match:
            scores.append(1)
        else:
            scores.append(0)
    end = time.time()
    
    score = sum(scores) / len(scores)
    duration = end - start
    save_benchmark("alignment_free", duration, score)


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
    best_score = 0.0
    best_match = None
    
    #add your code here to find the best kmer match
    for k in database_kmers:
        if len(minimizer.intersection(library_minimizers[k][1])) >= minimum_minimizer_overlap:
            score = JaccardIndex(query_kmers, database_kmers[k])
            if score > best_score:
                best_score = score
                best_match = k

    return best_score, best_match


def run_minimizers(database, queries, query_type, BEST_K=15):
    num_window_sizes = 10
    window_sizes = [15 * i + 20 for i in range(num_window_sizes)]
    
    file_path = f"saved_results/minimizers_{query_type}.json"
    if not os.path.isfile(file_path):
        overall_scores = []

        best_m = None
        best_score = None

        # Calculate local align ground truth
        local_alignment_results = get_local_alignment_results(database, queries, query_type)

        # Calculate kmers for queries and database based on best_K from previous part
        database_kmers = ConvertLibaryToKmerSets(database, BEST_K)
        queries_kmers = ConvertLibaryToKmerSets(queries, BEST_K)

        for window_size in window_sizes:
            # print(f"Testing agreement for window size of length {window_size}")
            database_minimizers = calculate_database_minimizers(window_size, BEST_K, database)
            queries_minimizers = calculate_database_minimizers(window_size, BEST_K, queries)
            scores = []
            for k, sequence in queries.items():
                _, minimizer_best_match = MinimizerMatch(sequence, queries_minimizers[k][1], database_minimizers, queries_kmers[k], database_kmers)
                alignment_best_match = local_alignment_results[k]
                if minimizer_best_match == alignment_best_match:
                    scores.append(1)
                else:
                    scores.append(0)
            
            score = sum(scores) / len(scores)
            overall_scores.append(score)
            if best_score is None or score > best_score:
                best_m = window_size
                best_score = score
    
        with open(file_path, 'w') as f:
            results = {
                "overall_scores": overall_scores,
                "best_m": best_m
            }

            json.dump(results, f)

    else:
        with open(file_path) as f:
            results = json.load(f)
            overall_scores = results["overall_scores"]
            best_m = results["best_m"]


    plt.figure()
    plt.plot(window_sizes, overall_scores, marker='o')
    plt.xlabel('Window size')
    plt.ylabel('Agreement')
    plt.title('Window size vs Alignment Agreement')
    plt.ylim(0, 1)
    plt.savefig(f"graphs/agreement_against_window_size_{query_type}.png", bbox_inches='tight', dpi=300)
    plt.close()

    return best_m


def benchmark_minimizers(database, queries, query_type, best_K, best_m):
    finished = True
    for min_matches in [1,2,4,6]:
        if not check_benchmark(f"minimizers_{min_matches}"):
            finished = False
    if finished:
        return
    
    local_alignment_results = get_local_alignment_results(database, queries, query_type)

    for min_matches in [1,2,4,6]:
        start = time.time()
        database_kmers = ConvertLibaryToKmerSets(database, best_K)
        queries_kmers = ConvertLibaryToKmerSets(queries, best_K)
        database_minimizers = calculate_database_minimizers(best_m, best_K, database)
        queries_minimizers = calculate_database_minimizers(best_m, best_K, queries)
        scores = []
        for k, sequence in queries.items():
            _, minimizer_best_match = MinimizerMatch(sequence, queries_minimizers[k][1], database_minimizers, queries_kmers[k], database_kmers, minimum_minimizer_overlap=min_matches)
            alignment_best_match = local_alignment_results[k]
            if minimizer_best_match == alignment_best_match:
                scores.append(1)
            else:
                scores.append(0)
        score = sum(scores) / len(scores)
        end = time.time()
        duration = end - start

        save_benchmark(f"minimizers_{min_matches}", duration, score)


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


def compare_sequence_with_database(query_file, database):
    my_query_sequence = ""
    
    # constructing the sequence from the query_file
    try:
        with open(query_file, 'r') as e:
            next(e) 
            for line in e:
                my_query_sequence += line.strip()
    
    except FileNotFoundError:
        print("Not found, error")
        my_query_sequence = "" 

    else:
        print(f"Loaded query sequence, length: {len(my_query_sequence)})")

        K = 15 
        # convert to kmer set
        print("converting library to kmer set")
        database_kmers = ConvertLibaryToKmerSets(database, K)
        
        query_kmer_lib = ConvertLibaryToKmerSets({"my_query": my_query_sequence}, K)
        my_query_kmer_set = query_kmer_lib["my_query"]
        
        print("Matching k-mers...")
        (kmer_score, kmer_match) = KmerMatch(my_query_kmer_set, database_kmers)
        
        print("\n--- Results ---")
        print(f"Best K-mer Match: {kmer_match}")
        print(f"Score: {kmer_score}")
        

def plot_benchmark_graph():
    benchmark = {}
    with open("saved_results/benchmark.json") as f:
        benchmark = json.load(f)

    names = []
    durations = []
    accuracies = []
    for name, val in benchmark.items():
        duration, acc = val
        names.append(name)
        durations.append(duration)
        accuracies.append(acc)

    if len(durations) == 0:
        print("No benchmark entries to plot.")
        return

    plt.figure()
    plt.scatter(durations, accuracies)
    for i, label in enumerate(names):
        plt.annotate(label, (durations[i], accuracies[i]), xytext=(5, 2), textcoords='offset points', fontsize=8)

    plt.xlabel("Duration (s)")
    plt.ylabel("Max agreement")
    plt.ylim(0, 1)
    plt.title("Benchmark: Max Agreement vs Duration")

    out_path = "graphs/benchmark_plot.png"
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == "__main__":
    # Task 1/2
    fn = "bacterial_16s_genes.fa"
    database, queries = Load16SFastA(fn, fraction=0.5, database_size=10, query_size=10)

    print("Loaded %d 16s database sequences." % len(database))
    print("Loaded %d 16s query sequences." % len(queries))

    # Task 3/4
    print("Running alignment free for regular")
    best_K = run_alignment_free(database, queries, "regular")
    print("Benchmarking alignment free")
    benchmark_alignment_free(database, queries, "regular", best_K)

    # Task 5/6
    print("Running minimizers for regular")
    best_m = run_minimizers(database, queries, "regular", BEST_K=best_K)
    print("Benchmarking minimizer free")
    benchmark_minimizers(database, queries, "regular", best_K, best_m)

    # Task 7
    print("Running for illumina mutation")
    illumina_queries = illumina_mutation(queries)
    best_K_illumina = run_alignment_free(database, illumina_queries, "illumina")
    best_m_illumina = run_minimizers(database, illumina_queries, "illumina", BEST_K=best_K_illumina)

    # Task 7
    print("Running for nanopore mutation")
    nanopore_queries = nanopore_mutation(queries)

    # Task 10 comparing sequence to database
    database, _ = Load16SFastA(fn, fraction=1.0, database_size=20486, query_size=0)

    # file to read from
    query_fasta_file = 'Kangas0346_21_R1-16S-rRNA-seqR.fasta'
    compare_sequence_with_database(query_fasta_file, database)
    