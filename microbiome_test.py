# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import alignment
import random
import matplotlib.pyplot as plt
import seaborn as sns
import time
import json
import os
import pandas as pd
import numpy as np
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
    for i, (k, sequence) in enumerate(queries.items()):
        if i % 5 == 0:
            print(f"Finished aligning {i}/{len(queries)} queries")
        _, alignment_best_match = AlignmentMatch(sequence, database)
        agreement_results[k] = alignment_best_match
    print(f"Finished aligning {len(queries)}/{len(queries)} queries")
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
    plt.title('Agreement of alignment free sequence matching against K-mer length')
    plt.ylim(0, 0.65)
    plt.xticks(thresholds)
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
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

    window_kmer_count = window_size - k + 1
    num_windows = len(kmers) - window_kmer_count + 1
    
    minimizers = set()
    for wstart in range(num_windows):
        window_kmers = kmers[wstart : wstart + window_kmer_count]
        minimizer = min(window_kmers)
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
    with open("benchmark.json") as f:
        benchmark = json.load(f)

    names = []
    durations = []
    accuracies = []
    for name, val in benchmark.items():
        duration, acc = float(val[0]), float(val[1])
        names.append(name)
        durations.append(duration)
        accuracies.append(acc)

    if len(durations) == 0:
        print("No benchmark entries to plot.")
        return

    pos = np.arange(len(names))
    width = 0.35

    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Blue bars: durations (log scale)
    bars1 = ax1.bar(pos - width/2, durations, width, color='tab:green', label='Total runtime')
    ax1.set_yscale('log')
    ax1.set_ylabel("Total runtime (s, log scale)")
    ax1.set_xticks(pos)
    ax1.set_ylim(min(durations) * 0.5, max(durations) * 5)
    ax1.set_xticklabels(names, rotation=30, ha='right')

    ax2 = ax1.twinx()
    bars2 = ax2.bar(pos + width/2, accuracies, width, color='tab:red', label='Max agreement')
    ax2.set_ylim(0, 1.2)
    ax2.set_ylabel("Max agreement")

    # Combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper right')

    # Annotate bars
    for i, d in enumerate(durations):
        ax1.text(pos[i] - width/2, d * 1.05, f"{d:.2f}s", ha='center', va='bottom', fontsize=8)
    for i, a in enumerate(accuracies):
        ax2.text(pos[i] + width/2, a + 0.02, f"{a:.2f}", ha='center', va='bottom', fontsize=8)

    plt.title("Benchmark runtimes and max agreements for different methods")
    plt.grid(True, which='both', ls='--', lw=0.5, alpha=0.3)

    os.makedirs("graphs", exist_ok=True)
    out_path = "graphs/benchmark_plot.png"
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()
    print("Saved benchmark plot to", out_path)


def get_results(name):
    with open(f"saved_results/{name}.json", 'r') as f:
        results = json.load(f)
    return results


def plot_mutations_alignment_free():
    regular = get_results("alignment_free_regular")["overall_scores"]
    illumina = get_results("alignment_free_illumina")["overall_scores"]
    nanopore = get_results("alignment_free_nanopore")["overall_scores"]
    thresholds = [2 * i + 1 for i in range(len(regular))]

    data = pd.DataFrame({
        'x': thresholds,
        'No mutation': regular,
        'Illumina': illumina,
        'Nanopore': nanopore
    })

    plt.figure()
    sns.lineplot(data=data, x='x', y='No mutation', label='No mutation', marker='o')
    sns.lineplot(data=data, x='x', y='Nanopore', label='Nanopore', marker='o')
    sns.lineplot(data=data, x='x', y='Illumina', label='Illumina', marker='o')
    plt.xlabel('K-mer length')
    plt.ylabel('Agreement')
    plt.title('Agreement of alignment free sequence matching against k-mer length with different mutations')
    plt.ylim(0, 0.65)
    plt.xticks(thresholds)
    plt.savefig(f"graphs/mutations_afsm.png", bbox_inches='tight', dpi=300)
    plt.close()

    
def plot_mutations_minimizers():
    regular = get_results("minimizers_regular")["overall_scores"]
    illumina = get_results("minimizers_illumina")["overall_scores"]
    nanopore = get_results("minimizers_nanopore")["overall_scores"]
    window_sizes = [15 * i + 20 for i in range(10)]

    data = pd.DataFrame({
        'x': window_sizes,
        'No mutation': regular,
        'Illumina': illumina,
        'Nanopore': nanopore
    })


    plt.figure()
    sns.lineplot(data=data, x='x', y='No mutation', label='No mutation', marker='o')
    sns.lineplot(data=data, x='x', y='Nanopore', label='Nanopore', marker='o')
    sns.lineplot(data=data, x='x', y='Illumina', label='Illumina', marker='o')
    plt.xlabel('Window Size')
    plt.ylabel('Agreement')
    plt.title('Agreement of minimizers matching against window size with different mutations')
    plt.ylim(0, 0.8)
    plt.xticks(window_sizes)
    plt.savefig(f"graphs/mutations_minimizers.png", bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == "__main__":
    # Task 1/2
    fn = "bacterial_16s_genes.fa"
    database, queries = Load16SFastA(fn, fraction=0.5)

    print("Loaded %d 16s database sequences." % len(database))
    print("Loaded %d 16s query sequences." % len(queries))

    # Task 3/4
    print("\nRunning alignment free for regular")
    best_K = run_alignment_free(database, queries, "regular")
    print("\nBenchmarking alignment free")
    # benchmark_alignment_free(database, queries, "regular", best_K)

    # Task 5/6
    print("\nRunning minimizers for regular")
    best_m = run_minimizers(database, queries, "regular", BEST_K=best_K)
    print(f"\nBenchmarking minimizer free with best k {best_K}, best_m {best_m}")
    # benchmark_minimizers(database, queries, "regular", best_K, best_m)

    # Task 7
    print("\nRunning for illumina mutation")
    illumina_queries = illumina_mutation(queries)
    best_K_illumina = run_alignment_free(database, illumina_queries, "illumina")
    best_m_illumina = run_minimizers(database, illumina_queries, "illumina", BEST_K=best_K_illumina)

    # Task 7
    print("\nRunning for nanopore mutation")
    nanopore_queries = nanopore_mutation(queries)
    best_K_nanopore = run_alignment_free(database, nanopore_queries, "nanopore")
    best_m_nanopore = run_minimizers(database, nanopore_queries, "nanopore", BEST_K=best_K_nanopore)

    plot_mutations_alignment_free()
    plot_mutations_minimizers()
    plot_benchmark_graph()

    # Task 10 comparing sequence to database
    # print("\nComparing sequence to database")
    # database, _ = Load16SFastA(fn, fraction=1.0, database_size=20486, query_size=0)

    # # file to read from
    # query_fasta_file = 'Kangas0346_21_R1-16S-rRNA-seqR.fasta'
    # compare_sequence_with_database(query_fasta_file, database)
    