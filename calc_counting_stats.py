import os
import sys
from itertools import combinations

"""
Usage: python calc_counting_stats.py path/to/kmer/counts/dir
"""

def make_one_set(file_path):
    kmers = set({})
    counts = []

    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            kmers.add(columns[0])
            counts.append(int(columns[1]))

    return kmers, counts


def calc_stats(directory_path):
    all_kmer_sets = []
    all_counts_list = []
    for fname in os.listdir(directory_path):
        if fname.endswith('.txt'):
            fpath = os.path.join(directory_path, fname)
            kmers, counts = make_one_set(fpath)
            all_kmer_sets.append(kmers)
            all_counts_list.append(counts)

    calc_count_stats(all_counts_list)
    calc_kmer_stats(all_kmer_sets)


def calc_count_stats(all_counts_list):
    aves = []
    for counts_list in all_counts_list:
        ave = sum(counts_list) / len(counts_list)
        aves.append(ave)

    overall_ave = sum(aves) / len(aves)
    print("average kmer count: ", overall_ave)


def calc_kmer_stats(all_kmer_sets):
    # find the average pairwise intersection length
    intersection_lens = []
    
    for kmer_set_1, kmer_set_2 in combinations(all_kmer_sets, 2):
        intersection = kmer_set_1 & kmer_set_2
        intersection_lens.append(len(intersection))
    
    avg_intersection_len = sum(intersection_lens) / len(intersection_lens)
    print("average pairwise intersection length: ", avg_intersection_len)

    # find the length of the intersection between all kmer sets
    intersection = all_kmer_sets[0]
    intersection_increments = []
    for kmer_set in all_kmer_sets[1:]:
        intersection &= kmer_set
        intersection_increments.append(len(intersection))
    
    print("number of kmers all files have in common: ", len(intersection))
    with open('intersection_increments.txt', 'w') as f:
        for item in intersection_increments:
            f.write(f"{item}\n")

    # find the average length of the kmer sets
    tot_len = sum(len(kmer_set) for kmer_set in all_kmer_sets)
    ave_len = tot_len / len(all_kmer_sets)

    print("average number of kmers: ", ave_len)

    unique_kmers = set().union(*all_kmer_sets)
    print("number of feature matrix columns: ", len(unique_kmers))
    

if __name__ == "__main__":
    dir_path = sys.argv[1]

    if not os.path.isdir(dir_path):
        print(f"not valid directory")
        sys.exit(1)

    calc_stats(dir_path)
