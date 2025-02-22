import mmh3
import screed
import argparse

"""
Usage:
python /path/to/adapted_sourmash.py ksize scaled /path/to/input_file.fasta /path/to/output_file.txt
"""


parser = argparse.ArgumentParser(description="apply FracMinHashing to a DNA sequence and retain the kmers")
parser.add_argument("ksize", type=int, help="kmer length (an int)")
parser.add_argument("scaled", type=int, help="sampling rate")
parser.add_argument("file_path", type=str, help="path to fasta file")
parser.add_argument("output_file_path", type=str, help="path to the output txt file")
args = parser.parse_args()

ksize = args.ksize # 31 default in sourmash
scaled = args.scaled # 1000 default in sourmash

MAX_HASH = 2**64
keep_below = MAX_HASH / scaled
# print(keep_below)

"""
get all kmers - directly from sourmash
"""
def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i : i + ksize]
        kmers.append(kmer)

    return kmers


"""
read in file - directly from sourmash
"""
def read_kmers_from_file(filename):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    return all_kmers


"""
hash the kmers - adapted (barely) from sourmash
"""
def hash_kmer(kmer):
    # calculate the reverse complement
    rc_kmer = screed.rc(kmer)

    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer

    # calculate murmurhash using a hash seed of 42
    hash = mmh3.hash64(canonical_kmer, 42)[0]
    if hash < 0:
        hash += 2**64

     # done
    return canonical_kmer, hash


"""
filter the kmers - adapted from sourmash
"""
def subsample_kmers(kmers):
    keep = {}
    for kmer in kmers:
        canonical_kmer, hash_val = hash_kmer(kmer)
        if hash_val < keep_below:
            if canonical_kmer in keep:
                keep[canonical_kmer] += 1
            else:
                keep[canonical_kmer] = 1
    return keep


kmers = read_kmers_from_file(args.file_path)
keep = subsample_kmers(kmers)
# print(keep)
# print("total kmers:", len(kmers), "kept kmers", len(keep))
with open(args.output_file_path, "w") as output_file:
    for key, value in keep.items():
        output_file.write(f"{key} {value}\n")