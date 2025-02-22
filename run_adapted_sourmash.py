import os
import subprocess
import argparse

"""
Usage:
python /path/to/run_adapted_sourmash.py /path/to/fasta/files/directory ksize scaled
"""

def run_adapted_sourmash(ksize, scaled, file_path, output_file_path):
    command = ["python", "adapted_sourmash.py", str(ksize), str(scaled), file_path, output_file_path] # put the two files in the same dir
    subprocess.run(command, check=True)

def parse_args():
    parser = argparse.ArgumentParser(description="run adapted_sourmash on all fasta files in a directory.")
    parser.add_argument("input_dir", type=str, help="directory containing the .fasta files.")
    parser.add_argument("ksize", type=int, help="kmer length (an int)")
    parser.add_argument("scaled", type=int, help="sampling rate")
    return parser.parse_args()

def main():
    args = parse_args()
    for filename in os.listdir(args.input_dir):
        if filename.endswith(".fasta"):
            fpath = os.path.join(args.input_dir, filename)
            out_fpath = os.path.join(args.input_dir, f"{os.path.splitext(filename)[0]}.txt")

            run_adapted_sourmash(args.ksize, args.scaled, fpath, out_fpath)
            print(f"finished signature for: {out_fpath}")

if __name__ == "__main__":
    main()
