import os
import csv
import sys

def make_feature_matrix(directory):
    """
    makes feature matrix from a directory
    """
    all_kmers = set() # will be columns of mat
    file_data = []
    file_names = []
    
    """
    go through each file and compile a list of kmers present across all of them.
    at the same time, make the file contents faster to access (dict)
    """
    for fname in os.listdir(directory):
        if fname.endswith(".txt"):
            fpath = os.path.join(directory, fname)
            file_names.append(fname)
            file_content = {}
            
            with open(fpath, 'r') as f:
                for line in f:
                    """
                    because each file is in the format
                    kmer count
                    kmer count
                    kmer count
                    ...
                    """
                    parts = line.split()
                    if len(parts) == 2:
                        kmer, count = parts
                        count = int(count)
                        file_content[kmer] = count
                        all_kmers.add(kmer)
            
            file_data.append(file_content)
    
    # sort the kmers alphabetically
    all_kmers = sorted(all_kmers)
    
    matrix = []
    for file_content in file_data:
        row = []
        for kmer in all_kmers:
            row.append(file_content.get(kmer, 0))  # defaults to 0 if the kmer isn't present
        matrix.append(row)
    print("finished making unfiltered matrix")

    # remove singletons (only one row has a nonzero value for that kmer)
    cols_to_remove = set({})
    for j in range(len(all_kmers)):
        nonzero_count = sum(1 for i in range(len(matrix)) if matrix[i][j] > 0)
        if nonzero_count == 1:
            cols_to_remove.add(j)
    print("finished finding columns to remove")

    # remove the singleton columns from the matrix and update all_kmers
    matrix_filtered = []
    all_kmers_filtered = []
    for i in range(len(matrix)):
        new_row = [matrix[i][j] for j in range(len(all_kmers)) if j not in cols_to_remove]
        matrix_filtered.append(new_row)

    all_kmers_filtered = [all_kmers[j] for j in range(len(all_kmers)) if j not in cols_to_remove]
    print("finished filtering")
    
    # write files
    output_csv = os.path.join(directory, "feature_matrix_not_normalized.csv")
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i, row in enumerate(matrix_filtered):
            writer.writerow(row)
    
    out_fnames = os.path.join(directory, "row_fnames.txt")
    with open(out_fnames, 'w') as row_fnames_file:
        for fname in file_names:
            row_fnames_file.write(fname + "\n")
    
    out_columns = os.path.join(directory, "column_kmers.txt")
    with open(out_columns, 'w') as column_kmers_file:
        for kmer in all_kmers_filtered:
            column_kmers_file.write(kmer + "\n")
    
    print("finished making feature matrix (unnormalized)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python script.py /path/to/adapted_sourmash_output_dir")
        sys.exit(1)
    
    directory = sys.argv[1]
    
    if not os.path.isdir(directory):
        print(f"director dne")
        sys.exit(1)
    
    make_feature_matrix(directory)
