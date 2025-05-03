import pandas as pd

eval_data_kmers = "path_to_diabimmune_column_kmers.txt"
training_data_kmers = "path_to_training_data_column_kmers.txt"
eval_data_fm = "path_to_diabimmune_feature_matrix.csv"

# get the diabimmune kmers
with open(eval_data_kmers, 'r') as f:
    column_names = [line.strip() for line in f]

# read in the diabimmune feature matrix
df = pd.read_csv(eval_data_fm, header=None)
df.columns = column_names

# print(df.head())

# get the training kmers
with open(training_data_kmers, 'r') as f:
    training_kmers = [line.strip() for line in f]

df_adjusted = df.copy()
df_adjusted_colset = set(df_adjusted.columns)
training_kmers_set = set(training_kmers)

# add missing columns to the diabimmune feature matrix 
# (kmers that are present in training_kmers but not in df)
# so that the feature matrix maintains the same order as in training
missing_cols = [kmer for kmer in training_kmers_set if kmer not in df_adjusted_colset]
if missing_cols:
    zeros_df = pd.DataFrame(0, index=df_adjusted.index, columns=missing_cols)
    df_adjusted = pd.concat([df_adjusted, zeros_df], axis=1)

df_adjusted_colset = set(df_adjusted.columns)

# drop extra columns 
# (kmers that are present in df but not in training_kmers)
# because the model did not learn anything about these so they'll be
# meaningless and also mess up the kmer order
extra_cols = [col for col in df_adjusted_colset if col not in training_kmers_set]
df_adjusted.drop(columns=extra_cols, inplace=True)

# reorder the columns to match training_kmers
df_adjusted = df_adjusted[training_kmers]


df_adjusted.to_csv('diabimmune_feature_matrix_output_path.csv', index=False, header=False) #TODO: replace path