# gut-microbiome-dl
## In progress!

## File Guide:
* ```model/```
  * ```autoencoder_final.ipynb``` - final autoencoder architecture. Use this file if you want to train a model on a feature matrix, and then save the model as a .pth file.
  * ```get_embeddings_from_autoencoder.ipynb``` - Use this file if you want to load an existing model (with the correct architecture) and run data (in the form of a feature matrix) through it to get the ebeddings. Saves the embeddings as a .csv file.
* ```scripts/```
  * ```feature_matrix_scripts/```
    * ```adapted_sourmash.py``` - counts k-mers in a mash sketch for a single fasta file. outputs them as a text file, where each row is in the format "kmer #".
    * ```run_adapted_sourmash.py``` - runs ```adapted_sourmash.py``` on a whole directory. The two files must be in the same directory.
    * ```aggregate_adapted_sourmash_results.py``` - makes a feature matrix from the kmer counts. This feature matrix has dropped singletons and is not normalized. It will be normalized when loaded into either file in ```model/```.
    * ```prep_eval_feature_matrix.py``` - makes sure that the feature matrix created for any data used in evaluation has the correct k-mers as columns. Drops k-mers that were not present in the training data and adds empty columns for k-mers present in the training data but not evaluation data.
    * ```calc_counting_stats.py``` - an extra file that takes in a directory of k-mer counting files (generated from ```run_adapted_sourmash.py```) and outputs some statistics about them.
  * ```analysis_visualization_scripts/```
    * ```visualize_training_embeddings.py``` - script for visualizing the training data embeddings. This outputs the hierarchical clustering plots, with samples colored by hardcoded categories.
    * ```visualize_diabimmune_embedding_data.py``` - script for visualizing the diabimmune data from the saved model (doesn't generalize to any trained model or any evaluation data, as the node numbers and metadata fields are hardcoded). This outputs the scatterplots for node vs participant age (and significance info) and the stripplots for node activations colored by abx exposure (and significance info).
* ```data/```
  * ```column_kmers.txt``` - the k-mers used as the columns of the feature matrix used to train the model (and will also be the column names, in order, of any evaluation feature matrix). These are in the same order as the columns of the feature matrix. This file is needed for running ```prep_eval_feature_matrix.py```. This file was created by ```aggregate_adapted_sourmash_results.py```.
  * ```row_fnames_training.txt``` - the row names (file names) of the feature matrix used to train the model, in order. This file will be needed to run ```visualize_training_embeddings.py```, but you may need to remove the ".txt"s first. This file was created by ```aggregate_adapted_sourmash_results.py```.
  * ```row_fnames_diabimmune.txt``` - the row names (file names) of the feature matrix made from the diabimmune data, in order. This file will be needed to run ```visualize_diabimmune_embedding_data.py```. This file was created by ```aggregate_adapted_sourmash_results.py```.
  * ```embeddings_training.csv``` - the embeddings of the training samples, in the same order as the training feature matrix (and therefore as ```row_fnames_training.txt```). This file was created by ```get_embeddings_from_autoencoder.ipynb```.
  * ```embeddings_diabimmune.csv``` - the embeddings of the diabimmune samples, in the same order as the diabimmune feature matrix (and therefore as ```row_fnames_diabimmune.txt```). This file was created by ```get_embeddings_from_autoencoder.ipynb```.
  * ```diabimmune_t1d_wgs_metadata.csv``` - the diabimmune study metadata, taken from the website
  * ```diabimmune_metadata_and_embeddings_merged.csv``` - the result of merging the diabimmune metadata with the embeddings, where each embedding is assigned to the correct participant. This file will be needed to run ```visualize_diabimmune_embedding_data.py```.

Note that paths in some files will need to be updated based on where you are storing the data - see TODOs within the code files.

## General Workflows:
### If you want to train the model and then evaluate:
### If you want to evaluate an existing model on new data:
