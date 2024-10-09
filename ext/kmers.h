//
// Created by miken on 10/5/2018.
//

#ifndef EXT_KMERS_H
#define EXT_KMERS_H


// Struct definition:


// Utility functions:
char UPPERCASE(char c1);
char LOWERCASE(char c1);
int k_to_ull_stride(int k);
int ull_array_compare(unsigned long long *w1, unsigned long long *w2, int pos);
void inplace_reverse(char * str);
char nucleotide_binary_lookup(char nuc);
long matrix_index(long i, long j, long w);
unsigned long long kmer32_binary_RC(unsigned long long kmer, int k);
unsigned long long kmer_ull_array_moveright32(unsigned long long kmer_p1, unsigned long long kmer_p2, int move_len);


// Helper functions:
unsigned long long dna_to_binary(char *seq, int k);
void dna_to_binary_ull_array(char *seq, int k, unsigned long long *w);
void binary_to_dna(unsigned long long binseq, int seqlen, char *seq);
void binary_ull_array_to_dna(unsigned long long *w, int k, char *seq);
void binary_ull_array_pop_append_nucleotide(unsigned long long *w, int k, char N, unsigned long long *w_new);
void binary_ull_array_reverse_complement(unsigned long long *w, int k, unsigned long long *w_rc);

//void kmer_ull_array_min_of_self_RC(unsigned long long *w, int k);
//void kmer_ull_array_replace_with_RC(unsigned long long *w, int k);
// Main functions:
unsigned long long *dna_sequence_to_kmer_ull_array_list(char *seq, size_t seqlen, int k);
void dna_sequence_to_ull_array_list_provided(char *seq, size_t seqlen, int k, unsigned long long *w);
unsigned long ull_array_list_set_to_min_RC(unsigned long long *w, int k, size_t array_len);

/*
 * DEBUGGING FUNCTIONS:
 * */
void DEBUG_print_int_matrix(int *i_mat, int nrows, int ncols);
void DEBUG_print_long_matrix(long *l_mat, int nrows, int ncols);
void DEBUG_print_w_row(unsigned long long *w, int row, int stride);
void test_move_matrix_elements(unsigned long long *w, int k);

#endif //EXT_KMERS_H
