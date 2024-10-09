//
// Created by miken on 10/5/2018.
//
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<string.h>
#include "kmers.h"
#include "myutils.h"

// Technically this is the maximum stride, so the max-K this module supports is 32*MAX_K = 512
#define MAX_K 16

/*
 * Function: matrix_index(i,j,w)
 * ------------------------
 * Takes coordinates i,j and a matrix-width w and returns a linear row-major coordinate
 * coresponding to (i,j) in the (H X W)-matrix (for any H).
 * */
long matrix_index(long i, long j, long w) {
    return (i*w+j);
}

/*
 * Functions: UPPERCASE(c1) and LOWERCASE(c1)
 * ------------------------------------------
 * Each one takes a single 'char' value and converts it to the appropriate case. If it is already in
 * the right case, the result is unchanged.
 *      UPPERCASE: if c1 > 90 (90=Z), then subtract 32
 *      LOWERCASE: if c1 >=65 (65=A) AND c1<=90 (90=Z), then add 32
 * */
char UPPERCASE(char c1) {
    return (c1 > 90 ? c1 - 32 : c1);
}
char LOWERCASE(char c1) {
    return (c1 >= 65 && c1 <= 90 ? c1 + 32 : c1 );
}
int k_to_ull_stride(int k) {
    return (int)ceil(k/32.);
}

/*
 * Function: ull_array_compare( w1, w2, pos)
 * -----------------------------------------
 * Comparison function for two ULL-arrays representing a k-mer. Essentially an 'is_greater_than' function, so
 * returns 1 if w1>w2, -1 if w1<w2 and 0 if they are equal. This particular function takes the argument 'pos' which
 * should be (stride-1) when called on two ULL-arrays because they are stored in little endian order when they are
 * multiple ULLs long. So when calling this function externally, call it with the highest position value first, and
 * the comparisons proceed down recursively to pos=0 in the case of equality, and return a result if not.
 * */
int ull_array_compare(unsigned long long *w1, unsigned long long *w2, int pos) {
    int diff;
    diff = w1[pos] < w2[pos] ? -1 : (w1[pos] > w2[pos] ? 1 : 0);
    // printf("pos=%d, w1[pos]=%llu, w2[pos]=%llu, diff=%d\n", pos, w1[pos], w2[pos], diff);
    return (pos==0 ? diff : (diff==0 ? ull_array_compare(w1, w2, pos-1) : diff));
}

/*
 * Function: inplace_reverse
 * -------------------------
 * Utility function to take a string pointer and reverse the string in-place.
 * */
void inplace_reverse(char * str)
{
    if (str)
    {
        char * end = str + strlen(str) - 1;

        // swap the values in the two given variables
        // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a,b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
    } while (0)

        // walk inwards from both ends of the string,
        // swapping until we get to the middle
        while (str < end)
        {
            XOR_SWAP(*str, *end);
            str++;
            end--;
        }
#   undef XOR_SWAP
    }
}

/*
 * Function: nucleotide_binary_lookup
 * ----------------------------------
 * Converts the nucleotide character to a binary value where:
 *      A = 0x00;  C = 0x01;  G = 0x02;  T = 0x03;
 *
 *  NOTE: if the input is not one of 'A','C','G','T','a','c','g', or 't' the function
 *        will return 0x00 by default.
 * */
char nucleotide_binary_lookup(char nuc) {
    switch (UPPERCASE(nuc)) {
        case 'A':
            return 0x00;
        case 'C':
            return 0x01;
        case 'G':
            return 0x02;
        case 'T':
            return 0x03;
        case 'U':
            return 0x03;
        default:
            return 0x00;
    }
}

/*
 * Function: kmer32_binary_RC
 * --------------------------
 * Basic binary RC operation on a kmer where k >=32 and k is stored in a single ULL.
 * */
unsigned long long kmer32_binary_RC(unsigned long long kmer, int k) {
    unsigned long long notKmer, rc;
    notKmer = ~kmer;
    rc = 0x00;
    int i;
    for (i=0; i<k; i++) {
        rc = (rc << 2) | (notKmer & 0x03);
        notKmer >>= 2;
    }
    return rc;
}

/*
 * Function: kmer_ull_array_moveright32(kmer_p1, kmer_p2, move_len)
 * ----------------------------------------------------------------
 * Takes the first two ULLs of 'kmer' (which correspond to the right-most and second to right-most 32-mers in the
 * sequence) and shifts the positions right by 'move_len', which must be less than or equal to 32. Effectively this
 * takes 'kmer_p2' and moves it right by 'move_len', then appends 'move_len' characters (i.e. 2-bits) from the right
 * side of kmer_p1 to the left of it. So the result is (in python notation):
 *      kmer_p2[-move_len:] + kmer_p1[:(32-move_len)]
 *
 * **This function is used notably in the RC operation.
 * */
unsigned long long kmer_ull_array_moveright32(unsigned long long kmer_p1, unsigned long long kmer_p2, int move_len) {
    unsigned long long res;
    res = (kmer_p1 >> (move_len * 2)) | (kmer_p2 << ((32-move_len)*2));
    return res;
}

// **********************************************************************
// *            HELPER FUNCTIONS:
// **********************************************************************

/*
 * Function: dna_to_binary
 * -----------------------
 * Converts a string representing a DNA sequence to an unsigned long long integer representing the first
 * 'k' characters of that sequence. We must have k<=32 and this is not tested for.
 * */
unsigned long long dna_to_binary(char *seq, int k) {
    int i;
    unsigned long long res;
    res = 0x00;
    for (i=0; i<k; i++) {
        res <<= 2;
        res |= nucleotide_binary_lookup(seq[i]);
    }
    return res;
}

/*
 * Function: dna_to_binary_ull_array
 * -----------------------
 * Converts a sequence to a binary k-mer representation for arbitrary size K, using an array of
 * unsigned long long, which must be provided and must be long enough. (No bounds checking). Here
 * the sequence 'seq' can be longer than 'k', but we are only going to return the ULL_array for
 * the first 'k'-bases of it.
 * */
void dna_to_binary_ull_array(char *seq, int k, unsigned long long *w) {
    int stride, k_remainder, pos;
    stride = k_to_ull_stride(k);
    k_remainder = k - (stride-1)*32;

    pos = 0;
    w[stride - pos - 1] = dna_to_binary(seq, k_remainder);
    for (pos=1; pos<stride; pos++) {
        w[stride - pos - 1] = dna_to_binary(&seq[k_remainder + (pos-1)*32], 32);
    }
}

/*
 * Function: binary_to_dna
 * -----------------------
 * Converts the k-mer <value> to an ascii sequence of (uppercase) nucleotide letters. Note:
 *          1) the variable 'seq' must be allocated to a string of length <seqlen>.
 *          2) seqlen must be no larger than 32.
 * */
void binary_to_dna(unsigned long long binseq, int seqlen, char *seq) {
    char nucs[] = {'A', 'C', 'G', 'T'};
    int i;
    for (i=seqlen-1; i>=0; i--) {
        seq[i] = nucs[binseq & 0x03];
        binseq >>= 2;
    }
}

/*
 * Takese a ULL array plus a k and populates a string called 'seq' which must be length
 * at least (k+1).
 * */
void binary_ull_array_to_dna(unsigned long long *w, int k, char *seq) {
    int stride, k_remainder, i;
    stride = k_to_ull_stride(k);
    k_remainder = k - (stride-1)*32;
    // do the leftmost byte with the remainder:
    binary_to_dna(w[stride-1], k_remainder, &seq[0]);
    // do the rest:
    for (i=stride-2; i>=0; i--) {
        binary_to_dna(w[i], 32, &seq[k_remainder + 32 * i]);
    }
    seq[k]='\0';
}


/*
 * Function: binary_ull_array_pop_append_nucleotide
 * ------------------------------------------------
 * Takes a ULL array 'w', a 'k' value and a new nucleotide and fills the ULL array 'w_new'
 * with the new ULL for the kmer representing the same sequence with the first base removed
 * and the new base appended to the end.
 * */
void binary_ull_array_pop_append_nucleotide(unsigned long long *w, int k, char N, unsigned long long *w_new) {
    int stride, k_remainder, curr_pos;
    unsigned long long remainder_mask;
    stride = k_to_ull_stride(k);
    k_remainder = k - (stride-1)*32;
    remainder_mask = rightmask(k_remainder*2);
    curr_pos = stride-1; // starting at the rightmost ULL
    // If it's only one ULL, simple operation:
    if (k<=32) {
        w_new[0] = ((w[0]<<2) | nucleotide_binary_lookup(N) ) & remainder_mask;
        return;
    }
    // Otherwise, do the remainder first with the masK:
    w_new[curr_pos] = ((w[curr_pos]<<2) | (w[curr_pos-1]>>62)) & remainder_mask;
    curr_pos--;
    // .... then the rest:
    while (curr_pos>0) {
        w_new[curr_pos] = ((w[curr_pos]<<2) | (w[curr_pos-1]>>62));
        curr_pos--;
    }
    w_new[curr_pos] = ((w[curr_pos]<<2) | nucleotide_binary_lookup(N));
}

/*
 * Function: binary_ull_array_reverse_complement
 * ---------------------------------------------
 * Takes the ULL array 'w' and length k and populates a ULL array named 'w_rc' with the reverse
 * complement also in ULL array form, without going back to the character space.
 * */
void binary_ull_array_reverse_complement(unsigned long long *w, int k, unsigned long long *w_rc) {
    int stride; stride=k_to_ull_stride(k);
    int i, k_remainder;
    unsigned long long temp, mymask;
    k_remainder = k-(stride-1)*32;
    mymask = rightmask(k_remainder*2);
    // For the first stride-1 bytes, use the right-move operation:
    for (i=0; i<stride-1; i++) {
        temp = kmer_ull_array_moveright32(w[stride-2-i], w[stride-1-i], k_remainder);
        w_rc[i] = kmer32_binary_RC(temp, 32);
    }
    // For the last one, get the RC of the right (k_remainder*2) bits:
    w_rc[stride-1] = kmer32_binary_RC(w[0] & mymask, k_remainder);
}


// **********************************************************************
// *            MAIN FUNCTIONS:
// **********************************************************************

/*
 * Function: dna_sequence_to_kmer_ull_array_list
 * ---------------------------------------------
 * Takes a sequence and returns an array that is [len(seq)-k+1]*stride long and has a list
 * of the kmers occuring in the sequence in order.
 * */
unsigned long long *dna_sequence_to_kmer_ull_array_list(char *seq, size_t seqlen, int k) {
    int stride; //, seqlen;
    stride = k_to_ull_stride(k);
    unsigned long long *w = (unsigned long long *) malloc((seqlen-k+1)*stride*sizeof(unsigned long long));

//    dna_to_binary_ull_array(&seq[0], k, &w[0]);
//    for (pos=k; pos<seqlen; pos++) {
//        binary_ull_array_pop_append_nucleotide(&w[(pos-k)*stride],k,seq[pos],&w[(pos-k+1)*stride]);
//    }
    dna_sequence_to_ull_array_list_provided(seq, seqlen, k, w);
    return w;
}

/*
 * Function: dna_sequence_to_kmer_ull_array_list_provided
 * ---------------------------------------------
 * Takes a sequence and **populates** an array that is [len(seq)-k+1]*stride long and has a list
 * of the kmers occuring in the sequence in order. In this case that array is provided rather than
 * created from scratch (as it is in the function above). The function above creates W, calls this
 * one, then returns W.
 * */
void dna_sequence_to_ull_array_list_provided(char *seq, size_t seqlen, int k, unsigned long long *w) {
    int pos, stride; //, seqlen;
    stride = k_to_ull_stride(k);
    dna_to_binary_ull_array(&seq[0], k, &w[0]);
    for (pos=k; pos<seqlen; pos++) {
        binary_ull_array_pop_append_nucleotide(&w[(pos-k)*stride],k,seq[pos],&w[(pos-k+1)*stride]);
    }
}

/*
 * Function: ull_array_list_set_to_min_RC
 * --------------------------------------
 * Takes an array of binary_kmer_ull_arrays along with 'k' and goes down the list, setting every entry
 * to be the minimimum of itself and it's reverse complement.
 * */
unsigned long ull_array_list_set_to_min_RC(unsigned long long *w, int k, size_t array_len) {
    int stride; stride=k_to_ull_stride(k);
    unsigned long i, num_swaps; unsigned char j;
    unsigned long long temp_rc[MAX_K];
    zero_out_ull(temp_rc,MAX_K);
    num_swaps = 0;

    // Go down the array and replace them one by one:
    for (i=0; i<array_len; i++) {
        binary_ull_array_reverse_complement(&w[stride*i], k, &temp_rc[0]);
        // if the comparison is > 0, meaning that w is GREATER than w_rc, so replace with w_rc...
        if (ull_array_compare(&w[stride*i], &temp_rc[0], stride-1) > 0) {
            // replace w[stride * i] with values from temp_rc
            for (j=0; j<stride; j++) w[stride*i + j] = temp_rc[j];
            num_swaps++;
        }
        zero_out_ull(temp_rc,stride);  //not positive this is necessary but we're re-using the variable...
    }
    return num_swaps;
}


// **********************************************************************
// *            DEBUGGING FUNCTIONS:
// **********************************************************************

/*
 * DEBUGGING FUNCTIONS:
 * */
void DEBUG_print_int_matrix(int *i_mat, int nrows, int ncols) {
    printf("**** DEBUGGING: printing integer matrix size (%d x %d):\n", nrows, ncols);
    int i, j;
    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++) {
            printf("%4d ", i_mat[i*ncols + j]);
            if (j > 0 && (j % 20)==0) printf("\n  ...");
        }
        printf("\n");
    }
}

void DEBUG_print_long_matrix(long *l_mat, int nrows, int ncols) {
    printf("**** DEBUGGING: printing LONG matrix size (%d x %d):\n", nrows, ncols);
    int i, j;
    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++) {
            printf("%4ld ", l_mat[i*ncols + j]);
            if (j > 0 && (j % 20)==0) printf("\n  ...");
        }
        printf("\n");
    }
    printf("\n");
}

void DEBUG_print_w_row(unsigned long long *w, int row, int stride) {
    int i,pos; pos = row*stride;
    printf("(%llu",w[pos]);
    for (i=1; i<stride; i++) {printf(",%llu", w[pos+i]);}
    printf(")\n");
}
// reverses order of first k elements of w:
void test_move_matrix_elements(unsigned long long *w, int k) {
    int i;
    unsigned long long temp;
    for (i=0; i<k-i-1; i++) {
        temp = w[i];
        w[i] = w[k-i-1];
        w[k-i-1] = temp;
    }
}