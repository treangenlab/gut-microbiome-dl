//
// Created by miken on 3/22/2024.
//

#ifndef PYTHON_EXTENSION_MODULE_EXAMPLE_MYUTILS_H
#define PYTHON_EXTENSION_MODULE_EXAMPLE_MYUTILS_H
#include <stdio.h>

/*
 * Function: leftmask
 * -------------------
 * returns a 64-bit unsigned int with the rightmost 'n_bits' bits set to 1 and the
 * rest set to 0.
 * */
unsigned long long rightmask(unsigned int n_bits) {
    unsigned long long mask = ~0x00;
    int i;
    for (i=0; i<64-n_bits; i++) {
        mask >>= 1;
    }
    return mask;
}

unsigned long long leftmask(unsigned int n_bits) {
    unsigned long long mask = ~0x00;
    int i;
    for (i=0; i<64-n_bits; i++) {
        mask <<= 1;
    }
    return mask;
}

void zero_out_ull(unsigned long long *w, unsigned int w_len) {
    unsigned int i;
    for (i=0; i<w_len; i++) w[i]=0x00;
}

#endif //PYTHON_EXTENSION_MODULE_EXAMPLE_MYUTILS_H
