from phylogeny_utilities.utilities import *
import numpy as np
import kmers

separate_tests = lambda : print('\n------------------------------------------------\n')

tseq = "ATCAGAAGAATCATCGGTCTGGTTTTCATTCTCCTGTCTGGTAAGCGAAGGATTTCCGGGATAAGAAATATATTCTACAAATTTTCTGACGTCACTGATGTCATCTGCATCTTTTGCATATCCCACTGCCATATACAGGCTATTTCCCTCTGAATAAATTGTCTGATAGCAGGCATCGTATGTGTCATTGCATTTGGCAATCCC"

w=np.zeros(10,dtype=np.uint64)
kmers.seq_to_kmer_binary(tseq,34,w)
print(w)
w_shouldbe='''
[5197386428564066165,                   3,                   0,
                   0,                   0,                   0,
                   0,                   0,                   0,
                   0], dtype=uint64)
'''
print(f'should be: {w_shouldbe}')

separate_tests()

w0_rc = kmers.test_ull_binary_rc(5197386428564066165, 32)
s1 = binary_to_dna(w0_rc, 32)
s2 = rc(tseq[2:34])
print(f's1: {s1}\ns2: {s2}')
print(f's1 == s2: {s1==s2}  (should be True)')


