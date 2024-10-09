//
// Created by miken on 3/22/2024.
//
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "kmers.h"

#define DEBUG_ON 0

/* Docstrings */
//#include "_nwops_docstr.h"

/* Available functions */
//static PyObject *nwops_testPrintStringSizes(PyObject *self, PyObject *args);
static PyObject *kmers_getVersionString(PyObject *self, PyObject *args);
static PyObject *kmers_test_InplaceNumpyMove(PyObject *self, PyObject *args);
static PyObject *kmers_test_binaryRC(PyObject *self, PyObject *args);
static PyObject *kmers_test_kmer_moveright(PyObject *self, PyObject *args);
//static PyObject *kmers_test_kmerull_inplace_rc(PyObject *self, PyObject *args);
static PyObject *kmers_test_kmer_ull_to_dna(PyObject *self, PyObject *args);
static PyObject *kmers_test_kmer_ull_array_compare(PyObject *self, PyObject *args);
static PyObject *kmers_seqToKmer(PyObject *self, PyObject *args);
static PyObject *kmers_seq_to_kmers_ull_array(PyObject *self, PyObject *args);
static PyObject *kmers_seqiter_to_kmers_ull_array(PyObject *self, PyObject *args);
static PyObject *kmers_ull_array_batch_set_min_RC(PyObject *self, PyObject *args);

/* Docstrings */
static char version_str[] = "Version 0.1 (beta)";
static char get_version_docstring[] =
        "Just prints the version number of this module. For easy testing that import is ok.                  \n";
static char module_docstring[] =
        "Provides some methods to work with kmers                                                            \n";
static char seq_to_kmer_binary_docstring[] =
        "Takes a sequence, a value of k, and a numpy uint64 array and fills the numpy array with the binary  \n"
        "representation of the k-mer. Works for arbitrary size k.                                            \n"
        "";
static char seq_to_kmers_ull_array_docstring[] =
        "Takes a sequence and a value of k and converts it to a numpy array of tuples of UINT64s where the   \n"
        "size of the tuple is given by stride, which is roundup(k/32). Optionally also accepts a numpy       \n"
        "UINT64 array of sufficient size to hold all the kmers, which is then populated. In either case the  \n"
        "final array is returned by the function.                                                            \n"
        "";
static char seqiter_to_kmers_ull_array_docstring[] =
        "Similar to the function 'seq_to_kmer_ull_array' but in this case it accepts an iterator of          \n"
        "rather than a single sequence. The resulting kmer array will include all sequences combined.        \n"
        "";
static char ull_arr_batch_min_RC_docstring[] =
        "Function: kmer_ull_array_batch_set_min_RC(w_array, k, array_length)                                 \n"
        "Takes a completed array of kmers in binary ULL format and sets every entry to be the minimum of     \n"
        "itself and its reverse-complement. Since we are doing a lot of counting unique k-mers this cuts down\n"
        "on the total kmer counts and is a more parsimonious representation.                                 \n"
        "";


/* Module specification */
static PyMethodDef module_methods[] = {
        {"get_version", kmers_getVersionString, METH_VARARGS, get_version_docstring},
        {"seq_to_kmer_binary", kmers_seqToKmer, METH_VARARGS, seq_to_kmer_binary_docstring},
        {"seq_to_kmer_ull_array", kmers_seq_to_kmers_ull_array, METH_VARARGS, seq_to_kmers_ull_array_docstring},
        {"seqiter_to_kmer_ull_array", kmers_seqiter_to_kmers_ull_array, METH_VARARGS, seqiter_to_kmers_ull_array_docstring},
        {"kmer_ull_array_batch_set_min_RC", kmers_ull_array_batch_set_min_RC, METH_VARARGS, ull_arr_batch_min_RC_docstring},
        // Testing functions w/ no docstring for now...
        {"test_inplace_numpy_replace", kmers_test_InplaceNumpyMove, METH_VARARGS, seq_to_kmer_binary_docstring},
        {"test_ull_binary_rc", kmers_test_binaryRC, METH_VARARGS, NULL},
        {"test_ull_move_right", kmers_test_kmer_moveright, METH_VARARGS, NULL},
        // {"test_ull_full_inplace_rc", kmers_test_kmerull_inplace_rc, METH_VARARGS, NULL},
        {"test_ull_array_to_dna", kmers_test_kmer_ull_to_dna, METH_VARARGS, NULL},
        {"test_ull_array_compare", kmers_test_kmer_ull_array_compare, METH_VARARGS, NULL},
        {NULL, NULL, 0, NULL}
};

/* Initialize the module */
static struct PyModuleDef kmers =
{
    PyModuleDef_HEAD_INIT,
    "kmers", /* name of module */
    module_docstring, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_kmers(void)
{
    PyObject *module = PyModule_Create(&kmers);

    /* Load `numpy` functionality. */
    import_array();
    return module;
}

static PyObject *kmers_getVersionString(PyObject *self, PyObject *args)
{
    PyObject *ret = Py_BuildValue("s", version_str);
    return ret;
}

static PyObject *kmers_test_InplaceNumpyMove(PyObject *self, PyObject *args)
{
    PyObject *w = NULL;
    int k;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Oi",&w, &k))
        return NULL;

    unsigned long long *w_ptr;
    PyObject *w_array;
    // Allocate W
    if (w)
    {
        /* Interpret the input objects as numpy arrays. */
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long*)PyArray_DATA(w_array); //Get pointer to data as C-type.
    }

    test_move_matrix_elements(w_ptr, k);

    // Clean up W refs
    if (w) {
        Py_XDECREF(w_array);
    }

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *kmers_test_binaryRC(PyObject *self, PyObject *args)
{
    unsigned long long my_kmer, my_rc;
    int k;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Ki",&my_kmer, &k))
        return NULL;

    my_rc = kmer32_binary_RC(my_kmer, k);

    // Build the output tuple
    PyObject *ret;
    ret = Py_BuildValue("K", my_rc);
    return ret;
}

static PyObject *kmers_test_kmer_moveright(PyObject *self, PyObject *args)
{
    unsigned long long my_kmer_p1, my_kmer_p2, my_mr;
    int k;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "KKi",&my_kmer_p1, &my_kmer_p2, &k))
        return NULL;

    my_mr = kmer_ull_array_moveright32(my_kmer_p1, my_kmer_p2, k);

    // Build the output tuple
    PyObject *ret;
    ret = Py_BuildValue("K", my_mr);
    return ret;
}

static PyObject *kmers_test_kmer_ull_to_dna(PyObject *self, PyObject *args)
{
    PyObject *w = NULL;
    int k;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Oi",&w, &k))
        return NULL;

    unsigned long long *w_ptr;
    PyObject *w_array;
    // Allocate W
    if (w)
    {
        /* Interpret the input objects as numpy arrays. */
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long*)PyArray_DATA(w_array); //Get pointer to data as C-type.
    }

    char *seq = (char *)malloc(k+1);
    seq[k]='\0';
    binary_ull_array_to_dna(w_ptr, k, seq);

    // Build the output tuple
    PyObject *ret;
    ret = Py_BuildValue("s", seq);
    return ret;
}

static PyObject *kmers_test_kmer_ull_array_compare(PyObject *self, PyObject *args)
{
    // Module function: test_ull_array_compare(w1, w2, k)
    PyObject *w1 = NULL;
    PyObject *w2 = NULL;
    int k, diffcomp, mystride;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOi", &w1, &w2, &k))
        return NULL;

    unsigned long long *w1_ptr;
    unsigned long long *w2_ptr;
    PyObject *w1_array;
    PyObject *w2_array;
    // Allocate W
    if (w1)
    {
        /* Interpret the input objects as numpy arrays. */
        w1_array = PyArray_FROM_OTF(w1, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w1_array == NULL) {
            Py_XDECREF(w1_array); // Kill and throw exception on failure.
            return NULL;
        }
        w1_ptr = (unsigned long long*)PyArray_DATA(w1_array); //Get pointer to data as C-type.
    }
    if (w2)
    {
        /* Interpret the input objects as numpy arrays. */
        w2_array = PyArray_FROM_OTF(w2, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w2_array == NULL) {
            Py_XDECREF(w2_array); // Kill and throw exception on failure.
            return NULL;
        }
        w2_ptr = (unsigned long long*)PyArray_DATA(w2_array); //Get pointer to data as C-type.
    }
    mystride = k_to_ull_stride(k);
    diffcomp = ull_array_compare(w1_ptr,w2_ptr, mystride-1);

    // Build the output tuple
    PyObject *ret;
    ret = Py_BuildValue("i", diffcomp);
    return ret;
}

//static PyObject *kmers_test_kmerull_inplace_rc(PyObject *self, PyObject *args)
//{
//    PyObject *w = NULL;
//    int k;
//    /* Parse the input tuple */
//    if (!PyArg_ParseTuple(args, "Oi",&w, &k))
//        return NULL;
//
//    unsigned long long *w_ptr;
//    PyObject *w_array;
//    // Allocate W
//    if (w)
//    {
//        /* Interpret the input objects as numpy arrays. */
//        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
//        if (w_array == NULL) {
//            Py_XDECREF(w_array); // Kill and throw exception on failure.
//            return NULL;
//        }
//        w_ptr = (unsigned long long*)PyArray_DATA(w_array); //Get pointer to data as C-type.
//    }
//
//    kmer_ull_array_replace_with_RC(w_ptr, k);
//
//    // Clean up W refs
//    if (w) {
//        Py_XDECREF(w_array);
//    }
//
//    /* Build the output tuple */
//    Py_INCREF(Py_None);
//    return Py_None;
//}

static PyObject *kmers_seq_to_kmers_ull_array(PyObject *self, PyObject *args)
{
    // MODULE FUNCTION: seq_to_kmer_ull_array(seq (str), k (int), OPTIONAL w_array (np.ndarray, dtype = np.uint64))
    // Takes a string and a k value and creates the numpy array containing all the consecutive k-mers in the
    //   sequence for the given 'k'. K-mer values are stored in binary format as ULL_arrays with the stride of the
    //   ULL depending on k. The array 'W' can either be provided in advance as a properly sized numpy array or it
    //   can be omitted and it will be allocated and sized. In either case the finished W array is returned with
    //   the function.

    PyObject *w = NULL; // numpy array to store the result in.
    int k, stride, i;
    char *seq;
    size_t l_seq;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#i|O", &seq, &l_seq, &k, &w))
        return NULL;

    // Get the W-array representing all the kmers:
    stride = k_to_ull_stride(k);
    unsigned long long *w_ptr;
    PyObject *w_array;
    if (w) {
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long *)PyArray_DATA(w_array);
        dna_sequence_to_ull_array_list_provided(seq, l_seq, k, w_ptr);
    } else {
        w_ptr = dna_sequence_to_kmer_ull_array_list(seq, l_seq, k);
        // Convert w_ptr to a Numpy array:
        npy_intp outdims[2];
        outdims[0] = l_seq - k + 1;
        outdims[1] = stride;
        w = PyArray_SimpleNewFromData(2, outdims, NPY_UINT64, (void *) w_ptr);
    }

    if (!w) {Py_XDECREF(w); return NULL;}

    /* Build the output tuple */
    PyObject *ret;
    ret = Py_BuildValue("O", w);
    return ret;
}

static PyObject *kmers_seqiter_to_kmers_ull_array(PyObject *self, PyObject *args)
{
    // MODULE FUNCTION: seqiter_to_kmer_ull_array(sequence_iterator, k, w_array, OPTIONAL: use_min_RC = False)
    // -------------------------------------------------------------------------
    // Takes a sequence iterator (i.e. a list or a dict.values()) and a k, plus a numpy
    //  uint64 array that is assumed to be large enough to record the kmers for all of
    //  the sequences in the iterator. If 'use_min_RC' is True, each k-mer value will be
    //  stored as the minimum of itself and its reverse-compliment.
    PyObject *w = NULL; // numpy array to store the result in.
    PyObject *str_iter_obj;
    int k, stride, i, use_min_RC;
    unsigned long w_pos, tot_w_len;
    char *seq;
    size_t l_seq;
    use_min_RC = 0;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "iOO|p", &k, &str_iter_obj, &w, &use_min_RC))
        return NULL;

    stride = k_to_ull_stride(k);
    // Get the provided W matrix:
    unsigned long long *w_ptr;
    PyObject *w_array;
    if (w) {
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long *) PyArray_DATA(w_array);
    }
    // Get the string iterator:
    PyObject *str_iter = PyObject_GetIter(str_iter_obj);
    PyObject *str_item;
    if (str_iter == NULL) {
        Py_XDECREF(str_iter); return NULL;
    }

    // iterate through and populate W with each one:
    w_pos = 0; i = 0; tot_w_len = 0;
    size_t str_size;
    if (DEBUG_ON) printf("k=%d, stride=%d, use-min-RC=%d\n", k, stride, use_min_RC);
    while ((str_item = PyIter_Next(str_iter))) {
        seq = PyUnicode_AsUTF8AndSize(str_item, &str_size);
        l_seq = strlen(seq);
        dna_sequence_to_ull_array_list_provided(seq, str_size, k, &w_ptr[w_pos]);
        w_pos = w_pos + (l_seq - k + 1)*stride;
        tot_w_len = tot_w_len + l_seq - k + 1;
        if (DEBUG_ON) printf("i=%d, l_seq=%lu, str_size=%lu, w_pos=%lu, tot_w_len=%lu\n", i, l_seq, str_size, w_pos, tot_w_len);
        i++;
        Py_DECREF(str_item);
    }

    // If selected to use the minimum of self and RC, run the batch correction at the end:
    if (use_min_RC) {
        printf("running final min-RC swap on results.\n");
        unsigned long n_swaps;
        n_swaps = ull_array_list_set_to_min_RC(w_ptr, k, tot_w_len);
    }

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *kmers_seqToKmer(PyObject *self, PyObject *args)
{
    // MODULE FUNCTION: seq_to_kmer_binary(seq, k, w_array).
    // Takes a string and a k value, plus a numpy uint64 array to store the result, and fills the numpy array
    // with the integers representing the k-mer. Assumes the w_array is properly sized and if it's bigger than
    // needed, only fill the first 'stride' positions of it. Likewise does not care if 'seq' is longer than 'k',
    // only uses the first K-bases.

    PyObject *w = NULL; // numpy array to store the result in.
    int k;
    char *seq;
    int l_seq;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#iO",&seq, &l_seq, &k, &w))
        return NULL;

    unsigned long long *w_ptr;

    PyObject *w_array;

    // Allocate W
    if (w)
    {
        /* Interpret the input objects as numpy arrays. */
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long*)PyArray_DATA(w_array); //Get pointer to data as C-type.
    }

    // Call the external C function.
    dna_to_binary_ull_array(seq, k, w_ptr);
    // Clean up W matrices:
    if (w) {
        Py_XDECREF(w_array);
    }

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *kmers_ull_array_batch_set_min_RC(PyObject *self, PyObject *args) {
    // MODULE FUNCTION: kmer_ull_array_batch_set_min_RC(w_array, k, w_length)
    // Takes the W array, k and the legnth of W (in k-mers) and goes through it one at
    //  a time substituting in the RC of each k-mer if the RC is the lower value. This
    //  is the canonical representation.
    PyObject *w = NULL; // numpy array to store the result in.
    int k, stride;
    unsigned long num_swaps;
    size_t array_size;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Oii", &w, &k, &array_size))
        return NULL;

    // Get the W-array representing all the kmers:
    stride = k_to_ull_stride(k);
    unsigned long long *w_ptr;
    PyObject *w_array;
    if (w) {
        w_array = PyArray_FROM_OTF(w, NPY_UINT64, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (w_array == NULL) {
            Py_XDECREF(w_array); // Kill and throw exception on failure.
            return NULL;
        }
        w_ptr = (unsigned long long *)PyArray_DATA(w_array);
    }

    num_swaps = ull_array_list_set_to_min_RC(w_ptr, k, array_size);

    if (!w) {Py_XDECREF(w); return NULL;}

    /* Build the output tuple */
    PyObject *ret;
    ret = Py_BuildValue("i", num_swaps);
    return ret;
}

static PyObject *kmers_seqiter_depthiter_to_kmers_ull(PyObject *self, PyObject *args) {

}
