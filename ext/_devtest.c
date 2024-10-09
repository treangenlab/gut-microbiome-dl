//
// Created by miken on 3/22/2024.
//
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "devtest.h"
#include "myutils.h"

#define DEBUG_ON 1

/* Docstrings */
//#include "_nwops_docstr.h"

/* Available functions */
static PyObject *devtest_getVersionString(PyObject *self, PyObject *args);
static PyObject *devtest_testFunction(PyObject *self, PyObject *args);
//static PyObject *nwops_nwalign_tuple(PyObject *self, PyObject *args);
//static PyObject *nwops_nwalign_score(PyObject *self, PyObject *args);
//static PyObject *nwops_swalign_score_custom_costs(PyObject *self, PyObject *args);

/* Docstrings */
static char version_str[] = "Version 0.1 (beta)";
static char get_version_docstring[] =
        "Just prints the version number of this module. For easy testing that import is ok.                  \n";
static char module_docstring[] =
        "Misc testing functions for python C extension modules.                                              \n";

/* Module specification */
static PyMethodDef module_methods[] = {
        {"get_version", devtest_getVersionString, METH_VARARGS, get_version_docstring},
        {"test_function", devtest_testFunction, METH_VARARGS, NULL},
        {NULL, NULL, 0, NULL}
};

/* Initialize the module */
static struct PyModuleDef devtest =
{
    PyModuleDef_HEAD_INIT,
    "devtest", /* name of module */
    module_docstring, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_devtest(void)
{
    PyObject *module = PyModule_Create(&devtest);

    /* Load `numpy` functionality. */
    import_array();
    return module;
}

static PyObject *devtest_getVersionString(PyObject *self, PyObject *args)
{
    PyObject *ret = Py_BuildValue("s", version_str);
    return ret;
}

static PyObject *devtest_testFunction(PyObject *self, PyObject *args)
{
    int k, do_right;
    unsigned long long mask;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "ii",&k,&do_right))
        return NULL;

    if (do_right) {
        mask = rightmask(k);
    } else {
        mask = leftmask(k);
    }


    PyObject *ret = Py_BuildValue("K", mask);
    return ret;
}
