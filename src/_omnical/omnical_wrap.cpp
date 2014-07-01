#include "include/omnical_wrap.h"

PyObject* phase_wrap(PyObject *self, PyObject *args){
    float a, b;

    if (!PyArg_ParseTuple(args, "ff", &a, &b))
        return NULL;
    float result = phase(a, b);
    //cout << a << endl; cout.flush();
    return Py_BuildValue("f", result);
}

PyObject* norm_wrap(PyObject *self, PyObject *args){
    PyArrayObject *in_array;
    float capturedata;
    //PyObject      *out_array;
    //NpyIter *in_iter;
    //NpyIter *out_iter;
    //NpyIter_IterNextFunc *in_iternext;
    //NpyIter_IterNextFunc *out_iternext;

    /*  parse single numpy array argument */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &in_array))
        return NULL;
    capturedata = ((float *)PyArray_BYTES(in_array))[0];

    return Py_BuildValue("f", capturedata);
}

// Module methods
static PyMethodDef omnical_methods[] = {
    {"phase", (PyCFunction)phase_wrap, METH_VARARGS,
        "Return the phase of a + bi."},
    {"norm", (PyCFunction)norm_wrap, METH_VARARGS,
        "Return the norm of input array."},
    {NULL, NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC init_omnical(void) {
    PyObject* m;
    m = Py_InitModule3("_omnical", omnical_methods,
    "Wrapper for Omnical redundant calibration code.");
    import_array();
}
