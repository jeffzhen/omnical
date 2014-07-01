#include "include/omnical_wrap.h"

PyObject * test(PyObject *self, PyObject *args) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) return NULL;
    return Py_BuildValue("i", i+j);
}

// Module methods
static PyMethodDef omnical_methods[] = {
    {"test", (PyCFunction)test, METH_VARARGS,
        "test(i,j)\nReturn i+j"},
    {NULL, NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC init_omnical(void) {
    PyObject* m;
    m = Py_InitModule3("_omnical", omnical_methods,
    "Wrapper for Omnical redundant calibration code.");
}
