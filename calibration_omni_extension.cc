#include <Python.h>
#include <calibration_omni.h>
#include <numpy/arrayobject.h>
#include <math.h>
PyObject*
PHASE(PyObject *self, PyObject *args){
    float a, b;

    if (!PyArg_ParseTuple(args, "ff", &a, &b))
        return NULL;
    float result = phase(a, b);
    //cout << a << endl; cout.flush();
    return Py_BuildValue("f", result);
}

PyObject*
NORM(PyObject *self, PyObject *args){
    PyArrayObject *in_array;
    //PyObject      *out_array;
    //NpyIter *in_iter;
    //NpyIter *out_iter;
    //NpyIter_IterNextFunc *in_iternext;
    //NpyIter_IterNextFunc *out_iternext;

    /*  parse single numpy array argument */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &in_array))
        return NULL;
    float * capturedata = (float *)PyArray_BYTES(in_array);

    return Py_BuildValue("f", *capturedata);
}

static PyMethodDef
calibration_omni_extension_methods[] = {
        {"phase",         PHASE,      METH_VARARGS,
         "Return the phase of a + bi."},
         {"norm",         NORM,      METH_VARARGS,
         "Return the norm of input array."},
        {NULL,          NULL}           /* sentinel */
};

extern "C" void
initcalibration_omni_extension(void)
{
        PyImport_AddModule("calibration_omni_extension");
        Py_InitModule("calibration_omni_extension", calibration_omni_extension_methods);
}
