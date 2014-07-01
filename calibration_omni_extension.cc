#include <Python.h>
#include <calibration_omni.h>

static PyObject*
phase(PyObject *self, PyObject *args){
    const char *command;
    float a, b;

    if (!PyArg_ParseTuple(args, "ff", &a, &b))
        return NULL;
    float result = phase(a, b);
    cout << a << endl; cout.flush();
    return Py_BuildValue("f", a);
}

static PyMethodDef
calibration_omni_extension_methods[] = {
        {"phase",         phase,      METH_VARARGS,
         "Return the phase of a + bi."},
        {NULL,          NULL}           /* sentinel */
};

extern "C" void
initcalibration_omni_extension(void)
{
        PyImport_AddModule("calibration_omni_extension");
        Py_InitModule("calibration_omni_extension", calibration_omni_extension_methods);
}
