#include <Python.h>
#include <calibration_omni.h>

static PyObject* phase(PyObject *self, PyObject *args){
    const char *command;
    int a, b;

    if (!PyArg_ParseTuple(args, "ff", &command, &a, &b))
        return NULL;
    int result = phase(a, b);
    return Py_BuildValue("f", result);
}

static PyMethodDef calibration_omni_extension_methods[] = {
        {"phase",         phase,      METH_VARARGS,
         "Return the phase of a + bi."},
        {NULL,          NULL}           /* sentinel */
};

void
initcalibration_omni_extension(void)
{
        PyImport_AddModule("calibration_omni_extension");
        Py_InitModule("calibration_omni_extension", calibration_omni_extension_methods);
}
