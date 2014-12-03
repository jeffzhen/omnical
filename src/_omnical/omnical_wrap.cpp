#include "include/omnical_wrap.h"
#include <math.h>
#define QUOTE(a) # a
#define uint unsigned int
#define CHK_NULL(a) \
    if (a == NULL) { \
        PyErr_Format(PyExc_MemoryError, "Failed to allocate %s", QUOTE(a)); \
        return NULL; }

/*____                           _                    _
 / ___|_ __ ___  _   _ _ __   __| |_      _____  _ __| | __
| |  _| '__/ _ \| | | | '_ \ / _` \ \ /\ / / _ \| '__| |/ /
| |_| | | | (_) | |_| | | | | (_| |\ V  V / (_) | |  |   <
 \____|_|  \___/ \__,_|_| |_|\__,_| \_/\_/ \___/|_|  |_|\_\
*/
// Python object that holds instance of redundantinfo
typedef struct {
    PyObject_HEAD
    redundantinfo info;
} RedInfoObject;

// Deallocate memory when Python object is deleted
static void RedInfoObject_dealloc(RedInfoObject* self) {
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object and redundantinfo (__new__)
static PyObject *RedInfoObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    RedInfoObject *self;
    self = (RedInfoObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int RedInfoObject_init(RedInfoObject *self) {
    return 0;
}


/*___          _ ___        __                    _            _
|  _ \ ___  __| |_ _|_ __  / _| ___     __ _  ___| |_ ___  ___| |_
| |_) / _ \/ _` || || '_ \| |_ / _ \   / _` |/ _ \ __/ __|/ _ \ __|
|  _ <  __/ (_| || || | | |  _| (_) | | (_| |  __/ |_\__ \  __/ |_
|_| \_\___|\__,_|___|_| |_|_|  \___/   \__, |\___|\__|___/\___|\__|
                                       |___/                       */

// RedundantInfo.nAntenna
PyObject *RedInfoObject_get_nAntenna(RedInfoObject *self, void *closure) {
    return Py_BuildValue("l", self->info.nAntenna);
}

int RedInfoObject_set_nAntenna(RedInfoObject *self, PyObject *value, void *closure) {
    if (!PyInt_Check(value)) {
        PyErr_Format(PyExc_ValueError, "nAntenna must be an integer");
        return -1;
    }
    self->info.nAntenna = (int) PyInt_AsLong(value);
    return 0;
}

// RedundantInfo.nUBL
PyObject *RedInfoObject_get_nUBL(RedInfoObject *self, void *closure) {
    return Py_BuildValue("l", self->info.nUBL);
}

int RedInfoObject_set_nUBL(RedInfoObject *self, PyObject *value, void *closure) {
    if (!PyInt_Check(value)) {
        PyErr_Format(PyExc_ValueError, "nUBL must be an integer");
        return -1;
    }
    self->info.nUBL = (int) PyInt_AsLong(value);
    return 0;
}

// RedundantInfo.nBaseline
PyObject *RedInfoObject_get_nBaseline(RedInfoObject *self, void *closure) {
    return Py_BuildValue("l", self->info.nBaseline);
}

int RedInfoObject_set_nBaseline(RedInfoObject *self, PyObject *value, void *closure) {
    if (!PyInt_Check(value)) {
        PyErr_Format(PyExc_ValueError, "nBaseline must be an integer");
        return -1;
    }
    self->info.nBaseline = (int) PyInt_AsLong(value);
    return 0;
}

// RedundantInfo.nCross
PyObject *RedInfoObject_get_nCross(RedInfoObject *self, void *closure) {
    return Py_BuildValue("l", self->info.nCross);
}

int RedInfoObject_set_nCross(RedInfoObject *self, PyObject *value, void *closure) {
    if (!PyInt_Check(value)) {
        PyErr_Format(PyExc_ValueError, "nCross must be an integer");
        return -1;
    }
    self->info.nCross = (int) PyInt_AsLong(value);
    return 0;
}

// RedundantInfo.subsetant (1D integer array)
PyObject *RedInfoObject_get_subsetant(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.subsetant.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.subsetant[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_subsetant(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "subsetant must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "subsetant must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.subsetant.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.subsetant[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.antloc (1D integer array)
PyObject *RedInfoObject_get_antloc(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.antloc.size(), self->info.antloc[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.antloc[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_antloc(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "antloc must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "antloc must be a 2D array of floats");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.antloc.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.antloc[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.antloc[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

// RedundantInfo.subsetbl (1D integer array)
PyObject *RedInfoObject_get_subsetbl(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.subsetbl.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.subsetbl[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_subsetbl(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "subsetbl must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "subsetbl must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.subsetbl.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.subsetbl[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.ubl (1D integer array)
PyObject *RedInfoObject_get_ubl(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.ubl.size(), self->info.ubl[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.ubl[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_ubl(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "ubl must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "ubl must be a 2D array of floats");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.ubl.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.ubl[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.ubl[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

// RedundantInfo.bltoubl (1D integer array)
PyObject *RedInfoObject_get_bltoubl(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.bltoubl.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.bltoubl[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_bltoubl(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "bltoubl must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "bltoubl must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.bltoubl.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.bltoubl[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.reversed (1D integer array)
PyObject *RedInfoObject_get_reversed(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.reversed.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.reversed[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_reversed(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "reversed must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "reversed must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.reversed.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.reversed[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.reversedauto (1D integer array)
PyObject *RedInfoObject_get_reversedauto(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.reversedauto.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.reversedauto[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_reversedauto(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "reversedauto must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "reversedauto must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.reversedauto.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.reversedauto[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.autoindex (1D integer array)
PyObject *RedInfoObject_get_autoindex(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.autoindex.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.autoindex[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_autoindex(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "autoindex must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "autoindex must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.autoindex.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.autoindex[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.crossindex (1D integer array)
PyObject *RedInfoObject_get_crossindex(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.crossindex.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.crossindex[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_crossindex(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "crossindex must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "crossindex must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.crossindex.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.crossindex[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.bl2d
PyObject *RedInfoObject_get_bl2d(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.bl2d.size(), self->info.bl2d[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((int *) PyArray_GETPTR2(rv,i,j))[0] = self->info.bl2d[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_bl2d(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1,dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "bl2d must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "bl2d must be a 2D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.bl2d.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.bl2d[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.bl2d[i][j] = ((int *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

// RedundantInfo.ublcount (1D integer array)
PyObject *RedInfoObject_get_ublcount(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[1] = {self->info.ublcount.size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(1, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
        ((int *) PyArray_GETPTR1(rv,i))[0] = self->info.ublcount[i];
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_ublcount(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "ublcount must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 1 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "ublcount must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    self->info.ublcount.resize(dim1);
    for (int i=0; i < dim1; i++) {
        self->info.ublcount[i] = ((int *) PyArray_GETPTR1(v,i))[0];
    }
    return 0;
}

// RedundantInfo.ublindex
PyObject *RedInfoObject_get_ublindex(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.crossindex.size(), 5};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    int cnter = 0;
    for (unsigned int i=0; i < self->info.ublindex.size(); i++) {
      for (unsigned int j=0; j < self->info.ublindex[i].size(); j++) {
            ((int *) PyArray_GETPTR2(rv,cnter,0))[0] = i;
            ((int *) PyArray_GETPTR2(rv,cnter,1))[0] = j;
            ((int *) PyArray_GETPTR2(rv,cnter,2))[0] = self->info.ublindex[i][j][0];
            ((int *) PyArray_GETPTR2(rv,cnter,3))[0] = self->info.ublindex[i][j][1];
            ((int *) PyArray_GETPTR2(rv,cnter,4))[0] = self->info.ublindex[i][j][2];
            cnter ++;
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_ublindex(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "ublindex must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT || PyArray_DIM(v,1) != 5) {
        PyErr_Format(PyExc_ValueError, "ublindex must be a 2D array of ints consisting of (ublindex, index in this ubl, a1, a2, bl)");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    //dim2 = PyArray_DIM(v,1);//always5
    vector<int> dummy (3,0);


    for (int n=0; n < dim1; n++) {
        uint i = ((uint *) PyArray_GETPTR2(v,n,0))[0];
        uint j = ((uint *) PyArray_GETPTR2(v,n,1))[0];
        dummy[0] = ((int *) PyArray_GETPTR2(v,n,2))[0];
        dummy[1] = ((int *) PyArray_GETPTR2(v,n,3))[0];
        dummy[2] = ((int *) PyArray_GETPTR2(v,n,4))[0];
        if (self->info.ublindex.size() < i + 1){
            self->info.ublindex.resize(i + 1);
        }
        if(self->info.ublindex[i].size() < j + 1){
            self->info.ublindex[i].resize(j + 1);
        }
        self->info.ublindex[i][j] = dummy;
    }
    //self->info.ublindex.resize(dim1);
    //for (int i=0; i < dim1; i++) {
      //self->info.ublindex[i].resize(dim2);
      //for (int j=0; j < dim2; j++) {
        //self->info.ublindex[i][j].resize(dim3);
        //for (int k=0; k < dim3; k++) {
            //self->info.ublindex[i][j][k] = ((int *) PyArray_GETPTR3(v,i,j,k))[0];
        //}
      //}
    //}
    return 0;
}

// RedundantInfo.bl1dmatrix (1D integer array)
PyObject *RedInfoObject_get_bl1dmatrix(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.bl1dmatrix.size(), self->info.bl1dmatrix[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((int *) PyArray_GETPTR2(rv,i,j))[0] = self->info.bl1dmatrix[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_bl1dmatrix(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "bl1dmatrix must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT) {
        PyErr_Format(PyExc_ValueError, "bl1dmatrix must be a 1D array of ints");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.bl1dmatrix.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.bl1dmatrix[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.bl1dmatrix[i][j] = ((int *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

// RedundantInfo.degenM (1D integer array)
PyObject *RedInfoObject_get_degenM(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.degenM.size(), self->info.degenM[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.degenM[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_degenM(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "degenM must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "degenM must be a 2D array of floats");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.degenM.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.degenM[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.degenM[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

////// RedundantInfo.A (1D integer array)
////PyObject *RedInfoObject_get_A(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.A.size(), self->info.A[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((int *) PyArray_GETPTR2(rv,i,j))[0] = self->info.A[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_A(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "A must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT) {
        ////PyErr_Format(PyExc_ValueError, "A must be a 2D array of ints");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.A.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.A[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.A[i][j] = ((int *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.B (1D integer array)
////PyObject *RedInfoObject_get_B(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.B.size(), self->info.B[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((int *) PyArray_GETPTR2(rv,i,j))[0] = self->info.B[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_B(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "B must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT) {
        ////PyErr_Format(PyExc_ValueError, "B must be a 2D array of ints");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.B.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.B[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.B[i][j] = ((int *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

// RedundantInfo.Atsparse (1D integer array)
PyObject *RedInfoObject_get_Atsparse(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {3 * (self->info.crossindex.size()), 3};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    int cnter = 0;
    for (unsigned int i=0; i < self->info.Atsparse.size(); i++) {
      for (unsigned int j=0; j < self->info.Atsparse[i].size(); j++) {
            ((int *) PyArray_GETPTR2(rv,cnter,0))[0] = i;
            ((int *) PyArray_GETPTR2(rv,cnter,0))[1] = self->info.Atsparse[i][j];
            ((int *) PyArray_GETPTR2(rv,cnter,0))[2] = 1;
            cnter++;
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_Atsparse(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "Atsparse must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT || PyArray_DIM(v,1) != 3) {
        PyErr_Format(PyExc_ValueError, "Atsparse must be a 2D array of ints consisting of list of (i, j, value)");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    //dim2 = PyArray_DIM(v,1);//always 3

    for (int n=0; n < dim1; n++) {
        uint i = ((uint *) PyArray_GETPTR2(v,n,0))[0];
        uint j = ((uint *) PyArray_GETPTR2(v,n,1))[0];
        ////uint k = ((uint *) PyArray_GETPTR2(v,n,2))[0];//always 1

        if (self->info.Atsparse.size() < i + 1){
            self->info.Atsparse.resize(i + 1);
        }

        self->info.Atsparse[i].push_back(j);
    }
    return 0;
}

// RedundantInfo.Btsparse (2D integer array IO and 3D vector in C++)
PyObject *RedInfoObject_get_Btsparse(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {3 * (self->info.crossindex.size()), 3};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_INT);
    int cnter = 0;
    for (unsigned int i=0; i < self->info.Btsparse.size(); i++) {
      for (unsigned int j=0; j < self->info.Btsparse[i].size(); j++) {
            ((int *) PyArray_GETPTR2(rv,cnter,0))[0] = i;
            ((int *) PyArray_GETPTR2(rv,cnter,0))[1] = self->info.Btsparse[i][j][0];
            ((int *) PyArray_GETPTR2(rv,cnter,0))[2] = self->info.Btsparse[i][j][1];
            cnter++;
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_Btsparse(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "Btsparse must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_INT || PyArray_DIM(v,1) != 3) {
        PyErr_Format(PyExc_ValueError, "Btsparse must be a 2D array of ints consisting of list of (i, j, value).");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);

    vector<int> dummy(2,0);
    for (int n=0; n < dim1; n++) {
        uint i = ((uint *) PyArray_GETPTR2(v,n,0))[0];
        dummy[0] = ((int *) PyArray_GETPTR2(v,n,1))[0];
        dummy[1] = ((int *) PyArray_GETPTR2(v,n,2))[0];
        if (self->info.Btsparse.size() < i + 1){
            self->info.Btsparse.resize(i + 1);
        }
        self->info.Btsparse[i].push_back(dummy);
    }
    return 0;
}

// RedundantInfo.AtAi (1D integer array)
PyObject *RedInfoObject_get_AtAi(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.AtAi.size(), self->info.AtAi[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.AtAi[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_AtAi(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "AtAi must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "AtAi must be a 2D array of floats");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.AtAi.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.AtAi[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.AtAi[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

// RedundantInfo.BtBi (1D integer array)
PyObject *RedInfoObject_get_BtBi(RedInfoObject *self, void *closure) {
    PyArrayObject *rv;
    npy_intp data_dims[2] = {self->info.BtBi.size(), self->info.BtBi[0].size()};
    rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    for (int i=0; i < data_dims[0]; i++) {
      for (int j=0; j < data_dims[1]; j++) {
        ((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.BtBi[i][j];
      }
    }
    return PyArray_Return(rv);
}

int RedInfoObject_set_BtBi(RedInfoObject *self, PyObject *value, void *closure) {
    PyArrayObject *v;
    npy_intp dim1, dim2;
    if (!PyArray_Check(value)) {
        PyErr_Format(PyExc_ValueError, "BtBi must be a numpy array");
        return -1;
    }
    v = (PyArrayObject *) value;
    if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "BtBi must be a 2D array of floats");
        return -1;
    }
    dim1 = PyArray_DIM(v,0);
    dim2 = PyArray_DIM(v,1);
    self->info.BtBi.resize(dim1);
    for (int i=0; i < dim1; i++) {
      self->info.BtBi[i].resize(dim2);
      for (int j=0; j < dim2; j++) {
        self->info.BtBi[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      }
    }
    return 0;
}

////// RedundantInfo.AtAiAt (1D integer array)
////PyObject *RedInfoObject_get_AtAiAt(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.AtAiAt.size(), self->info.AtAiAt[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.AtAiAt[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_AtAiAt(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "AtAiAt must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "AtAiAt must be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.AtAiAt.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.AtAiAt[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.AtAiAt[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.BtBiBt (1D integer array)
////PyObject *RedInfoObject_get_BtBiBt(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.BtBiBt.size(), self->info.BtBiBt[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.BtBiBt[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_BtBiBt(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "BtBiBt must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "BtBiBt must be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.BtBiBt.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.BtBiBt[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.BtBiBt[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.PA (1D integer array)
////PyObject *RedInfoObject_get_PA(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.PA.size(), self->info.PA[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.PA[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_PA(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "PA must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "PA must be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.PA.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.PA[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.PA[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.PB (1D integer array)
////PyObject *RedInfoObject_get_PB(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.PB.size(), self->info.PB[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.PB[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_PB(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "PB must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "PB must be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.PB.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.PB[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.PB[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.ImPA (1D integer array)
////PyObject *RedInfoObject_get_ImPA(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.ImPA.size(), self->info.ImPA[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.ImPA[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_ImPA(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "ImPA must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "ImPA be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.ImPA.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.ImPA[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.ImPA[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

////// RedundantInfo.ImPB (1D integer array)
////PyObject *RedInfoObject_get_ImPB(RedInfoObject *self, void *closure) {
    ////PyArrayObject *rv;
    ////npy_intp data_dims[2] = {self->info.ImPB.size(), self->info.ImPB[0].size()};
    ////rv = (PyArrayObject *) PyArray_SimpleNew(2, data_dims, PyArray_FLOAT);
    ////for (int i=0; i < data_dims[0]; i++) {
      ////for (int j=0; j < data_dims[1]; j++) {
        ////((float *) PyArray_GETPTR2(rv,i,j))[0] = self->info.ImPB[i][j];
      ////}
    ////}
    ////return PyArray_Return(rv);
////}

////int RedInfoObject_set_ImPB(RedInfoObject *self, PyObject *value, void *closure) {
    ////PyArrayObject *v;
    ////npy_intp dim1, dim2;
    ////if (!PyArray_Check(value)) {
        ////PyErr_Format(PyExc_ValueError, "ImPB must be a numpy array");
        ////return -1;
    ////}
    ////v = (PyArrayObject *) value;
    ////if (PyArray_NDIM(v) != 2 || PyArray_TYPE(v) != PyArray_FLOAT) {
        ////PyErr_Format(PyExc_ValueError, "ImPB must be a 2D array of floats");
        ////return -1;
    ////}
    ////dim1 = PyArray_DIM(v,0);
    ////dim2 = PyArray_DIM(v,1);
    ////self->info.ImPB.resize(dim1);
    ////for (int i=0; i < dim1; i++) {
      ////self->info.ImPB[i].resize(dim2);
      ////for (int j=0; j < dim2; j++) {
        ////self->info.ImPB[i][j] = ((float *) PyArray_GETPTR2(v,i,j))[0];
      ////}
    ////}
    ////return 0;
////}

/*___          _ ___        __                        _   _               _
|  _ \ ___  __| |_ _|_ __  / _| ___    _ __ ___   ___| |_| |__   ___   __| |
| |_) / _ \/ _` || || '_ \| |_ / _ \  | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` |
|  _ <  __/ (_| || || | | |  _| (_) | | | | | | |  __/ |_| | | | (_) | (_| |
|_| \_\___|\__,_|___|_| |_|_|  \___/  |_| |_| |_|\___|\__|_| |_|\___/ \__,_|*/

PyObject* RedInfoObject_readredundantinfo(RedInfoObject *self, PyObject *args){
    char *filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;
    readredundantinfo(filename, &(self->info));

    Py_INCREF(Py_None);
    return Py_None;
}


/*___          _ ___        __        ___  _     _           _
|  _ \ ___  __| |_ _|_ __  / _| ___  / _ \| |__ (_) ___  ___| |_
| |_) / _ \/ _` || || '_ \| |_ / _ \| | | | '_ \| |/ _ \/ __| __|
|  _ <  __/ (_| || || | | |  _| (_) | |_| | |_) | |  __/ (__| |_
|_| \_\___|\__,_|___|_| |_|_|  \___/ \___/|_.__// |\___|\___|\__|
                                              |__/       */

static PyMethodDef RedInfoObject_methods[] = {
    {"readredundantinfo", (PyCFunction)RedInfoObject_readredundantinfo, METH_VARARGS,
        "readredundantinfo(filename)\nRead data in from specified filename."},
    {NULL}  // Sentinel
};

static PyGetSetDef RedInfoObject_getseters[] = {
    {"nAntenna", (getter)RedInfoObject_get_nAntenna, (setter)RedInfoObject_set_nAntenna, "nAntenna", NULL},
    {"nUBL", (getter)RedInfoObject_get_nUBL, (setter)RedInfoObject_set_nUBL, "nUBL", NULL},
    {"nBaseline", (getter)RedInfoObject_get_nBaseline, (setter)RedInfoObject_set_nBaseline, "nBaseline", NULL},
    {"nCross", (getter)RedInfoObject_get_nCross, (setter)RedInfoObject_set_nCross, "nCross", NULL},
    {"subsetant", (getter)RedInfoObject_get_subsetant, (setter)RedInfoObject_set_subsetant, "subsetant", NULL},
    {"antloc", (getter)RedInfoObject_get_antloc, (setter)RedInfoObject_set_antloc, "antloc", NULL},
    {"subsetbl", (getter)RedInfoObject_get_subsetbl, (setter)RedInfoObject_set_subsetbl, "subsetbl", NULL},
    {"ubl", (getter)RedInfoObject_get_ubl, (setter)RedInfoObject_set_ubl, "ubl", NULL},
    {"bltoubl", (getter)RedInfoObject_get_bltoubl, (setter)RedInfoObject_set_bltoubl, "bltoubl", NULL},
    {"reversed", (getter)RedInfoObject_get_reversed, (setter)RedInfoObject_set_reversed, "reversed", NULL},
    {"reversedauto", (getter)RedInfoObject_get_reversedauto, (setter)RedInfoObject_set_reversedauto, "reversedauto", NULL},
    {"autoindex", (getter)RedInfoObject_get_autoindex, (setter)RedInfoObject_set_autoindex, "autoindex", NULL},
    {"crossindex", (getter)RedInfoObject_get_crossindex, (setter)RedInfoObject_set_crossindex, "crossindex", NULL},
    {"bl2d", (getter)RedInfoObject_get_bl2d, (setter)RedInfoObject_set_bl2d, "bl2d", NULL},
    {"ublcount", (getter)RedInfoObject_get_ublcount, (setter)RedInfoObject_set_ublcount, "ublcount", NULL},
    {"ublindex", (getter)RedInfoObject_get_ublindex, (setter)RedInfoObject_set_ublindex, "ublindex", NULL},
    {"bl1dmatrix", (getter)RedInfoObject_get_bl1dmatrix, (setter)RedInfoObject_set_bl1dmatrix, "bl1dmatrix", NULL},
    {"degenM", (getter)RedInfoObject_get_degenM, (setter)RedInfoObject_set_degenM, "degenM", NULL},
    //{"A", (getter)RedInfoObject_get_A, (setter)RedInfoObject_set_A, "A", NULL},
    //{"B", (getter)RedInfoObject_get_B, (setter)RedInfoObject_set_B, "B", NULL},
    {"Atsparse", (getter)RedInfoObject_get_Atsparse, (setter)RedInfoObject_set_Atsparse, "Atsparse", NULL},
    {"Btsparse", (getter)RedInfoObject_get_Btsparse, (setter)RedInfoObject_set_Btsparse, "Btsparse", NULL},
    {"AtAi", (getter)RedInfoObject_get_AtAi, (setter)RedInfoObject_set_AtAi, "AtAi", NULL},
    {"BtBi", (getter)RedInfoObject_get_BtBi, (setter)RedInfoObject_set_BtBi, "BtBi", NULL},
    ////{"AtAiAt", (getter)RedInfoObject_get_AtAiAt, (setter)RedInfoObject_set_AtAiAt, "AtAiAt", NULL},
    ////{"BtBiBt", (getter)RedInfoObject_get_BtBiBt, (setter)RedInfoObject_set_BtBiBt, "BtBiBt", NULL},
    ////{"PA", (getter)RedInfoObject_get_PA, (setter)RedInfoObject_set_PA, "PA", NULL},
    ////{"PB", (getter)RedInfoObject_get_PB, (setter)RedInfoObject_set_PB, "PB", NULL},
    ////{"ImPA", (getter)RedInfoObject_get_ImPA, (setter)RedInfoObject_set_ImPA, "ImPA", NULL},
    ////{"ImPB", (getter)RedInfoObject_get_ImPB, (setter)RedInfoObject_set_ImPB, "ImPB", NULL},
    {NULL}  /* Sentinel */
};

PyTypeObject RedInfoType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "RedundantInfo", /*tp_name*/
    sizeof(RedInfoObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)RedInfoObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "This class provides a basic interface to omnical's redundantinfo.  RedundantInfo()",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    RedInfoObject_methods,     /* tp_methods */
    NULL,                    /* tp_members */
    RedInfoObject_getseters,     /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)RedInfoObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    RedInfoObject_new,       /* tp_new */
};


/*_  __           _       _                       _   _               _
|  \/  | ___   __| |_   _| | ___   _ __ ___   ___| |_| |__   ___   __| |___
| |\/| |/ _ \ / _` | | | | |/ _ \ | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
| |  | | (_) | (_| | |_| | |  __/ | | | | | |  __/ |_| | | | (_) | (_| \__ \
|_|  |_|\___/ \__,_|\__,_|_|\___| |_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/ */

PyObject *redcal_wrap(PyObject *self, PyObject *args, PyObject *kwds) {
    int uselogcal = 1, removedegen = 1, maxiter = 10, dummy = 0, computeUBLFit = 1, trust_period = 1;
    float stepsize=.3, conv=.01;
    npy_intp dims[3] = {0, 0, 0}; // time, fq, bl
    npy_intp nint, nfreq, nbls;
    RedInfoObject *redinfo;
    PyArrayObject *data, *additivein, *calpar; // input arrays
    PyArrayObject *additiveout=NULL; // output arrays
    PyObject *rv;
    static char *kwlist[] = {"data", "calpar", "info", "additivein", "uselogcal", "removedegen", "maxiter", "stepsize", "conv", "computeUBLFit", "trust_period", "dummy"};
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"O!O!O!O!|iiiffiii", kwlist,
            &PyArray_Type, &data, &PyArray_Type, &calpar, &RedInfoType, &redinfo, &PyArray_Type, &additivein,
            &uselogcal, &removedegen, &maxiter, &stepsize, &conv, &computeUBLFit, &trust_period, &dummy))
        return NULL;
    // check shape and type of data
    if (PyArray_NDIM(data) != 3 || PyArray_TYPE(data) != PyArray_CFLOAT) {
        PyErr_Format(PyExc_ValueError, "data must be a (nint,nfreq,nbls) array of complex floats");
        return NULL;
    }
    nint = PyArray_DIM(data,0);
    nfreq = PyArray_DIM(data,1);
    nbls = PyArray_DIM(data,2);
    vector<vector<float> > data_v(nbls, vector<float>(2, 0));
    vector<float> calpar_v(3 + 2*(redinfo->info.nUBL + redinfo->info.nAntenna), 0);
    vector<vector<float> >additivein_v(nbls, vector<float>(2, 0));
    vector<vector<float> >additiveout_v(nbls, vector<float>(2, 0));
    // check that dims of additivein and data match
    if (PyArray_NDIM(additivein) != 3 || PyArray_TYPE(additivein) != PyArray_CFLOAT
            || PyArray_DIM(additivein,0) != nint || PyArray_DIM(additivein,1) != nfreq || PyArray_DIM(additivein,2) != nbls) {
        PyErr_Format(PyExc_ValueError, "additivein must be of the same type and shape as data");
        return NULL;
    }
    if (PyArray_NDIM(calpar) != 3 || PyArray_TYPE(calpar) != PyArray_FLOAT
            || PyArray_DIM(calpar,0) != nint || PyArray_DIM(calpar,1) != nfreq || (uint)PyArray_DIM(calpar,2) != calpar_v.size()) {
        PyErr_Format(PyExc_ValueError, "calpar is expected to be a 3D numpy array of float32 with the first 2 dimensions identical to those of data");
        return NULL;
    }

    // allocate output additiveout array
    dims[0] = nint;
    dims[1] = nfreq;
    dims[2] = nbls;
    additiveout = (PyArrayObject *) PyArray_SimpleNew(3, dims, PyArray_CFLOAT);
    CHK_NULL(additiveout);
    calmemmodule module;////memory module to be reused in order to avoid redeclaring all sorts of long vectors
    initcalmodule(&module, &(redinfo->info));

    for (int t = 0; t < nint; t++){
        for (int f = 0; f < nfreq; f++){
            // copy from input arrays
            for (int b = 0; b < nbls; b++) {
                data_v[b][0] = ((float *) PyArray_GETPTR3(data,t,f,b))[0];
                data_v[b][1] = ((float *) PyArray_GETPTR3(data,t,f,b))[1];
                additivein_v[b][0] = ((float *) PyArray_GETPTR3(additivein,t,f,b))[0];
                additivein_v[b][1] = ((float *) PyArray_GETPTR3(additivein,t,f,b))[1];
            }

            if (uselogcal) {
                logcaladd(
                    &data_v, //(vector<vector<float> > *) PyArray_GETPTR3(data,t,f,0),
                    &additivein_v, //(vector<vector<float> > *) PyArray_GETPTR3(additivein,t,f,0),
                    &(redinfo->info),
                    &calpar_v, //(vector<float> *) PyArray_GETPTR3(calpar,t,f,0),
                    &additiveout_v, //(vector<vector<float> > *) PyArray_GETPTR3(additiveout,t,f,0),
                    1,
                    &module
                );
                //lincal(
                    //&data_v, //(vector<vector<float> > *) PyArray_GETPTR3(data,t,f,0),
                    //&additivein_v, //(vector<vector<float> > *) PyArray_GETPTR3(additivein,t,f,0),
                    //&(redinfo->info),
                    //&calpar_v, //(vector<float> *) PyArray_GETPTR3(calpar,t,f,0),
                    //&additiveout_v, //(vector<vector<float> > *) PyArray_GETPTR3(additiveout,t,f,0),
                    //0,
                    //&module,
                    //conv,
                    //maxiter,
                    //stepsize
                //);
            } else {
                if (t % trust_period == 0){//whether to start from logcal calpar or the result of revious lincal result
                    for (unsigned int n = 0; n < calpar_v.size(); n ++){
                        calpar_v[n] = ((float *) PyArray_GETPTR2(calpar, t, f))[n];
                    }
                } else {
                    for (unsigned int n = 0; n < calpar_v.size(); n ++){
                        calpar_v[n] = ((float *) PyArray_GETPTR2(calpar, t - 1, f))[n];
                    }
                }
                lincal(
                    &data_v, //(vector<vector<float> > *) PyArray_GETPTR3(data,t,f,0),
                    &additivein_v, //(vector<vector<float> > *) PyArray_GETPTR3(additivein,t,f,0),
                    &(redinfo->info),
                    &calpar_v, //(vector<float> *) PyArray_GETPTR3(calpar,t,f,0),
                    &additiveout_v, //(vector<vector<float> > *) PyArray_GETPTR3(additiveout,t,f,0),
                    computeUBLFit,
                    &module,
                    conv,
                    maxiter,
                    stepsize
                );
            }
            //if (removedegen) removeDegen((vector<float> *) PyArray_GETPTR3(calpar,t,f,0), &(redinfo->info), &module);
            if (removedegen) removeDegen(&calpar_v, &(redinfo->info), &module);
            // copy to output arrays
            for (int b = 0; b < nbls; b++) {
                ((float *) PyArray_GETPTR3(additiveout,t,f,b))[0] = additiveout_v[b][0];
                ((float *) PyArray_GETPTR3(additiveout,t,f,b))[1] = additiveout_v[b][1];
            }
            for (unsigned int b = 0; b < calpar_v.size(); b++) {
                ((float *) PyArray_GETPTR3(calpar,t,f,b))[0] = calpar_v[b];
            }
        }
    }

    // Make sure to DECREF when using Py_BuildValue() with objects!!
    ////rv = Py_BuildValue("(OO)", PyArray_Return(calpar), PyArray_Return(additiveout));
    rv = Py_BuildValue("O", PyArray_Return(additiveout));
    Py_DECREF(additiveout);
    return rv;
}

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
    float *capturedata;
    //PyObject      *out_array;
    //NpyIter *in_iter;
    //NpyIter *out_iter;
    //NpyIter_IterNextFunc *in_iternext;
    //NpyIter_IterNextFunc *out_iternext;

    /*  parse single numpy array argument */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &in_array))
        return NULL;
    Py_INCREF(in_array);
    capturedata = ((float *)PyArray_BYTES(in_array));
    vector<int> v(capturedata, capturedata + sizeof capturedata / sizeof capturedata[0]);
    return Py_BuildValue("ffffff", capturedata[0], capturedata[1], capturedata[2], capturedata[3], capturedata[4], capturedata[5]);
}


PyObject *unwrap_phase(PyObject *self, PyObject *args, PyObject *kwds) {
    npy_intp dims[1] = {0}; // fq
    npy_intp nfreq;
    int dummy = 0;
    PyArrayObject *data; // input arrays
    PyArrayObject *out=NULL; // output arrays

    PyObject *rv;
    static char *kwlist[] = {"data", "dummy"};
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"O!|i", kwlist,
            &PyArray_Type, &data, &dummy))
        return NULL;
    // check shape and type of data
    if (PyArray_NDIM(data) != 1 || PyArray_TYPE(data) != PyArray_FLOAT) {
        PyErr_Format(PyExc_ValueError, "data must be a 1D array of floats");
        return NULL;
    }
    nfreq = PyArray_DIM(data,0);
    dims[0] = nfreq;
    out = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_FLOAT);
    CHK_NULL(out);

    ((float *) PyArray_GETPTR1(out,0))[0] = ((float *) PyArray_GETPTR1(data,0))[0];
    float p = 3.141592653589793238462643383279502884197169399;
    float tp = 2 * p;

    for (int f = 1; f < nfreq; f++){
        ((float *) PyArray_GETPTR1(out,f))[0] = ((float *) PyArray_GETPTR1(data,f))[0];
        while (((float *) PyArray_GETPTR1(out,f))[0] < ((float *) PyArray_GETPTR1(out,f - 1))[0] - p){
            ((float *) PyArray_GETPTR1(out,f))[0] += tp;
        }
        while (((float *) PyArray_GETPTR1(out,f))[0] > ((float *) PyArray_GETPTR1(out,f - 1))[0] + p){
            ((float *) PyArray_GETPTR1(out,f))[0] -= tp;
        }
    }
    rv = Py_BuildValue("O", PyArray_Return(out));
    Py_DECREF(out);
    return rv;
}

// Module methods
static PyMethodDef omnical_methods[] = {
    {"phase", (PyCFunction)phase_wrap, METH_VARARGS,
        "Return the phase of a + bi."},
    {"norm", (PyCFunction)norm_wrap, METH_VARARGS,
        "Return the norm of input array."},
    {"redcal", (PyCFunction)redcal_wrap, METH_VARARGS | METH_KEYWORDS,
        "redcal(data,calpar,info,additivein,uselogcal=1,removedegen=0,maxiter=20,stepsize=.3,conv=.001)\nRun redundant calibration on data (3D array of complex floats)."},
    {"unwrap_phase", (PyCFunction)unwrap_phase, METH_VARARGS | METH_KEYWORDS,
        "unwrap_phase(calpar)\nUnwrap phase in radians."},
    ////{"omnical_old", (PyCFunction)cal_wrap_old, METH_VARARGS,
        ////"omnical outdated version that relies on hard disk I/O."},
    ////{"omnical", (PyCFunction)cal_wrap, METH_VARARGS,
        ////"omnical."},
    {NULL, NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC init_omnical(void) {
    PyObject* m;
    RedInfoType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&RedInfoType) < 0) return;
    m = Py_InitModule3("_omnical", omnical_methods,
    "Wrapper for Omnical redundant calibration code.");
    Py_INCREF(&RedInfoType);
    PyModule_AddObject(m, "RedundantInfo", (PyObject *)&RedInfoType);
    import_array();
}
