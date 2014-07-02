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

PyObject* cal_wrap(PyObject *self, PyObject *args){
    string FILENAME = "omnical_wrap.cpp";
    string METHODNAME = "cal_wrap";

    clock_t tStart = clock();

    char* visin_char;
    char* infopath_char;


    int nInt, nFreq, nAnt;
    bool removedegen, removeadd;
    int additivePeriod;
    char* additivePeriodstr_char;
    bool use_logcal;
    float converge_percent;
    int max_iter, step_size;

    if (!PyArg_ParseTuple(args, "ssiiiiisifii", &visin_char, &infopath_char, &nInt, &nFreq, &nAnt, &removedegen, &removeadd, &additivePeriodstr_char, &use_logcal, &converge_percent, &max_iter, &step_size))
        return NULL;

    string visin(visin_char);
    string infopath(infopath_char);
    additivePeriod = atoi(additivePeriodstr_char);
    string additivePeriodstr(additivePeriodstr_char);

    string calparout;
    if (removeadd){
        calparout = visin + "_add" + additivePeriodstr + ".omnical";
    }else{
        calparout = visin + ".omnical";
    }

    cout << "##" << FILENAME << "##" << METHODNAME << ": Starting " << visin << " " << nInt << " by " << nFreq << endl;
    cout << "##" << FILENAME << "##" << METHODNAME << ": Reading redundant baseline information and pre-computed matrices:" << endl;//generated from 16. additive noise investigation _from_17.nb
    redundantinfo info;
    readredundantinfo(infopath, &info);

    cout << "Good antenna count: " << info.nAntenna << ". UBL count: " << info.nUBL << "." << endl;
    cout.flush();


    int nBaseline = nAnt * (nAnt + 1) / 2;
    int nCross = nAnt * (nAnt - 1) / 2;

    clock_t tStartVar = clock();
    ////allocate big memories for calibration operations
    cout << "##" << FILENAME << "##" << METHODNAME << ": Allocating big memories for calibration operations...";
    cout.flush();
    vector<vector<vector<vector<float> > > > rawdata(1, vector<vector<vector<float> > >(nInt, vector<vector<float> >(nFreq, vector<float>(2 * nBaseline, 0))));
    vector<vector<vector<vector<float> > > > data(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    vector<vector<vector<float> > > calpar(nInt, vector<vector<float> >(nFreq, vector<float>(3 + 2*(info.nUBL + info.nAntenna), 0)));
    vector<vector<float> >additiveplaceholder(data[0][0].size(), vector<float>(data[0][0][0].size(), 0));
    vector<vector<float> >additiveplaceholder2(data[0][0].size(), vector<float>(data[0][0][0].size(), 0));
    //vector<vector<vector<vector<float> > > > additivein(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    vector<vector<vector<vector<float> > > > additiveout(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    calmemmodule module;////memory module to be reused in order to avoid redeclaring all sorts of long vectors
    initcalmodule(&module, &info);
    cout << "Done." << endl;
    cout.flush();




    ////////////Start calibration///////////

    readBinaryVisibilityLarge((visin).c_str(), &rawdata, 1, nInt, nFreq, nBaseline);

    cout << "##" << FILENAME << "##" << METHODNAME << ": Loading good visibilities...";
    cout.flush();
    loadGoodVisibilities(&rawdata, &data, &info, 0);
    cout << "Done." << endl;
    cout.flush();
    //printvv(&(data[5][50]));
    //return 0;
    //logcaladd(&(data[5][50]), &(additiveplaceholder), &info, &(calpar[5][50]), &(additiveplaceholder), 1, &module);
    //lincal(&(data[5][50]), &(additiveplaceholder2), &info, &(calpar[5][50]), &module, 0.01, 10, 0.3);
    //printv(&(calpar[5][50]), 0,10);
    //return 0;
    cout << "##" << FILENAME << "##" << METHODNAME << ": Calibrating...";
    cout.flush();
    for (int t = 0; t < data.size(); t++){
        for (int f = 0; f < data[0].size(); f++){
            if(use_logcal){
                logcaladd(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveplaceholder2), 1, &module);
                lincal(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveout[t][f]), 0, &module, converge_percent, max_iter, step_size);
            } else{
                lincal(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveout[t][f]), 1, &module, converge_percent, max_iter, step_size);
            }
            //if (f==50) cout << calpar[t][f][0] << " " << calpar[t][f][1] << " " << calpar[t][f][2] << endl;
            if((not removeadd) and removedegen) removeDegen(&(calpar[t][f]), &info, &module);
        }
    }
    if(removeadd){
        runAverage(&additiveout, 0, additivePeriod);
        for (int t = 0; t < data.size(); t++){
            for (int f = 0; f < data[0].size(); f++){
                lincal(&(data[t][f]), &(additiveout[t][f]), &info, &(calpar[t][f]), &(additiveplaceholder2), 0, &module, converge_percent, max_iter, step_size);
                if(removedegen) removeDegen(&(calpar[t][f]), &info, &module);
            }
        }
    }
    cout << "Done." << endl;
    cout.flush();


    cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing results...";
    cout.flush();
    if(removeadd){
        stringstream ss;
        ss << additivePeriod;
        //string str = ss.str();
        outputDataLarge(&additiveout, (visin + ".omniadd" + ss.str()).c_str());
    }
    outputCalparSP(&calpar, calparout, false, info.nAntenna);
    cout << "##" << FILENAME << "##" << METHODNAME << ": Done. ";
    cout.flush();
    printf("Calibration time taken: %.2fs; ", (double)(clock() - tStartVar)/CLOCKS_PER_SEC);
    printf("Total time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* cal_wrap_old(PyObject *self, PyObject *args){
    string FILENAME = "omnical_wrap.cpp";
    string METHODNAME = "cal_wrap";

    clock_t tStart = clock();

    char* visin_char;
    char* infopath_char;


    int nInt, nFreq, nAnt;
    bool removedegen, removeadd;
    int additivePeriod;
    char* additivePeriodstr_char;
    bool use_logcal;
    float converge_percent;
    int max_iter, step_size;

    if (!PyArg_ParseTuple(args, "ssiiiiisifii", &visin_char, &infopath_char, &nInt, &nFreq, &nAnt, &removedegen, &removeadd, &additivePeriodstr_char, &use_logcal, &converge_percent, &max_iter, &step_size))
        return NULL;

    string visin(visin_char);
    string infopath(infopath_char);
    additivePeriod = atoi(additivePeriodstr_char);
    string additivePeriodstr(additivePeriodstr_char);

    string calparout;
    if (removeadd){
        calparout = visin + "_add" + additivePeriodstr + ".omnical";
    }else{
        calparout = visin + ".omnical";
    }

    cout << "##" << FILENAME << "##" << METHODNAME << ": Starting " << visin << " " << nInt << " by " << nFreq << endl;
    cout << "##" << FILENAME << "##" << METHODNAME << ": Reading redundant baseline information and pre-computed matrices:" << endl;//generated from 16. additive noise investigation _from_17.nb
    redundantinfo info;
    readredundantinfo(infopath, &info);

    cout << "Good antenna count: " << info.nAntenna << ". UBL count: " << info.nUBL << "." << endl;
    cout.flush();


    int nBaseline = nAnt * (nAnt + 1) / 2;
    int nCross = nAnt * (nAnt - 1) / 2;

    clock_t tStartVar = clock();
    ////allocate big memories for calibration operations
    cout << "##" << FILENAME << "##" << METHODNAME << ": Allocating big memories for calibration operations...";
    cout.flush();
    vector<vector<vector<vector<float> > > > rawdata(1, vector<vector<vector<float> > >(nInt, vector<vector<float> >(nFreq, vector<float>(2 * nBaseline, 0))));
    vector<vector<vector<vector<float> > > > data(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    vector<vector<vector<float> > > calpar(nInt, vector<vector<float> >(nFreq, vector<float>(3 + 2*(info.nUBL + info.nAntenna), 0)));
    vector<vector<float> >additiveplaceholder(data[0][0].size(), vector<float>(data[0][0][0].size(), 0));
    vector<vector<float> >additiveplaceholder2(data[0][0].size(), vector<float>(data[0][0][0].size(), 0));
    //vector<vector<vector<vector<float> > > > additivein(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    vector<vector<vector<vector<float> > > > additiveout(nInt, vector<vector<vector<float> > >(nFreq, vector<vector<float> >(info.subsetbl.size(), vector<float>(2, 0))));
    calmemmodule module;////memory module to be reused in order to avoid redeclaring all sorts of long vectors
    initcalmodule(&module, &info);
    cout << "Done." << endl;
    cout.flush();




    ////////////Start calibration///////////

    readBinaryVisibilityLarge((visin).c_str(), &rawdata, 1, nInt, nFreq, nBaseline);

    cout << "##" << FILENAME << "##" << METHODNAME << ": Loading good visibilities...";
    cout.flush();
    loadGoodVisibilities(&rawdata, &data, &info, 0);
    cout << "Done." << endl;
    cout.flush();
    //printvv(&(data[5][50]));
    //return 0;
    //logcaladd(&(data[5][50]), &(additiveplaceholder), &info, &(calpar[5][50]), &(additiveplaceholder), 1, &module);
    //lincal(&(data[5][50]), &(additiveplaceholder2), &info, &(calpar[5][50]), &module, 0.01, 10, 0.3);
    //printv(&(calpar[5][50]), 0,10);
    //return 0;
    cout << "##" << FILENAME << "##" << METHODNAME << ": Calibrating...";
    cout.flush();
    for (int t = 0; t < data.size(); t++){
        for (int f = 0; f < data[0].size(); f++){
            if(use_logcal){
                logcaladd(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveplaceholder2), 1, &module);
                lincal(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveout[t][f]), 0, &module, converge_percent, max_iter, step_size);
            } else{
                lincal(&(data[t][f]), &(additiveplaceholder), &info, &(calpar[t][f]), &(additiveout[t][f]), 1, &module, converge_percent, max_iter, step_size);
            }
            //if (f==50) cout << calpar[t][f][0] << " " << calpar[t][f][1] << " " << calpar[t][f][2] << endl;
            if((not removeadd) and removedegen) removeDegen(&(calpar[t][f]), &info, &module);
        }
    }
    if(removeadd){
        runAverage(&additiveout, 0, additivePeriod);
        for (int t = 0; t < data.size(); t++){
            for (int f = 0; f < data[0].size(); f++){
                lincal(&(data[t][f]), &(additiveout[t][f]), &info, &(calpar[t][f]), &(additiveplaceholder2), 0, &module, converge_percent, max_iter, step_size);
                if(removedegen) removeDegen(&(calpar[t][f]), &info, &module);
            }
        }
    }
    cout << "Done." << endl;
    cout.flush();


    cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing results...";
    cout.flush();
    if(removeadd){
        stringstream ss;
        ss << additivePeriod;
        //string str = ss.str();
        outputDataLarge(&additiveout, (visin + ".omniadd" + ss.str()).c_str());
    }
    outputCalparSP(&calpar, calparout, false, info.nAntenna);
    cout << "##" << FILENAME << "##" << METHODNAME << ": Done. ";
    cout.flush();
    printf("Calibration time taken: %.2fs; ", (double)(clock() - tStartVar)/CLOCKS_PER_SEC);
    printf("Total time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    Py_INCREF(Py_None);
    return Py_None;
}

// Module methods
static PyMethodDef omnical_methods[] = {
    {"phase", (PyCFunction)phase_wrap, METH_VARARGS,
        "Return the phase of a + bi."},
    {"norm", (PyCFunction)norm_wrap, METH_VARARGS,
        "Return the norm of input array."},
    {"omnical_old", (PyCFunction)cal_wrap_old, METH_VARARGS,
        "omnical outdated version that relies on hard disk I/O."},
    {"omnical", (PyCFunction)cal_wrap, METH_VARARGS,
        "omnical outdated version that relies on hard disk I/O."},
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
