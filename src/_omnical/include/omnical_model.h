#ifndef OMNICAL_MODEL_H
#define OMNICAL_MODEL_H

#include <stdint.h>
#include <stdio.h>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include <numeric>
#define uint unsigned int
using namespace std;

///////////////////////////////////////
//command line and Python interaction//
///////////////////////////////////////

string ftostr(float f);//float to string
string itostr(int i, int len);//int to string of specified length len
vector<float> strtovf(string in);


///////////////////////////////////////
//Helper methods///////////////////////
///////////////////////////////////////
vector<float> tp2xyz (vector<float> thephi);
vector<float> tp2xyz (float t, float p);
vector<float> xyz2tp (vector<float> xyz);
vector<float> xyz2tp (float x, float y, float z);

vector<float> tp2rd (vector<float> thephi);
vector<float> tp2rd (float t, float p);
vector<float> rd2tp (vector<float> rd);
vector<float> rd2tp (float r, float d);
vector<float> tp2aa (vector<float> thephi);//alt-az
vector<float> tp2aa (float t, float p);//alt-az
vector<float> aa2tp (vector<float> aa);//alt-az
vector<float> aa2tp (float alt, float az);



void matrixDotV(vector<vector<float> > * A, vector<float> * x, vector<float> * b);

////void iqDemod(vector<vector<vector<vector<vector<float> > > > > *data, vector<vector<vector<vector<vector<float> > > > > *data_out, int nIntegrations, int nFrequencies, int nAnt);

////void iqDemodLarge(vector<vector<vector<vector<float> > > > *data, vector<vector<vector<vector<float> > > > *data_out, int nIntegrations, int nFrequencies, int nAnt);//For large set of data (more than ~200 time slices), it takes more than 24G of memory to allocate a large enough 5D vector. So we consolidate the last two dimensions (baseline and re/im) together

int gc(int a, int b);//Get coordinate for large vectors when the last two coordinates are consolidated into one dimension. 2a+b

/*void iqDemodBinary(string data, string data_out, int nIntegrations, int nFrequencies, int nAnt);//Used for large odfs for which there's not enough memory to take in all the data, so we have to read and iq and write without saving too much into memory*/

float square(float x);

int get1DBL(int i, int j, int nAntenna);

vector<int> get2DBL(int bl, int nAntenna);//bl counts cross corrs AND auto corrs

vector<int> get2DBLCross(int bl, int nAntenna);//bl only counts cross corrs

bool contains(vector<vector<float> > * UBL, vector<float> bl);//automatically checks for the opposite direction

int indexUBL(vector<vector<float> > * UBL, vector<float> bl);//give the index of a baseline inside the unique baseline list; automatically checks for the opposite direction, the opposite direction will give -index

bool contains_int(vector<int> * list, int j);

float amp(vector<float> * x);

float amp(float x, float y);

float phase(float re, float im);

float phase(vector<float> * c);

float norm(vector<vector<float> > * v);

vector<float> conjugate (vector<float> x);

void addPhase(vector<float> * x, float phi);

float vectorDot(vector<float>* v1, vector<float>* v2);//v1.v2

vector<float> matrixDotV(vector<vector<float> >* m, vector<float>* v);//m.v

vector<vector<float> > rotationMatrix(float x, float y, float z);//approximation for a rotation matrix rotating around x,y,z axis, {{1, -z, y}, {z, 1, -x}, {-y, x, 1}}

vector<vector<float> > rotationMatrixZ(float z);

///////////////PHASE CALIBRATE STUFF///////////////////
/////////////////////////////////////////////

vector<float> getBL(int i, int j, vector<vector<float> > *antloc);

int countUBL(vector<vector<float> > *antloc );

int lookupAnt(float x, float y, vector<vector<float> > antloc);

float phaseWrap (float x, float offset = -atan2(0,-1));//Wrap phase to be on (offset, offset+2pi]



///////////////CHI SQUARE STUFF///////////////////
/////////////////////////////////////////////
void computeUBL(vector<vector<float> > * antloc, vector<vector<float> > * listUBL);


vector<float> modelToMeasurement(vector<float> *modelCor, float ampcal, float phasecal);

void computeUBLcor(vector<vector<float> >* calibratedData, vector<vector<float> > *listUBL, vector<vector<float> > *UBLcor, vector<bool> *goodAnt);//average each group of calibrated redundant baselines to get the estimate for that ubl, only useful when directly applying a set of calpars instead of using logcal to fit for them.


vector<float> getModel(int i, int j, vector<vector<float> > *antloc, vector<vector<float> > *listUBL, vector<vector<float> > *UBLcor);

vector<vector<float> > ReverseEngineer(vector<float> * ampcalpar, vector<float> * phasecalpar, vector<vector<float> > * UBLcor, vector<vector<float> > * antloc, vector<vector<float> > * listUBL);

void ReverseEngineer(vector<vector<float> >* output, vector<float> * calpar, int numAntenna, vector<vector<float> > * UBLindex);

float chiSq(vector<vector<float> > * dataf, vector<vector<float> > * sdevf, vector<vector<float> > * antloc, vector<float> * ampcalpar, vector<float> * phasecalpar, vector<vector<float> > * UBLcor, int numAntenna, vector<vector<float> > * listUBL);

bool fillChiSq(vector<vector<float> >* dataf, vector<vector<float> >* sdevf, vector<float>* calpar, int numAntenna, vector<int>* UBLindex, vector<bool>* goodAnt);

///////////////POINT SOURCE STUFF///////////////////
/////////////////////////////////////////////
float median (vector<float> list); //Not using pointer because the process involves sorting which will modify the list, which may be bad

float medianAngle (vector<float> *list); //Not using pointer because the process involves sorting which will modify the list, which may be bad

float mean (vector<float> *v, int start = -1, int end = -1);// take mean from start to end indices of vector v. 0 indexed

vector<float> stdev(vector<float> *v);//returns {mean, sample standard deviation}. Created by Hrant and modified by Jeff

float meanAngle (vector<float> *list);

vector<float> stdevAngle(vector<float> *v);

//float mode (vector<float> list);//Didn't finish!!!! Decided to write medianAngle() instead.

void substractComplexPhase(vector<float> *a, vector<float> *b, float angle);

bool invert(vector<vector<int> > * AtNinvAori, vector<vector<double> > * AtNinvAinv );
bool invert(vector<vector<float> > * AtNinvAori, vector<vector<double> > * AtNinvAinv );
///////////////REDUNDANT BASELINE CALIBRATION STUFF///////////////////
/////////////////////////////////////////////

void vecmatmul(vector<vector<double> > * Afitting, vector<float> * v, vector<float> * ampfit);
void vecmatmul(vector<vector<float> > * Afitting, vector<float> * v, vector<float> * ampfit);
void runAverage1d(vector<float> *in, vector<float> *out, uint w);//compute running average with running length 2w+1. The first and last w elements are averaged with less elements.
void runAverage(vector<vector<vector<vector<float> > > > *in, int dimension, int w);//compute running average with running length 2w+1 along dimension. The first and last w elements are averaged with less elements. Input array is modified!
#endif
