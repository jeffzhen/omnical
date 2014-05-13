#ifndef CALIBRATION_OMNI_H
#define CALIBRATION_OMNI_H

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
using namespace std;

struct redundantinfo{
	int nAntenna;//number of good antennas among all (64) antennas, same as the length of subsetant
	int nUBL;//number of unique baselines
	int nBaseline;
	int nCross;
	vector<int> subsetant;//the index of good antennas in all (64) antennas, currently not used
	vector<vector<float> > antloc;//3d antloc for each good antenna. strongly recommend using idealized antloc rather than literal ones
	vector<int> subsetbl;//the index of good baselines (auto included) in all baselines
	vector<vector<float> > ubl;//unique baseline vectors
	vector<int> bltoubl;//cross bl number to ubl index
	vector<int> reversed;//cross only bl if reversed -1, otherwise 1
	vector<int> reversedauto;//auto included bl if reversed -1, otherwise 1, currently not used
	vector<int> autoindex;//index of auto bls among good bls, currently not used
	vector<int> crossindex;//index of cross bls among good bls
	vector<vector<int> > bl2d;//from 1d bl to a pair of antenna numbers, (0,0), (0,1) (0,2) etc
	vector<int> ublcount;//for each ubl, the number of good cross bls corresponding to it
	vector<vector<vector<int> > > ublindex;//for each ubl, the vector<int> contains (ant1, ant2, autobl)
	vector<vector<int> > bl1dmatrix;//a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index counting auto corr
	vector<vector<float> > degenM;//degenM.(phase calibrations as a vector of nAnt) will generate a vector of length (nAnt + nUBL) which, when added to the existing calibration parameters' phases, will remove the linear phase field 
	vector<vector<int> > A;//A matrix for logcal amplitude
	vector<vector<int> > B;//B matrix for logcal phase
	vector<vector<int> > Atsparse;//At matrix for logcal amplitude, in sparse form, all entries are one
	vector<vector<vector<int> > > Btsparse;//Bt matrix for logcal phase, entries can be -1 or one so there's one more dimension than Atsparse
	vector<vector<float> > AtAi;//(AtA)^-1
	vector<vector<float> > BtBi;//(BtB)^-1
	vector<vector<float> > AtAiAt;//(AtA)^-1At
	vector<vector<float> > BtBiBt;//(BtB)^-1Bt
	vector<vector<float> > PA;//A(AtA)^-1At
	vector<vector<float> > PB;//B(BtB)^-1Bt
	vector<vector<float> > ImPA;//I-PA
	vector<vector<float> > ImPB;//I-PB

};

struct calmemmodule{//temporary memory modules for logcaladditive and lincal
	vector<float> amp1;
	vector<float> amp2;
	vector<float> amp3;
	vector<float> pha1;
	vector<float> pha2;
	vector<float> pha3;
	vector<float> x1;
	vector<float> x2;
	vector<float> x3;
	vector<float> x4;
	vector<vector<float> > g0;
	vector<vector<float> > g1;
	vector<vector<float> > g2;
	vector<vector<float> > g3;
	vector<vector<float> > adata1;
	vector<vector<float> > adata2;
	vector<vector<float> > adata3;
	vector<vector<float> > cdata1;
	vector<vector<float> > cdata2;
	vector<vector<float> > cdata3;
	vector<vector<float> > ubl1;//nubl by 2
	vector<vector<float> > ubl2;
	vector<vector<float> > ubl3;
	vector<vector<float> > ubl0;
	vector<vector<float> > ublgrp1;//regrouped baseline according to ubl, second dimension of unequal length of ubl count
	vector<vector<float> > ublgrp2;
	vector<vector<vector<float> > > ubl2dgrp1;//regrouped baseline according to ubl, second dimension of unequal length of ubl count
	vector<vector<vector<float> > > ubl2dgrp2;
};

bool readredundantinfo(string filepath, redundantinfo* info);
void initcalmodule(calmemmodule* module, redundantinfo* info);

///////////////////////////////////////
//command line and Python interaction//
///////////////////////////////////////
void printv(vector<float> *v, int mini = 0, int maxi = 200);
void printvv(vector<vector<float> > *v, int mini = 0, int maxi = 200);
void printvv(vector<vector<int> > *v, int mini = 0, int maxi = 200);
void printv(vector<double> *v, int mini = 0, int maxi = 200);
void printv(vector<int> *v, int mini = 0, int maxi = 200);

string exec(string input_command);//return result of stdout of a system command such as "python BLAH_BLAH". Won't catch stderr. Direct copy from internet http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c

string ftostr(float f);//float to string

string itostr(int i, int len);//int to string of specified length len

vector<float> strtovf(string in);

//!!!!!!pyEphem always deals with UTC, i.e., Coordinated Universal Time, and never the time in any particular time zone. 

vector<float> pySunPosForks(string date, string time, float time_offset);//date and time in UTC such as ($python sunpos.py '2012/5/24' '12:22:56') and returns altitude(degree from horizon to north) and az (degree from north towards east)

void pyStarPosForks(string date, string time, float time_offset, vector<vector<float> >* results);//date and time in UTC such as ($python sunpos.py '2012/5/24' '12:22:56') and fills in altitude(degree from horizon to north) and az (degree from north towards east) for sun, sagA, casA, cygA, crab, and virgo A

void pySatPosForks(string date, string time, float time_offset, vector<vector<float> >* results);

void pySatPos(string date, string time, float time_offset, float lat, float lon, float ele, vector<vector<float> >* results);//lat lon in degrees, ele in meters

void pyRequ2hor(string date, string time, float time_offset, float lat, float lon, vector<vector<float> >* results);//lat lon in degrees. Using this matrix in mathematica generates results that agree with pyephem output within 0.02 degrees on the few sources I checked at Forks

vector<float> pySunPos(string date, string time, float time_offset, float lat, float lon, float ele);


vector<string> pyTimeShift(string date, string time, float time_offset);//date and time in UTC such as 2012/5/24 12:22:56 and offset in seconds and returns [y/m/d, h:m:s])

float pyTimeDif(string time1, string time2);//returns time1-time2 in seconds

string pyOdfDate(string date, string time); //This script takes in UTC date and time such as '2012/5/24' '12:22:56' and prints the date format required by omniviewer Date('May 27, 2012,  20:31:15.0 UTC')

///////////////////////////////////////
//File I/O interaction/////////////////
///////////////////////////////////////
int cmdCountLine(const char* inputfilename);

string cmdAbsPath(string relativePath);

vector<string> cmdLs(string options);// options such as "-A ../test.odf". will return all items in their absolute path form **remember to use -d option for directories such as .odf**

void cmdMove(string a, string b);

void cmdCopy(string a, string b);

void readBinaryVisibility(const char* inputfilename, vector<vector<vector<vector<vector<float> > > > > * data, int nPol, int nIntegrations, int nFrequencies, int nBase);//Modified from Devon Rosner's OmniODF_IQ_ODF.cc code

void readBinaryVisibilityLarge(const char* inputfilename, vector<vector<vector<vector<float> > > > * data, int nPol, int nIntegrations, int nFrequencies, int nBase);//Modified from Devon Rosner's OmniODF_IQ_ODF.cc code

void readBinaryVisibilityLargeConj(const char* inputfilename, vector<vector<vector<vector<float> > > > * data, int nPol, int nIntegrations, int nFrequencies, int nBase);//Modified from Devon Rosner's OmniODF_IQ_ODF.cc code

void breakLarge(vector<float> *largeslice, vector<vector<float> > * smallslice);// breaks up the frequency slice in large format (1D of length 2*nBaseline) into small format(2D of nBaseline by re/im) 

void padSmall(vector<vector<float> > *smallslice, vector<float> * largeslice);// pad the frequency slice in small format(2D of nBaseline by re/im) into large format (1D of length 2*nBaseline)   

void readBinaryCalparSP(const char* inputfilename, vector<vector<vector<float> > > * data, int nIntegrations, int nFrequencies, int nAnt, int nUBL);//read binary calpar single polarization, assumes chisqx3, log10(ampcalpar) x nant, phasecalpar degrees x nant, ubl fits in r/i r/i; //turn degree into rad when reading

vector<float> readAscii(const char* inputfilename, int count = -1, bool verbose = false);

void readVisibility(const char* inputfilename, vector<vector<vector<float> > > * data, int nFreq, int nBL);

bool readCalparAscii(const char* inputfilename, vector<float> *chisq, vector<vector<float> > *ampcalpar, vector<vector<float> > *phasecalpar, vector<vector<vector<float> > > * UBLcalpar, int nFreq, int nAnt, int nUBL);

void readAntloc(const char* inputfilename, vector<vector<float> > * antloc, vector<vector<float> > * cablelen, int nAnt);//cablelen structured as [x/y][nAnt]

void readSunpos(const char* inputfilename, vector<vector<float> > * sunpos);//read sunpos.dat from x5 odf (which is in alt/az, and return a list of pairs in k vector (x y z) or (S E U), note that this is technically -k vector since k vector should point inwards.

void outputChiSqAscii(vector<float>  * data); //takes chi square for one time slice each freq each antenna and write a set of calibration parameters under current directory

void outputPhaseCalParAscii(vector<vector<float> >  * phasecalpar, int numUBL, string output_name);

void outputCalParAscii(vector<float> * chisq, vector<vector<float> > * ampcalpar, vector<vector<float> >  * phasecalpar, vector<vector<vector<float> > > * UBLcalpar, string output_name); //takes phase corrections for one time slice each freq each antenna and write a set of calibration parameters under current directory

void outputAscii(vector<float>* data, string output_name, int type = 1, bool cmdoutput = true);//0 for overwrite, 1 for append
void outputAscii(vector<vector<float> >  * data, string output_name, int type = 1, bool cmdoutput = true);//0 for overwrite, 1 for append
void outputAscii(vector<vector<vector<float> > >  * data, string output_name, int type = 1, bool cmdoutput = true);//0 for overwrite, 1 for append

void outputDataAscii(vector<vector<vector<float> > >  * data, string output_name);

bool outputCalpar(vector<vector<vector<vector<float> > > > * data, string outputfilename, bool in_degree = true, int nAnt = -1);// outputs binary calpar file. in_degree means if phase calpar is in degree, default true
bool outputCalparSP(vector<vector<vector<float> > > * data, string outputfilename, bool in_degree = true, int nAnt = -1);// outputs binary calpar file. in_degree means if phase calpar is in degree, default true

bool outputData(vector<vector<vector<vector<vector<float> > > > > * data, string outputfilename);// outputs binary file Modified from Devon Rosner's OmniODF_IQ_ODF.cc code

bool outputDataLarge(vector<vector<vector<vector<float> > > > * data, string outputfilename);// outputs binary file Modified from Devon Rosner's OmniODF_IQ_ODF.cc code

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

void iqDemod(vector<vector<vector<vector<vector<float> > > > > *data, vector<vector<vector<vector<vector<float> > > > > *data_out, int nIntegrations, int nFrequencies, int nAnt);

void iqDemodLarge(vector<vector<vector<vector<float> > > > *data, vector<vector<vector<vector<float> > > > *data_out, int nIntegrations, int nFrequencies, int nAnt);//For large set of data (more than ~200 time slices), it takes more than 24G of memory to allocate a large enough 5D vector. So we consolidate the last two dimensions (baseline and re/im) together

int gc(int a, int b);//Get coordinate for large vectors when the last two coordinates are consolidated into one dimension. 2a+b

/*void iqDemodBinary(string data, string data_out, int nIntegrations, int nFrequencies, int nAnt);//Used for large odfs for which there's not enough memory to take in all the data, so we have to read and iq and write without saving too much into memory*/

string strReplace(string input, string a, string b);//replace string a contained in input by string b

string getFileName(string fileNameWithPath);//get the filename in a long path/path/filename

string extFileName(string fileName, string ext); //extend a file name by a string, such as extend way.cool.odf with shit to get way.coolshit.odf

vector<string> parseLines(string bigLine);//can be used to parse result returned by exec() which is a big string containing multiple lines

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

void mergeCalibrationPar (vector<float> * ampcalpar1, vector<float> * ampcalpar2, vector<float> * ampcalparM, vector<float> * phasecalpar1, vector<float> * phasecalpar2, vector<float> * phasecalparM);//Only deals with ampcalpar and phasecalpar, Chisq should be the chisq of the second set, so are the UBLs

void addPhase(vector<float> * x, float phi);

float vectorDot(vector<float>* v1, vector<float>* v2);//v1.v2

vector<float> matrixDotV(vector<vector<float> >* m, vector<float>* v);//m.v

vector<vector<float> > rotationMatrix(float x, float y, float z);//approximation for a rotation matrix rotating around x,y,z axis, {{1, -z, y}, {z, 1, -x}, {-y, x, 1}}

vector<vector<float> > rotationMatrixZ(float z);

bool createAmatrix(vector<vector<int> > *receiverAmatrix, vector<vector<float> > *antloc);

bool createBmatrix(vector<vector<int> > *receiverBmatrix, vector<vector<float> > *antloc);//convention for reversed baselines are to flip the two ant calpars signs, and expect user to flip the corresponding data's sign as well

bool findReversedBaselines(vector<int> *receiverList, vector<vector<float> > *antloc);

///////////////PHASE CALIBRATE STUFF///////////////////
/////////////////////////////////////////////
float getCableCorrection(float length, float freq);//measured phase - phaseCor = ideal phase


float getFreq(int f, int nFrequency, float startFreq, float endFreq);//To be consistent with OmniViewer frequency interpretation!!!

vector<float> getBL(int i, int j, vector<vector<float> > *antloc);

int countUBL(vector<vector<float> > *antloc );

int lookupAnt(float x, float y, vector<vector<float> > antloc);

float phaseWrap (float x, float offset = -atan2(0,-1));//Wrap phase to be on (offset, offset+2pi]


void phaseCalibrate120(vector<float>* calpar120, vector<float>* calpar16, int nAnt, vector<bool>* badAnt);//find the median solution of 16 antenna calpars from 120 visibility calpars

vector<float> phaseCalibrate(vector<vector<float> > *dataf, string pol, float freq, vector<vector<float> > *antloc, vector<vector<float> > *cablelen, int startingAnt1, int startingAnt2, int startingAnt3, int nAntenna);



///////////////CHI SQUARE STUFF///////////////////
/////////////////////////////////////////////
void computeUBL(vector<vector<float> > * antloc, vector<vector<float> > * listUBL);


vector<float> modelToMeasurement(vector<float> *modelCor, float ampcal, float phasecal);

vector<float> correctMeasurement(vector<float> *measuredCor, float ampcal1, float phasecal1, float ampcal2, float phasecal2);//assumes <x1* , x2>, correct a specific visibility OPTIM

void correctMeasurementPhase(vector<float> *measuredCor, float phasecal1, float phasecal2);

void correctMeasurementMatrix(vector<vector<float> > *correctedData, vector<float> *ampcalpar, vector<float> *phasecalpar);//correct a specific frequency OPTIM

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

void pointSourceCalAccordingTo(int referenceAnt, vector<vector<float> > *data, vector<vector<float> > *ampcalparArray, vector<vector<float> > *phasecalparArray);

void pointSourceCal(vector<vector<float> > *data, vector<float> *ampcalpar, vector<float> *phasecalpar, vector<vector<float> > *UBLcalpar);

void substractComplexPhase(vector<float> *a, vector<float> *b, float angle);

void rotateCalpar(vector<float> *originalPhase, vector<float> *rotatedPhase, vector<vector<float> > *originalUBLcor, vector<vector<float> > *rotatedUBLcor, vector<vector<float> > *antloc, vector<vector<float> > * listUBL, float theta, float phi, float freq);//theta in [0, PI/2] (rotating z axis down), phi [0, 2PI), freq in MHz, ONLY WORKS FOR ROTATING POINT SOURCE STRAIGHT ABOVE *TO* THETA AND PHI 

bool invert(vector<vector<int> > * AtNinvAori, vector<vector<double> > * AtNinvAinv );
bool invert(vector<vector<float> > * AtNinvAori, vector<vector<double> > * AtNinvAinv );
///////////////REDUNDANT BASELINE CALIBRATION STUFF///////////////////
/////////////////////////////////////////////

void vecmatmul(vector<vector<double> > * Afitting, vector<float> * v, vector<float> * ampfit);
void vecmatmul(vector<vector<float> > * Afitting, vector<float> * v, vector<float> * ampfit);
void logcaladd(vector<vector<float> >* data, vector<vector<float> >* additivein, redundantinfo* info, vector<float>* calpar, vector<vector<float> >* additiveout, int command, calmemmodule* module);//option 0 only computes additive term and chi^2 but does not compute calibration parameters, option 1 does compute everything
void lincal(vector<vector<float> >* data, vector<vector<float> >* additivein, redundantinfo* info, vector<float>* calpar, vector<vector<float> >* additiveout, int command, calmemmodule* module, float convergethresh, int maxiter, float stepsize);//if command is 1, compute the ubl estimates given data and calpars, rather than read ubl estimates from input; additive term will only be updated if lincal can achieve a lower chi^2
void loadGoodVisibilities(vector<vector<vector<vector<float> > > > * rawdata, vector<vector<vector<vector<float> > > >* receiver, redundantinfo* info, int xy);
void removeDegen(vector<float> *calpar, redundantinfo * info, calmemmodule* module);//forces the calibration parameters to have average 1 amp, and no shifting the image in phase. Note: 1) If you have not absolute calibrated the data, there's no point in running this, because this can only ensure that the calpars don't screw up already absolute calibrated data. 2) the antloc and ubl stored in redundant info must be computed from idealized antloc, otherwise the error in antloc from perfect redundancy will roll into this result, in an unknown fashion.
void runAverage1d(vector<float> *in, vector<float> *out, int w);//compute running average with running length 2w+1. The first and last w elements are averaged with less elements.
void runAverage(vector<vector<vector<vector<float> > > > *in, int dimension, int w);//compute running average with running length 2w+1 along dimension. The first and last w elements are averaged with less elements. Input array is modified!
#endif
