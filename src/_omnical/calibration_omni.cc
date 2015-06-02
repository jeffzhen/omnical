//Jeff 2012-07-19
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
#include <numeric>
#include <algorithm>
#include <functional>
#include <numeric>
#include "include/calibration_omni.h"
#include <algorithm>
#define uint unsigned int
using namespace std;
const string FILENAME = "calibration_omni.cc";
const float SPEEDC = 299.792458;
const float PI = atan2(0, -1);
const bool DEBUG = false;
const float UBLPRECISION = pow(10, -3);
const float MIN_NONE_ZERO = pow(10, -10);
const float MAX_NONE_INF = pow(10, 10);
const float MAX_10_POW = 20; //single precision float max is 3.4*10^38, so I'm limiting the power to be 20.
const float MAX_POW_2 = pow(10, 10); //limiting the base of ^2 to be 10^10
const float X4_LONGITUDE = -69.987182;
const float X4_LATITUDE = 45.297728;
const float X4_ELEVATION = 171;
const float X4_TIMESHIFT = 28957;//seconds to add to raw data header time to get correct UTC
const float DEF_LONGITUDE = -69.987182;
const float DEF_LATITUDE = 45.297728;
const float DEF_ELEVATION = 171;
const int NUM_OBJECTS = 30;//Number of satellites we have in tracked_bodies_X4.tle;

void initcalmodule(calmemmodule* module, redundantinfo* info){
	int nant = info->nAntenna;
	int nbl = info->nBaseline;
	int nubl = info->nUBL;
	int ncross = info->crossindex.size();
	(module->amp1).resize(ncross);
	(module->amp2).resize(ncross);
	(module->amp3).resize(ncross);
	(module->pha1).resize(ncross);
	(module->pha2).resize(ncross);
	(module->pha3).resize(ncross);
	(module->x1).resize(nubl + nant);
	(module->x2).resize(nubl + nant);
	(module->x3).resize(nubl + nant);
	(module->x4).resize(nubl + nant);

	(module->g1).resize(nant);
	for (int i = 0; i < nant; i++){
		(module->g1)[i].resize(2);
	}
	(module->g0) = (module->g1);
	(module->g2) = (module->g1);
	(module->g3) = (module->g1);

	(module->adata1).resize(nbl);
	for (int i = 0; i < nbl; i++){
		(module->adata1)[i].resize(2);
	}
	(module->adata2) = (module->adata1);
	(module->adata3) = (module->adata1);

	(module->cdata1).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(module->cdata1)[i].resize(2);
	}
	(module->cdata2) = (module->cdata1);
	(module->cdata3) = (module->cdata1);


	(module->ubl1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ubl1)[i].resize(2);
	}
	(module->ubl0) = (module->ubl1);
	(module->ubl2) = (module->ubl1);
	(module->ubl3) = (module->ubl1);

	(module->ublgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ublgrp1)[i].resize(info->ublcount[i]);
	}
	(module->ublgrp2) = (module->ublgrp1);

	(module->ubl2dgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ubl2dgrp1)[i].resize(info->ublcount[i]);
		for (int j = 0; j < info->ublcount[i]; j ++){
			(module->ubl2dgrp1)[i][j].resize(2);
		}
	}
	(module->ubl2dgrp2) = (module->ubl2dgrp1);
	return;
}

///////////////////////////////////////
//command line and Python interaction//
///////////////////////////////////////

string ftostr(float f){
	ostringstream buffer;
	buffer << f;
	return buffer.str();
}

string itostr(int i, uint len){
	ostringstream buffer;
	buffer << abs(i);
	string raw = buffer.str();//unpadded int
	string output;
	if ( i >= 0) {
		output = "";
	} else {
		output = "-";
	}
	for (uint n = 0; n < len - raw.size(); n ++) output += "0";
	output += raw;
	return output;
};

vector<float> strtovf(string in){
	//cout << "DEBUG " << in << endl;
	stringstream stream(in);
	vector<float> output;
	float tmpf;
	while (stream.good()){
		stream >> tmpf;
		//cout << tmpf << endl;
		output.push_back(tmpf);
	}
	//output.pop_back(); //sometimes get a 0 in end..tricky...
	return output;
}

void breakLarge(vector<float> *largeslice, vector<vector<float> > * smallslice){// breaks up the frequency slice in large format (1D of length 2*nBaseline) into small format(2D of nBaseline by re/im)
	string METHODNAME = "breakLarge";
	if ( largeslice->size() != 2*(smallslice->size()) or (smallslice->at(0)).size() != 2 ){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array is initialized at (nBaseline, 2) = (" << smallslice->size() << " by " << (smallslice->at(0)).size() << "), where as the large slice is specified as (" << largeslice->size() << "). Exiting!!" << endl;
		return;
	}
	for (unsigned int i = 0; i < largeslice->size(); i++){
		(smallslice->at(floor(i/2)))[i%2] = largeslice->at(i);
	}
	return;
}

void padSmall(vector<vector<float> > * smallslice, vector<float> * largeslice){// pad the frequency slice in small format(2D of nBaseline by re/im) into large format (1D of length 2*nBaseline)
	string METHODNAME = "padSmall";
	if ( largeslice->size() != 2*(smallslice->size()) or (smallslice->at(0)).size() != 2 ){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The large array is initialized at (nBaseline, 2) = (" << smallslice->size() << " by " << (smallslice->at(0)).size() << "), where as the small slice is specified as (" << largeslice->size() << "). Exiting!!" << endl;
		return;
	}
	for (unsigned int i = 0; i < largeslice->size(); i++){
		largeslice->at(i) = (smallslice->at(floor(i/2)))[i%2];
	}
	return;
}

vector<float> tp2xyz (vector<float> thephi){
	vector<float> xyz(3,0);
	xyz[0] = sin(thephi[0])*cos(thephi[1]);
	xyz[1] = sin(thephi[0])*sin(thephi[1]);
	xyz[2] = cos(thephi[0]);
	return xyz;
}
vector<float> tp2xyz (float t, float p){
	vector<float> xyz(3,0);
	xyz[0] = sin(t)*cos(p);
	xyz[1] = sin(t)*sin(p);
	xyz[2] = cos(t);
	return xyz;
}

vector<float> xyz2tp (vector<float> xyz){
	string METHODNAME = "xyz2tp";
	vector<float> thephi(2,0);
	float r = sqrt(xyz[2]*xyz[2]+xyz[1]*xyz[1]+xyz[0]*xyz[0]);
	if( r == 0){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": !!FATAL error!! input x,y,z are all 0!" << endl;
		return thephi;
	}
	thephi[0] = acos(xyz[2] / r);

	thephi[1] = atan2(xyz[1],xyz[0]);
	return thephi;
}

vector<float> xyz2tp (float x, float y, float z){
	string METHODNAME = "xyz2tp";
	vector<float> thephi(2,0);
	float r = sqrt(x*x + y*y + z*z);
	if( r == 0){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": !!FATAL error!! input x,y,z are all 0!" << endl;
		return thephi;
	}
	thephi[0] = acos(z / r);

	thephi[1] = atan2(y,x);
	return thephi;
}

vector<float> tp2rd (vector<float> thephi){
	vector<float> rd(2,0);
	rd[0] = thephi[1];
	rd[1] = PI/2 - thephi[0];
	return rd;
}
vector<float> tp2rd (float t, float p){
	vector<float> rd(2,0);
	rd[0] = p;
	rd[1] = PI/2 - t;
	return rd;
}
vector<float> rd2tp (vector<float> rd){
	vector<float> tp(2,0);
	tp[1] = rd[0];
	tp[0] = PI/2 - rd[1];
	return tp;
}
vector<float> rd2tp (float r, float d){
	vector<float> tp(2,0);
	tp[1] = r;
	tp[0] = PI/2 - d;
	return tp;
}
vector<float> tp2aa (vector<float> thephi){//alt-az
	vector<float> aa(2,0);
	aa[0] = PI/2 - thephi[0];
	aa[1] = PI - thephi[1];
	return aa;
}
vector<float> tp2aa (float t, float p){//alt-az
	vector<float> aa(2,0);
	aa[0] = PI/2 - t;
	aa[1] = PI - p;
	return aa;
}
vector<float> aa2tp (vector<float> aa){//alt-az
	vector<float> tp(2,0);
	tp[1] = PI - aa[1];
	tp[0] = PI/2 - aa[0];
	return tp;
}
vector<float> aa2tp (float alt, float az){
	vector<float> tp(2,0);
	tp[1] = PI - az;
	tp[0] = PI/2 - alt;
	return tp;
}

void matrixDotV(vector<vector<float> > * A, vector<float> * b, vector<float> * x){
	int i, j;
	double sum;
	int n = min(A->size(),x->size());
	int m = min(A->at(0).size(),b->size());
	for(i = 0; i < n; i++){
		sum = 0.0;
		for(j = 0; j < m; j++){
			sum = sum + (A->at(i))[j] * (x->at(j));
		}
		(x->at(i)) = sum;
	}
	return;
}

//void iqDemod(vector<vector<vector<vector<vector<float> > > > > *data, vector<vector<vector<vector<vector<float> > > > > *data_out, int nIntegrations, int nFrequencies, int nAnt){
	//string METHODNAME = "iqDemod";
	//int nChannels = nAnt * 4; //a factor of 2 from consolidating x and y polarizations, and another factor of 2 from consolidating iq
	//int n_xxi = nAnt * (nAnt + 1)/2;

	//if ( data->size() != 1 or data_out->size() != 4 or (data->at(0)).size() != nIntegrations or (data_out->at(0)).size() != nIntegrations or (data->at(0))[0].size() != nFrequencies or (data_out->at(0))[0].size() != 2 * nFrequencies or (data->at(0))[0][0].size() != nChannels * ( nChannels + 1 ) / 2  or (data_out->at(0))[0][0].size() != nAnt * ( nAnt + 1 ) / 2 ){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input array and IQ array are initialized at (p, t, f, bl) = (" << data->size() << ", " << (data->at(0)).size() << ", " <<  (data->at(0))[0].size()  << ", " <<  (data->at(0))[0][0].size() << ") and (" << data_out->size() << ", " << (data_out->at(0)).size() << ", " <<  (data_out->at(0))[0].size()  << ", " <<  (data_out->at(0))[0][0].size() << "), where as the parameters are specified as (t, f, ant) = (" << nIntegrations << ", "  << nFrequencies << ", " << nAnt << "). Exiting!!" << endl;
		//return;
	//}
	//vector<vector<float> > *freq_slice;
	//int prevk, k1i, k1q, prevk1i, prevk1q, k2xi, k2xq, k2yi, k2yq, prevk2xi, prevk2yi, bl;
	//float a1xx_re, a1xx_im, a2xx_re, a2xx_im, a3xx_re, a3xx_im, a1xy_re, a1xy_im, a2xy_re, a2xy_im, a3xy_re, a3xy_im, a1yx_re, a1yx_im, a2yx_re, a2yx_im, a3yx_re, a3yx_im, a1yy_re, a1yy_im, a2yy_re, a2yy_im, a3yy_re, a3yy_im;
	//int c2nchan1 = 2 * nChannels - 1; //frequently used constant
	//for (int t = 0; t < nIntegrations; t++){
		////cout << t << endl;
		//for (int f = 0; f < nFrequencies; f++){
			//freq_slice = &((data->at(0))[t][f]);
			////loop for xx and xy
			//for (int k1 = 0; k1 < nAnt; k1++){
				//prevk = (2 * nAnt - k1 - 1) * k1 / 2;
				//k1i = 2*k1;
				//k1q = k1i + 2 * nAnt;
				//prevk1i = (c2nchan1 - k1i)*k1i/2;
				//prevk1q = (c2nchan1 - k1q)*k1q/2;
				//for (int k2 = k1; k2 < nAnt; k2++){
					//k2xi = 2 * k2;
					//k2xq = k2xi + 2 * nAnt;
					//k2yi = k2xi + 1;
					//k2yq = k2xq + 1;
					//prevk2xi = (c2nchan1 - k2xi) * k2xi / 2;
					//prevk2yi = (c2nchan1-k2yi) * k2yi / 2;
					//// performing complex arithmetic: 0 index --> real
					//// 1 index --> imag
					//a1xx_re = freq_slice->at(prevk1i+k2xi)[0] + freq_slice->at(prevk1q+k2xq)[0];
					//a1xx_im = freq_slice->at(prevk1i+k2xi)[1] + freq_slice->at(prevk1q+k2xq)[1];
					//a2xx_re = freq_slice->at(prevk1i+k2xq)[0] - freq_slice->at(prevk2xi+k1q)[0];
					//a2xx_im = freq_slice->at(prevk1i+k2xq)[1] + freq_slice->at(prevk2xi+k1q)[1];
					//a3xx_re = -1 * a2xx_im;
					//a3xx_im = a2xx_re;
					//a1xy_re = freq_slice->at(prevk1i+k2yi)[0] + freq_slice->at(prevk1q+k2yq)[0];
					//a1xy_im = freq_slice->at(prevk1i+k2yi)[1] + freq_slice->at(prevk1q+k2yq)[1];
					//a2xy_re = freq_slice->at(prevk1i+k2yq)[0] - freq_slice->at(prevk2yi+k1q)[0];
					//a2xy_im = freq_slice->at(prevk1i+k2yq)[1] + freq_slice->at(prevk2yi+k1q)[1];
					//a3xy_re = -1 * a2xy_im;
					//a3xy_im = a2xy_re;

					////writing to output matrix
					//bl = prevk + k2;
					//if (f == 0){
						//(data_out->at(0))[t][2*nFrequencies-1][bl][0] = ( a1xx_re + a3xx_re);
						//(data_out->at(0))[t][2*nFrequencies-1][bl][1] = -1*( a1xx_im + a3xx_im);
						//(data_out->at(1))[t][2*nFrequencies-1][bl][0] = (a1xy_re + a3xy_re);
						//(data_out->at(1))[t][2*nFrequencies-1][bl][1] = -1*(a1xy_im + a3xy_im);
					//}

					//(data_out->at(0))[t][nFrequencies-1+f][bl][0] = ( a1xx_re + a3xx_re);
					//(data_out->at(0))[t][nFrequencies-1+f][bl][1] = -1*( a1xx_im + a3xx_im);
					//(data_out->at(0))[t][nFrequencies-1-f][bl][0] = a1xx_re - a3xx_re;
					//(data_out->at(0))[t][nFrequencies-1-f][bl][1] = a1xx_im - a3xx_im;
					//(data_out->at(1))[t][nFrequencies-1+f][bl][0] = (a1xy_re + a3xy_re);
					//(data_out->at(1))[t][nFrequencies-1+f][bl][1] = -1*(a1xy_im + a3xy_im);
					//(data_out->at(1))[t][nFrequencies-1-f][bl][0] = a1xy_re - a3xy_re;
					//(data_out->at(1))[t][nFrequencies-1-f][bl][1] = a1xy_im - a3xy_im;
				//}
			//}
				////loop for yy and yx
				////computational difference: k1i = 2*k1 (+ 1)
			//for (int k1=0; k1 < nAnt; k1++){
				//prevk = (2*nAnt-k1-1)*k1/2;
				//k1i = 2*k1 + 1;
				//k1q = k1i + 2 * nAnt;
				//prevk1i = (c2nchan1 - k1i)*k1i/2;
				//prevk1q = (c2nchan1 - k1q)*k1q/2;
				//for (int k2=k1; k2 < nAnt; k2++){
					//k2xi = 2*k2;
					//k2xq = k2xi + 2*nAnt;
					//k2yi = k2xi + 1;
					//k2yq = k2xq + 1;
					//prevk2xi = (c2nchan1-k2xi)*k2xi/2;
					//prevk2yi = (c2nchan1-k2yi)*k2yi/2;
					//// performing complex arithmetic: 0 index --> real
					//// 1 index --> imag
					//a1yx_re = freq_slice->at(prevk1i+k2xi)[0] + freq_slice->at(prevk1q+k2xq)[0];
					//a1yx_im = freq_slice->at(prevk1i+k2xi)[1] + freq_slice->at(prevk1q+k2xq)[1];
					//a2yx_re = freq_slice->at(prevk1i+k2xq)[0] - freq_slice->at(prevk2xi+k1q)[0];
					//a2yx_im = freq_slice->at(prevk1i+k2xq)[1] + freq_slice->at(prevk2xi+k1q)[1];
					//a3yx_re = -1 * a2yx_im;
					//a3yx_im = a2yx_re;
					//a1yy_re = freq_slice->at(prevk1i+k2yi)[0] + freq_slice->at(prevk1q+k2yq)[0];
					//a1yy_im = freq_slice->at(prevk1i+k2yi)[1] + freq_slice->at(prevk1q+k2yq)[1];
					//a2yy_re = freq_slice->at(prevk1i+k2yq)[0] - freq_slice->at(prevk2yi+k1q)[0];
					//a2yy_im = freq_slice->at(prevk1i+k2yq)[1] + freq_slice->at(prevk2yi+k1q)[1];
					//a3yy_re = -1 * a2yy_im;
					//a3yy_im = a2yy_re;

					////writing to output matrix
					//bl = prevk + k2;
					//if (f == 0){
						//(data_out->at(2))[t][2*nFrequencies-1][bl][0] = ( a1yx_re + a3yx_re);
						//(data_out->at(2))[t][2*nFrequencies-1][bl][1] = -1*( a1yx_im + a3yx_im);
						//(data_out->at(3))[t][2*nFrequencies-1][bl][0] = (a1yy_re + a3yy_re);
						//(data_out->at(3))[t][2*nFrequencies-1][bl][1] = -1*(a1yy_im + a3yy_im);
					//}
					//(data_out->at(2))[t][nFrequencies-1+f][bl][0] = (a1yx_re + a3yx_re);
					//(data_out->at(2))[t][nFrequencies-1+f][bl][1] = -1*(a1yx_im + a3yx_im);
					//(data_out->at(2))[t][nFrequencies-1-f][bl][0] = a1yx_re - a3yx_re;
					//(data_out->at(2))[t][nFrequencies-1-f][bl][1] = a1yx_im - a3yx_im;
					//(data_out->at(3))[t][nFrequencies-1+f][bl][0] = (a1yy_re + a3yy_re);
					//(data_out->at(3))[t][nFrequencies-1+f][bl][1] = -1*(a1yy_im + a3yy_im);
					//(data_out->at(3))[t][nFrequencies-1-f][bl][0] = a1yy_re - a3yy_re;
					//(data_out->at(3))[t][nFrequencies-1-f][bl][1] = a1yy_im - a3yy_im;
				//}
			//}
		//}
	//}
	//return;
//}

//void iqDemodLarge(vector<vector<vector<vector<float> > > > *data, vector<vector<vector<vector<float> > > > *data_out, int nIntegrations, int nFrequencies, int nAnt){
	//string METHODNAME = "iqDemodLarge";
	//int nChannels = nAnt * 4; //a factor of 2 from consolidating x and y polarizations, and another factor of 2 from consolidating iq
	//int n_xxi = nAnt * (nAnt + 1)/2;

	//if ( data->size() != 1 or data_out->size() != 4 or (data->at(0)).size() != nIntegrations or (data_out->at(0)).size() != nIntegrations or (data->at(0))[0].size() != nFrequencies or (data_out->at(0))[0].size() != 2 * nFrequencies or (data->at(0))[0][0].size() != nChannels * ( nChannels + 1 ) or (data_out->at(0))[0][0].size() != nAnt * ( nAnt + 1 ) ){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input array and IQ array are initialized at (p, t, f, bl) = (" << data->size() << ", " << (data->at(0)).size() << ", " <<  (data->at(0))[0].size()  << ", " <<  (data->at(0))[0][0].size() << ") and (" << data_out->size() << ", " << (data_out->at(0)).size() << ", " <<  (data_out->at(0))[0].size()  << ", " <<  (data_out->at(0))[0][0].size() << "), where as the parameters are specified as (t, f_in, f_out, bl_in, bl_out) = (" << nIntegrations << ", "  << nFrequencies << ", "  << 2 * nFrequencies << ", " << nChannels * ( nChannels + 1 ) << ", " << nAnt * ( nAnt + 1 ) << "). Exiting!!" << endl;
		//return;
	//}
	//vector<float> *freq_slice;
	//int prevk, k1i, k1q, prevk1i, prevk1q, k2xi, k2xq, k2yi, k2yq, prevk2xi, prevk2yi, bl;
	//float a1xx_re, a1xx_im, a2xx_re, a2xx_im, a3xx_re, a3xx_im, a1xy_re, a1xy_im, a2xy_re, a2xy_im, a3xy_re, a3xy_im, a1yx_re, a1yx_im, a2yx_re, a2yx_im, a3yx_re, a3yx_im, a1yy_re, a1yy_im, a2yy_re, a2yy_im, a3yy_re, a3yy_im;
	//int c2nchan1 = 2 * nChannels - 1; //frequently used constant
	//for (int t = 0; t < nIntegrations; t++){
		////cout << t << endl;
		//for (int f = 0; f < nFrequencies; f++){
			//freq_slice = &((data->at(0))[t][f]);
			////loop for xx and xy
			//for (int k1 = 0; k1 < nAnt; k1++){
				//prevk = (2 * nAnt - k1 - 1) * k1 / 2;
				//k1i = 2*k1;
				//k1q = k1i + 2 * nAnt;
				//prevk1i = (c2nchan1 - k1i)*k1i/2;
				//prevk1q = (c2nchan1 - k1q)*k1q/2;
				//for (int k2 = k1; k2 < nAnt; k2++){
					//k2xi = 2 * k2;
					//k2xq = k2xi + 2 * nAnt;
					//k2yi = k2xi + 1;
					//k2yq = k2xq + 1;
					//prevk2xi = (c2nchan1 - k2xi) * k2xi / 2;
					//prevk2yi = (c2nchan1-k2yi) * k2yi / 2;
					//// performing complex arithmetic: 0 index --> real
					//// 1 index --> imag
					//a1xx_re = freq_slice->at(gc(prevk1i+k2xi, 0)) + freq_slice->at(gc(prevk1q+k2xq, 0));
					//a1xx_im = freq_slice->at(gc(prevk1i+k2xi, 1)) + freq_slice->at(gc(prevk1q+k2xq, 1));
					//a2xx_re = freq_slice->at(gc(prevk1i+k2xq, 0)) - freq_slice->at(gc(prevk2xi+k1q, 0));
					//a2xx_im = freq_slice->at(gc(prevk1i+k2xq, 1)) + freq_slice->at(gc(prevk2xi+k1q, 1));
					//a3xx_re = -1 * a2xx_im;
					//a3xx_im = a2xx_re;
					//a1xy_re = freq_slice->at(gc(prevk1i+k2yi, 0)) + freq_slice->at(gc(prevk1q+k2yq, 0));
					//a1xy_im = freq_slice->at(gc(prevk1i+k2yi, 1)) + freq_slice->at(gc(prevk1q+k2yq, 1));
					//a2xy_re = freq_slice->at(gc(prevk1i+k2yq, 0)) - freq_slice->at(gc(prevk2yi+k1q, 0));
					//a2xy_im = freq_slice->at(gc(prevk1i+k2yq, 1)) + freq_slice->at(gc(prevk2yi+k1q, 1));
					//a3xy_re = -1 * a2xy_im;
					//a3xy_im = a2xy_re;

					////writing to output matrix
					//bl = prevk + k2;
					//if (f == 0){
						//(data_out->at(0))[t][2*nFrequencies-1][gc(bl, 0)] = ( a1xx_re + a3xx_re);
						//(data_out->at(0))[t][2*nFrequencies-1][gc(bl, 1)] = -1*( a1xx_im + a3xx_im);
						//(data_out->at(1))[t][2*nFrequencies-1][gc(bl, 0)] = (a1xy_re + a3xy_re);
						//(data_out->at(1))[t][2*nFrequencies-1][gc(bl, 1)] = -1*(a1xy_im + a3xy_im);
					//}

					//(data_out->at(0))[t][nFrequencies-1+f][gc(bl, 0)] = ( a1xx_re + a3xx_re);
					//(data_out->at(0))[t][nFrequencies-1+f][gc(bl, 1)] = -1*( a1xx_im + a3xx_im);
					//(data_out->at(0))[t][nFrequencies-1-f][gc(bl, 0)] = a1xx_re - a3xx_re;
					//(data_out->at(0))[t][nFrequencies-1-f][gc(bl, 1)] = a1xx_im - a3xx_im;
					//(data_out->at(1))[t][nFrequencies-1+f][gc(bl, 0)] = (a1xy_re + a3xy_re);
					//(data_out->at(1))[t][nFrequencies-1+f][gc(bl, 1)] = -1*(a1xy_im + a3xy_im);
					//(data_out->at(1))[t][nFrequencies-1-f][gc(bl, 0)] = a1xy_re - a3xy_re;
					//(data_out->at(1))[t][nFrequencies-1-f][gc(bl, 1)] = a1xy_im - a3xy_im;
				//}
			//}
				////loop for yy and yx
				////computational difference: k1i = 2*k1 (+ 1)
			//for (int k1=0; k1 < nAnt; k1++){
				//prevk = (2*nAnt-k1-1)*k1/2;
				//k1i = 2*k1 + 1;
				//k1q = k1i + 2 * nAnt;
				//prevk1i = (c2nchan1 - k1i)*k1i/2;
				//prevk1q = (c2nchan1 - k1q)*k1q/2;
				//for (int k2=k1; k2 < nAnt; k2++){
					//k2xi = 2*k2;
					//k2xq = k2xi + 2*nAnt;
					//k2yi = k2xi + 1;
					//k2yq = k2xq + 1;
					//prevk2xi = (c2nchan1-k2xi)*k2xi/2;
					//prevk2yi = (c2nchan1-k2yi)*k2yi/2;
					//// performing complex arithmetic: 0 index --> real
					//// 1 index --> imag
					//a1yx_re = freq_slice->at(gc(prevk1i+k2xi, 0)) + freq_slice->at(gc(prevk1q+k2xq, 0));
					//a1yx_im = freq_slice->at(gc(prevk1i+k2xi, 1)) + freq_slice->at(gc(prevk1q+k2xq, 1));
					//a2yx_re = freq_slice->at(gc(prevk1i+k2xq, 0)) - freq_slice->at(gc(prevk2xi+k1q, 0));
					//a2yx_im = freq_slice->at(gc(prevk1i+k2xq, 1)) + freq_slice->at(gc(prevk2xi+k1q, 1));
					//a3yx_re = -1 * a2yx_im;
					//a3yx_im = a2yx_re;
					//a1yy_re = freq_slice->at(gc(prevk1i+k2yi, 0)) + freq_slice->at(gc(prevk1q+k2yq, 0));
					//a1yy_im = freq_slice->at(gc(prevk1i+k2yi, 1)) + freq_slice->at(gc(prevk1q+k2yq, 1));
					//a2yy_re = freq_slice->at(gc(prevk1i+k2yq, 0)) - freq_slice->at(gc(prevk2yi+k1q, 0));
					//a2yy_im = freq_slice->at(gc(prevk1i+k2yq, 1)) + freq_slice->at(gc(prevk2yi+k1q, 1));
					//a3yy_re = -1 * a2yy_im;
					//a3yy_im = a2yy_re;

					////writing to output matrix
					//bl = prevk + k2;
					//if (f == 0){
						//(data_out->at(2))[t][2*nFrequencies-1][gc(bl, 0)] = ( a1yx_re + a3yx_re);
						//(data_out->at(2))[t][2*nFrequencies-1][gc(bl, 1)] = -1*( a1yx_im + a3yx_im);
						//(data_out->at(3))[t][2*nFrequencies-1][gc(bl, 0)] = (a1yy_re + a3yy_re);
						//(data_out->at(3))[t][2*nFrequencies-1][gc(bl, 1)] = -1*(a1yy_im + a3yy_im);
					//}
					//(data_out->at(2))[t][nFrequencies-1+f][gc(bl, 0)] = (a1yx_re + a3yx_re);
					//(data_out->at(2))[t][nFrequencies-1+f][gc(bl, 1)] = -1*(a1yx_im + a3yx_im);
					//(data_out->at(2))[t][nFrequencies-1-f][gc(bl, 0)] = a1yx_re - a3yx_re;
					//(data_out->at(2))[t][nFrequencies-1-f][gc(bl, 1)] = a1yx_im - a3yx_im;
					//(data_out->at(3))[t][nFrequencies-1+f][gc(bl, 0)] = (a1yy_re + a3yy_re);
					//(data_out->at(3))[t][nFrequencies-1+f][gc(bl, 1)] = -1*(a1yy_im + a3yy_im);
					//(data_out->at(3))[t][nFrequencies-1-f][gc(bl, 0)] = a1yy_re - a3yy_re;
					//(data_out->at(3))[t][nFrequencies-1-f][gc(bl, 1)] = a1yy_im - a3yy_im;
				//}
			//}
		//}
	//}
	//return;
//}

int gc(int a, int b){
	return 2 * a + b;
}

string getFileName(string fileNameWithPath){//get the filename in a long path/path/filename
	string output = fileNameWithPath;
	size_t found;
	found = output.rfind( "/" );
	if (found!=string::npos and found!=output.size() - 1){
		output.erase( 0, found + 1 );
	};
	return output;
}

string strReplace(string input, string a, string b){
	string output = input;
	size_t found;
	found = output.find( a );
	while (found!=string::npos){
		output.replace( found, a.size(), b );
		found = output.find( a );
	};
	return output;
}

string extFileName(string fileName, string ext){//extend a file name by a string, such as extend way.cool.odf with shit to get way.coolshit.odf
	string output = fileName;
	size_t found;
	found = output.rfind( "." );
	if (found!=string::npos){
		output.insert( found, ext );
	};
	return output;
}

vector<string> parseLines(string bigLine){
	string line;
	vector<string> lines;
	stringstream ssls (bigLine);

	while (getline( ssls, line )){
		lines.push_back(line);
	}
	return lines;
}

float square(float x){
	return pow( max(min(x, MAX_POW_2), -MAX_POW_2), 2);
}


int get1DBL(int i, int j, int nAntenna){//0 indexed
	int output;
	if (i <= j) {
		output = ( ( 2 * nAntenna - 1 - i ) * i / 2 + j );
	} else {
		output = ( ( 2 * nAntenna - 1 - j ) * j / 2 + i );
	}
	return output;
}


vector<int> get2DBL(int bl, int nAntenna){//bl counts cross corrs AND auto corrs
	if(bl < nAntenna){
		vector<int> v(2);
		v[0] = 0;
		v[1] = bl;
		return v;
	} else{
		vector<int> v;
		v = get2DBL(bl-nAntenna, nAntenna-1);
		v[0] = v[0] + 1;
		v[1] = v[1] + 1;
		return v;
	}
	vector<int> v(2, -1);
	return v;
}

vector<int> get2DBLCross(int bl, int nAntenna){//bl only counts cross corrs
	if(bl < nAntenna - 1){
		vector<int> v(2);
		v[0] = 0;
		v[1] = bl + 1;
		return v;
	} else{
		vector<int> v;
		v = get2DBLCross(bl - nAntenna + 1, nAntenna - 1);
		v[0] = v[0] + 1;
		v[1] = v[1] + 1;
		return v;
	}
	vector<int> v(2, -1);
	return v;
}

bool contains(vector<vector<float> > * UBL, vector<float> bl){//automatically checks for the opposite direction
	for (unsigned int i = 0; i < UBL->size(); i++){
		if ( ( fabs((&(UBL->at(i)))->at(0) - bl[0]) < UBLPRECISION && fabs((&(UBL->at(i)))->at(1) - bl[1]) < UBLPRECISION ) or ( fabs((&(UBL->at(i)))->at(0) + bl[0]) < UBLPRECISION && fabs((&(UBL->at(i)))->at(1) + bl[1]) < UBLPRECISION ) ){
			return true;
		}
	}
	return false;
}

int indexUBL(vector<vector<float> > * UBL, vector<float> bl){//give the 1-indexed index of a baseline inside the unique baseline list; the opposite direction will give -index
	for (unsigned int i = 0; i < UBL->size(); i++){
		if ( fabs((&(UBL->at(i)))->at(0) - bl[0]) < UBLPRECISION && fabs((&(UBL->at(i)))->at(1) - bl[1]) < UBLPRECISION ) {
			return 1+i;
		} else if ( fabs((&(UBL->at(i)))->at(0) + bl[0]) < UBLPRECISION && fabs((&(UBL->at(i)))->at(1) + bl[1]) < UBLPRECISION ){
			return -1-i;
		}
	}
	return 0;
}


float amp(vector<float> * x){
	return sqrt( square(x->at(0))  + square(x->at(1)) );
}

float amp(float x, float y){
	return sqrt( square(x)  + square(y) );
}

float phase(float re, float im){
	/*if (re == 0 and im == 0){
		return 0;
	}*/
	return atan2(im, re);
}



float norm(vector<vector<float> > * v){
	float res = 0;
	for (unsigned int i = 0; i < v->size(); i++){
		for (unsigned int j = 0; j < v->at(i).size(); j++){
			res += pow(v->at(i)[j], 2);
		}
	}
	return pow(res, 0.5);
}

float phase(vector<float> * c){
	return atan2(c->at(1), c->at(0));
};

vector<float> conjugate (vector<float> x){
	vector<float> y = x;
	y[1] = -x[1];
	return y;
}

bool contains_int(vector<int> * list, int j){
	for (uint i = 0; i < list->size(); i ++){
		if ( list->at(i) == j ){
			return true;
		}
	}
	return false;
}

void addPhase(vector<float> * x, float phi){
	float am = amp(x);
	float ph = phase(x->at(0), x->at(1));
	ph = phaseWrap(ph + phi);
	x->at(0) = am * cos(ph);
	x->at(1) = am * sin(ph);
	return;
}


float vectorDot(vector<float>* v1, vector<float>* v2){//v1.v2
	string METHODNAME = "vectorDot";
	if ( v1->size() != v2->size() ){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL INPUT MISMATCH! lengths of vectors are " << v1->size() << " and " << v2->size() << ". 0 returned!!!" << endl;
		return 0;
	}
	double sum = 0;
	for (unsigned int i = 0; i < v1->size(); i ++){
		sum += (v1->at(i) * v2->at(i));
	}
	return float(sum);
}


vector<float> matrixDotV(vector<vector<float> >* m, vector<float>* v){//m.v
	string METHODNAME = "maxtrixDotV";
	vector<float> u(m->size(), 0);
	if ( (m->at(0)).size() != v->size() ){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL INPUT MISMATCH! Dimensions of matrix and of vector are " << m->size() << "x" << (m->at(0)).size() << " and " << v->size() << ". 0 vector  returned!!!" << endl;
		return u;
	}
	for (unsigned int i = 0; i < m->size(); i ++){
		u[i] = vectorDot( &(m->at(i)), v);
	}
	return u;
}

vector<vector<float> > rotationMatrix(float x, float y, float z){//approximation for a rotation matrix rotating around x,y,z axis, {{1, -z, y}, {z, 1, -x}, {-y, x, 1}}
	vector<vector<float> > r(3, vector<float>(3,1));
	r[0][1] =-z;
	r[0][2] = y;
	r[1][0] = z;
	r[1][2] =-x;
	r[2][0] =-y;
	r[2][1] = x;
	return r;
}

vector<vector<float> > rotationMatrixZ(float z){/*approximation for a rotation matrix rotating around x,y,z axis, Ry[roty_] := {{Cos[roty], 0, Sin[roty]}, {0, 1, 0}, {-Sin[roty], 0,
    Cos[roty]}};
Rx[roty_] := {{1, 0, 0}, {0, Cos[roty], -Sin[roty]}, {0, Sin[roty],
    Cos[roty]}};
Rz[rotz_] := {{Cos[rotz], -Sin[rotz], 0}, {Sin[rotz], Cos[rotz],
    0}, {0, 0, 1}};*/
	vector<vector<float> > r(3, vector<float>(3, 0));
	r[0][0] = cos(z);
	r[0][1] =-sin(z);
	r[1][0] = sin(z);
	r[1][1] = cos(z);
	r[2][2] = 1;

	return r;
}


void mergeCalibrationPar (vector<float> * ampcalpar1, vector<float> * ampcalpar2, vector<float> * ampcalparM, vector<float> * phasecalpar1, vector<float> * phasecalpar2, vector<float> * phasecalparM)//Only deals with ampcalpar and phasecalpar, Chisq should be the chisq of the second set, so are the UBLs
{
	for (uint i = 0; i < ampcalparM->size(); i++){
		ampcalparM->at(i) = ampcalpar1->at(i) + ampcalpar2->at(i);
		phasecalparM->at(i) = phaseWrap( phasecalpar1->at(i) + phasecalpar2->at(i) );
	}
}

float getCableCorrection(float length, float freq){//measured phase - phaseCor = ideal phase
	return (length * freq / SPEEDC) * 2 * PI;
}

float getFreq(int f, int nFrequency, float startFreq, float endFreq){//To be consistent with OmniViewer frequency interpretation!!!
	return startFreq + (endFreq - startFreq) * f / ( nFrequency - 1 );
}

vector<float> getBL(int i, int j, vector<vector<float> > *antloc){
	vector<float> bl(2,0);
	bl[0] = (&(antloc->at(j)))->at(0) - (&(antloc->at(i)))->at(0);
	bl[1] = (&(antloc->at(j)))->at(1) - (&(antloc->at(i)))->at(1);
	return bl;
}

int countUBL(vector<vector<float> > *antloc ){
	vector<float> bl;
	vector<vector<float> > UBL;
	for (unsigned int i = 0; i < antloc->size() - 1; i++){
		for (unsigned int j = i + 1; j < antloc->size(); j++){
			bl = getBL(i, j, antloc);
			if (!contains(&UBL, bl)) {
				UBL.push_back(bl);
			}
		}
	}
	return UBL.size();
}

int lookupAnt(float x, float y, vector<vector<float> > antloc){
	for (unsigned int i = 0; i < antloc.size(); i++){
		if ( x == antloc[i][0] && y == antloc[i][1]){
			return i;
		}
	}
	return -1;
}

float phaseWrap (float x, float offset/*default -pi*/){
	while ( x <= offset ){
		x = x + 2 * PI;
	}
	while ( x > offset + 2 * PI ){
		x = x - 2 * PI;
	}
	return x;
}

void phaseCalibrate120(vector<float>* calpar120, vector<float>* calpar16, uint nAnt, vector<bool>* badAnt){//find the median solution of 16 antenna calpars from 120 visibility calpars
	string METHODNAME = "phaseCalibrate";
	if (calpar120->size() != nAnt * (nAnt - 1) / 2) {
		cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: input calpar120 length is " << calpar120->size() << ", not the expected length of " << nAnt * (nAnt - 1) / 2 << "! Abort!" << endl;
		return;
	}
	if (calpar16->size() != nAnt) {
		cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: input calpar16 length is " << calpar16->size() << ", not the expected length of nAnt " << nAnt << "! Abort!" << endl;
		return;
	}
	vector<vector<float> > calpar(nAnt, vector<float>(nAnt, 0));
	int counter = 0;
	for (uint i = 0; i < nAnt; i++){
		for (uint j = i + 1; j < nAnt; j++){
			calpar[i][j] = calpar120->at(counter);
			counter++;
		}
	}
	for (uint i = 0; i < nAnt; i++){
		if(!(badAnt->at(i))){
			for (uint j = 0; j < i; j++){
				calpar[i][j] = -calpar[j][i];
			}
		}
	}
/*
	for (int i = 0; i < nAnt; i++){
		for (int j = 0; j < nAnt; j++){//reversed order such that the first element won't get zeroed out before substracted from the rest
			cout << calpar[i][j] << " ";//DEBUG
		}
		cout << endl;	//DEBUG
	}
	cout << endl << "**********" << endl;	//DEBUG
*/
	for (uint i = 0; i < nAnt; i++){
		if(!(badAnt->at(i))){
			for (int j = nAnt - 1; j > -1; j = j - 1){//reversed order such that the first element won't get zeroed out before substracted from the rest
				calpar[i][j] = phaseWrap(calpar[i][j] - calpar[i][0]);
			}
		}

	}

	for (int i = nAnt -1 ; i > -1; i = i - 1){
		if(badAnt->at(i)){
			calpar.erase(calpar.begin() + i);
		}

	}

/*
	for (int i = 0; i < nAnt; i++){
		for (int j = 0; j < nAnt; j++){//reversed order such that the first element won't get zeroed out before substracted from the rest
			cout << calpar[i][j] << " ";//DEBUG
		}
		cout << endl;	//DEBUG
	}
	cout << endl << "**********" << endl;	//DEBUG
*/
	vector<vector<float> > calparT (calpar[0].size(), vector<float>(calpar.size(), 0));
	for (unsigned int i = 0; i < calparT.size(); i++){
		for (unsigned int j = 0; j < calparT[0].size(); j++){
			calparT[i][j] = calpar[j][i];
		}
	}
	for (uint i = 0; i < nAnt; i++){
		calpar16->at(i) = medianAngle(&(calparT[i]));
//		cout << calpar16->at(i) << " ";//DEBUG
	}
//	cout << endl << "*************************" << endl;	//DEBUG
	return;
}

vector<float> phaseCalibrate(vector<vector<float> > *dataf, string pol, float freq, vector<vector<float> > *antloc, vector<vector<float> > *cablelen, int startingAnt1, int startingAnt2, int startingAnt3, uint nAntenna){
	string METHODNAME = "phaseCalibrate";
	vector<float> output( nAntenna, 0 );
	vector<bool> calibrated( nAntenna, false);//denote whether each antenna has been calibrated
	vector<int> calibratedList;//a growing list of calibrated antenna number
	vector<float> bl(2);
	vector<vector<float> > knownUBL;
	vector<float> knownUBLphase;
	vector<float> cableleni( nAntenna );
	vector<float> cablelenj( nAntenna );
	if ( pol == "xx"){
		cableleni = cablelen->at(0);
		cablelenj = cablelen->at(0);
	}else if ( pol == "xy"){
		cableleni = cablelen->at(0);
		cablelenj = cablelen->at(1);
		cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: xy polarization not yet supported! Abort!" << endl;
		return output;
	}else if ( pol == "yx"){
		cableleni = cablelen->at(1);
		cablelenj = cablelen->at(0);
		cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: yx polarization not yet supported! Abort!" << endl;
		return output;
	}else if ( pol == "yy"){
		cableleni = cablelen->at(1);
		cablelenj = cablelen->at(1);
	}


	float phase12 = phase((dataf->at(get1DBL(startingAnt1, startingAnt2, nAntenna)))[0], (dataf->at(get1DBL(startingAnt1, startingAnt2, nAntenna)))[1]);
	float phaseCor2;//measured phase - phaseCor = ideal phase
	if (startingAnt1 < startingAnt2){
		 phaseCor2 = getCableCorrection( (cablelenj[startingAnt2] - cableleni[startingAnt1]) , freq);
		 phase12 = phaseWrap(phase12 - phaseCor2);//<x*, x> as expected by adrian's fortran codes
	} else {
		 phaseCor2 = getCableCorrection( (cableleni[startingAnt2] - cablelenj[startingAnt1]) , freq);
		 phase12 = phaseWrap( - ( phase12 + phaseCor2 ) );
	}
	float phase13 = phase((dataf->at(get1DBL(startingAnt1, startingAnt3, nAntenna)))[0], (dataf->at(get1DBL(startingAnt1, startingAnt3, nAntenna)))[1]);
	float phaseCor3;
	if (startingAnt1 < startingAnt3){
		 phaseCor3 = getCableCorrection( (cablelenj[startingAnt3] - cableleni[startingAnt1]) , freq);
		 phase13 = phaseWrap( phase13 - phaseCor3 );
	} else {
		 phaseCor3 = getCableCorrection( (cableleni[startingAnt3] - cablelenj[startingAnt1]) , freq);
		 phase13 = phaseWrap( - ( phase13 + phaseCor3 ) );
	}
	output[startingAnt2] = phaseCor2;
	output[startingAnt3] = phaseCor3;
	calibrated[startingAnt1] = true;
	calibrated[startingAnt2] = true;
	calibrated[startingAnt3] = true;
	calibratedList.push_back(startingAnt1);
	calibratedList.push_back(startingAnt2);
	calibratedList.push_back(startingAnt3);
	knownUBL.push_back(getBL(startingAnt1, startingAnt2, antloc));
	knownUBLphase.push_back(phase12);
	knownUBL.push_back(getBL(startingAnt1, startingAnt3, antloc));
	knownUBLphase.push_back(phase13);
	knownUBL.push_back(getBL(startingAnt2, startingAnt3, antloc));
	knownUBLphase.push_back( phaseWrap( phase13 - phase12 ));
	knownUBL.push_back(getBL(startingAnt2, startingAnt1, antloc));
	knownUBLphase.push_back(-phase12);
	knownUBL.push_back(getBL(startingAnt3, startingAnt1, antloc));
	knownUBLphase.push_back(-phase13);
	knownUBL.push_back(getBL(startingAnt3, startingAnt2, antloc));
	knownUBLphase.push_back( phaseWrap( - phase13 + phase12 ));
	int numCorrected = 1;
	while ( numCorrected > 0 ){
		numCorrected = 0;
		for (uint i = 0; i < nAntenna; i++ ){
			for (uint jj = 0; jj < calibratedList.size(); jj++){
				if (!calibrated[i]){
					uint j = calibratedList[jj];
					vector<float> ijBL = getBL(i, j, antloc);
					float ijphase = phase((dataf->at(get1DBL(i, j, nAntenna)))[0], (dataf->at(get1DBL(i, j, nAntenna)))[1]);
					for (uint k = 0; k < knownUBL.size(); k++){
						if (ijBL[0] == knownUBL[k][0] && ijBL[1] == knownUBL[k][1]){
							if ( i < j ){
								output[i] = phaseWrap( knownUBLphase[k] - ijphase + output[j]);
							} else if ( i > j ){
								output[i] = phaseWrap( knownUBLphase[k] + ijphase + output[j]);
							} else{
								cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL CONFUSION IN BASELINE SEEKING PROCESS! i = j = " << i << endl;
							}
							if(DEBUG){
								cout << "##" << FILENAME << "##" << METHODNAME << ": Calibrating antenna #" << i << " with known antenna #" << j << " whose correction is " << output[j] << ". Known phase difference is " << knownUBLphase[k] << " and in-data phase difference is " << ijphase << ". Thus phase correction parameter for antenna #" << i << " is " << output[i] << "." << endl;
							}
							calibrated[i] = true;
							calibratedList.push_back(i);
							numCorrected ++;
							break;
						}
					}
				}
			}
		}
	}
	if (calibratedList.size() < nAntenna) {
		cout << "##" << FILENAME << "##" << METHODNAME << " WARNING!!!: Only calibrated " << calibratedList.size() << " antennas at frequency " << freq << "MHz!!! " << endl;
	}
	/*for (int i = 0; i < nAntenna; i ++ ){
		cout << output[i] << " ";
	}
	cout << endl;*/
	return output;
}

void computeUBL(vector<vector<float> > * antloc, vector<vector<float> > * listUBL){
	int numAntenna = antloc->size();
	vector<float> baseline(2);
	vector<vector<float> > UBLtmp;
	for (int i = 0; i < numAntenna; i++){
		for (int j = i + 1; j < numAntenna; j++){
			baseline = getBL(i, j, antloc);
			if ( ! (contains(&UBLtmp, baseline)) ) {
				UBLtmp.push_back(baseline);
			}
		}
	}
	for (unsigned int i = 0; i < UBLtmp.size(); i++){
		listUBL->at(i) = UBLtmp[i];
	}
	return;
}

vector<float> modelToMeasurement(vector<float> *modelCor, float ampcal, float phasecal){
	string METHODNAME = "modelToMeasurement";
	vector<float> measurement (2, 0.0);
	float modelAmp = sqrt( modelCor->at(0) * modelCor->at(0) + modelCor->at(1) * modelCor->at(1) );
	float modelPhase = phase( modelCor->at(0), modelCor->at(1) );
	float measurementAmp = modelAmp * pow( 10, min(MAX_10_POW,ampcal) );
	float measurementPhase = phaseWrap( modelPhase + phasecal );

	//cout <<  "#!#" << FILENAME << "#!#" << METHODNAME << ": " << "modelAmp: " << modelAmp << " modelPhase: " << modelPhase << " measurementAmp: " << measurementAmp << " measurementPhase: " << measurementPhase << " ampcal: " << ampcal << " phasecal: " << phasecal << endl;
	measurement[0] = measurementAmp * cos(measurementPhase);
	measurement[1] = measurementAmp * sin(measurementPhase);
	return measurement;
}

vector<float> correctMeasurement(vector<float> *measuredCor, float ampcal1, float phasecal1, float ampcal2, float phasecal2){//assumes <x1* , x2>
	vector<float> corrected (2, 0.0);
	float a = amp(measuredCor) / max( MIN_NONE_ZERO, (float)pow( 10, min(MAX_10_POW, (ampcal1 + ampcal2)) ) );
	float p = phaseWrap( phase(measuredCor->at(0), measuredCor->at(1)) + phasecal1 - phasecal2 );
	corrected[0] = a * cos(p);
	corrected[1] = a * sin(p);
	return corrected;
};

void correctMeasurementPhase(vector<float> *measuredCor, float phasecal1, float phasecal2){//assumes <x1* , x2>
	float a = amp(measuredCor);
	float p = phaseWrap( phase(measuredCor->at(0), measuredCor->at(1)) + phasecal1 - phasecal2 );
	measuredCor->at(0) = a * cos(p);
	measuredCor->at(1) = a * sin(p);
	return;
};

void correctMeasurementMatrix(vector<vector<float> > *correctedData, vector<float> *ampcalpar, vector<float> *phasecalpar){
	int nAntenna = ampcalpar->size();
	for ( int i = 0; i < nAntenna; i++){
		for ( int j = i; j < nAntenna; j++){
			correctedData->at(get1DBL(i, j, nAntenna)) = correctMeasurement(&(correctedData->at(get1DBL(i, j, nAntenna))), ampcalpar->at(i), phasecalpar->at(i), ampcalpar->at(j), phasecalpar->at(j));
		}
	}
}

void computeUBLcor(vector<vector<float> >* calibratedData, vector<int> *UBLindex, vector<vector<float> > *UBLcor, vector<bool> *goodAnt){//average each group of calibrated redundant baselines to get the estimate for that ubl, only useful when directly applying a set of calpars instead of using logcal to fit for them.
}

vector<float> getModel(int i, int j, vector<vector<float> > *antloc, vector<vector<float> > *listUBL, vector<vector<float> > *UBLcor){
	vector<float> baseline = getBL(i, j, antloc);
	for (uint k = 0; k < listUBL->size(); k++){
		if ( baseline[0] == (&(listUBL->at(k)))->at(0) && baseline[1] == (&(listUBL->at(k)))->at(1) ){
			return UBLcor->at(k);
		} else	if ( baseline[0] == -(&(listUBL->at(k)))->at(0) && baseline[1] == -(&(listUBL->at(k)))->at(1) ){
			return conjugate( UBLcor->at(k) );
		}
	}
	return vector<float>(UBLcor->at(0).size(), 0);
}

vector<vector<float> > ReverseEngineer(vector<float> * ampcalpar, vector<float> * phasecalpar, vector<vector<float> > * UBLcor, vector<vector<float> > * antloc, vector<vector<float> > * listUBL){
	string METHODNAME = "ReverseEngineer";
	int numAntenna = antloc->size();
	int numCrosscor = numAntenna * ( numAntenna - 1 ) / 2;
	vector<float> complexdummy (2, 0.0);//Just an awkward way to initialize output
	vector<vector<float> > output (numCrosscor, complexdummy);

	int cnter = 0;
	for (int i = 0; i < numAntenna; i++){
		for (int j = i + 1; j < numAntenna; j++){
			vector<float> cor(2, 0.0);
			cor = getModel(i, j, antloc, listUBL, UBLcor);
			output[cnter] = modelToMeasurement( &cor, (ampcalpar->at(i) + ampcalpar->at(j)), (phasecalpar->at(j) - phasecalpar->at(i)) );
			cnter ++;
		}
	}
	return output;
}

void ReverseEngineer(vector<vector<float> >* output, vector<float> * calpar, int numAntenna, vector<int> * UBLindex){
	string METHODNAME = "ReverseEngineer";
	//int numCrosscor = numAntenna * ( numAntenna - 1 ) / 2;
	int cnter = 0;
	vector<float> cor(2, 0.0);
	for (int i = 0; i < numAntenna; i++){
		for (int j = i + 1; j < numAntenna; j++){
			int ubl = fabs(UBLindex->at(cnter)) - 1;
			if(UBLindex->at(cnter) > 0){
				cor[0] = calpar->at(3 + 2 * numAntenna + 2 * ubl);
				cor[1] = calpar->at(3 + 2 * numAntenna + 2 * ubl + 1);
			} else{
				cor[0] = calpar->at(3 + 2 * numAntenna + 2 * ubl);
				cor[1] = -calpar->at(3 + 2 * numAntenna + 2 * ubl + 1);
			}
			output->at(cnter) = modelToMeasurement( &cor, (calpar->at(3 + i) + calpar->at(3 + j)), (calpar->at(3 + numAntenna + j) - calpar->at(3 + numAntenna + i)) );
			cnter ++;
		}
	}
	return;
}

float chiSq(vector<vector<float> > * dataf, vector<vector<float> > * sdevf, vector<vector<float> > * antloc, vector<float> * ampcalpar, vector<float> * phasecalpar, vector<vector<float> > * UBLcor, int numAntenna, vector<vector<float> > * listUBL){
	string METHODNAME = "chiSq";
	uint numCrosscor = numAntenna * ( numAntenna - 1 ) / 2;
	uint numAutocor = numAntenna * ( numAntenna + 1 ) / 2;
	if ( dataf->size() != numAutocor) {
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << dataf->size() << ", not consistent with expected " << numAutocor << " from specified " << numAntenna << " antenna!" << endl;
		return 0.0;
	}
	if ( dataf->size() != sdevf->size()) {
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << dataf->size() << ", not consistent with length of standard deviation input " << sdevf->size() << "!" << endl;
		return 0.0;
	}
	float output = 0.0;
	vector<float> bl(2);

	vector<vector<float> > Ax(numCrosscor, bl);
	vector<vector<float> > y(numCrosscor, bl);
	vector<vector<float> > N(numCrosscor, bl);// Here N is technicaly n, noise in units of sdev, not the N covariant matrix

	int cnter = 0;
	for ( int i = 0; i < numAntenna; i++){
		for ( int j = i + 1; j < numAntenna; j++){
			int index = get1DBL(i, j, numAntenna);
			y[cnter] = dataf->at(index);
			N[cnter] = sdevf->at(index);
			cnter ++;
		}
	}

	Ax = ReverseEngineer(ampcalpar, phasecalpar, UBLcor, antloc, listUBL);

	for (uint i = 0; i < numCrosscor; i ++){
		output = output + square(( Ax[i][0] - y[i][0] ) / max(N[i][0], MIN_NONE_ZERO)) + square( ( Ax[i][1] - y[i][1] ) / max(N[i][1], MIN_NONE_ZERO));
	}

	return output;
}

bool fillChiSq(vector<vector<float> >* dataf, vector<vector<float> >* sdevf, vector<float>* calpar, int numAntenna, vector<int>* UBLindex, vector<bool>* goodAnt){
	string METHODNAME = "fillChiSq";
	uint numCrosscor = numAntenna * ( numAntenna - 1 ) / 2;
	uint numAutocor = numAntenna * ( numAntenna + 1 ) / 2;
	if ( dataf->size() != numAutocor) {
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << dataf->size() << ", not consistent with expected " << numAutocor << " from specified " << numAntenna << " antenna!" << endl;
		return false;
	}
	if ( dataf->size() != sdevf->size()) {
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << dataf->size() << ", not consistent with length of standard deviation input " << sdevf->size() << "!" << endl;
		return false;
	}
	float output = 0.0;
	vector<float> bl(2);

	vector<vector<float> > Ax(numCrosscor, bl);
	ReverseEngineer(&Ax, calpar, numAntenna, UBLindex);

	for (uint i = 0; i < numCrosscor; i ++){
		vector<int> a = get2DBLCross(i, numAntenna);
		if(goodAnt->at(a[0]) and goodAnt->at(a[1])){
			int blauto = get1DBL(a[0], a[1], numAntenna);
			//cout << a[0] << " " << a[1] << " (" << Ax[i][0] << "," << Ax[i][1] << ")" << "(" << dataf->at(blauto)[0] << "," << dataf->at(blauto)[1] << ") " << square(( Ax[i][0] - dataf->at(blauto)[0] ) / max(sdevf->at(blauto)[0], MIN_NONE_ZERO)) + square( ( Ax[i][1] - dataf->at(blauto)[1] ) / max(sdevf->at(blauto)[1], MIN_NONE_ZERO)) << endl;
			output = output + square(( Ax[i][0] - dataf->at(blauto)[0] ) / max(sdevf->at(blauto)[0], MIN_NONE_ZERO)) + square( ( Ax[i][1] - dataf->at(blauto)[1] ) / max(sdevf->at(blauto)[1], MIN_NONE_ZERO));
		}
	}
	calpar->at(0) = output;
	return true;
}




float median (vector<float> list){
	int l = list.size();
	if (l == 0) return 0;
	sort(list.begin(), list.end());
	int index = floor( l / 2 );
	//cout << list[0] << " " << list[list.size() - 1] << endl;
	/*cout << "DBG: Input list: ";
	for ( int i = 0; i < list.size(); i++){
		cout << list[i] << " ";
	}
	cout << endl << "median: " << list[index] << endl;*/
	if(l % 2 == 1) return list[index];
	else return (list[index] + list[index - 1]) / 2;
}

float medianAngle (vector<float> *list){
	string METHODNAME = "medianAngle";
	//cout << "#!#" << FILENAME << "#!#" << METHODNAME << " DBG ";
	vector<float> xList(list->size());
	vector<float> yList(list->size());
	for (unsigned int i = 0; i < list->size(); i++){
		//cout << list[i] << " ";
		xList[i] = cos(list->at(i));
		yList[i] = sin(list->at(i));
	}
	//cout << " median is " << atan2(median(yList), median(xList)) << endl;
	return atan2(median(yList), median(xList));
}

float mean (vector<float> *v, int start, int end){// take mean from start to end indices of vector v. 0 indexed
	string METHODNAME = "mean";
	int size = (int) (v->size());//size is originally unsigned

	if (size <= 0){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << " !!WARNING!! mean of an empty array requested!!";
		return 0;
	}


	if (end > size - 1 or start > size - 1){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << " !!WARNING!! start/end index requested at "<< start << "/" << end << " is out of array length of " << v->size() << "!!";
	}
	int a,b;
	if (start < 0 or start > size - 1) a = 0; else a = start;
	if (end < 0 or end > size - 1) b = size - 1; else b = end;
	float sum = accumulate(v->begin() + a, v->begin() + b, 0.0);
	//cout <<  start << " " << end << " " << a << " " << b << " " << sum << " " << endl;
	//cout << sum << endl;
	return sum / (b - a + 1);
}

vector<float> stdev(vector<float> *v){//returns {mean, sample standard deviation}. Created by Hrant and modified by Jeff
	string METHODNAME = "stdev";
	vector <float> result(2,0);
	int n = v-> size();
	if ( n <= 1){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << " !!WARNING!! standard deviation of an empty or unit array requested!!";
		return result;
	}
	float m = mean(v);
	//cout << m << endl;
	result[0] = m;
	float var = 0;
	for (int i = 0; i < n; i ++){
		var += square( v->at(i) - m );
	}
	/*standard deviation*/result[1] = sqrt( var / (n - 1) );
	return result;
}

float meanAngle (vector<float> *list){
	string METHODNAME = "meanAngle";
	//cout << "#!#" << FILENAME << "#!#" << METHODNAME << " DBG ";
	vector<float> xList(list->size());
	vector<float> yList(list->size());
	for (unsigned int i = 0; i < list->size(); i++){
		//cout << list[i] << " ";
		xList[i] = cos(list->at(i));
		yList[i] = sin(list->at(i));
	}
	//cout << " median is " << atan2(median(yList), median(xList)) << endl;
	return atan2(mean(&yList), mean(&xList));
}

vector<float> stdevAngle(vector<float> *v){
	string METHODNAME = "stdevAngle";
	vector <float> result(2,0);
	int n = v-> size();
	if ( n <= 1){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << " !!WARNING!! standard deviation of an empty or unit array requested!!";
		return result;
	}
	float m = meanAngle(v);
	//cout << m << endl;
	result[0] = m;
	float var = 0;
	for (int i = 0; i < n; i ++){
		var += square( phaseWrap(v->at(i) - m) );
	}
	/*standard deviation*/result[1] = sqrt( var / n );
	return result;
}

///////////////MAJOR STUFF///////////////////
/////////////////////////////////////////////
void pointSourceCalAccordingTo(uint referenceAnt, vector<vector<float> > *data, vector<vector<float> > *ampcalparArray, vector<vector<float> > *phasecalparArray){//referenceAnt is 0 indexed, ampcalparArray has first dimension represent each antenna's calibration parameter, the second dimention represent the same parameter computed using different reference antennas. This method fills up one column of ampcalparArray (which is essentially a square matrix) at a time
	string METHODNAME = "pointSourceCalAccordingTo";
	uint nAnt = ampcalparArray->size();
	if ( data->size() != ( nAnt * (nAnt + 1) / 2 ) ){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << data->size() << ", not consistent with expected number of crosscorrelations " << nAnt * (nAnt + 1) / 2 << " computed from " << nAnt << " antennas!" << endl;
		return;
	}

	//Extracting amplitudes off all cross-correlations related to the reference ant
	vector<float> amps (nAnt - 1);
	for ( uint i = 0; i < referenceAnt; i ++){
		amps[i] = amp( &(data->at(get1DBL(i, referenceAnt, nAnt))) );
	}
	for ( uint i = referenceAnt + 1; i < nAnt; i ++){
		amps[i - 1] = amp( &(data->at(get1DBL(i, referenceAnt, nAnt))) );
	}
	float standardAmp = median(amps);
	standardAmp = max(standardAmp , MIN_NONE_ZERO);

	float refphase = phase( (data->at(referenceAnt))[0], -(data->at(referenceAnt))[1] );//reference phase, < x_ref*, x_0>. If we were not to demand x_0's phase correction to be 0, this would have been x_0's phase correction for < x_ref*, x_0> to have 0 phase. However, we substract this from all subsequent phase corrections to demand phase correction for x_0 is 0.
	for (uint i = 0; i < nAnt; i++){
		(&(ampcalparArray->at(i)))->at(referenceAnt) = log10( max( amp(&(data->at(get1DBL(referenceAnt, i, nAnt)))) / standardAmp, MIN_NONE_ZERO ) );
		//cout << amp( &(data->at(get1DBL(referenceAnt, i, nAnt))) ) << " " << standardAmp << endl;

		float rawphase = phase( (data->at(get1DBL(referenceAnt, i, nAnt)))[0], (data->at(get1DBL(referenceAnt, i, nAnt)))[1] );
		if ( i < referenceAnt ){
			rawphase = - rawphase;
		}
		(&(phasecalparArray->at(i)))->at(referenceAnt) = phaseWrap(rawphase - refphase);
		//cout << "#!#" << METHODNAME << " DBG: refAnt: #" << referenceAnt << " currentAnt: #" << i << " refphase:" << refphase << " rawphase: " << rawphase << " raw data: " << (data->at(get1DBL(referenceAnt, i, nAnt)))[0] << " " << (data->at(get1DBL(referenceAnt, i, nAnt)))[1] << " calpar: " << phaseWrap(rawphase - refphase) << endl;
	}
}


void pointSourceCal(vector<vector<float> > *data, vector<float> *ampcalpar, vector<float> *phasecalpar, vector<vector<float> > *UBLcalpar){
	string METHODNAME = "pointSourceCal";
	uint nAnt = ampcalpar->size();
	if ( data->size() != ( nAnt * (nAnt + 1) / 2 ) ){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of data is " << data->size() << ", not consistent with expected number of crosscorrelations " << nAnt * (nAnt + 1) / 2 << " computed from " << nAnt << " antennas!" << endl;
		return;
	}

	vector<float> dummy (nAnt);
	vector<vector<float> > ampcalparArray (nAnt, dummy);//contains an aray of ampcalpar parameters from different ampcalpars based on different reference antennas
	vector<vector<float> > phasecalparArray (nAnt, dummy);

	for (uint i = 0; i < nAnt; i++){
		pointSourceCalAccordingTo(i, data, &ampcalparArray, &phasecalparArray);
	}
	for (uint i = 0; i < nAnt; i++){
		ampcalpar->at(i) = median(ampcalparArray[i]);
		phasecalpar->at(i) = medianAngle(&(phasecalparArray[i]));
		//if ( i == 14 ){
		//	for ( int j = 0; j < nAnt; j ++){
		//		cout << ampcalparArray[i][j] << " " << endl;
		//	}
		//}
	}
	vector<float> autocor(nAnt);
	for (uint i = 0; i < nAnt; i++){
		autocor[i] = (data->at(get1DBL(i, i, nAnt)))[0];
	}
	float autoMedian = median(autocor);
	for (unsigned int i = 0; i < UBLcalpar->size(); i++){
		(UBLcalpar->at(i))[0] = autoMedian;
		(UBLcalpar->at(i))[1] = 0;
	}
}

void substractComplexPhase(vector<float> *a, vector<float> *b, float angle){
	float amptmp = amp(a);
	float phasetmp = phase(a->at(0), a->at(1));
	phasetmp = phaseWrap( phasetmp - angle );
	b->at(0) = amptmp * cos(phasetmp);
	b->at(1) = amptmp * sin(phasetmp);
	return;
}

void rotateCalpar(vector<float> *originalPhase, vector<float> *rotatedPhase, vector<vector<float> > *originalUBLcor, vector<vector<float> > *rotatedUBLcor, vector<vector<float> > *antloc, vector<vector<float> > * listUBL, float theta, float phi, float freq){//theta in [0, PI/2] (rotating z axis down), phi [0, 2PI), freq in MHz,  ONLY WORKS FOR ROTATING POINT SOURCE STRAIGHT ABOVE *TO* THETA AND PHI
	string METHODNAME = "rotatePhasecalpar";
	float k = 2 * PI / SPEEDC * freq;
	float sint = sin(theta);
	////float cost = cos(theta);
	float sinp = sin(phi);
	float cosp = cos(phi);

	//vector<float> rotation (originalPhase->size(), 0.0);
	if ( originalPhase->size() != rotatedPhase->size() or originalPhase->size() != antloc->size()){
		cout << "#!#" << FILENAME << "#!#" << METHODNAME << ": !!!!FATAL ERROR!!!! Length of original phasecalpar is " << originalPhase->size() << ", not consistent with length of rotated phasecalpar " << rotatedPhase->size() << " or that implied by antloc " << antloc->size() << endl;
		return;
	}

	for (uint i = 0; i < originalPhase->size(); i ++){
		//rotation[i] = k * (sint * cosp * ((&(antloc->at(i)))->at(0)) + sint * sinp * ((&(antloc->at(i)))->at(1)));
		rotatedPhase->at(i) = phaseWrap( originalPhase->at(i) + k * (sint * sinp * ((&(antloc->at(i)))->at(0)) - sint * cosp * ((&(antloc->at(i)))->at(1))) );
	}

	for (uint i = 0; i < originalUBLcor->size(); i ++){
		//rotation[i] = k * (sint * cosp * ((&(antloc->at(i)))->at(0)) + sint * sinp * ((&(antloc->at(i)))->at(1)));
		substractComplexPhase(
			&(originalUBLcor->at(i)),
			&(rotatedUBLcor->at(i)),
			phaseWrap(
				k * (sint * sinp * ((&(listUBL->at(i)))->at(0)) - sint * cosp * ((&(listUBL->at(i)))->at(1)))
				- phase((originalUBLcor->at(i))[0], (originalUBLcor->at(i))[1])
			)
		);
	}
}


//Logcal functions

bool invert(vector<vector<int> > * AtNinvAori, vector<vector<double> > * AtNinvAinv ){//Gauss–Jordan elimination
	string METHODNAME = "invert";
	//WARNING: PASSING OF ARGUMENT VALUES NOT TESTED !!!!!!!!!!!!

	/* 	subroutine invert(np,n,A,B)
	! B = A^{-1}, where A is an arbitraty non-singular matrix.
	! In the process, A gets destroyed (A=I on exit).
	! Note that A need NOT be symmetric.
	! This routine is designed for small matrices where you dont care about
	! speed or numerical stability.
	! Written & tested by Max Tegmark 000506
	*/

	uint i=0, j=0, k=0, n = AtNinvAori->size();//todo check size
	vector<vector<float> > AtNinvAqaz(AtNinvAori->size(), vector<float>(AtNinvAori->at(0).size(),0));
	vector<vector<float> > *AtNinvA = &AtNinvAqaz;
	for( i = 0; i < AtNinvAori->size(); i++){
		for( j = 0; j < AtNinvAori->at(0).size(); j++){
			AtNinvA->at(i)[j] = AtNinvAori->at(i)[j];
		}
	}
	double r;
	//clock_t t1,t2;
	//t1=clock();


	for( i = 0; i < n; i++){// set matrix to one
		for (j = 0; j < n; j++){
			(AtNinvAinv->at(i))[j]=0.0;
		}
		(AtNinvAinv->at(i))[i]=1.0;
	}
	//t2=clock();
	//if(TIME) cout << ((float)t2-(float)t1) / CLOCKS_PER_SEC << "sec ";
	for ( i = 0; i < n; i++){ // find inverse by making A(i,i)=1

		if(isnan(AtNinvA->at(i)[i])){
			printf("%s: FATAL ERROR: input matrix AtNinvA of size %i by %i has NaN on diagonal! ABORT!\n", METHODNAME.c_str(), n, n);
			return false;
		}

		if(fabs((AtNinvA->at(i))[i]) < MIN_NONE_ZERO){
			if (i == n - 1){
				return false;
			}
			for(j = i+1; j < n; j++){//find a row to add to row i to make its diagonal non zero
				if(fabs((AtNinvA->at(j))[i]) >= MIN_NONE_ZERO and !isnan((AtNinvA->at(j))[i])){
					for(k = 0; k < n; k++){
						(AtNinvA->at(i))[k]=(AtNinvA->at(i))[k] + (AtNinvA->at(j))[k];
						(AtNinvAinv->at(i))[k]= (AtNinvAinv->at(i))[k] + (AtNinvAinv->at(j))[k];
					}
					break;
				}
				if (j == n - 1){
					return false;
				}
			}
		}
		r = (AtNinvA->at(i))[i];
		for (j = 0; j < n; j++){
			(AtNinvA->at(i))[j]=(AtNinvA->at(i))[j]/r;
			(AtNinvAinv->at(i))[j]= (AtNinvAinv->at(i))[j]/r;
			//if(j == i and (AtNinvAinv->at(i))[j]<0) {
				//printf("%s: FATAL ERROR: inverse matrix AtNinvAinv of size %i by %i has negative number on diagonal! ABORT!\n", METHODNAME.c_str(), n, n);
				//return false;
			//}
		}
		// Zero remaining elements A(*,i)
		for (k = 0; k < n; k++){
			if( k != i ){
				r = (AtNinvA->at(k))[i];
				for(j = 0; j < n; j++){
					(AtNinvA->at(k))[j] = (AtNinvA->at(k))[j] - r*(AtNinvA->at(i))[j];
					(AtNinvAinv->at(k))[j] = (AtNinvAinv->at(k))[j] - r*(AtNinvAinv->at(i))[j];
				}
			}
		}
	}
	//t2=clock();
	//if(TIME) cout << ((float)t2-(float)t1) / CLOCKS_PER_SEC << "sec ";
	return true;
}

bool invert(vector<vector<float> > * AtNinvAori, vector<vector<double> > * AtNinvAinv ){//Gauss–Jordan elimination
	string METHODNAME = "invert";
	//WARNING: PASSING OF ARGUMENT VALUES NOT TESTED !!!!!!!!!!!!

	/* 	subroutine invert(np,n,A,B)
	! B = A^{-1}, where A is an arbitraty non-singular matrix.
	! In the process, A gets destroyed (A=I on exit).
	! Note that A need NOT be symmetric.
	! This routine is designed for small matrices where you dont care about
	! speed or numerical stability.
	! Written & tested by Max Tegmark 000506
	*/

	uint i=0, j=0, k=0, n = AtNinvAori->size();//todo check size
	vector<vector<float> > AtNinvAqaz(AtNinvAori->size(), vector<float>(AtNinvAori->at(0).size(),0));
	vector<vector<float> > *AtNinvA = &AtNinvAqaz;
	for( i = 0; i < AtNinvAori->size(); i++){
		for( j = 0; j < AtNinvAori->at(0).size(); j++){
			AtNinvA->at(i)[j] = AtNinvAori->at(i)[j];
		}
	}
	double r;
	//clock_t t1,t2;
	//t1=clock();


	for( i = 0; i < n; i++){// set matrix to one
		for (j = 0; j < n; j++){
			(AtNinvAinv->at(i))[j]=0.0;
		}
		(AtNinvAinv->at(i))[i]=1.0;
	}
	//t2=clock();
	//if(TIME) cout << ((float)t2-(float)t1) / CLOCKS_PER_SEC << "sec ";
	for ( i = 0; i < n; i++){ // find inverse by making A(i,i)=1

		if(isnan(AtNinvA->at(i)[i])){
			printf("%s: FATAL ERROR: input matrix AtNinvA of size %i by %i has NaN on diagonal! ABORT!\n", METHODNAME.c_str(), n, n);
			return false;
		}

		if(fabs((AtNinvA->at(i))[i]) < MIN_NONE_ZERO){
			if (i == n - 1){
				return false;
			}
			for(j = i+1; j < n; j++){//find a row to add to row i to make its diagonal non zero
				if(fabs((AtNinvA->at(j))[i]) >= MIN_NONE_ZERO and !isnan((AtNinvA->at(j))[i])){
					for(k = 0; k < n; k++){
						(AtNinvA->at(i))[k]=(AtNinvA->at(i))[k] + (AtNinvA->at(j))[k];
						(AtNinvAinv->at(i))[k]= (AtNinvAinv->at(i))[k] + (AtNinvAinv->at(j))[k];
					}
					break;
				}
				if (j == n - 1){
					return false;
				}
			}
		}
		r = (AtNinvA->at(i))[i];
		for (j = 0; j < n; j++){
			(AtNinvA->at(i))[j]=(AtNinvA->at(i))[j]/r;
			(AtNinvAinv->at(i))[j]= (AtNinvAinv->at(i))[j]/r;
			//if(j == i and (AtNinvAinv->at(i))[j]<0) {
				//printf("%s: FATAL ERROR: inverse matrix AtNinvAinv of size %i by %i has negative number on diagonal! ABORT!\n", METHODNAME.c_str(), n, n);
				//return false;
			//}
		}
		// Zero remaining elements A(*,i)
		for (k = 0; k < n; k++){
			if( k != i ){
				r = (AtNinvA->at(k))[i];
				for(j = 0; j < n; j++){
					(AtNinvA->at(k))[j] = (AtNinvA->at(k))[j] - r*(AtNinvA->at(i))[j];
					(AtNinvAinv->at(k))[j] = (AtNinvAinv->at(k))[j] - r*(AtNinvAinv->at(i))[j];
				}
			}
		}
	}
	//t2=clock();
	//if(TIME) cout << ((float)t2-(float)t1) / CLOCKS_PER_SEC << "sec ";
	return true;
}



///////////////REDUNDANT BASELINE CALIBRATION STUFF///////////////////
/////////////////////////////////////////////


/******************************************************/
/******************************************************/
void vecmatmul(vector<vector<double> > * Afitting, vector<float> * v, vector<float> * ampfit){
	int i, j;
	double sum;
	int n = Afitting->size();//todo size check
	int m = v->size();
	for(i=0; i < n; i++){
		sum = 0.0;
		for(j = 0; j < m; j++){
			sum = sum + (Afitting->at(i))[j] * (v->at(j));
		}
		(ampfit->at(i)) = sum;
	}
	return;
}

void vecmatmul(vector<vector<float> > * Afitting, vector<float> * v, vector<float> * ampfit){
	int i, j;
	double sum;
	int n = Afitting->size();//todo size check
	int m = v->size();
	for(i=0; i < n; i++){
		sum = 0.0;
		for(j = 0; j < m; j++){
			sum = sum + (Afitting->at(i))[j] * (v->at(j));
		}
		(ampfit->at(i)) = sum;
	}
	return;
}

void vecmatmul(vector<vector<int> > * A, vector<float> * v, vector<float> * yfit){
	int i, j;
	double sum;
	int n = A->size();//todo size check
	int m = v->size();
	for(i=0; i < n; i++){
		sum = 0.0;
		for(j = 0; j < m; j++){
			sum = sum + (A->at(i))[j] * (v->at(j));
		}
		(yfit->at(i)) = sum;
	}
	return;
}


/******************************************************/
/******************************************************/


void logcaladd(vector<vector<float> >* data, vector<vector<float> >* additivein, redundantinfo* info, vector<float>* calpar, vector<vector<float> >* additiveout, int command, calmemmodule* module){
	int nant = info->nAntenna;
	int nubl = info->nUBL;
	//int nbl = info->nBaseline;
	int ncross = info->crossindex.size();
	////read in amp and args
	for (int b = 0; b < ncross; b++){

		if ((data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0] == 0) and (data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1] == 0)){//got 0, quit
			for(int i = 3; i < 3 + 2 * nant + 2 * nubl; i++){
				calpar->at(i) = 0;
			}
			calpar->at(1) = INFINITY;
			return;
		}

		module->amp1[b] = log10(amp(data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0], data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1]));
		module->pha1[b] = phase(data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0], data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1]) * info->reversed[b];
	}
	////rewrap args
	for(int i = 0; i < nubl; i ++){
		for (uint j = 0; j < (module->ublgrp1)[i].size(); j ++){
			(module->ublgrp1)[i][j] = module->pha1[info->ublindex[i][j][2]];
		}
	}

	for (int i = 0; i < nubl; i++){
		(module->ubl1)[i][1] = medianAngle(&((module->ublgrp1)[i]));
	}

	for (int b = 0; b < ncross; b++) {
		module->pha1[b] = phaseWrap(module->pha1[b], (module->ubl1)[info->bltoubl[b]][1] - PI);
	}

	fill(module->x3.begin(), module->x3.end(), 0);////At.y
	for (unsigned int i = 0; i < info->Atsparse.size(); i++){
		for (unsigned int j = 0; j < info->Atsparse[i].size(); j++){
			module->x3[i] += module->amp1[info->Atsparse[i][j]];
		}
	}
	fill(module->x4.begin(), module->x4.end(), 0);////Bt.y
	for (unsigned int i = 0; i < info->Btsparse.size(); i++){
		for (unsigned int j = 0; j < info->Btsparse[i].size(); j++){
			module->x4[i] += module->pha1[info->Btsparse[i][j][0]] * info->Btsparse[i][j][1];
		}
	}
	vecmatmul(&(info->AtAi), &(module->x3), &(module->x1));
	vecmatmul(&(info->BtBi), &(module->x4), &(module->x2));
	//vecmatmul(&(info->AtAiAt), &(module->amp1), &(module->x1));////This is actually slower than seperate multiplications
	//vecmatmul(&(info->BtBiBt), &(module->pha1), &(module->x2));


	for(int b = 0; b < ncross; b++) {
		float amp = pow(10, module->x1[nant + info->bltoubl[b]] + module->x1[info->bl2d[info->crossindex[b]][0]] + module->x1[info->bl2d[info->crossindex[b]][1]]);
		float phase =  module->x2[nant + info->bltoubl[b]] * info->reversed[b] - module->x2[info->bl2d[info->crossindex[b]][0]] + module->x2[info->bl2d[info->crossindex[b]][1]];
		additiveout->at(info->crossindex[b])[0] = data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0] - amp * cos(phase);
		additiveout->at(info->crossindex[b])[1] = data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1] - amp * sin(phase);
	}
	if(command == 0){////compute additive term only
		calpar->at(1) = pow(norm(additiveout), 2);
		//cout << norm(additiveout) << endl;
		return;
	} else if(command == 1){////compute full set of calpars
		for(int a = 0; a < nant; a++){
			calpar->at(3 + a) = module->x1[a];
			calpar->at(3 + nant + a) = module->x2[a];
		}
		for(int u = 0; u < nubl; u++){
			calpar->at(3 + 2 * nant + 2 * u) = pow(10, module->x1[nant + u]) * cos(module->x2[nant + u]);
			calpar->at(3 + 2 * nant + 2 * u + 1) = pow(10, module->x1[nant + u]) * sin(module->x2[nant + u]);
		}
		calpar->at(1) = pow(norm(additiveout), 2);
	}
	return;
}

vector<float> minimizecomplex(vector<vector<float> >* a, vector<vector<float> >* b){
	vector<float> sum1(2, 0);
	for (uint i =0; i < a->size(); i++){
		sum1[0] += a->at(i)[0] * b->at(i)[0] + a->at(i)[1] * b->at(i)[1];
		sum1[1] += a->at(i)[1] * b->at(i)[0] - a->at(i)[0] * b->at(i)[1];
	}
	float sum2 = pow(norm(b), 2);
	sum1[0] = sum1[0] / sum2;
	sum1[1] = sum1[1] / sum2;
	return sum1;
}

void lincal(vector<vector<float> >* data, vector<vector<float> >* additivein, redundantinfo* info, vector<float>* calpar, vector<vector<float> >* additiveout, int command, calmemmodule* module, float convergethresh, int maxiter, float stepsize){
	//cout << "lincal DBG" << info->ublindex[(info->nUBL)-1][0][0] << " " << info->ublindex[(info->nUBL)-1][0][1] << " " << info->ublindex[(info->nUBL)-1][0][2] <<endl<<flush;
	//int DBGg1 = info->ublindex[(info->nUBL)-1][0][0];
	//int DBGg2 = info->ublindex[(info->nUBL)-1][0][1];
	//int DBGbl = info->ublindex[(info->nUBL)-1][0][2];

	////initialize data and g0 ubl0
	for (unsigned int b = 0; b < (module->cdata1).size(); b++){
		module->cdata1[b][0] = data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0];
		module->cdata1[b][1] = data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1];
	}
	float amptmp;
	unsigned int cbl;
	float stepsize2 = 1 - stepsize;
	for (int a = 0; a < info->nAntenna; a++){
		amptmp = pow(10, calpar->at(3 + a));
		module->g0[a][0] = amptmp * cos(calpar->at(3 + info->nAntenna + a));
		module->g0[a][1] = amptmp * sin(calpar->at(3 + info->nAntenna + a));
	}
	if (command != 1){
		for (int u = 0; u < info->nUBL; u++){
			module->ubl0[u][0] = calpar->at(3 + 2 * info->nAntenna + 2 * u);
			module->ubl0[u][1] = calpar->at(3 + 2 * info->nAntenna + 2 * u + 1);
		}
	} else{//if command is 1, compute the ubl estimates given data and calpars, rather than read ubl estimates from input
		for (int u = 0; u < info->nUBL; u++){
			for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
				cbl = info->ublindex[u][i][2];
				module->ubl2dgrp1[u][i][0] = module->cdata1[cbl][0];
				module->ubl2dgrp1[u][i][1] = module->cdata1[cbl][1] * info->reversed[cbl];
				module->ubl2dgrp2[u][i][0] = module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][0] + module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][1];
				module->ubl2dgrp2[u][i][1] = (module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][1] - module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][0]) * info->reversed[cbl];
			}

			module->ubl0[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
		}
	}

	float gre, gim, starting_chisq, chisq, chisq2, delta;
	int a1, a2;
	chisq = 0;
	for (unsigned int b = 0; b < (module->cdata2).size(); b++){
		a1 = info->bl2d[info->crossindex[b]][0];
		a2 = info->bl2d[info->crossindex[b]][1];
		gre = module->g0[a1][0] * module->g0[a2][0] + module->g0[a1][1] * module->g0[a2][1];
		gim = module->g0[a1][0] * module->g0[a2][1] - module->g0[a1][1] * module->g0[a2][0];
		module->cdata2[b][0] = gre * module->ubl0[info->bltoubl[b]][0] - gim * module->ubl0[info->bltoubl[b]][1] * info->reversed[b];
		module->cdata2[b][1] = gre * module->ubl0[info->bltoubl[b]][1] * info->reversed[b] + gim * module->ubl0[info->bltoubl[b]][0];
		delta = (pow(module->cdata2[b][0] - module->cdata1[b][0], 2) + pow(module->cdata2[b][1] - module->cdata1[b][1], 2));
		chisq += delta;
		//if (delta != 0){
			//cout << delta << " " << module->cdata2[b][0]-1 << " " << module->cdata2[b][1] << " " << module->ubl0[info->bltoubl[b]][0]-1 << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
		//}
		//cout << gre << " " << gim << " " << module->ubl0[info->bltoubl[b]][0] << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
	}
	starting_chisq = chisq;
	//cout << "lincal DBG v " << module->cdata1[DBGbl][0] << " " <<  module->cdata1[DBGbl][1] << endl<<flush;
	//cout << "lincal DBG c0 g0 g0 " << module->ubl0[info->nUBL - 1][0] << " " <<  module->ubl0[info->nUBL -1][1]  << " " << module->g0[DBGg1][0] << " " <<  module->g0[DBGg1][1]  << " " << module->g0[DBGg2][0] << " " <<  module->g0[DBGg2][1] << endl<<flush;
	//cout << "lincal DBG c0g0g0 "  << module->cdata2[DBGbl][0] << " " << module->cdata2[DBGbl][1] << endl<<flush;

	////start iterations
	int iter = 0;
	float componentchange = 100;
	while(iter < maxiter and componentchange > convergethresh){
		iter++;
		//cout << "iteration #" << iter << endl; cout.flush();
		////calpar g

		for (unsigned int a3 = 0; a3 < module->g3.size(); a3++){////g3 will be containing the final dg, g1, g2 will contain a and b as in the cost function LAMBDA = ||a + b*g||^2
			for (unsigned int a = 0; a < module->g3.size(); a++){
				cbl = info->bl1dmatrix[a3][a];
				if (cbl < 0 or cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){//badbl or ubl has only 1 bl
					module->g1[a] = vector<float>(2,0);
					module->g2[a] = vector<float>(2,0);
				}else if(info->bl2d[info->crossindex[cbl]][1] == a3){
					module->g1[a] = module->cdata1[cbl];
					module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * info->reversed[cbl]);
					module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * info->reversed[cbl] - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				}else{
					module->g1[a][0] = module->cdata1[cbl][0];
					module->g1[a][1] = -module->cdata1[cbl][1];////vij needs to be conjugated
					module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]));////Mi-j needs to be conjugated
					module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]) - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				}
			}
			//(module->g1)[a3] = vector<float>(2,0);
			//(module->g2)[a3] = (module->g1)[a3];
			//for (unsigned int a = a3 + 1; a < module->g3.size(); a++){
				//cbl = info->bl1dmatrix[a3][a];
				//if (cbl < 0 or cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){//badbl or ubl has only 1 bl
					//module->g1[a] = vector<float>(2,0);
					//module->g2[a] = vector<float>(2,0);
				//}else{
					//module->g1[a][0] = module->cdata1[cbl][0];
					//module->g1[a][1] = -module->cdata1[cbl][1];////vij needs to be conjugated
					//module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]));////Mi-j needs to be conjugated
					//module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]) - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				//}
			//}
			module->g3[a3] = minimizecomplex(&(module->g1), &(module->g2));
		}

		////ubl M
		for (int u = 0; u < info->nUBL; u++){
			for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
				cbl = info->ublindex[u][i][2];
				module->ubl2dgrp1[u][i][0] = module->cdata1[cbl][0];
				module->ubl2dgrp1[u][i][1] = module->cdata1[cbl][1] * info->reversed[cbl];
				module->ubl2dgrp2[u][i][0] = module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][0] + module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][1];
				module->ubl2dgrp2[u][i][1] = (module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][1] - module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][0]) * info->reversed[cbl];
			}

			module->ubl3[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
		}


		////Update g and ubl, do not update single-bl bls since they are not reversible. Will reverse this step later is chisq increased
		//float fraction;
		for (unsigned int a = 0; a < module->g3.size(); a++){
			module->g0[a][0] = stepsize2 * module->g0[a][0] + stepsize * module->g3[a][0];
			module->g0[a][1] = stepsize2 * module->g0[a][1] + stepsize * module->g3[a][1];

		}
		for (unsigned int u = 0; u < module->ubl3.size(); u++){
			if ((info->ublcount)[u] > 1){
				module->ubl0[u][0] = stepsize2 * module->ubl0[u][0] + stepsize * module->ubl3[u][0];
				module->ubl0[u][1] = stepsize2 * module->ubl0[u][1] + stepsize * module->ubl3[u][1];
			}
			//else{
				////make sure there's no error on unique baselines with only 1 baseline
				//for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
					//cbl = info->ublindex[u][i][2];
					//module->ubl2dgrp1[u][i][0] = module->cdata1[cbl][0];
					//module->ubl2dgrp1[u][i][1] = module->cdata1[cbl][1] * info->reversed[cbl];
					//module->ubl2dgrp2[u][i][0] = module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][0] + module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][1];
					//module->ubl2dgrp2[u][i][1] = (module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][1] - module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][0]) * info->reversed[cbl];
				//}

				//module->ubl3[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
				//module->ubl0[u][0] = module->ubl3[u][0];
				//module->ubl0[u][1] = module->ubl3[u][1];
			//}
		}

		//compute chisq and decide convergence
		chisq2 = 0;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			if ((info->ublcount)[info->bltoubl[b]] > 1){//automatically use 0 for single-bl ubls, their actaul values are not updated yet
				a1 = info->bl2d[info->crossindex[b]][0];
				a2 = info->bl2d[info->crossindex[b]][1];
				gre = module->g0[a1][0] * module->g0[a2][0] + module->g0[a1][1] * module->g0[a2][1];
				gim = module->g0[a1][0] * module->g0[a2][1] - module->g0[a1][1] * module->g0[a2][0];
				module->cdata2[b][0] = gre * module->ubl0[info->bltoubl[b]][0] - gim * module->ubl0[info->bltoubl[b]][1] * info->reversed[b];
				module->cdata2[b][1] = gre * module->ubl0[info->bltoubl[b]][1] * info->reversed[b] + gim * module->ubl0[info->bltoubl[b]][0];
				chisq2 += (pow(module->cdata2[b][0] - module->cdata1[b][0], 2) + pow(module->cdata2[b][1] - module->cdata1[b][1], 2));
				//cout << gre << " " << gim << " " << module->ubl0[info->bltoubl[b]][0] << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
			}
		}
		componentchange = (chisq - chisq2) / chisq;

		if (componentchange > 0){//if improved, keep g0 and ubl0 updates, and update single-bl ubls and chisq
			chisq = chisq2;
			for (unsigned int u = 0; u < module->ubl3.size(); u++){
			//make sure there's no error on unique baselines with only 1 baseline
				for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
					cbl = info->ublindex[u][i][2];
					module->ubl2dgrp1[u][i][0] = module->cdata1[cbl][0];
					module->ubl2dgrp1[u][i][1] = module->cdata1[cbl][1] * info->reversed[cbl];
					module->ubl2dgrp2[u][i][0] = module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][0] + module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][1];
					module->ubl2dgrp2[u][i][1] = (module->g0[info->ublindex[u][i][0]][0] * module->g0[info->ublindex[u][i][1]][1] - module->g0[info->ublindex[u][i][0]][1] * module->g0[info->ublindex[u][i][1]][0]) * info->reversed[cbl];
				}

				module->ubl3[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
				module->ubl0[u][0] = module->ubl3[u][0];
				module->ubl0[u][1] = module->ubl3[u][1];
			}
		} else {//reverse g0 and ubl0 changes
			iter--;
			for (unsigned int a = 0; a < module->g3.size(); a++){
				module->g0[a][0] = (module->g0[a][0] - stepsize * module->g3[a][0]) / stepsize2;
				module->g0[a][1] = (module->g0[a][1] - stepsize * module->g3[a][1]) / stepsize2;

			}
			for (unsigned int u = 0; u < module->ubl3.size(); u++){
				if ((info->ublcount)[u] > 1){
					module->ubl0[u][0] = (module->ubl0[u][0] - stepsize * module->ubl3[u][0]) / stepsize2;
					module->ubl0[u][1] = (module->ubl0[u][1] - stepsize * module->ubl3[u][1]) / stepsize2;
				}
			}
		}
	}


	////update calpar and additive term
	if(componentchange > 0 or iter > 1){
		for (unsigned int a = 0; a < module->g0.size(); a++){
			calpar->at(3 + a) = log10(amp(&(module->g0[a])));
			calpar->at(3 + info->nAntenna + a) = phase(&(module->g0[a]));
		}
		int tmp = 3 + 2 * info->nAntenna;
		for (unsigned int u = 0; u < module->ubl0.size(); u++){
			calpar->at(tmp + 2 * u) = module->ubl0[u][0];
			calpar->at(tmp + 2 * u + 1) = module->ubl0[u][1];
		}

		calpar->at(0) += iter;
		calpar->at(2) = chisq;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			additiveout->at(info->crossindex[b])[0] = module->cdata1[b][0] - module->cdata2[b][0];
			additiveout->at(info->crossindex[b])[1] = module->cdata1[b][1] - module->cdata2[b][1];
		}
	}else{////if chisq didnt decrease, keep everything untouched
		calpar->at(0) += 0;
		calpar->at(2) = starting_chisq;
	}
	//cout << "lincal DBG v "  << module->cdata1[DBGbl][0] << " " << module->cdata1[DBGbl][1] << endl<<flush;
	//cout << "lincal DBG c0g0g0 "  << module->cdata2[DBGbl][0] << " " << module->cdata2[DBGbl][1] << endl<<flush;
	return;
}


void gaincal(vector<vector<float> >* data, vector<vector<float> >* additivein, redundantinfo* info, vector<float>* calpar, vector<vector<float> >* additiveout, calmemmodule* module, float convergethresh, int maxiter, float stepsize){
	//cout << "lincal DBG" << info->ublindex[(info->nUBL)-1][0][0] << " " << info->ublindex[(info->nUBL)-1][0][1] << " " << info->ublindex[(info->nUBL)-1][0][2] <<endl<<flush;
	//int DBGg1 = info->ublindex[(info->nUBL)-1][0][0];
	//int DBGg2 = info->ublindex[(info->nUBL)-1][0][1];
	//int DBGbl = info->ublindex[(info->nUBL)-1][0][2];

	////initialize data and g0 ubl0
	for (unsigned int b = 0; b < (module->cdata1).size(); b++){
		module->cdata1[b][0] = data->at(info->crossindex[b])[0] - additivein->at(info->crossindex[b])[0];
		module->cdata1[b][1] = data->at(info->crossindex[b])[1] - additivein->at(info->crossindex[b])[1];
	}
	float amptmp;
	unsigned int cbl;
	float stepsize2 = 1 - stepsize;
	for (int a = 0; a < info->nAntenna; a++){
		amptmp = pow(10, calpar->at(3 + a));
		module->g0[a][0] = amptmp * cos(calpar->at(3 + info->nAntenna + a));
		module->g0[a][1] = amptmp * sin(calpar->at(3 + info->nAntenna + a));
	}

	for (int u = 0; u < info->nUBL; u++){
		module->ubl0[u][0] = 1;
		module->ubl0[u][1] = 0;
	}


	float gre, gim, starting_chisq, chisq, chisq2, delta;
	int a1, a2;
	chisq = 0;
	for (unsigned int b = 0; b < (module->cdata2).size(); b++){
		a1 = info->bl2d[info->crossindex[b]][0];
		a2 = info->bl2d[info->crossindex[b]][1];
		gre = module->g0[a1][0] * module->g0[a2][0] + module->g0[a1][1] * module->g0[a2][1];
		gim = module->g0[a1][0] * module->g0[a2][1] - module->g0[a1][1] * module->g0[a2][0];
		module->cdata2[b][0] = gre * module->ubl0[info->bltoubl[b]][0] - gim * module->ubl0[info->bltoubl[b]][1] * info->reversed[b];
		module->cdata2[b][1] = gre * module->ubl0[info->bltoubl[b]][1] * info->reversed[b] + gim * module->ubl0[info->bltoubl[b]][0];
		delta = (pow(module->cdata2[b][0] - module->cdata1[b][0], 2) + pow(module->cdata2[b][1] - module->cdata1[b][1], 2));
		chisq += delta;
		//if (delta != 0){
			//cout << delta << " " << module->cdata2[b][0]-1 << " " << module->cdata2[b][1] << " " << module->ubl0[info->bltoubl[b]][0]-1 << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
		//}
		//cout << gre << " " << gim << " " << module->ubl0[info->bltoubl[b]][0] << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
	}
	starting_chisq = chisq;
	//cout << "lincal DBG v " << module->cdata1[DBGbl][0] << " " <<  module->cdata1[DBGbl][1] << endl<<flush;
	//cout << "lincal DBG c0 g0 g0 " << module->ubl0[info->nUBL - 1][0] << " " <<  module->ubl0[info->nUBL -1][1]  << " " << module->g0[DBGg1][0] << " " <<  module->g0[DBGg1][1]  << " " << module->g0[DBGg2][0] << " " <<  module->g0[DBGg2][1] << endl<<flush;
	//cout << "lincal DBG c0g0g0 "  << module->cdata2[DBGbl][0] << " " << module->cdata2[DBGbl][1] << endl<<flush;

	////start iterations
	int iter = 0;
	float componentchange = 100;
	while(iter < maxiter and componentchange > convergethresh){
		iter++;
		//cout << "iteration #" << iter << endl; cout.flush();
		////calpar g

		for (unsigned int a3 = 0; a3 < module->g3.size(); a3++){////g3 will be containing the final dg, g1, g2 will contain a and b as in the cost function LAMBDA = ||a + b*g||^2
			for (unsigned int a = 0; a < module->g3.size(); a++){
				cbl = info->bl1dmatrix[a3][a];
				if (cbl < 0 or cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){//badbl or ubl has only 1 bl
					module->g1[a] = vector<float>(2,0);
					module->g2[a] = vector<float>(2,0);
				}else if(info->bl2d[info->crossindex[cbl]][1] == a3){
					module->g1[a] = module->cdata1[cbl];
					module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * info->reversed[cbl]);
					module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * info->reversed[cbl] - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				}else{
					module->g1[a][0] = module->cdata1[cbl][0];
					module->g1[a][1] = -module->cdata1[cbl][1];////vij needs to be conjugated
					module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]));////Mi-j needs to be conjugated
					module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]) - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				}
			}
			//(module->g1)[a3] = vector<float>(2,0);
			//(module->g2)[a3] = (module->g1)[a3];
			//for (unsigned int a = a3 + 1; a < module->g3.size(); a++){
				//cbl = info->bl1dmatrix[a3][a];
				//if (cbl < 0 or cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){//badbl or ubl has only 1 bl
					//module->g1[a] = vector<float>(2,0);
					//module->g2[a] = vector<float>(2,0);
				//}else{
					//module->g1[a][0] = module->cdata1[cbl][0];
					//module->g1[a][1] = -module->cdata1[cbl][1];////vij needs to be conjugated
					//module->g2[a][0] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][0] + module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]));////Mi-j needs to be conjugated
					//module->g2[a][1] = (module->g0[a][0] * module->ubl0[info->bltoubl[cbl]][1] * (-info->reversed[cbl]) - module->g0[a][1] * module->ubl0[info->bltoubl[cbl]][0]);
				//}
			//}
			module->g3[a3] = minimizecomplex(&(module->g1), &(module->g2));
		}


		////Update g and ubl, do not update single-bl bls since they are not reversible. Will reverse this step later is chisq increased
		//float fraction;
		for (unsigned int a = 0; a < module->g3.size(); a++){
			module->g0[a][0] = stepsize2 * module->g0[a][0] + stepsize * module->g3[a][0];
			module->g0[a][1] = stepsize2 * module->g0[a][1] + stepsize * module->g3[a][1];

		}

		//compute chisq and decide convergence
		chisq2 = 0;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			if ((info->ublcount)[info->bltoubl[b]] > 1){//automatically use 0 for single-bl ubls, their actaul values are not updated yet
				a1 = info->bl2d[info->crossindex[b]][0];
				a2 = info->bl2d[info->crossindex[b]][1];
				gre = module->g0[a1][0] * module->g0[a2][0] + module->g0[a1][1] * module->g0[a2][1];
				gim = module->g0[a1][0] * module->g0[a2][1] - module->g0[a1][1] * module->g0[a2][0];
				module->cdata2[b][0] = gre * module->ubl0[info->bltoubl[b]][0] - gim * module->ubl0[info->bltoubl[b]][1] * info->reversed[b];
				module->cdata2[b][1] = gre * module->ubl0[info->bltoubl[b]][1] * info->reversed[b] + gim * module->ubl0[info->bltoubl[b]][0];
				chisq2 += (pow(module->cdata2[b][0] - module->cdata1[b][0], 2) + pow(module->cdata2[b][1] - module->cdata1[b][1], 2));
				//cout << gre << " " << gim << " " << module->ubl0[info->bltoubl[b]][0] << " " << module->ubl0[info->bltoubl[b]][1] * info->reversed[b] << " " <<  a1 << " " <<  a2 << " " <<  b << " " << info->reversed[b] << endl;
			}
		}
		componentchange = (chisq - chisq2) / chisq;

		if (componentchange > 0){//if improved, keep g0 and ubl0 updates, and update single-bl ubls and chisq
			chisq = chisq2;
		} else {//reverse g0 and ubl0 changes
			iter--;
			for (unsigned int a = 0; a < module->g3.size(); a++){
				module->g0[a][0] = (module->g0[a][0] - stepsize * module->g3[a][0]) / stepsize2;
				module->g0[a][1] = (module->g0[a][1] - stepsize * module->g3[a][1]) / stepsize2;
			}
		}
	}


	////update calpar and additive term
	if(componentchange > 0 or iter > 1){
		for (unsigned int a = 0; a < module->g0.size(); a++){
			calpar->at(3 + a) = log10(amp(&(module->g0[a])));
			calpar->at(3 + info->nAntenna + a) = phase(&(module->g0[a]));
		}

		calpar->at(0) += iter;
		calpar->at(2) = chisq;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			additiveout->at(info->crossindex[b])[0] = module->cdata1[b][0] - module->cdata2[b][0];
			additiveout->at(info->crossindex[b])[1] = module->cdata1[b][1] - module->cdata2[b][1];
		}
	}else{////if chisq didnt decrease, keep everything untouched
		calpar->at(0) += 0;
		calpar->at(2) = starting_chisq;
	}
	//cout << "lincal DBG v "  << module->cdata1[DBGbl][0] << " " << module->cdata1[DBGbl][1] << endl<<flush;
	//cout << "lincal DBG c0g0g0 "  << module->cdata2[DBGbl][0] << " " << module->cdata2[DBGbl][1] << endl<<flush;
	return;
}

void loadGoodVisibilities(vector<vector<vector<vector<float> > > > * rawdata, vector<vector<vector<vector<float> > > >* receiver, redundantinfo* info, int xy){////0 for xx 3 for yy
	for (unsigned int t = 0; t < receiver->size(); t++){
		for (unsigned int f = 0; f < receiver->at(0).size(); f++){
			for (unsigned int bl = 0; bl < receiver->at(0)[0].size(); bl++){
				receiver->at(t)[f][bl][0] = rawdata->at(xy)[t][f][2 * info->subsetbl[bl]];
				receiver->at(t)[f][bl][1] = rawdata->at(xy)[t][f][2 * info->subsetbl[bl] + 1];
			}
		}
	}
	return;
}


void removeDegen(vector<float> *calpar, redundantinfo * info, calmemmodule* module){//forces the calibration parameters to have average 1 amp, and no shifting the image in phase. Note: 1) If you have not absolute calibrated the data, there's no point in running this, because this can only ensure that the calpars don't screw up already absolute calibrated data. 2) the antloc and ubl stored in redundant info must be computed from idealized antloc, otherwise the error in antloc from perfect redundancy will roll into this result, in an unknown fashion.
	////load data
	vector<float> pha1(info->nAntenna, 0);
	for (int a = 0 ; a < info->nAntenna; a ++){
		pha1[a] = calpar->at(3 + info->nAntenna + a);
	}
	for (int u = 0 ; u < info->nUBL; u ++){
		module->ubl1[u][0] = amp(calpar->at(3 + 2 * info->nAntenna + 2 * u), calpar->at(3 + 2 * info->nAntenna + 2 * u + 1));
		module->ubl1[u][1] = phase(calpar->at(3 + 2 * info->nAntenna + 2 * u), calpar->at(3 + 2 * info->nAntenna + 2 * u + 1));
	}

	////compute amp delta
	float ampfactor = 0;//average |g|, divide ant calpar by this, multiply ublfit by square this
	for (int a = 0 ; a < info->nAntenna; a ++){
		ampfactor += pow(10, calpar->at(3 + a));
	}
	ampfactor = ampfactor / info->nAntenna;
	//cout << ampfactor << endl;

	////compute phase delta
	vecmatmul(&(info->degenM), &(pha1), &(module->x1));//x1: add ant calpar and ubl fit by this

	////correct ant calpar
	for (int a = 0 ; a < info->nAntenna; a ++){
		calpar->at(3 + a) = calpar->at(3 + a) - log10(ampfactor);
		calpar->at(3 + info->nAntenna + a) = phaseWrap(calpar->at(3 + info->nAntenna + a) + module->x1[a]);
	}

	////correct ublfit
	for (int u = 0 ; u < info->nUBL; u ++){
		module->ubl2[u][0] = module->ubl1[u][0] * ampfactor * ampfactor;
		module->ubl2[u][1] = module->ubl1[u][1] + module->x1[info->nAntenna + u];
		calpar->at(3 + 2 * info->nAntenna + 2 * u) = module->ubl2[u][0] * cos(module->ubl2[u][1]);
		calpar->at(3 + 2 * info->nAntenna + 2 * u + 1) = module->ubl2[u][0] * sin(module->ubl2[u][1]);
	}
	return;
}

void runAverage1d(vector<float> *in, vector<float> *out, uint w){//compute running average with running length 2w+1. The first and last w elements are averaged with less elements.
	string METHODNAME = "runAverage1d";
	if(in->size() != out->size()){
		printf("#!!#%s#!!#%s: FATAL ERROR: input and output arrays have different dimensions: %lu vs %lu. ABORT!\n", FILENAME.c_str(), METHODNAME.c_str(), in->size(), out->size());
		return;
	}
	uint l = in->size();
	double sum = 0;
	uint n = 0;//number of items in sum
	for (uint i = 0; i < min(w, l); i++){
		sum += in->at(i);
		n++;
	}
	for (uint i = 0; i < out->size(); i++){
		if(i + w < l){
			sum += in->at(i + w);
			n++;
		}
		if(i - w -1 >= 0){
			sum += -(in->at(i - w - 1));
			n += -1;
		}
		out->at(i) = float(sum / n);
	}
}
void runAverage(vector<vector<vector<vector<float> > > > *in, int dimension, int w){//compute running average along dimension with running length 2w+1. The first and last w elements are averaged with less elements. Input array is modified!
	string METHODNAME = "runAverage";
	//if(in->size() != out->size() or in->at(0).size() != out->at(0).size()){
		//printf("#!!#%s#!!#%s: FATAL ERROR: input and output arrays have different dimensions: (%i, %i) vs (%i, %i). ABORT!\n", FILENAME, METHODNAME.c_str(), in->size(), in->at(0).size(), out->size(), out->at(0).size());
		//return;
	//}
	vector<float> dummy, dummy2;
	switch(dimension){
		case 0:
			dummy = vector<float>(in->size(), 0);
			break;
		case 1:
			dummy = vector<float>(in->at(0).size(), 0);
			break;
		case 2:
			dummy = vector<float>(in->at(0)[0].size(), 0);
			break;
		case 3:
			dummy = vector<float>(in->at(0)[0][0].size(), 0);
			break;
		default:
			printf("#!!#%s#!!#%s: FATAL ERROR: input array does not contain dimension %i. ABORT!\n", FILENAME.c_str(), METHODNAME.c_str(), dimension);
			return;
			break;
	}
	dummy2 = dummy;


	for (unsigned int t = 0; t < in->size(); t++){
		for (unsigned int f = 0; f < in->at(0).size(); f++){
			for (unsigned int b = 0; b < in->at(0)[0].size(); b++){
				for (unsigned int ri = 0; ri < in->at(0)[0][0].size(); ri++){
					switch(dimension){
						case 0:
							for (unsigned int i = 0; i < dummy.size(); i ++){
								dummy[i] = in->at(i)[f][b][ri];
							}
								runAverage1d(&dummy, &dummy2, w);
							for (unsigned int i = 0; i < dummy.size(); i ++){
								in->at(i)[f][b][ri] = dummy2[i];
							}
							break;
						case 1:
							for (unsigned int i = 0; i < dummy.size(); i ++){
								dummy[i] = in->at(t)[i][b][ri];
							}
								runAverage1d(&dummy, &dummy2, w);
							for (unsigned int i = 0; i < dummy.size(); i ++){
								in->at(t)[i][b][ri] = dummy2[i];
							}
							break;
						case 2:
							for (unsigned int i = 0; i < dummy.size(); i ++){
								dummy[i] = in->at(t)[f][i][ri];
							}
								runAverage1d(&dummy, &dummy2, w);
							for (unsigned int i = 0; i < dummy.size(); i ++){
								in->at(t)[f][i][ri] = dummy2[i];
							}
							break;
						case 3:
							for (unsigned int i = 0; i < dummy.size(); i ++){
								dummy[i] = in->at(t)[f][b][i];
							}
								runAverage1d(&dummy, &dummy2, w);
							for (unsigned int i = 0; i < dummy.size(); i ++){
								in->at(t)[f][b][i] = dummy2[i];
							}
							break;
					}
					if (dimension == 3) break;
				}
				if (dimension == 2) break;
			}
			if (dimension == 1) break;
		}
		if (dimension == 0) break;
	}

}

