//saved from shanalogcalfast and modified for x5. Original fortran code by Adrian Liu and tranlated to C++ by Shana Tribiano. Improvements made by Jeff Zheng.
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
#include "include/calibration_omni.h"
#include <iomanip>

using namespace std;
const string FILENAME = "omnical.cc";
const float PI = atan2(0,-1);
const float SPEEDC = 299.792458;
const bool DBG = false;
const bool TIME = true;
const float MAX_NONE_INF = pow(10, 10);
const float MIN_NONE_ZERO = pow(10, -10);
const double MAX_NONE_INFD = pow(10, 10);
const float MAX_POW_10 = 10;


/******************************************************/
/******************************************************/

/******************************************************/
/******************************************************/


int main(int argc, char *argv[]){
	string METHODNAME = "main";
	if (argc != 13){
		//cout << argc << endl;
		cout << "##" << FILENAME << "##" << METHODNAME << ": FALTAL ERROR: Incorrect input format! Expecting data path, info path, time count, frequency count, total number of antenna contained in the data, whether or not to remove degeneracy(1 for yes), and whether or not to remove additive term(1 for yes), with removal period, and whether to use logcal before lincal." << endl;
		return 0;
	}
	clock_t tStart = clock();
	string visin(argv[1]);
	string infopath(argv[2]);
	int nInt = atoi(argv[3]);
	int nFreq = atoi(argv[4]);
	int nAnt = atoi(argv[5]);
	bool removedegen = atoi(argv[6]);
	bool removeadd = atoi(argv[7]);
	int additivePeriod = atoi(argv[8]);
	string additivePeriodstr = argv[8];
	bool use_logcal = atoi(argv[9]);
	float converge_percent = atof(argv[10]);
	int max_iter  = atoi(argv[11]);
	float step_size = atof(argv[12]);
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

	return 0;
}
