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

bool readredundantinfo(string filepath, redundantinfo* info){//data file generated from 16. additive noise investigation _from_17.nb
	string METHODNAME = "readredundantinfo";
	uint parser = 0;
	vector<float> rawinfo = readAscii(filepath.c_str());
	if(rawinfo.size() < 3){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": length inconsistency when parsing " << filepath << "! Make sure file exists." << endl;
		return false;
	}

	info->nAntenna = rawinfo[parser];//number of good antennas among all (64) antennas, same as the length of subsetant
	parser++;

	info->nUBL = rawinfo[parser];//number of unique baselines
	parser++;

	int nbl = rawinfo[parser];//number f good baselines, auto included
	parser++;
	int ncross = nbl - info->nAntenna;
	info->nBaseline = nbl;
	info->nCross = ncross;


	(info->subsetant).resize(info->nAntenna);
	for (int i = 0; i < info->nAntenna; i++){
		(info->subsetant)[i] = rawinfo[parser];//the index of good antennas in all (64) antennas
		parser++;
	}
	(info->antloc).resize(info->nAntenna);
	for (int i = 0; i < info->nAntenna; i++){
		(info->antloc)[i].resize(3);
		for (int j = 0; j < 3; j++){
			(info->antloc)[i][j] = rawinfo[parser];//antloc vectors
			parser++;
		}
	}


	(info->subsetbl).resize(nbl);
	for (int i = 0; i < nbl; i++){
		(info->subsetbl)[i] = rawinfo[parser];//the index of good baselines (auto included) in all baselines
		parser++;
	}
	(info->ubl).resize(info->nUBL);
	for (int i = 0; i < info->nUBL; i++){
		(info->ubl)[i].resize(3);
		for (int j = 0; j < 3; j++){
			(info->ubl)[i][j] = rawinfo[parser];//unique baseline vectors
			parser++;
		}
	}
	(info->bltoubl).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->bltoubl)[i] = rawinfo[parser];//cross bl number to ubl index
		parser++;
	}
	(info->reversed).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->reversed)[i] = rawinfo[parser];//cross only bl if reversed -1, otherwise 1
		parser++;
	}

	(info->reversedauto).resize(nbl);
	for (int i = 0; i < nbl; i++){
		(info->reversedauto)[i] = rawinfo[parser];//the index of good baselines (auto included) in all baselines
		parser++;
	}
	(info->autoindex).resize(info->nAntenna);
	for (int i = 0; i < info->nAntenna; i++){
		(info->autoindex)[i] = rawinfo[parser];//index of auto bls among good bls
		parser++;
	}
	(info->crossindex).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->crossindex)[i] = rawinfo[parser];//index of cross bls among good bls
		parser++;
	}
	(info->bl2d).resize(nbl);
	for (int i = 0; i < nbl; i++){
		(info->bl2d)[i].resize(2);
		for (int j = 0; j < 2; j++){
			(info->bl2d)[i][j] = rawinfo[parser];//from 1d bl index to a pair of antenna numbers
			parser++;
		}
	}


	(info->ublcount).resize(info->nUBL);
	for (int i = 0; i < info->nUBL; i++){
		(info->ublcount)[i] = rawinfo[parser];//for each ubl, the number of good cross bls corresponding to it
		parser++;
	}
	(info->ublindex).resize(info->nUBL);
	for (int i = 0; i < info->nUBL; i++){
		(info->ublindex)[i].resize((info->ublcount)[i]);
		for (int j = 0; j < (info->ublcount)[i]; j++){
			(info->ublindex)[i][j].resize(3);
			for (int k = 0; k < 3; k++){
				(info->ublindex)[i][j][k] = rawinfo[parser];//for each ubl, the vector<int> contains (ant1, ant2, crossbl)
				parser++;
			}
		}
	}

	(info->bl1dmatrix).resize(info->nAntenna);
	for (int i = 0; i < info->nAntenna; i++){
		(info->bl1dmatrix)[i].resize(info->nAntenna);
		for (int j = 0; j < info->nAntenna; j++){
			(info->bl1dmatrix)[i][j] = rawinfo[parser];//a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
			parser++;
		}
	}

	//matrices
	(info->degenM).resize(info->nAntenna + info->nUBL);
	for (int i = 0; i < info->nAntenna + info->nUBL; i++){
		(info->degenM)[i].resize(info->nAntenna);
		for (int j = 0; j < info->nAntenna; j++){
			(info->degenM)[i][j] = rawinfo[parser];//A matrix for logcal amplitude
			parser++;
		}
	}

	(info->A).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->A)[i].resize(info->nUBL + info->nAntenna);
		for (int j = 0; j < info->nUBL + info->nAntenna; j++){
			(info->A)[i][j] = rawinfo[parser];//A matrix for logcal amplitude
			parser++;
		}
	}
	(info->B).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->B)[i].resize(info->nUBL + info->nAntenna);
		for (int j = 0; j < info->nUBL + info->nAntenna; j++){
			(info->B)[i][j] = rawinfo[parser];//B matrix for logcal amplitude
			parser++;
		}
	}

	////The sparse matrices are treated a little differently because they are not rectangular
	(info->Atsparse).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->Atsparse)[i].resize(0);
	}
	for (int i = 0; i < 3 * ncross; i++){
		(info->Atsparse)[rawinfo[parser]].push_back(rawinfo[parser + 1]);////parser + 2 is always 1 for At so omitted
		parser+=3;
	}

	(info->Btsparse).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->Btsparse)[i].resize(0);
	}

	vector<int> dummybl(2,0);
	for (int i = 0; i < 3 * ncross; i++){
		dummybl[0] = rawinfo[parser + 1];
		dummybl[1] = rawinfo[parser + 2];
		(info->Btsparse)[rawinfo[parser]].push_back(dummybl);
		parser+=3;
	}

	(info->AtAi).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->AtAi)[i].resize(info->nUBL + info->nAntenna);
		for (int j = 0; j < info->nUBL + info->nAntenna; j++){
			(info->AtAi)[i][j] = rawinfo[parser];//(AtA)^-1
			parser++;
		}
	}
	(info->BtBi).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->BtBi)[i].resize(info->nUBL + info->nAntenna);
		for (int j = 0; j < info->nUBL + info->nAntenna; j++){
			(info->BtBi)[i][j] = rawinfo[parser];//(BtB)^-1
			parser++;
		}
	}

	(info->AtAiAt).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->AtAiAt)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->AtAiAt)[i][j] = rawinfo[parser];//(AtA)^-1At
			parser++;
		}
	}
	(info->BtBiBt).resize(info->nUBL + info->nAntenna);
	for (int i = 0; i < info->nUBL + info->nAntenna; i++){
		(info->BtBiBt)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->BtBiBt)[i][j] = rawinfo[parser];//(BtB)^-1Bt
			parser++;
		}
	}

	(info->PA).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->PA)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->PA)[i][j] = rawinfo[parser];//A(AtA)^-1At
			parser++;
		}
	}
	(info->PB).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->PB)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->PB)[i][j] = rawinfo[parser];//B(BtB)^-1Bt
			parser++;
		}
	}
	(info->ImPA).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->ImPA)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->ImPA)[i][j] = rawinfo[parser];//I-PA
			parser++;
		}
	}
	(info->ImPB).resize(ncross);
	for (int i = 0; i < ncross; i++){
		(info->ImPB)[i].resize(ncross);
		for (int j = 0; j < ncross; j++){
			(info->ImPB)[i][j] = rawinfo[parser];//I-PB
			parser++;
		}
	}


	if(parser != rawinfo.size()){
		cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": length inconsistency when parsing " << filepath << "! Expecting " << parser << " elements in file but read " <<  rawinfo.size() << "!" << endl;
		return false;
	}

	return true;
}

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
void printv(vector<float> *v, int mini, int maxi){
	for (int i = max(0,mini); i < min(int(v->size()),maxi); i++) printf("%10.5f", v->at(i));
	if (maxi == 200 and v->size() > 200) printf("...more entries not shown");
	cout << endl;
	return;
}

void printvv(vector<vector<float> > *v, int mini, int maxi){
	for (int i = max(0,mini); i < min(int(v->size()),maxi); i++){
		cout << "(";
		for (uint j = 0; j < (v->at(i)).size(); j++){
			printf("%10.5f,", (v->at(i))[j]);
		}
		cout << ") ";
	}
	if (maxi == 200 and v->size() > 200) printf("...more entries not shown");
	cout << endl;
	return;
}

void printvv(vector<vector<int> > *v, int mini, int maxi){
	for (int i = max(0,mini); i < min(int(v->size()),maxi); i++){
		cout << "(";
		for (uint j = 0; j < (v->at(i)).size(); j++){
			printf("%10i,", (v->at(i))[j]);
		}
		cout << ") ";
	}
	if (maxi == 200 and v->size() > 200) printf("...more entries not shown");
	cout << endl;
	return;
}

void printv(vector<double> *v, int mini, int maxi){
	for (int i = max(0,mini); i < min(int(v->size()),maxi); i++) printf("%10.5f", v->at(i));
	if (maxi == 200 and v->size() > 200) printf("...more entries not shown");
	cout << endl;
	return;
}
void printv(vector<int> *v, int mini, int maxi){
	for (int i = max(0,mini); i < min(int(v->size()),maxi); i++) printf("%5.1i", v->at(i));
	if (maxi == 200 and v->size() > 200) printf("...more entries not shown");
	cout << endl;
	return;
}

string exec(string input_command) {//return result of stdout of a system command such as "python BLAH_BLAH". Won't catch stderr. Direct copy from internet http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
	//cout << " test-2 " + input_command;
	char* cmd = (char*)input_command.c_str();
	//cout << " test-1 " + input_command;
	//cout << cmd << endl;
	FILE* pipe = popen(cmd, "r");
	if (!pipe) return "ERROR";
	char buffer[1024];
	std::string result = "";
	//cout << " test0 " + input_command;
	while(!feof(pipe)) {
		if(fgets(buffer, 1024, pipe) != NULL){
			//cout << buffer << " ";
			result += buffer;
		}
	}
	pclose(pipe);
	if ( result.size() > 0 )
		result.erase(result.end()-1,result.end());
	//cout << " test1 " + input_command;
	return result;
}


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

////vector<float> pySunPosForks(string date, string time, float time_offset){//date and time in UTC such as ($python sunpos.py '2012/5/24' '12:22:56') and returns altitude(degree from horizon to north) and az (degree from north towards east)

	////string system_command = "python sunPos.py ";
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////return strtovf(exec(system_command));

////}

////void pyStarPosForks(string date, string time, float time_offset, vector<vector<float> >* results){
	////string METHODNAME = "pyStarPosForks";
	////string system_command = "python starPos.py ";
	////uint objects = 6;
	////if (results->size() != objects or (results->at(0)).size() != 2){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! the receiver array is not initialized as " << objects << " by 2. Exiting!!" << endl;
		////return;
	////}
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////vector<float> pos = strtovf(exec(system_command));
	////if (pos.size() < objects * 2){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << system_command << "! Most likely the fixed_bodies_X4.db file is not places under the same folder as starPos.py. Exiting!!" << endl;
		////return;
	////}
	////for (int i = 0; i < results->size(); i++){
		////(results->at(i))[0] = pos[2 * i];
		////(results->at(i))[1] = pos[2 * i + 1];
	////}
	////return;

////}

////void pySatPosForks(string date, string time, float time_offset, vector<vector<float> >* results){
	////string METHODNAME = "pySatPosForks";
	////string system_command = "python satPos.py ";
	////if (results->size() != NUM_OBJECTS or (results->at(0)).size() != 3){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! the receiver array is not initialized as " << NUM_OBJECTS << " by 3. Exiting!!" << endl;
		////return;
	////}
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////vector<float> pos = strtovf(exec(system_command));
	////if (pos.size() < NUM_OBJECTS * 3){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << system_command << "! Only grabbed " << pos.size() << " numbers from python. " << exec(system_command) << " Most likely the fixed_bodies_X4.tle file is not places under the same folder as starPos.py. Exiting!!" << endl;
		////return;
	////}
	////for (int i = 0; i < results->size(); i++){
		////(results->at(i))[0] = pos[3 * i];
		////(results->at(i))[1] = pos[3 * i + 1];
		////(results->at(i))[2] = pos[3 * i + 2];
	////}
	////return;

////}


////void pySatPos(string date, string time, float time_offset, float lat, float lon, float ele, vector<vector<float> >* results){//lat lon in degrees, ele in meters
	////string METHODNAME = "pySatPos";
	////string system_command = "python satPos.py ";
	////if (results->size() != NUM_OBJECTS or (results->at(0)).size() != 3){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! the receiver array is not initialized as " << NUM_OBJECTS << " by 3. Exiting!!" << endl;
		////return;
	////}
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset << " " << lat << " " << lon << " " << ele;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////vector<float> pos = strtovf(exec(system_command));
	////if (pos.size() < NUM_OBJECTS * 3){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << system_command << "! Only grabbed " << pos.size() << " numbers from python. " << exec(system_command) << " Most likely the fixed_bodies_X4.tle file is not places under the same folder as starPos.py. Exiting!!" << endl;
		////return;
	////}
	////for (int i = 0; i < results->size(); i++){
		////(results->at(i))[0] = pos[3 * i];
		////(results->at(i))[1] = pos[3 * i + 1];
		////(results->at(i))[2] = pos[3 * i + 2];
	////}
	////return;

////}

////void pyRequ2hor(string date, string time, float time_offset, float lat, float lon, vector<vector<float> >* results){//lat lon in degrees. Using this matrix in mathematica generates results that agree with pyephem output within 0.02 degrees on the few sources I checked at Forks
	////string METHODNAME = "pyRequtohor";
	////string system_command = "python rotationmatrix.py ";
	////if (results->size() != 3 or (results->at(0)).size() != 3){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! the receiver array is not initialized as 3 by 3. Exiting!!" << endl;
		////return;
	////}
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset << " " << lat << " " << lon;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////vector<float> pos = strtovf(exec(system_command));
	////if (pos.size() < 9){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << system_command << "! Only grabbed " << pos.size() << " numbers from python. " << exec(system_command) << " Most likely the C++ code is not places under the same folder as rotationmatrix.py. Exiting!!" << endl;
		////return;
	////}
	////for (int i = 0; i < 3; i++){
		////for (int j = 0; j < 3; j++){
			////(results->at(j))[i] = pos[3 * i + j];
		////}
	////}
	////return;


////}


////vector<float> pySunPos(string date, string time, float time_offset, float lat, float lon, float ele){//date and time in UTC such as ($python sunpos.py '2012/5/24' '12:22:56') and returns altitude(degree from horizon to north) and az (degree from north towards east)

	////string system_command = "python sunPos.py ";
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset << " " << lat << " " << lon << " " << ele;
	////system_command += buffer.str();
	//////cout << system_command << endl;
	////return strtovf(exec(system_command));

////}

////vector<string> pyTimeShift(string date, string time, float time_offset){
	////string system_command = "python timeShift.py ";
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////system_command += " ";
	////ostringstream buffer;
	////buffer << time_offset;
	////system_command += buffer.str();
	//////cout << "DEBUG " << system_command << endl;
	////vector<float> newTimes = strtovf(exec(system_command));
	//////cout << "DEBUG " << newTimes.size() << endl;
	////vector<string> returnString(2);
	////stringstream ss1, ss2;//create a stringstream
   	////ss1 << (int)(newTimes[0]) << "/" << (int)(newTimes[1]) << "/" << (int)(newTimes[2]) ;//add number to the stream
  	////returnString[0] = ss1.str();
   	////ss2 << (int)(newTimes[3]) << ":" << (int)(newTimes[4]) << ":" << (int)(newTimes[5]) ;//add number to the stream
  	////returnString[1] = ss2.str();
	////return returnString;
////}

////float pyTimeDif(string time1, string time2){//returns time1-time2 in seconds
	////string system_command = "python timeDif.py \"";
	////system_command += time1;
	////system_command += "\" \"";
	////system_command += time2 + "\"";
	////return atof(exec(system_command).c_str());
////}

////string pyOdfDate(string date, string time){
	////string system_command = "python odfDate.py ";
	////system_command += date;
	////system_command += " ";
	////system_command += time;
	////return exec(system_command);
////}

////int cmdCountLine(const char* inputfilename){
	////string command(inputfilename);
	////command = "wc " + command;
	////vector<float> results = strtovf(exec(command));
	////return (int)(results[0]);
////}

////string cmdAbsPath(string relativePath){
	////string command = "readlink -f " + relativePath;
	////return exec(command);
////}

////vector<string> cmdLs(string options){
	////string ls = exec("ls " + options);
	//////cout << ls << endl;
	////string line;
	////vector<string> lines;
	////stringstream ssls (ls);

	////while (getline( ssls, line )){
		////lines.push_back(line);
	////}
	////return lines;
////}

////void cmdMove(string a, string b){
	////string command = "mv " + a + " " + b;
	////exec(command);
	////return;
////}

////void cmdCopy(string a, string b){
	////string command = "cp -R " + a + " " + b;
	////exec(command);
	////return;
////}



//void readBinaryVisibility(const char* inputfilename, vector<vector<vector<vector<vector<float> > > > >* data, uint nPol, uint nIntegrations, uint nFrequencies, uint nBase){
	//string METHODNAME = "readBinaryVisibility";
	//if ( data->size() != nPol or (data->at(0)).size() != nIntegrations or (data->at(0))[0].size() != nFrequencies or (data->at(0))[0][0].size() != nBase){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << inputfilename << "! The receiver array is initialized at (pol, t, f, bl) = (" << data->size() << ", " << (data->at(0)).size() << ", " << (data->at(0))[0].size() << ", " << (data->at(0))[0][0].size() << "), where as the parameters are specified as (" << nPol << ", " << nIntegrations << ", " << nFrequencies << ", " << nBase << "). Exiting!!" << endl;
		//return;
	//}
	//FILE * file;
	//float * fmemblock;
	//unsigned long  len;
	//size_t result;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Reading File " << inputfilename << endl;

	//// Open file
	//file = fopen(inputfilename,"rb");
	//if (file == NULL) {
		//cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Reading " << inputfilename << " FAILED! Check if the file exists!" << endl;
		//return;
	//};
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Getting Length... ";
	//// Get length
	//fseek (file, 0, SEEK_END);
	//len = ftell (file);//Number of Bytes, where each float takes 4 Bytes (32bits), so number of floats is len / 4
	//rewind(file);
	//cout << len / 4 << " floats." << endl;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Importing Data...";cout.flush();
	//fmemblock = (float*) malloc (sizeof(float)*len / 4);
	//result = fread(fmemblock, sizeof(float), len/4, file);
	//// Close file and clear buffer
	//fclose(file);
	//cout << " Done!" << endl;

	////Transfer data to vectors
	//if (nIntegrations != len / 4 / nPol/ nFrequencies / nBase / 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The binary data seems to contain " << len / 4 / nPol / nFrequencies / nBase / 2 << " time slices, where as the parameter is specified as " << nIntegrations << " time slices. Exiting!!" << endl;
		//return;
	//};

	//int n  = 0;
	//for(int p = 0; p < nPol; p++){
		//string pol;
		//switch (p){
  		//case 0:
   			//pol = "xx";
   			//break;
		//case 1:
			//pol = "xy";
			//break;
		//case 2:
			//pol = "yx";
			//break;
		//case 3:
			//pol = "yy";
			//break;

		//default:
			//pol = "xx";
		//}
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Parsing polarization: " << pol << "\n";
		//for(int t = 0; t < nIntegrations; t++){
			//for(int f=0; f < nFrequencies; f++){
				//for(int v = 0; v < nBase; v++){
					//for(int c = 0; c < 2; c++){
						//(data->at(p))[t][f][v][c] = fmemblock[n];
						//n++;
					//}
				//}
			//}
		//}
	//}
	////printf("%e \n",data[0][0][1][0]);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Done Reading Vector \n";
	////delete[] memblock;
	//free (fmemblock);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Deleted memblock \n";
//}

//void readBinaryVisibilityLarge(const char* inputfilename, vector<vector<vector<vector<float> > > > * data, int nPol, int nIntegrations, int nFrequencies, int nBase){
	//string METHODNAME = "readBinaryVisibilityLarge";
	//if ( data->size() != nPol or (data->at(0)).size() != nIntegrations or (data->at(0))[0].size() != nFrequencies or (data->at(0))[0][0].size() != nBase * 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << inputfilename << "! The receiver array is initialized at (pol, t, f, bl) = (" << data->size() << ", " << (data->at(0)).size() << ", " << (data->at(0))[0].size() << ", " << (data->at(0))[0][0].size() / 2 << "), where as the parameters are specified as (" << nPol << ", " << nIntegrations << ", " << nFrequencies << ", " << nBase << "). Exiting!!" << endl;
		//return;
	//}
	//FILE * file;
	//float * fmemblock;
	//unsigned long  len;
	//size_t result;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Reading File " << inputfilename;

	//// Open file
	//file = fopen(inputfilename,"rb");
	//if (file == NULL) {
		//cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Reading " << inputfilename << " FAILED! Check if the file exists!" << endl;
		//return;
	//};
	//cout << ": Getting Length... ";
	//// Get length
	//fseek (file, 0, SEEK_END);
	//len = ftell (file);//Number of Bytes, where each float takes 4 Bytes (32bits), so number of floats is len / 4
	//rewind(file);
	//cout << len / 4 << " 32 bit floats. ";
	//cout << "Importing data into memory..."; cout.flush();
	//fmemblock = (float*) malloc (sizeof(float) * len / 4);
	//result = fread(fmemblock, sizeof(float), len / 4, file);
	//// Close file and clear buffer
	//fclose(file);


	////Transfer data to vectors
	////cout << len / 4 << " " << nPol << " " <<  nFrequencies << " " <<  nBase << endl;
	//if (nIntegrations != len / 4 / nPol/ nFrequencies / nBase / 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The binary data seems to contain " << len / 4 / nPol / nFrequencies / nBase / 2 << " time slices, where as the parameter is specified as " << nIntegrations << " time slices. Exiting!!" << endl;
		//return;
	//};

	//int n  = 0;
	//for(int p = 0; p < nPol; p++){
		//string pol;
		//switch (p){
  		//case 0:
   			//pol = "xx";
   			//break;
		//case 1:
			//pol = "xy";
			//break;
		//case 2:
			//pol = "yx";
			//break;
		//case 3:
			//pol = "yy";
			//break;

		//default:
			//pol = "xx";
		//}
		//if(nPol > 1) cout << " Parsing polarization: " << pol << "...";
		//for(int t = 0; t < nIntegrations; t++){
			//for(int f=0; f < nFrequencies; f++){
				//for(int v = 0; v < nBase; v++){
					//for(int c = 0; c < 2; c++){
						//(data->at(p))[t][f][gc(v, c)] = fmemblock[n];
						//n++;
					//}
				//}
			//}
		//}
	//}

	//free (fmemblock);
	//cout << " All done!" << endl;
//}


//void readBinaryVisibilityLargeConj(const char* inputfilename, vector<vector<vector<vector<float> > > > * data, int nPol, int nIntegrations, int nFrequencies, int nBase){
	//string METHODNAME = "readBinaryVisibilityLargeConj";
	//if ( data->size() != nPol or (data->at(0)).size() != nIntegrations or (data->at(0))[0].size() != nFrequencies or (data->at(0))[0][0].size() != nBase * 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << inputfilename << "! The receiver array is initialized at (pol, t, f, bl) = (" << data->size() << ", " << (data->at(0)).size() << ", " << (data->at(0))[0].size() << ", " << (data->at(0))[0][0].size() / 2 << "), where as the parameters are specified as (" << nPol << ", " << nIntegrations << ", " << nFrequencies << ", " << nBase << "). Exiting!!" << endl;
		//return;
	//}
	//FILE * file;
	//float * fmemblock;
	//unsigned long  len;
	//size_t result;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Reading File " << inputfilename << endl;

	//// Open file
	//file = fopen(inputfilename,"rb");
	//if (file == NULL) {
		//cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Reading " << inputfilename << " FAILED! Check if the file exists!" << endl;
		//return;
	//};
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Getting Length... ";
	//// Get length
	//fseek (file, 0, SEEK_END);
	//len = ftell (file);//Number of Bytes, where each float takes 4 Bytes (32bits), so number of floats is len / 4
	//rewind(file);
	//cout << len / 4 << " floats." << endl;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Importing data into memory...";
	//fmemblock = (float*) malloc (sizeof(float) * len / 4);
	//result = fread(fmemblock, sizeof(float), len / 4, file);
	//// Close file and clear buffer
	//fclose(file);
	//cout << " Done!" << endl;

	////Transfer data to vectors
	//if (nIntegrations != len / 4 / nPol/ nFrequencies / nBase / 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The binary data seems to contain " << len / 4 / nPol / nFrequencies / nBase / 2 << " time slices, where as the parameter is specified as " << nIntegrations << " time slices. Exiting!!" << endl;
		//return;
	//};

	//int n  = 0;
	//for(int p = 0; p < nPol; p++){
		//string pol;
		//switch (p){
  		//case 0:
   			//pol = "xx";
   			//break;
		//case 1:
			//pol = "xy";
			//break;
		//case 2:
			//pol = "yx";
			//break;
		//case 3:
			//pol = "yy";
			//break;

		//default:
			//pol = "xx";
		//}
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Parsing polarization: " << pol << "\n";
		////cout << data->size() << " " << data->at(0).size() << " " << data->at(0)[0].size() << " " << data->at(0)[0][0].size() << endl << fmemblock[len / 4 - 1]; cout.flush();
		//for(int t = 0; t < nIntegrations; t++){
			//for(int f = 0; f < nFrequencies; f++){
				//for(int v = 0; v < nBase; v++){
					////cout << p << " " << t << " " << f << " " << gc(v,0) << " " << n << " " << fmemblock[n]; cout.flush();
					//(data->at(p))[t][f][gc(v, 0)] = fmemblock[n];
					//n++;
					////cout << n << " " ; cout.flush();
					//(data->at(p))[t][f][gc(v, 1)] = -fmemblock[n];
					//n++;

				//}
			//}
		//}
	//}
	////printf("%e \n",data[0][0][1][0]);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Done Reading Vector \n";
	////delete[] memblock;
	//free (fmemblock);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Deleted memblock \n";
//}


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

//void readBinaryCalparSP(const char* inputfilename, vector<vector<vector<float> > > * data, int nIntegrations, int nFrequencies, int nAnt, int nUBL){//keep log10, turn degree into rad when reading
	//string METHODNAME = "readBinaryCalparSP";
	//if ( data->size() != nIntegrations or (data->at(0)).size() != nFrequencies or (data->at(0))[0].size() != 3 + 2 * (nAnt + nUBL)){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH for " << inputfilename << "! The receiver array is initialized at (t, f, 3+2(nAnt+nUBL)) = (" << data->size() << ", " << (data->at(0)).size() << ", " << (data->at(0))[0].size() << "), where as the parameters are specified as (" << nIntegrations << ", " << nFrequencies << ", " << 3 + 2 * (nAnt + nUBL) << "). Exiting!!" << endl;
		//return;
	//}
	//FILE * file;
	//float * fmemblock;
	//unsigned long  len;
	//size_t result;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Reading File " << inputfilename << endl;

	//// Open file
	//file = fopen(inputfilename,"rb");
	//if (file == NULL) {
		//cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Reading " << inputfilename << " FAILED! Check if the file exists!" << endl;
		//return;
	//};
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Getting Length... ";
	//// Get length
	//fseek (file, 0, SEEK_END);
	//len = ftell (file);//Number of Bytes, where each float takes 4 Bytes (32bits), so number of floats is len / 4
	//rewind(file);
	//cout << len / 4 << " floats." << endl;
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": Importing data into memory..."; cout.flush();
	//fmemblock = (float*) malloc (sizeof(float) * len / 4);
	//result = fread(fmemblock, sizeof(float), len / 4, file);
	//// Close file and clear buffer
	//fclose(file);
	//cout << " Done!" << endl;

	////Transfer data to vectors
	//if (nIntegrations != len / 4 / nFrequencies / (3 + 2 * (nAnt + nUBL))){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The binary data seems to contain " << len / 4 / nFrequencies / (3 + 2 * (nAnt + nUBL)) << " time slices, where as the parameter is specified as " << nIntegrations << " time slices. Exiting!!" << endl;
		//return;
	//};

	//int n  = 0;
	//for(int t = 0; t < nIntegrations; t++){
		//for(int f=0; f < nFrequencies; f++){
			//for(int chi = 0; chi < 3; chi++){
					//(data->at(t))[f][chi] = fmemblock[n];
					//n++;
			//}
			//for(int i = 3; i < 3 + nAnt; i++){//keep log10
					//(data->at(t))[f][i] = fmemblock[n];
					//n++;
			//}
			//for(int i = 3 + nAnt; i < 3 + 2 * nAnt; i++){//turn degree into rad
					//(data->at(t))[f][i] = fmemblock[n] * PI / 180;
					//n++;
			//}
			//for(int i = 3 + 2 * nAnt; i < 3 + 2 * nAnt + 2 * nUBL; i++){
					//(data->at(t))[f][i] = fmemblock[n];
					//n++;
			//}
		//}
	//}

	////printf("%e \n",data[0][0][1][0]);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Done Reading Vector \n";
	////delete[] memblock;
	//free (fmemblock);
	//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << " Deleted memblock \n";
	//return;
//}

vector<float> readAscii(const char* inputfilename, int count, bool verbose) {
	string METHODNAME = "readAscii";
	ifstream myReadFile;
	float tmpf;
	int cnter = 0;
	if (count <= 0){
		myReadFile.open(inputfilename);
		if (myReadFile.is_open()) {
		  if(verbose) cout << "The file is open" << endl;
		  while (!myReadFile.eof()) {
			 if(verbose) cout << "p";
			 myReadFile >> tmpf;
			 cnter++;
			 if(verbose and cnter%100 == 0) cout << cnter << " " << tmpf;
		  }
		  cnter--;
		  myReadFile.close();
		}
		else {
		  cout << "Error on file open " << inputfilename << endl;
		}
		if(verbose) cout << cnter << " numbers parsed in " << inputfilename << endl;
	}else{
		cnter = count;
	}

	vector<float> data(cnter);
	myReadFile.open(inputfilename);
	for (int i = 0; i < cnter; i++) {
	     if(myReadFile.eof()){
			cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: Ascii read failure: File seem to end at " << i << " while expecting " << cnter << " number!" << endl;
			break;
		 }
		 //if(verbose) cout << "r";
		 myReadFile >> data[i];
		 if (verbose and i < 5) cout << i << ": " << data[i] << "; ";
		 //cout << cnter << endl;
	}
	//if (DEBUG) cout << endl;
	myReadFile.close();
	if( DEBUG ){
		cout << "##" << FILENAME << "##" << METHODNAME << ": " << cnter << " numbers successfully read from " << inputfilename << endl;
	}
	if (cnter == 0){
		cout << "##" << FILENAME << "##" << METHODNAME << ": FATAL ERROR: Ascii read failure: No numbers read from " << inputfilename << endl;
	}
   if(verbose) cout << "First number:" << data[0] << endl;
   if(verbose) cout << "Last number:" << data[data.size() - 1] << endl;

   return data;
}

//void readVisibility(const char* inputfilename, vector<vector<vector<float> > > * data, int nFreq, int nBaseline){
	//string METHODNAME = "readVisibility";
	//if ( data->size() != nFreq or (data->at(0)).size() != nBaseline){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array is initialized at (f, bl) = (" << data->size() << ", " << (data->at(0)).size() << "), where as the parameters are specified as (" << nFreq << ", " << nBaseline << "). Exiting!!" << endl;
		//return;
	//}
	//int realNumberOfLine = cmdCountLine(inputfilename);
	//if (realNumberOfLine != nFreq){
		//cout << "#!!##!!#" << FILENAME << "#!!##!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input parameter specified number of frequency as " << nFreq << ", but there are " << realNumberOfLine << " lines in the data set " << inputfilename << ". Watch out!!" << endl;
	//}
	//vector<float> dataraw = readAscii(inputfilename);
	//int cnter = 0;
	//for (int f = 0; f < min(nFreq, realNumberOfLine); f++){
		//for (int bl = 0; bl < nBaseline; bl++){
			//for (int ri = 0; ri < 2; ri ++){
				////(&((&(data->at(f)))->at(bl)))->at(ri) = dataraw[cnter];
				//(data->at(f))[bl][ri] = dataraw[cnter];
				//cnter ++;
			//}
		//}
	//}
	//return;
//}

//bool readCalparAscii(const char* inputfilename, vector<float> *chisq, vector<vector<float> > *ampcalpar, vector<vector<float> > *phasecalpar, vector<vector<vector<float> > > * UBLcalpar, int nFreq, int nAnt, int nUBL){
	//string METHODNAME = "readCalpar";
	//if ( chisq->size() != nFreq or ampcalpar->size() != nFreq or phasecalpar->size() != nFreq or UBLcalpar->size() != nFreq ){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array is initialized at f = (" << chisq->size() << ", " << ampcalpar->size() << ", " << phasecalpar->size() << ", " << UBLcalpar->size() << "), where as the parameters are specified as (" << nFreq << ", " << nFreq << ", " << nFreq << ", " << nFreq << "). Exiting!!" << endl;
		//return false;
	//}
	//if ( (ampcalpar->at(0)).size() != nAnt or (phasecalpar->at(0)).size() != nAnt or (UBLcalpar->at(0)).size() != nUBL or (UBLcalpar->at(0))[0].size() != 2){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array is initialized at (ant, ant, (ubl, 2)) = (" << (ampcalpar->at(0)).size() << ", " << (phasecalpar->at(0)).size() << ", (" << (UBLcalpar->at(0)).size() << ", " << (UBLcalpar->at(0))[0].size() << ")), where as the parameters are specified as (" << nAnt << ", " << nAnt << ", (" << nUBL << ", " << 2 <<  ")). Exiting!!" << endl;
		//return false;
	//}
////	int realNumberOfLine = cmdCountLine(inputfilename);

	//vector<float> calparraw = readAscii(inputfilename);
	////cout << "rawcalpar size: " << calparraw.size() << endl;
	////for(int i = 0; i < calparraw.size(); i++) cout << calparraw[i] << " " <<
	////cout << endl;

	//int calparcnter = 0;
	//int nChisqConvention = 3; //number of entries for storing chi sq in the beginning of each freq in calpar files, currently following logcalpar convention, chisq, amp chisq, phase chisq, and most of the times we ignore the amp chisq and phase chisq
	//if (calparraw.size() != nFreq*(nChisqConvention + 2*nAnt + 2*nUBL)){
		//cout << "#!!##!!#" << FILENAME << "#!!##!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input parameter implied the calpar file to have " << nFreq*(nChisqConvention + 2*nAnt + 2*nUBL) << " entries, but there are " << calparraw.size() << " in the data set " << inputfilename << ". Exiting!!" << endl;
		//return false;
	//}
	//for (int f = 0; f < nFreq; f++){
		//chisq->at(f) = calparraw[calparcnter];
		//calparcnter = calparcnter + nChisqConvention;
		//for ( int i = 0; i < nAnt; i++){
			//(ampcalpar->at(f))[i] = calparraw[calparcnter];
			//calparcnter ++;
		//}
		//for ( int i = 0; i < nAnt; i++){
			//(phasecalpar->at(f))[i] = calparraw[calparcnter] * PI / 180;
			//calparcnter ++;
		//}
		//for ( int i = 0; i < nUBL; i++){
			//(UBLcalpar->at(f))[i][0] = calparraw[calparcnter];
			//(UBLcalpar->at(f))[i][1] = calparraw[calparcnter + 1];
			//calparcnter = calparcnter + 2;
		//}
	//}

	//return true;
//}

////void readAntloc(const char* inputfilename, vector<vector<float> > * antloc, vector<vector<float> > * cablelen, int nAnt){//cablelen structured as [x/y][nAnt]
	////string METHODNAME = "readAntloc";
	////if ( antloc->size() != nAnt or (cablelen->at(0)).size() != nAnt or (antloc->at(0)).size() != 3 or cablelen->size() != 2){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array for antloc is initialized at (" << antloc->size() << ", " << (antloc->at(0)).size() << "), where as the parameters are specified as (" << nAnt << ", " << 3 << "). The receiver array for cablelen is initialized at (" << cablelen->size() << ", " << (cablelen->at(0)).size() << "), where as the parameters are specified as (" << 2 << ", " << nAnt << "). Exiting!!" << endl;
		////return;
	////}
	////int realNumberOfLine = cmdCountLine(inputfilename);
	////if (realNumberOfLine != nAnt){
		////cout << "#!!##!!#" << FILENAME << "#!!##!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input parameter specified number of antenna as " << nAnt << ", but there are " << realNumberOfLine << " lines in the data set " << inputfilename << ". Watch out!!" << endl;
		////return;
	////}
	////vector<float> antlocraw = readAscii(inputfilename);
	////int jmax;
	////if (antlocraw.size() > 3*nAnt){
		////jmax = 5;
	////} else{
		////jmax = 3;
	////}

	////int cnter = 0;
	////for (int i = 0; i < min(nAnt, realNumberOfLine); i++){
		////for (int j = 0; j < jmax; j++){
			////if (j < 3) (antloc->at(i))[j] = antlocraw[cnter];
			////if (j == 3) cablelen->at(0)[i] = antlocraw[cnter];
			////if (j == 4) cablelen->at(1)[i] = antlocraw[cnter];
			////cnter ++;
		////}
	////}
	////return;
////}

////void readSunpos(const char* inputfilename, vector<vector<float> > * sunpos){//read sunpos.dat from x5 odf (which is in alt/az, and return a list of pairs in k vector (x y z) or (S E U), note that this is technically -k vector since k vector should point inwards.
	////string METHODNAME = "readSunpos";
	////if ( (sunpos->at(0)).size() != 3 ){
		////cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL I/O MISMATCH! The receiver array for sunpos is initialized at (" << sunpos->size() << ", " << (sunpos->at(0)).size() << "), where as the parameters are expected as (time, 3). Exiting!!" << endl;
		////return;
	////}
	////int realNumberOfLine = cmdCountLine(inputfilename);
	////if (realNumberOfLine != sunpos->size()){
		////cout << "#!!##!!#" << FILENAME << "#!!##!!#" << METHODNAME << ": FATAL I/O MISMATCH! The input parameter specified number of timeslice as " << sunpos->size() << ", but there are " << realNumberOfLine << " lines in the data set " << inputfilename << ". Watch out!!" << endl;
		////return;
	////}
	////vector<float> antlocraw = readAscii(inputfilename);

	////int cnter = 0;
	////float theta = 0, phi = 0;
	////for (int i = 0; i < min(int(sunpos->size()), realNumberOfLine); i++){
		////for (int j = 0; j < 2; j++){
			////if (j == 0) theta = PI/2 - antlocraw[cnter];
			////if (j == 1) phi = PI - antlocraw[cnter];
			////cnter ++;
		////}
		////sunpos->at(i) = tp2xyz(theta, phi);
	////}
	////return;
////}


////void outputChiSqAscii(vector<float>  * data) //takes chi square for one time slice each freq each antenna and write a set of calibration parameters under current directory
////{
 	 ////string METHODNAME = "outputChiSqAscii";
	 ////int numFreq = data->size();
	 ////char buffer[100];
	 ////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing Chi squares...";
	 ////cout.flush();
	 ////FILE * myFile;
	 ////sprintf(buffer,"./chiSquare.dat");
	 ////myFile = fopen(buffer,"w");
	 ////for(int f = 0; f < numFreq; f++){
		 ////fprintf(myFile, "%.5E", data->at(f));
	 	 ////fprintf(myFile, "\n");
	 ////}
	 ////fclose (myFile);
	 ////cout << "DONE!" << endl;
////}

////void outputPhaseCalParAscii(vector<vector<float> >  * data, int numUBL, string output_name) //takes phase corrections for one time slice each freq each antenna and write a set of calibration parameters under current directory
////{
 	 ////string METHODNAME = "outputPhaseCalParAscii";
	 ////int numFreq = data->size();
	 ////int numAnt = (data->at(0)).size();
	 ////char buffer[2 * output_name.size()];
	 ////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing " << output_name << " calibration parameters...";
	 ////cout.flush();
	 ////FILE * myFile;
	 ////sprintf(buffer, output_name.c_str());
	 ////myFile = fopen(buffer,"w");
	 ////for(int f = 0; f < numFreq; f++){
		 ////fprintf(myFile, "%.5E %.5E %.5E ", 0.0, 0.0, 0.0);
	 	 ////for(int b = 0; b < numAnt; b++){
			 ////fprintf(myFile, "%.5E ", 0.0);
		 ////}
		 ////for(int b = 0; b < numAnt; b++){
			 ////fprintf(myFile, "%.5E ", 180 * (&(data->at(f)))->at(b) / PI );//<========The only useful information, others are all dummy 0s for formatting convention. degree for ascii
		 ////}
		 ////for(int b = 0; b < numUBL; b++){
		 	////if ( b != ( numUBL - 1 ) ) {
				////fprintf(myFile, "%.5E %.5E ", 0.0, 0.0);
			////} else{
				////fprintf(myFile, "%.5E %.5E", 0.0, 0.0);
			////}
		 ////}
	 	 ////fprintf(myFile, "\n");
	 ////}
	 ////fclose (myFile);
	 ////cout << "DONE!" << endl;
////}


////void outputCalParAscii(vector<float> * chisq, vector<vector<float> > * ampcalpar, vector<vector<float> >  * phasecalpar, vector<vector<vector<float> > > * UBLcalpar, string output_name) //takes phase corrections for one time slice each freq each antenna and write a set of calibration parameters under current directory
////{
 	 ////string METHODNAME = "outputCalParAscii";
 	 ////int numUBL = (&(UBLcalpar->at(0)))->size();
	 ////int numFreq = ampcalpar->size();
	 ////int numAnt = (ampcalpar->at(0)).size();
	 ////char buffer[2 * output_name.size()];
	 ////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing "<< output_name << " calibration parameters...";
	 ////cout.flush();
	 ////FILE * myFile;
	 ////sprintf(buffer, output_name.c_str());
	 ////myFile = fopen(buffer,"w");
	 ////for(int f = 0; f < numFreq; f++){
		 ////fprintf(myFile, "%.5E %.5E %.5E ", chisq->at(f), 0.0, 0.0);
	 	 ////for(int b = 0; b < numAnt; b++){
			 ////fprintf(myFile, "%.5E ", (&(ampcalpar->at(f)))->at(b));
		 ////}
		 ////for(int b = 0; b < numAnt; b++){
			 ////fprintf(myFile, "%.5E ", 180 * (&(phasecalpar->at(f)))->at(b) / PI );//<========degree for ascii
		 ////}
		 ////for(int b = 0; b < numUBL; b++){
		 	////if ( b != ( numUBL - 1 ) ) {
				////fprintf(myFile, "%.5E %.5E ", ((&(UBLcalpar->at(f)))->at(b))[0], ((&(UBLcalpar->at(f)))->at(b))[1]);
			////} else{
				////fprintf(myFile, "%.5E %.5E", ((&(UBLcalpar->at(f)))->at(b))[0], ((&(UBLcalpar->at(f)))->at(b))[1]);
			////}
		 ////}
	 	 ////fprintf(myFile, "\n");
	 ////}
	 ////fclose (myFile);
	 ////cout << "DONE!" << endl;
////}

////void outputAscii(vector<float>* data, string output_name, int type, bool cmdoutput){//0 for overwrite, 1 for append
	////string METHODNAME = "outputAscii";
	////FILE * myFile;
	////if (type == 0){
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Writing " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "w");
	////}else{
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Appending to " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "a");
	////}
	////cout.flush();
	////for (unsigned int i = 0; i < data->size(); i++){
		////fprintf(myFile, "%.6E ", data->at(i));
	////}
	////fprintf(myFile, "\n");
	////fclose (myFile);
	////if (cmdoutput) cout << "DONE!" << endl;
	////cout.flush();
////}

////void outputAscii(vector<vector<float> >* data, string output_name, int type, bool cmdoutput){//0 for overwrite, 1 for append
	////string METHODNAME = "outputAscii";
	////FILE * myFile;
	////if (type == 0){
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Writing " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "w");
	////}else{
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Appending to " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "a");
	////}
	////cout.flush();
	////for (unsigned int i = 0; i < data->size(); i++){
		////for (unsigned int j = 0; j < data->at(i).size(); j++){
			////fprintf(myFile, "%.6E ", data->at(i)[j]);
		////}
		////fprintf(myFile, "\n");
	////}
	////fclose (myFile);
	////if (cmdoutput) cout << "DONE!" << endl;
	////cout.flush();
////}

////void outputAscii(vector<vector<vector<float> > >* data, string output_name, int type, bool cmdoutput){//0 for overwrite, 1 for append
	////string METHODNAME = "outputAscii";
	////FILE * myFile;
	////if (type == 0){
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Writing " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "w");
	////}else{
		////if (cmdoutput) cout << "##" << FILENAME << "##" << METHODNAME << ": Appending to " << output_name << "...";
		////myFile = fopen(output_name.c_str(), "a");
	////}
	////cout.flush();
	////for (unsigned int i = 0; i < data->size(); i++){
		////for (unsigned int j = 0; j < data->at(i).size(); j++){
			////for (unsigned int k = 0; k < data->at(i)[j].size(); k++){
				////fprintf(myFile, "%.6E ", data->at(i)[j][k]);
			////}
		////}
		////fprintf(myFile, "\n");
	////}
	////fclose (myFile);
	////if (cmdoutput) cout << "DONE!" << endl;
	////cout.flush();
////}

////void outputDataAscii(vector<vector<vector<float> > >  * data, string output_name){
	////string METHODNAME = "outputDataAscii";

	////int numBL = (&(data->at(0)))->size();
	////int numFreq = data->size();
	////char buffer[2 * output_name.size()];
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing "<< output_name << " visibilities data...";
	////cout.flush();
	////FILE * myFile;
	////sprintf(buffer, output_name.c_str());
	////myFile = fopen(buffer,"w");
	////for(int f = 0; f < numFreq; f++){
		 ////for(int b = 0; b < numBL; b++){
		 	////if ( b != ( numBL - 1 ) ) {
				////fprintf(myFile, "%.5E %.5E ", ((&(data->at(f)))->at(b))[0], ((&(data->at(f)))->at(b))[1]);
			////} else{
				////fprintf(myFile, "%.5E %.5E", ((&(data->at(f)))->at(b))[0], ((&(data->at(f)))->at(b))[1]);
			////}
		 ////}
	 	 ////fprintf(myFile, "\n");
	////}
	////fclose (myFile);
	////cout << "DONE!" << endl;
////}


////bool outputCalpar(vector<vector<vector<vector<float> > > > * data, string outputfilename, bool in_degree, int nAnt){// outputs binary calpar file. in_degree means if phase calpar is in degree, default true
	////string METHODNAME = "outputCalpar";

	////int numPol = data->size();
	////int numInt = (data->at(0)).size();
	////int numFreq = (data->at(0))[0].size();
	////int numBl = (data->at(0))[0][0].size();//NOT number of baselines but rather the number of entries for each frequency slice's calpar, usually 3 + 2*nant + 2*nUBL
	////if((not in_degree) and (nAnt <= 0 or nAnt >= (numBl - 3) / 2)){
		////cout << "##" << FILENAME << "##" << METHODNAME << ": ERROR: phasecalpar not in degree indicated, but number of antenna " << nAnt << " is incompatible with calpar array length of " << numBl << " Exiting!!" << endl;
		////return false;
	////}

	////cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing calpar: " << outputfilename << " of dimensions " << numPol << " by " <<  numInt << " by " <<  numFreq << " by " <<  numBl << endl;
	////union float_char{
		////float f;
		////unsigned char g[sizeof(float)];
	////};
	////float_char to_file;

	////FILE * myFile = fopen(outputfilename.c_str(), "wb");
	////if(myFile == NULL){
		////cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Outputing " << outputfilename << " FAILED! Check if the path directory exists!" << endl;
		////return false;
	////}
	////if(in_degree){
		////for(int p = 0; p < numPol; p++){
			////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing polarization: " << p << '\n';
			////for(int t = 0; t < numInt; t++){
				////for(int f = 0;f < numFreq; f++){
					////for(int b = 0; b < numBl; b++){
						////to_file.f = (data->at(p))[t][f][b];
						////fwrite(to_file.g, sizeof(float), 1, myFile);
						//////cout << p << " " << t << " " << f<< " " << b <<<< endl;
					////}
				////}
			////}
		////}
	////} else {
		////for(int p = 0; p < numPol; p++){
			////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing polarization: " << p << '\n';
			////for(int t = 0; t < numInt; t++){
				////for(int f = 0;f < numFreq; f++){
					////for(int b = 0; b < numBl; b++){
						////if(b >= 3 + nAnt and b < 3 + 2 * nAnt){
							////to_file.f = (data->at(p))[t][f][b] * 180 / PI;
						////}else{
							////to_file.f = (data->at(p))[t][f][b];
						////}
						////fwrite(to_file.g, sizeof(float), 1, myFile);
						//////cout << p << " " << t << " " << f<< " " << b <<<< endl;
					////}
				////}
			////}
		////}
	////}

	////fclose(myFile);
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Done outputing " << outputfilename << "." << endl;
	////return true;
////}

////bool outputCalparSP(vector<vector<vector<float> > > * data, string outputfilename, bool in_degree, int nAnt){// outputs binary calpar file in log10 and degree. in_degree means if phase calpar is in degree, default true
	////string METHODNAME = "outputCalparSP";

	////int numInt = data->size();
	////int numFreq = (data->at(0)).size();
	////int numBl = (data->at(0))[0].size();//NOT number of baselines but rather the number of entries for each frequency slice's calpar, usually 3 + 2*nant + 2*nUBL
	////if((not in_degree) and (nAnt <= 0 or nAnt >= (numBl - 3) / 2)){
		////cout << "##" << FILENAME << "##" << METHODNAME << ": ERROR: phasecalpar not in degree indicated, but number of antenna " << nAnt << " is incompatible with calpar array length of " << numBl << " Exiting!!" << endl;
		////return false;
	////}

	////cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing calpar: " << outputfilename << " of dimensions " <<  numInt << " by " <<  numFreq << " by " <<  numBl << endl;
	////union float_char{
		////float f;
		////unsigned char g[sizeof(float)];
	////};
	////float_char to_file;

	////FILE * myFile = fopen(outputfilename.c_str(), "wb");
	////if(myFile == NULL){
		////cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Outputing " << outputfilename << " FAILED! Check if the path directory exists!" << endl;
		////return false;
	////}
	////if(in_degree){
		////for(int t = 0; t < numInt; t++){
			////for(int f = 0; f < numFreq; f++){
				////for(int b = 0; b < numBl; b++){
					////to_file.f = (data->at(t))[f][b];
					////fwrite(to_file.g, sizeof(float), 1, myFile);
					//////cout << p << " " << t << " " << f<< " " << b <<<< endl;
				////}
			////}
		////}
	////} else {
		////for(int t = 0; t < numInt; t++){
			////for(int f = 0;f < numFreq; f++){
				////for(int b = 0; b < numBl; b++){
					////if(b >= 3 + nAnt and b < 3 + 2 * nAnt){
						////to_file.f = (data->at(t))[f][b] * 180 / PI;
					////}else{
						////to_file.f = (data->at(t))[f][b];
					////}
					////fwrite(to_file.g, sizeof(float), 1, myFile);
					//////cout << p << " " << t << " " << f<< " " << b <<<< endl;
				////}
			////}
		////}

	////}

	////fclose(myFile);
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Done outputing " << outputfilename << "." << endl;
	////return true;
////}


////bool outputData(vector<vector<vector<vector<vector<float> > > > > * data, string outputfilename){// outputs binary file. Can only be used on odf w/ ~<300 slices. It's recommended to use outputDataLarge under all circumstances
	////string METHODNAME = "outputData";
	////int numPol = data->size();
	////int numInt = (data->at(0)).size();
	////int numFreq = (data->at(0))[0].size();
	////int numBl = (data->at(0))[0][0].size();
	////int numRI = (data->at(0))[0][0][0].size();
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing " << outputfilename << " of dimensions " << numPol << " by " <<  numInt << " by " <<  numFreq << " by " <<  numBl << " by " << numRI << endl;
	////union float_char{
		////float f;
		////unsigned char g[sizeof(float)];
	////};
	////float_char to_file;

	////FILE * myFile = fopen(outputfilename.c_str(), "wb");
	////if(myFile == NULL){
		////cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Outputing " << outputfilename << " FAILED! Check if the path directory exists!" << endl;
		////return false;
	////}
	////for(int p = 0; p < numPol; p++){
		////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing polarization: " << p << '\n';
		////for(int t = 0; t < numInt; t++){
			////for(int f = 0;f < numFreq; f++){
				////for(int b = 0; b < numBl; b++){
					////for(int c = 0; c < numRI; c++){
						////to_file.f = (data->at(p))[t][f][b][c];
						////fwrite(to_file.g, sizeof(float), 1, myFile);
						//////cout << p << " " << t << " " << f<< " " << b << " " << c << endl;
					////}
				////}
			////}
		////}
	////}

	////fclose(myFile);
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Done outputing " << outputfilename << "." << endl;
	////return true;
////}

////bool outputDataLarge(vector<vector<vector<vector<float> > > > * data, string outputfilename){// outputs binary file. Can be used on odf w/ many slices. It's recommended to use outputDataLarge under all circumstances

	////string METHODNAME = "outputDataLarge";
	////int numPol = data->size();
	////int numInt = (data->at(0)).size();
	////int numFreq = (data->at(0))[0].size();
	////int numBl = (data->at(0))[0][0].size() / 2;
	////int numRI = 2;
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Outputing " << outputfilename << " of dimensions " << numPol << " by " <<  numInt << " by " <<  numFreq << " by " <<  numBl << " by " << numRI << endl;
	////union float_char{
		////float f;
		////unsigned char g[sizeof(float)];
	////};
	////float_char to_file;

	////FILE * myFile = fopen(outputfilename.c_str(), "wb");
	////if(myFile == NULL){
		////cout << "##" << FILENAME << "##" << METHODNAME << "!!!!FATAL I/O ERROR!!!!!!!!!: Outputing " << outputfilename << " FAILED! Check if the path directory exists!" << endl;
		////return false;
	////}
	////for(int p = 0; p < numPol; p++){
		//////cout << "##" << FILENAME << "##" << METHODNAME << ": Writing polarization: " << p << '\n';
		////for(int t = 0; t < numInt; t++){
			////for(int f = 0;f < numFreq; f++){
				////for(int b = 0; b < (data->at(0))[0][0].size() ; b++){

						////to_file.f = (data->at(p))[t][f][b];
						////fwrite(to_file.g, sizeof(float), 1, myFile);
						//////cout << p << " " << t << " " << f<< " " << b << " " << c << endl;

				////}
			////}
		////}
	////}

	////fclose(myFile);
	////cout << "##" << FILENAME << "##" << METHODNAME << ": Done outputing " << outputfilename << "." << endl;
	////return true;
////}


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


//bool createAmatrix(vector<vector<int> > *receiverAmatrix, vector<vector<float> > *antloc){
	//string METHODNAME = "createAmatrix";
	//int nAnt = antloc->size();
	//int nUBL = countUBL(antloc);
	//int nCrossBL = nAnt * (nAnt - 1) / 2;
	//vector<vector<float> > UBL(nUBL, vector<float> (3, 0));
	//computeUBL(antloc, &UBL);
	//if (receiverAmatrix->size() != (nCrossBL + 1) or (receiverAmatrix->at(0)).size() != (nAnt + nUBL)){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL INPUT MISMATCH! Dimensions of the receiver Amatrix are " << receiverAmatrix->size() << " by " << receiverAmatrix->size() << ", NOT expected " << (nCrossBL + 1) << " by " << (nAnt + nUBL) << " for " << nAnt << " antennae and " << nUBL << " unique baselines. Returned!!!" << endl;
		//return false;
	//}
	//for(int bl = 0; bl < nCrossBL; bl ++){
		//int ant1 = get2DBLCross(bl, nAnt)[0];
		//int ant2 = get2DBLCross(bl, nAnt)[1];
		//int UBLindex = abs(indexUBL(&UBL, getBL(ant1, ant2, antloc)))-1;
		//if(UBLindex < 0 or UBLindex >= nUBL){
			//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL ERROR! No corresponding unique baseline found for antenna pair" << ant1 << " " << ant2 << "! Amatrix Failed to generate!" << endl;
			//return false;
		//}
		//for(int col = 0; col < nAnt; col++){
			//if (col == ant1 or col == ant2){
				//(receiverAmatrix->at(bl))[col] = 1;
			//} else {
				//(receiverAmatrix->at(bl))[col] = 0;
			//}
		//}
		//for(int col = nAnt; col < nAnt + nUBL; col++){
			//if (col - nAnt == UBLindex){
				//(receiverAmatrix->at(bl))[col] = 1;
			//} else {
				//(receiverAmatrix->at(bl))[col] = 0;
			//}
		//}

	//}

	//int bl = nCrossBL;
	//for(int col = 0; col < nAnt; col++)	(receiverAmatrix->at(bl))[col] = 1;
	//for(int col = nAnt; col < nAnt + nUBL; col++) (receiverAmatrix->at(bl))[col] = 0;
	//return true;
//}

//bool createBmatrix(vector<vector<int> > *receiverBmatrix, vector<vector<float> > *antloc){
	//string METHODNAME = "createBmatrix";
	//int nAnt = antloc->size();
	//int nUBL = countUBL(antloc);
	//int nCrossBL = nAnt * (nAnt - 1) / 2;
	//vector<vector<float> > UBL(nUBL, vector<float> (3, 0));
	//computeUBL(antloc, &UBL);
	//if (receiverBmatrix->size() != (nCrossBL + 3) or (receiverBmatrix->at(0)).size() != (nAnt + nUBL)){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL INPUT MISMATCH! Dimensions of the receiver Amatrix are " << receiverBmatrix->size() << " by " << receiverBmatrix->size() << ", NOT expected " << (nCrossBL + 3) << " by " << (nAnt + nUBL) << " for " << nAnt << " antennae and " << nUBL << " unique baselines. Returned!!!" << endl;
		//return false;
	//}
	//for(int bl = 0; bl < nCrossBL; bl ++){
		//int ant1 = get2DBLCross(bl, nAnt)[0];
		//int ant2 = get2DBLCross(bl, nAnt)[1];
		//int UBLindex = indexUBL(&UBL, getBL(ant1, ant2, antloc));
		//if(UBLindex == 0){
			//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL ERROR! No corresponding unique baseline found for antenna pair" << ant1 << " " << ant2 << "! Bmatrix Failed to generate!" << endl;
			//return false;
		//}
		//int reversed = abs(UBLindex)/UBLindex;
		//for(int col = 0; col < nAnt; col++){
			//if (col == ant1){
				//(receiverBmatrix->at(bl))[col] = -reversed;
			//} else if (col == ant2){
				//(receiverBmatrix->at(bl))[col] = reversed;
			//} else {
				//(receiverBmatrix->at(bl))[col] = 0;
			//}
		//}
		//for(int col = nAnt; col < nAnt + nUBL; col++){
			//if (col - nAnt == abs(UBLindex)-1){
				//(receiverBmatrix->at(bl))[col] = 1;
			//} else {
				//(receiverBmatrix->at(bl))[col] = 0;
			//}
		//}

	//}

	//int bl = nCrossBL;
	//for(int col = 0; col < nAnt; col++)	(receiverBmatrix->at(bl))[col] = 1;
	//for(int col = nAnt; col < nAnt + nUBL; col++) (receiverBmatrix->at(bl))[col] = 0;
	//bl = nCrossBL + 1;
	//for(int col = 0; col < nAnt; col++)	(receiverBmatrix->at(bl))[col] = 0;
	//for(int col = nAnt; col < nAnt + nUBL; col++) (receiverBmatrix->at(bl))[col] = 0;
	//(receiverBmatrix->at(bl))[nAnt] = 1;
	//bl = nCrossBL + 2;
	//for(int col = 0; col < nAnt; col++)	(receiverBmatrix->at(bl))[col] = 0;
	//for(int col = nAnt; col < nAnt + nUBL; col++) (receiverBmatrix->at(bl))[col] = 0;
	//(receiverBmatrix->at(bl))[nAnt+1] = 1;
	//return true;
//}

//bool findReversedBaselines(vector<int> *receiverList, vector<vector<float> > *antloc){
	//string METHODNAME = "findReversedBaselines";
	//int nAnt = antloc->size();
	//int nUBL = countUBL(antloc);
	//int nCrossBL = nAnt * (nAnt - 1) / 2;
	//vector<vector<float> > UBL(nUBL, vector<float> (3, 0));
	//computeUBL(antloc, &UBL);
	//if (receiverList->size() != nCrossBL){
		//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL INPUT MISMATCH! Dimension of the receiver list is " << receiverList->size() << ", NOT expected " << nCrossBL << " for " << nAnt << " antennae. Returned!!!" << endl;
		//return false;
	//}
	//for(int bl = 0; bl < nCrossBL; bl++){
		//int ant1 = get2DBLCross(bl, nAnt)[0];
		//int ant2 = get2DBLCross(bl, nAnt)[1];
		//int UBLindex = indexUBL(&UBL, getBL(ant1, ant2, antloc));
		//if(UBLindex == 0){
			//cout << "#!!#" << FILENAME << "#!!#" << METHODNAME << ": FATAL ERROR! No corresponding unique baseline found for antenna pair" << ant1 << " " << ant2 << "! ReversedBaselines Failed to generate!" << endl;
			//return false;
		//}
		//if(UBLindex > 0){
			//receiverList->at(bl) = 0;
		//} else {
			//receiverList->at(bl) = 1;
		//}
	//}
	//return true;
//}

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

/*float mode (vector<float> list){//Didn't finish!!!!
	sort(list.begin(), list.end());
	int size = list.size();
	float min = list[0];
	float max = list[ size - 1 ];
	if ( max == min ) return max;
	int maxCount = 1;
	int maxDegeneracy = 0;


	float binSize =  ( max - min ) / size;
	vector<int> binCount(size);//divide the range of data into size number of bins, and this vector book keeps how many elements in each bin
	for ( int i = 0; i < size; i ++){
		int interval = (int)floor( ( list[i] - min ) / binSize ) % size;
		binCount[interval] ++;
	}

	for ( int i = 0; i < size; i ++){
		if ( binCount[i] > maxCount ) {
			maxCount = binCount[i];
			maxDegeneracy = 1;
		} else if (binCount[i] == maxCount){
			maxDegeneracy ++;
		}
	}
	int index = floor( list.size() / 2 );

	return list[index];
}*/

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
		calpar->at(0) = pow(norm(additiveout), 2);
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
