#include <iostream>
#include <fstream>
#include "LUT.h"

using namespace std;

LUT::LUT() : dim(2), lut_name("lut.dat"){
	size.push_back(50);		// Direction bin number
	size.push_back(100);	// Depth bin number
	LoadLUTData();
}

LUT::LUT(int dim, vector<int> dim_size, string LutDataFile) : dim(dim), lut_name(LutDataFile){
	size = dim_size;
	int totalsize = 1;
	for (int i = 0; i < dim; i++) totalsize *= dim_size[i];
	prob_table.reserve( totalsize );
	LoadLUTData(LutDataFile);
}

LUT::~LUT(){
	size.clear();
	prob_table.clear();
}




/*********************
* Load LookUpTable Data From DataFile
*********************/
void LUT::LoadLUTData(string DataFileName){
	cout << "Loading data from "<<DataFileName<<"......" << endl;
	ifstream inFile(DataFileName, ios::in | ios::binary);
	if (!inFile) {
		cout << "error while openning datafile"<< DataFileName << endl;
		return;
	}

	double tmpdata;
	while (inFile.read((char*)&tmpdata, sizeof(double))) {
		//cout << "datas: " << tmpdata << "	";
		prob_table.push_back(tmpdata);
	}
}

void LUT::PrintLUT(){
	for (int i = 0; i < prob_table.size(); i++) {
		cout << prob_table[i] << "\t";
	}
	cout << endl;
}
