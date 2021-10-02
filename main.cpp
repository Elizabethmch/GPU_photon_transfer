#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include "LUT.h"
#include "ManagePhotonAbsorption.h"
#include "MyCudaToolkit.h"
#include "RandomLUT.h"


using namespace std;

void PrintMeanProbOfDepth(LUT* lut, int depthid) {
	const vector<double>* tableptr = lut->getLUTAddress();
	auto sizevec = lut->getLUTSize();
	int anglebinsize = sizevec[0];

	double mean_prob = 0;
	for (int i = 0; i < anglebinsize; i++) {
		mean_prob += (*tableptr)[anglebinsize * depthid + i];
	}
	mean_prob /= anglebinsize;
	printf("Mean probability of depth %d is %f\n", depthid, mean_prob);
}

int main(int argc, char **argv) {
	//PrintRndmLUT();
	//PrintUniformLUT(100);
	//PrintUniformLUT(80);
	//PrintUniformLUT(30);

	vector<int> dim_size; dim_size.push_back(50); dim_size.push_back(100);
	//LUT* table = new LUT(2, dim_size, "lut100.dat");
	LUT* table = new LUT(2, dim_size, "lut80.dat");
	//LUT* table = new LUT(2, dim_size, "lut30.dat");
	//LUT* table = new LUT();
	
	//const vector<double>* test = table->getLUTAddress();
	//table->PrintLUT();


	//incident photon configuration
	vector<double> incident_depth;
	vector<long> incident_photon_num;
	//incident_depth.push_back(0.5);	incident_depth.push_back(1.7);	incident_depth.push_back(20.7);
	//incident_photon_num.push_back(1e7);	incident_photon_num.push_back(1e7);	incident_photon_num.push_back(1e7);
	for (int i = 0; i < 5; i++) {
		incident_depth.push_back(0.5);
		incident_photon_num.push_back(1e7);
	}

	double MAXDPT = 100., MINDPT = 0.;

	//double iStart, iElaps;
	//iStart = cpuSecond();
	ManagePhotonAbsorption* mpa = new ManagePhotonAbsorption(table, MAXDPT, MINDPT);
	//CHECK(cudaDeviceSynchronize());
	//iElaps = cpuSecond() - iStart;
	//printf("Time elapsed (main) %f ms\n", iElaps);

	//PrintMeanProbOfDepth(table, 0);
	//PrintMeanProbOfDepth(table, 1);
	//PrintMeanProbOfDepth(table, 20);


	//simulation start
	mpa->getAbsorbedPhotonNum(incident_depth, incident_photon_num);
	return 1;
}