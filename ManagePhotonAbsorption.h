#ifndef MANAGEPHOTONABSORPTION
#define MANAGEPHOTONABSORPTION
#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "LUT.h"
#include "MyCudaToolkit.h"


using std::vector;

/// <summary>
/// Simulate photon absorption
/// </summary>
class ManagePhotonAbsorption
{
public:
	ManagePhotonAbsorption(LUT* lut, double maxdepth, double mindepth, int blocksize = 512);
	~ManagePhotonAbsorption() {
		CHECK(cudaFree(states));
		CHECK(cudaFree(d_lut));
	}

	vector<long> getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, unsigned int rndmseed = 0);
	void PrintTimeConsume();


private:
	LUT* look_up_table;
	double max_depth;
	double min_depth;
	int threadBlockSize;
	int totalThreadNum;
	double* d_lut;	//look up table data for GPU
	curandStateXORWOW_t* states;
	double kernelStart;
	double kernelElaps;
	double transStart;
	double transElaps;
};
#endif

