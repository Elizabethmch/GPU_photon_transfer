#include "ManagePhotonAbsorption.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "MyCudaToolkit.h"

__global__ void SimulatePhotonAbsorption(long* absorb_num, int depth, double* lut) {

}


/**************************************************
* Calculate number of photon absorbed for given depth & incident photon number
**************************************************/
long ManagePhotonAbsorption::getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num) {
	
	vector<int> lut_size = look_up_table->getLUTSize();
	int lut_total_size = 1;
	for (int i = 0; i < lut_size.size(); i++)	lut_total_size *= lut_size[i];
	//Distribute depth array into given bins
	vector<int> depth_in_bin;

	//The last entry for lut_size is num of bins for depth
	int depthbinnum = lut_size[lut_size.size() - 1];
	for (int i = 0; i < depth.size(); i++) {
		int tmp =(int)( (depth[i] - min_depth) / (max_depth - min_depth) * depthbinnum );
		depth_in_bin.push_back(tmp);
	}

	//Look up table menmory allocation
	double* d_lut;
	const vector<double>* h_lut_ptr = look_up_table->getLUTAddress();
	CHECK( cudaMalloc((void**)&d_lut, lut_total_size * sizeof( double ) ) );
	CHECK( cudaMemcpy((void*)d_lut, &((*h_lut_ptr)[0]), lut_total_size * sizeof(double), cudaMemcpyHostToDevice));


	//Simulate photon absorption for each depth
	for (int dptid = 0; dptid < depth_in_bin.size(); dptid++) {
		//Photon-absorption state array
		long* d_absorb_num;
		long* h_absorb_num;
		CHECK(cudaMalloc((void**)&d_absorb_num, incident_photon_num[dptid] * sizeof(long)));

	}

	return 1;
}
