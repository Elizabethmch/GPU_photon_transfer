#include "ManagePhotonAbsorption.h"
#include <cuda_runtime.h>
#include <iostream>
#include <curand_kernel.h>
#include <cooperative_groups.h>
#include <vector>
#include "MyCudaToolkit.h"


using namespace std;
namespace cg = cooperative_groups;

__device__ unsigned int reduce_sum(long in, cg::thread_block cta)
{
	extern __shared__ long sdata[];
	// Perform first level of reduction:
	// - Write to shared memory
	unsigned int ltid = threadIdx.x;
	sdata[ltid] = in;
	if (in == -1) printf("error: countflag not initialzed\n");

	cg::sync(cta);
	// Do reduction in shared mem
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (ltid < s)
		{
			sdata[ltid] += sdata[ltid + s];
		}

		cg::sync(cta);
	}


	return sdata[0];
}

/// <summary>
/// kernel function for simulation
/// For every depth, do the loop #(photonnum/(blockdim*griddim)) times, calculate photon absorption, store results into parameter count;
/// Then continue the loop for next depth.
/// </summary>
/// <param name="count"></param>	Array for returning results
/// <param name="depthbin_ary"></param>	Incident photon depths (stored as bin id of the depth)
/// <param name="depthsize"></param>	Total number of photon sets
/// <param name="anglebinsize"></param>	Total bin number of angle
/// <param name="lut"></param>	Pointer to the lut data
/// <param name="photon_num"></param>	Photon number for each sets.
/// <param name="rndstates"></param>	Random number generater
/// <param name="seed"></param>	Random number seed
/// <returns></returns>
__global__ void SimulatePhotonAbsorption(long* count, int* depthbin_ary, int depthsize, int anglebinsize, double* lut, 
										long* photon_num, curandStateXORWOW_t* rndstates, unsigned int seed) {
	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();

	// Determine thread ID
	int bid = blockIdx.x;
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int local_tid = threadIdx.x;
	int step = blockDim.x * gridDim.x;

	// Initialise the RNG
	curand_init(seed, tid, 0, &rndstates[tid]);
	curandState localState = rndstates[tid];

	//begin simulation for each depth
	for (int depthid = 0; depthid < depthsize; depthid++) {
		unsigned int countflag = 0;
		for (unsigned photonid = tid; photonid < photon_num[depthid]; photonid += step) {
			int anglebin = (int)(anglebinsize * curand_uniform_double(&localState));
			double prob = lut[depthbin_ary[depthid] * anglebinsize + anglebin];
			double rndm = curand_uniform_double(&localState);
			if (rndm < prob)	countflag++;
		}
		countflag = reduce_sum(countflag, cta);
		if (threadIdx.x == 0) {
			count[bid + depthid * gridDim.x] = countflag;
		}
	}
}


/**************************************************
* Calculate number of photon absorbed for given depth & incident photon number
* ( If rndmseed==0, rndmseed=time(nullptr) )
**************************************************/
vector<long> ManagePhotonAbsorption::getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, unsigned int rndmseed) {
	unsigned int m_device = 0;
	struct cudaDeviceProp     deviceProperties;
	CHECK(cudaGetDeviceProperties(&deviceProperties, m_device));
	CHECK(cudaSetDevice(m_device));

	vector<int> lut_size = look_up_table->getLUTSize();
	//Distribute depth array into given depth bins
	vector<int> depth_in_bin;

	//The last entry for lut_size is num of bins for depth
	int depthbinnum = lut_size[0];
	for (int i = 0; i < depth.size(); i++) {
		int tmp = (int)((depth[i] - min_depth) / (max_depth - min_depth) * depthbinnum);
		depth_in_bin.push_back(tmp);
	}

	//Memory allocation in GPU for photonnum and depth_in_bin

	long* d_photonnum;
	CHECK(cudaMalloc((void**)&d_photonnum, depth.size() * sizeof(long)));
	CHECK(cudaMemcpy((void*)d_photonnum, (void*)&(incident_photon_num[0]), depth.size() * sizeof(long), cudaMemcpyHostToDevice));

	int* d_depth_in_bin;
	CHECK(cudaMalloc((void**)&d_depth_in_bin, depth.size() * sizeof(int)));
	CHECK(cudaMemcpy((void*)d_depth_in_bin, (void*)&(depth_in_bin[0]), depth.size() * sizeof(int), cudaMemcpyHostToDevice));

	
	//DepthCnt_ary: count absorbed photon for each depth	
	long* h_DepthCnt_ary, * d_DepthCnt_ary;

	//kernel function setup ( numSMs==30 )
	dim3 block, grid;
	block.x = threadBlockSize;
	grid.x = 0;
	grid.x = (incident_photon_num[0] - 1) / block.x + 1;
	unsigned int blocksPerSM = 10;
	unsigned int numSMs = deviceProperties.multiProcessorCount;
	while (grid.x > 2 * blocksPerSM * numSMs){
		grid.x >>= 1;
	}

	size_t count_array_size = grid.x * depth.size() * sizeof(long);
	CHECK(cudaMalloc((void**)&d_DepthCnt_ary, count_array_size));
	h_DepthCnt_ary = (long*)malloc(count_array_size);


	//Random number simulation
	if (totalThreadNum < grid.x * block.x) {
		totalThreadNum = grid.x * block.x;
		CHECK(cudaMalloc((void**)&states, sizeof(curandStateXORWOW_t) * block.x * grid.x));
	}


	if (rndmseed == 0) rndmseed = time(nullptr);
	SimulatePhotonAbsorption << <grid, block, block.x * sizeof(long) >> > (d_DepthCnt_ary, d_depth_in_bin, depth.size(), lut_size[1], d_lut, d_photonnum, states, rndmseed );
	cudaError_t cudaStatus = cudaGetLastError();
	CHECK(cudaStatus);


	//collect results
	CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, count_array_size, cudaMemcpyDeviceToHost));

	vector<long> absorbcnt;
	for (int i = 0; i < depth.size(); i++) {
		int tmpcnt = 0;
		for (int j = 0; j < grid.x; j++) {
			tmpcnt += h_DepthCnt_ary[i * grid.x + j];
		}
		absorbcnt.push_back(tmpcnt);
	}

	//Free storage
	CHECK(cudaFree(d_DepthCnt_ary));
	free(h_DepthCnt_ary);
	CHECK(cudaFree(d_depth_in_bin));

	return absorbcnt;
}

ManagePhotonAbsorption::ManagePhotonAbsorption(LUT* lut, double maxdepth, double mindepth, int blocksize) : look_up_table(lut),
max_depth(maxdepth), min_depth(mindepth), threadBlockSize(blocksize), totalThreadNum(0){
	const vector<double>* h_lut_ptr = look_up_table->getLUTAddress();
	vector<int> lut_size = look_up_table->getLUTSize();
	int lut_total_size = 1;
	for (int i = 0; i < lut_size.size(); i++)	lut_total_size *= lut_size[i];

	CHECK(cudaMalloc((void**)&d_lut, lut_total_size * sizeof(double)));
	CHECK(cudaMemcpy((void*)d_lut, &((*h_lut_ptr)[0]), lut_total_size * sizeof(double), cudaMemcpyHostToDevice));
}


