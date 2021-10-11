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
**************************************************/
vector<long> ManagePhotonAbsorption::getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, unsigned int rndmseed) {
	//cudaDeviceReset();
	//double iStart, iElaps;
	//iStart = cpuSecond();

	unsigned int m_device = 0;
	struct cudaDeviceProp     deviceProperties;
	CHECK(cudaGetDeviceProperties(&deviceProperties, m_device));
	CHECK(cudaSetDevice(m_device));

	vector<int> lut_size = look_up_table->getLUTSize();
	int lut_total_size = 1;
	for (int i = 0; i < lut_size.size(); i++)	lut_total_size *= lut_size[i];
	//Distribute depth array into given bins
	vector<int> depth_in_bin;

	//The last entry for lut_size is num of bins for depth
	int depthbinnum = lut_size[0];
	for (int i = 0; i < depth.size(); i++) {
		int tmp = (int)((depth[i] - min_depth) / (max_depth - min_depth) * depthbinnum);
		depth_in_bin.push_back(tmp);
		//printf("depth id: %d\n", tmp);
	}

	//Look up table menmory allocation
	double* d_lut;
	const vector<double>* h_lut_ptr = look_up_table->getLUTAddress();
	CHECK(cudaMalloc((void**)&d_lut, lut_total_size * sizeof(double)));
	CHECK(cudaMemcpy((void*)d_lut, &((*h_lut_ptr)[0]), lut_total_size * sizeof(double), cudaMemcpyHostToDevice));

	long* d_photonnum;
	CHECK(cudaMalloc((void**)&d_photonnum, depth.size() * sizeof(long)));
	CHECK(cudaMemcpy((void*)d_photonnum, (void*)&(incident_photon_num[0]), depth.size() * sizeof(long), cudaMemcpyHostToDevice));

	int* d_depth_in_bin;
	CHECK(cudaMalloc((void**)&d_depth_in_bin, depth.size() * sizeof(int)));
	CHECK(cudaMemcpy((void*)d_depth_in_bin, (void*)&(depth_in_bin[0]), depth.size() * sizeof(int), cudaMemcpyHostToDevice));


	/*
	* DepthCnt_ary: count absorbed photon for each depth
	* GridDivide: all photons are divided into several blocks. The division points (in grid id) are stored in this array
	* TailPhotonNum: photons in each depth are divided into an integer number of blocks. Photon num in the last block is stored in this array
	*/
	long* h_DepthCnt_ary, * d_DepthCnt_ary;
	//int* h_GridDivide, * d_GridDivide;
	//int* h_TailPhotonNum, * d_TailPhotonNum;
	//int* h_GridDivideId, * d_GridDibideId;
	//h_GridDivide = (int*)malloc((depth.size() + 1) * sizeof(int));
	//CHECK(cudaMalloc((void**)&d_GridDivide, (depth.size() + 1) * sizeof(int)));
	//h_TailPhotonNum = (int*)malloc(depth.size() * sizeof(int));
	//CHECK(cudaMalloc((void**)&d_TailPhotonNum, depth.size() * sizeof(int)));

	dim3 block, grid;
	block.x = threadBlockSize;
	grid.x = 0;
	grid.x = (1e7 - 1) / block.x + 1;
	unsigned int blocksPerSM = 10;
	unsigned int numSMs = deviceProperties.multiProcessorCount;
	while (grid.x > 2 * blocksPerSM * numSMs)
	{
		grid.x >>= 1;
	}
	//cout << "multiprocessor number: " << numSMs << endl;
	//printf("With grid size:%i, block size:%i\n", grid.x, block.x);

	//h_GridDivide[0] = 0;
	//for (int dptid = 0; dptid < depth.size(); dptid++) {
	//	//Allocate memory and size for gpu, grid and so on
	//	int tmp = (incident_photon_num[dptid] - 1) / block.x + 1;
	//	grid.x += tmp;
	//	h_GridDivide[dptid + 1] = grid.x;
	//	h_TailPhotonNum[dptid] = incident_photon_num[dptid] % block.x;
	//	if (h_TailPhotonNum[dptid] == 0) h_TailPhotonNum[dptid] = block.x;
	//	//printf("In dptid : % d: photon num: %d, block.x: %d, tail photon num : %d\n", dptid, incident_photon_num[dptid], block.x,  h_TailPhotonNum[dptid]);
	//}
	////for (int i = 0; i < depth.size() + 1; i++)	printf("Divide point %d: %d\n", i, h_GridDivide[i]);

	size_t count_array_size = grid.x * depth.size() * sizeof(long);
	CHECK(cudaMalloc((void**)&d_DepthCnt_ary, count_array_size));
	h_DepthCnt_ary = (long*)malloc(count_array_size);




	//cout << h_DepthCnt_ary[grid.x-1] << endl;
	//CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, grid.x * sizeof(long), cudaMemcpyDeviceToHost));
	//cout << h_DepthCnt_ary[grid.x-1] << endl;


	//CHECK(cudaMemcpy((void*)d_GridDivide, (void*)h_GridDivide, (depth.size() + 1) * sizeof(int), cudaMemcpyHostToDevice));
	//CHECK(cudaMemcpy((void*)d_TailPhotonNum, (void*)h_TailPhotonNum, depth.size() * sizeof(int), cudaMemcpyHostToDevice));

	//Random number simulation
	curandStateXORWOW_t* states;
	CHECK(cudaMalloc((void**)&states, sizeof(curandStateXORWOW_t) * block.x * grid.x));
	//size_t d_DepthCnt_ary_size = grid.x * sizeof(long);
	//cout << "grid.x: " << grid.x << " size of long: " << sizeof(long) << endl;
	////cout << "depth cnt ary size: " << d_DepthCnt_ary_size << endl;
	//printf("sizeof curandState:%d, sizeof curandStateXORWOW_t %d\n", sizeof(curandState), sizeof(curandStateXORWOW));
	//long long statesize = (long long)sizeof(curandStateXORWOW_t) * block.x * grid.x;
	//cout << "statesize: " << statesize << endl;
	//long a = sizeof(curandStateXORWOW_t);
	//size_t b = a * block.x * grid.x;
	//long c = b / 1024 / 1024;
	//size_t d = 4799692800;
	//printf("With grid size:%i, block size:%i\n", grid.x, block.x);
	//printf("size of 1 rnd states: %d bytes; total rnd states: %d MB\n", sizeof(curandStateXORWOW_t), sizeof(curandStateXORWOW_t) * block.x * grid.x/1024/1024);
	//cout << "a " << a << endl;
	//cout << "b " << b << endl;
	//cout << "c " << c << endl;
	//cout << "d " << d << endl;
	//cout << 97650 * 1024 * 48  << endl;
	//cout << 97650 * 1024 * 48 / 1024 / 1024 << endl;
	//cout << 97650 * 48 / 1024 << endl;
	//cout << h_DepthCnt_ary[grid.x-1] << endl;
	//CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, grid.x * sizeof(long), cudaMemcpyDeviceToHost));
	//cout << h_DepthCnt_ary[grid.x-1] << endl;

	//unsigned int rndmseed = 1234;
	//unsigned int rndmseed = time(nullptr);
	if (rndmseed == 0) rndmseed = time(nullptr);
	SimulatePhotonAbsorption << <grid, block, block.x * sizeof(long) >> > (d_DepthCnt_ary, d_depth_in_bin, depth.size(), lut_size[1], d_lut, d_photonnum, states, rndmseed );
	cudaError_t cudaStatus = cudaGetLastError();
	CHECK(cudaStatus);
	//printf("size of shared memory used per block: %i \n", block.x * sizeof(long));
	//printf("block number: %i\n", grid.x);
	//printf("Total shared memory used: %i kb\n", block.x * grid.x * sizeof(long) / 1024);
	//CHECK(cudaDeviceSynchronize());
	//printf("size of count array: %i kb\n", grid.x * sizeof(long) / 1024);
	CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, count_array_size, cudaMemcpyDeviceToHost));

	//collect information
	vector<long> absorbcnt;
	for (int i = 0; i < depth.size(); i++) {
		/*int tmpcnt = 0;
		for (int j = h_GridDivide[i]; j < h_GridDivide[i + 1]; j++) {
			tmpcnt += h_DepthCnt_ary[j];
		}
		cout << "In depth " << depth[i];
		cout << ", incident photon num:" << incident_photon_num[i] << ", absorbed photon num:" << tmpcnt << endl;
		absorbcnt.push_back(tmpcnt);*/
		int tmpcnt = 0;
		for (int j = 0; j < grid.x; j++) {
			tmpcnt += h_DepthCnt_ary[i * grid.x + j];
		}
		cout << "In depth " << depth[i] << ", depth bin num: "<< depth_in_bin[i];
		cout << ", incident photon num:" << incident_photon_num[i] << ", absorbed photon num:" << tmpcnt << endl;
		absorbcnt.push_back(tmpcnt);
	}


	//CHECK(cudaDeviceSynchronize());
	//iElaps = cpuSecond() - iStart;
	//printf("Time elapsed %f ms\n", iElaps);


	CHECK(cudaFree(d_DepthCnt_ary));	free(h_DepthCnt_ary);
	CHECK(cudaFree(states));
	//CHECK(cudaFree(d_GridDivide));		free(h_GridDivide);
	CHECK(cudaFree(d_depth_in_bin));
	CHECK(cudaFree(d_lut));
	//CHECK(cudaFree(d_TailPhotonNum));	free(h_TailPhotonNum);
	return absorbcnt;
}
