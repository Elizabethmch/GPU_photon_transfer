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
	//printf("1");
	sdata[ltid] = in;
	if (in == -1) printf("error: countflag not initialzed\n");

	cg::sync(cta);
	//printf("2");
		// Do reduction in shared mem
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (ltid < s)
		{
			//printf("%d\n",ltid + s);
			sdata[ltid] += sdata[ltid + s];
		}

		cg::sync(cta);
	}


	return sdata[0];
}

__global__ void SimulatePhotonAbsorption(long* count, int *depthbin_ary, int depthsize, int anglebinsize, double* lut, int* grid_divide_points, 
										 int* tail_photon_num, curandStateXORWOW_t* rndstates, unsigned long long seed) {
	cg::thread_block cta = cg::this_thread_block();
	int bid = blockIdx.x;
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int local_tid = threadIdx.x;
	int depthid = 0;
	long countflag = -1;
	//if (bid == gridDim.x - 1)	printf("enter griddim-1: %d \n", bid);

	for (int i = 0; i < depthsize; i++) {
		if (bid>grid_divide_points[i] && bid<grid_divide_points[i + 1]) {
			depthid = depthbin_ary[i];
			if (bid == grid_divide_points[i + 1] - 1 && local_tid >= tail_photon_num[i]) {
				countflag = 0;
			}
			break;
		}
	}
	//if (tid > size) return;
	//printf("1: %d, %d\n", count_flag, threadIdx.x);
	if (countflag != 0) {
		curand_init(seed, tid, 0, &rndstates[tid]);
		int anglebin = (int)(anglebinsize * curand_uniform_double(&rndstates[tid]));
		//int anglebin = 0;
		double prob = lut[depthid * anglebinsize + anglebin];
		double rndm = curand_uniform_double(&rndstates[tid]);
		//double rndm = 0;
		countflag = rndm < prob ? 1 : 0;
	}
//countflag = 1;

	countflag = reduce_sum(countflag, cta);
	//__syncthreads();

	if (threadIdx.x == 0) {
		count[bid] = countflag;
		//if(countflag != 1024) printf("bid: %d, count: %d\n", bid, countflag);
		//if (bid == gridDim.x - 1)	printf("in griddim-1: %d \n",countflag);
	}
}


/**************************************************
* Calculate number of photon absorbed for given depth & incident photon number
**************************************************/
vector<long> ManagePhotonAbsorption::getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, int threadBlockSize) {
	double iStart, iElaps;
	iStart = cpuSecond();
	
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
		//printf("depth id: %d\n", tmp);
	}

	//Look up table menmory allocation
	double* d_lut;
	const vector<double>* h_lut_ptr = look_up_table->getLUTAddress();
	CHECK( cudaMalloc((void**)&d_lut, lut_total_size * sizeof( double ) ) );
	CHECK( cudaMemcpy((void*)d_lut, &((*h_lut_ptr)[0]), lut_total_size * sizeof(double), cudaMemcpyHostToDevice));

	//long* d_photonnum;
	//CHECK(cudaMalloc((void**)&d_photonnum, depth.size() * sizeof(long)));
	//CHECK(cudaMemcpy((void*)d_photonnum, (void*)&(incident_photon_num[0]), depth.size() * sizeof(long), cudaMemcpyHostToDevice));

	int* d_depth_in_bin;
	CHECK(cudaMalloc((void**)&d_depth_in_bin, depth.size() * sizeof(int)));
	CHECK(cudaMemcpy((void*)d_depth_in_bin, (void*)&(depth_in_bin[0]), depth.size() * sizeof(int), cudaMemcpyHostToDevice));


	//DepthCnt_ary: count absorbed photon for each depth
	//GridDivide: all photons are divided into several blocks. The division points (in grid id) are stored in this array
	//TailPhotonNum: photons in each depth are divided into an integer number of blocks. Photon num in the last block is stored in this array
	long* h_DepthCnt_ary, * d_DepthCnt_ary;
	int* h_GridDivide, * d_GridDivide;
	int* h_TailPhotonNum, * d_TailPhotonNum;
	//int* h_GridDivideId, * d_GridDibideId;
	h_GridDivide = (int*)malloc((depth.size() + 1) *sizeof(int));
	CHECK(cudaMalloc((void**)&d_GridDivide, (depth.size() + 1) * sizeof(int)));
	h_TailPhotonNum = (int*)malloc(depth.size() * sizeof(int));
	CHECK(cudaMalloc((void**)&d_TailPhotonNum, depth.size() * sizeof(int)));

	dim3 block, grid;
	block.x = threadBlockSize;
	grid.x = 0;

	h_GridDivide[0] = 0;
	for (int dptid = 0; dptid < depth.size(); dptid++) {
		//Allocate memory and size for gpu, grid and so on
		int tmp = (incident_photon_num[dptid] - 1) / block.x + 1;
		grid.x += tmp;
		h_GridDivide[dptid+1] = grid.x;
		h_TailPhotonNum[dptid] = incident_photon_num[dptid] % block.x;
		//printf("In dptid : % d: photon num: %d, block.x: %d, tail photon num : %d\n", dptid, incident_photon_num[dptid], block.x,  h_TailPhotonNum[dptid]);
	}
	printf("With grid size:%i, block size:%i\n", grid.x, block.x);
	CHECK(cudaMalloc((void**)&d_DepthCnt_ary, grid.x * sizeof(long)));
	h_DepthCnt_ary = (long*)malloc(grid.x * sizeof(long));


	//cout << h_DepthCnt_ary[grid.x-1] << endl;
	//CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, grid.x * sizeof(long), cudaMemcpyDeviceToHost));
	//cout << h_DepthCnt_ary[grid.x-1] << endl;


	CHECK(cudaMemcpy((void*)d_GridDivide, (void*)h_GridDivide, (depth.size() + 1) * sizeof(int) , cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy((void*)d_TailPhotonNum, (void*)h_TailPhotonNum, depth.size() * sizeof(int) , cudaMemcpyHostToDevice));

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

	SimulatePhotonAbsorption << <grid, block, block.x * sizeof(long) >> > (d_DepthCnt_ary, d_depth_in_bin, depth.size(), lut_size[0], d_lut, d_GridDivide, d_TailPhotonNum, states, time(nullptr));
	//printf("size of shared memory used per block: %i \n", block.x * sizeof(long));
	//printf("block number: %i\n", grid.x);
	//printf("Total shared memory used: %i kb\n", block.x * grid.x * sizeof(long) / 1024);
	//CHECK(cudaDeviceSynchronize());
	//printf("size of count array: %i kb\n", grid.x * sizeof(long) / 1024);
	CHECK(cudaMemcpy(h_DepthCnt_ary, d_DepthCnt_ary, grid.x * sizeof(long), cudaMemcpyDeviceToHost));
	
	//collect information
	vector<long> absorbcnt;
	for (int i = 0; i < depth.size(); i++) {
		int tmpcnt = 0;
		for (int j = h_GridDivide[i]; j < h_GridDivide[i + 1]; j++) {
			tmpcnt += h_DepthCnt_ary[j];
		}
		cout << "In depth " << depth[i];
		cout << ", incident photon num:" << incident_photon_num[i] << ", absorbed photon num:" << tmpcnt << endl;
		absorbcnt.push_back(tmpcnt);
	}

	CHECK(cudaDeviceSynchronize());
	iElaps = cpuSecond() - iStart;
	printf("Time elapsed %f ms\n", iElaps);


	CHECK(cudaFree(d_DepthCnt_ary));	free(h_DepthCnt_ary);
	CHECK(cudaFree(states));
	CHECK(cudaFree(d_GridDivide));		free(h_GridDivide);
	CHECK(cudaFree(d_depth_in_bin));
	CHECK(cudaFree(d_lut));
	CHECK(cudaFree(d_TailPhotonNum));	free(h_TailPhotonNum);
	return absorbcnt;
}
