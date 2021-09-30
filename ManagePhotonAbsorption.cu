#include "ManagePhotonAbsorption.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cooperative_groups.h>
#include "MyCudaToolkit.h"

namespace cg = cooperative_groups;

__device__ unsigned int reduce_sum(long in, cg::thread_block cta)
{
	extern __shared__ long sdata[];


	// Perform first level of reduction:
	// - Write to shared memory
	unsigned int ltid = threadIdx.x;
	//printf("1");
	sdata[ltid] = in;

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

__global__ void SimulatePhotonAbsorption(long* absorb_num, int depth, double* lut) {
	cg::thread_block cta = cg::this_thread_block();
	int bid = blockIdx.x;
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid == 0) printf("flag %d processing...\n", flag);
	//if (tid > size) return;
	int count_flag = tid % 2;
	//printf("1: %d, %d\n", count_flag, threadIdx.x);

	if (tid >= size) count_flag = 0;
	count_flag = reduce_sum(count_flag, cta);
	//__syncthreads();

	if (threadIdx.x == 0) {
		//printf("2\n");
		//printf("%d\n", count_flag);
		printf("Check sync in GPU\n");
		count[bid] = count_flag;
	}
	if (tid == 0) printf("flag %d complete...\n", flag);
}


/**************************************************
* Calculate number of photon absorbed for given depth & incident photon number
**************************************************/
long ManagePhotonAbsorption::getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, int threadBlockSize) {
	
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
