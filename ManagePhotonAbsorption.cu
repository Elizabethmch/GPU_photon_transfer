#include "ManagePhotonAbsorption.h"
#include <cuda_runtime.h>
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

__global__ void SimulatePhotonAbsorption(long* count, int depth, double* lut, long* photon_num, int* grid_divide_points) {
	cg::thread_block cta = cg::this_thread_block();
	int bid = blockIdx.x;
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
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

	long* d_photonnum;
	CHECK(cudaMalloc((void**)d_photonnum, depth.size() * sizeof(long)));
	CHECK(cudaMemcpy((void*)d_photonnum, (void*)&(incident_photon_num[0]), depth.size() * sizeof(long), cudaMemcpyHostToDevice));

	long* h_DepthCnt_ary, * d_DepthCnt_ary;
	int *h_GridDivide, *d_GridDivide;
	//int* h_GridDivideId, * d_GridDibideId;
	h_GridDivide = (int*)malloc(depth.size() * sizeof(int));
	CHECK(cudaMalloc((void**)d_GridDivide, depth.size() * sizeof(int)));

	dim3 block, grid;
	block.x = threadBlockSize;
	grid.x = 0;
	//Simulate photon absorption for each depth
	for (int dptid = 0; dptid < depth.size(); dptid++) {
		//Allocate memory and size for gpu, grid and so on
		int tmp = (incident_photon_num[dptid] - 1) / block.x;
		grid.x += tmp;
		h_GridDivide[dptid] = grid.x;
	}
	CHECK( cudaMalloc((void**)&d_DepthCnt_ary, grid.x * sizeof(long)) );
	CHECK( cudaMemcpy( (void*)d_GridDivide, (void*)h_GridDivide, depth.size(), cudaMemcpyHostToDevice) );

	return 1;
}
