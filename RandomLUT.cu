#include <iostream>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <fstream>

#include "MyCudaToolkit.h"
#include "RandomLUT.h"

using namespace std;


__global__ void getLUT(double *lut, curandStateXORWOW_t* states, unsigned long long seed, size_t size) {
	unsigned int bid = blockIdx.x;
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid > size) return;
	curand_init(seed, tid, 0, &states[tid]);
	lut[tid] = curand_uniform_double(&states[tid]);
}

void PrintRndmLUT() {
	//Get LookUpTable Size
	const int dim = 2;
	int dim_size[dim] = { 100, 50 };
	int total_size = 1;
	cout << "Printing random LUT with dimension: ";
	for (int i = 0; i < dim; i++) {
		cout << dim_size[i];
		if (i < dim - 1)	cout << " * " ;
		total_size *= dim_size[i];
	}
	cout << endl;
	

	//Define Block Size
	dim3 block;
	dim3 grid;
	block.x = 1024;
	grid.x = (total_size - 1) / block.x + 1;


	//Allocate storage for LUT in host and device
	double iStart, iElaps;
	iStart = cpuSecond();
	double* h_lut = nullptr;
	double* d_lut = nullptr;
	CHECK( cudaMalloc((void**)&d_lut, total_size * sizeof(double)) );
	h_lut = (double*)malloc(total_size * sizeof(double));


	//Generate random numbers for lut
	curandStateXORWOW_t* states;
	CHECK( cudaMalloc( (void**)&states, sizeof(curandStateXORWOW_t) * total_size) ); 

	getLUT << <block, grid >> > (d_lut, states, time(nullptr), total_size);

	CHECK( cudaMemcpy(h_lut, d_lut, total_size * sizeof(double), cudaMemcpyDeviceToHost) );

	CHECK(cudaDeviceSynchronize());
	iElaps = cpuSecond() - iStart;
	printf("Time elapsed %f ms\n", iElaps);


	//Output The LUT into binary format
	//for (int i = 0; i < total_size; i++) {
	//	cout << h_lut[i] <<"\t";
	//}
	//cout << endl;
	ofstream outFile("lut.dat", ios::out | ios::binary);
	for (int i = 0; i < total_size; i++) {
		outFile.write((char*)&h_lut[i], sizeof(h_lut[i]));
	}
	outFile.close();
	cudaFree(d_lut);
	cudaFree(states);
	free(h_lut);

}

void PrintUniformLUT(int prob) {
	const int dim = 2;
	int dim_size[dim] = { 100, 50 };
	int total_size = 1;
	cout << "Printing uniform LUT with dimension: ";
	for (int i = 0; i < dim; i++) {
		cout << dim_size[i];
		if (i < dim - 1)	cout << " * ";
		total_size *= dim_size[i];
	}
	cout << endl;

	double* h_lut = (double*)malloc(sizeof(double) * total_size);


	ofstream outFile("lut"+to_string(prob)+".dat", ios::out | ios::binary);
	for (int i = 0; i < total_size; i++) {
		h_lut[i] = prob / 100.;
		outFile.write((char*)&h_lut[i], sizeof(h_lut[i]));
	}
	outFile.close();
	free(h_lut);

}