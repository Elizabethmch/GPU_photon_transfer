#include <iostream>
#include <iomanip>
#include <vector>
#include <cuda_runtime.h>
#include "LUT.h"
#include "ManagePhotonAbsorption.h"
#include "MyCudaToolkit.h"
#include "RandomLUT.h"
#include <helper_cuda.h>

using namespace std;

/// <summary>
/// Print the everage probability of the given depth from LUT
/// </summary>
/// <param name="lut"></param>	LookUpTable address
/// <param name="depthid"></param>	Depth number
void PrintMeanProbOfDepth(LUT* lut, int depthid) {
	const vector<double>* tableptr = lut->getLUTAddress();
	auto sizevec = lut->getLUTSize();
	int anglebinsize = sizevec[1];

	double mean_prob = 0;
	for (int i = 0; i < anglebinsize; i++) {
		int dataid = anglebinsize * depthid + i;
		mean_prob += (*tableptr)[dataid];
	}
	mean_prob /= anglebinsize;
	printf("Mean probability of depth %d is %f\n", depthid, mean_prob);
}

void showHelp(int argc, const char** argv)
{
	using std::cout;
	using std::endl;
	using std::left;
	using std::setw;

	if (argc > 0)
	{
		cout << endl << argv[0] << endl;
	}

	cout << endl << "Syntax:" << endl;
	cout << left;
	cout << "    " << setw(20) << "--lut=<lut.dat>" << "Specify look up table data file" << endl;
	cout << "    " << setw(20) << "--sim-num=<N>" << "Specify simulation number" << endl;
	cout << "    " << setw(20) << "--sets-num=<N>" << "Specify number of sets in each simulation" << endl;
	cout << "    " << setw(20) << "--photon-num=<N>" << "Specify photon number in each sets" << endl;
	cout << "    " << setw(20) << "--time-use" << "Printing time consumption for each sector" << endl;
	cout << "    " << setw(20) << "--verification" << "Enable verfication mode" << endl;
	cout << "    " << setw(20) << "--block-size=<N>" << "Specify number of threads per block" << endl;
	cout << "    " << setw(20) << "--seed=<N>" << "Specify the seed to use for the random number generator" << endl;
	cout << endl;

}

int main(int argc, char **argv) {
	double mainStart, mainElaps;
	mainStart = cpuSecond();


	using std::invalid_argument;
	using std::string;

	if (checkCmdLineFlag(argc, (const char**)argv, "help"))
	{
		printf("Displaying help on console\n");
		showHelp(argc, (const char**)argv);
		exit(EXIT_SUCCESS);
	}

	int sim_num = 1;
	int sets_num = 1;
	int photon_num = 1e7;
	int block_size = 512;
	unsigned int rndmseed = 0;
	bool time_use = false;
	bool verification_mode = false;
	vector<int> dim_size; dim_size.push_back(100); dim_size.push_back(50);
	LUT* table;
	//LUT* table = new LUT(2, dim_size, "lut100.dat");
	//LUT* table = new LUT(2, dim_size, "lut80.dat");
	//LUT* table = new LUT(2, dim_size, "lut30.dat");
	//LUT* table = new LUT();

	try
	{
		char* value;

		if (getCmdLineArgumentString(argc, (const char**) argv, "lut", &value))
		{
			table = new LUT(2, dim_size, value);
		}
		else
		{
			table = new LUT();
		}

		if (getCmdLineArgumentString(argc, (const char**) argv, "sim-num", &value))
		{
			sim_num = (int)atoi(value);
		}
		else
		{
			sim_num = 1;
		}

		if (getCmdLineArgumentString(argc, (const char**)argv, "sets-num", &value))
		{
			sets_num = (int)atoi(value);
		}
		else
		{
			sets_num = 1;
		}

		if (getCmdLineArgumentString(argc, (const char**)argv, "photon-num", &value))
		{
			photon_num = (int)atoi(value);
		}
		else
		{
			photon_num = 1e7;
		}

		if (checkCmdLineFlag(argc, (const char**)argv, "time_use"))
		{
			time_use = true;
		}

		if (checkCmdLineFlag(argc, (const char**)argv, "verification"))
		{
			printf("Verification mode enabled\n");
			verification_mode = true;
		}

		if (getCmdLineArgumentString(argc, (const char**)argv, "block-size", &value))
		{
			block_size = (int)atoi(value);
		}
		else
		{
			block_size = 512;
		}

		if (getCmdLineArgumentString(argc, (const char**)argv, "seed", &value))
		{
			rndmseed = (int)atoi(value);
		}
		else
		{
			rndmseed = 0;
		}
	}
	catch (invalid_argument& e)
	{
		printf("invalid command line argument (%s)\n", e.what());
		exit(EXIT_FAILURE);
	}

	//incident photon configuration
	vector<double> incident_depth;
	vector<long> incident_photon_num;
	double depth = 0.5;
	for (int i = 0; i < sets_num; i++) {
		incident_depth.push_back(depth);
		incident_photon_num.push_back(photon_num);
	}

	double MAXDPT = 100., MINDPT = 0.;

	ManagePhotonAbsorption* mpa = new ManagePhotonAbsorption(table, MAXDPT, MINDPT, block_size);



	//simulation start
	printf("depth vector size: %d, photon number in vec: %d\n", (int)incident_depth.size(), incident_photon_num[0]);
	vector<long> result;
	for (int i = 0; i < sim_num; i++) {
		result = mpa->getAbsorbedPhotonNum(incident_depth, incident_photon_num, rndmseed);
		double tmpresult = 0;
		for (int j = 0; j < incident_depth.size(); j++) {
			tmpresult += result[j];
		}
		cout << "Average absorbed photon number during simulation " << i << ": " << tmpresult / incident_depth.size() << endl;
	}

	if (time_use == true) {
		mpa->PrintTimeConsume();
	}

	//verification mode
	if (verification_mode == true) {
		
		printf("absorbed photon num:\n");
		for (int i = 0; i < incident_depth.size(); i++) {
			int tmpdepthid = (int)((incident_depth[i] - MINDPT) / (MAXDPT - MINDPT) * 100);
			PrintMeanProbOfDepth(table, tmpdepthid);
			printf("In set %d, depthid %d, absorbed photon num: %d\n", i, tmpdepthid, result[i]);
		}
	}
	
	mainElaps = cpuSecond() - mainStart;
	printf("Total time consumption for whole program: %f s\n", mainElaps);
	return 1;
}