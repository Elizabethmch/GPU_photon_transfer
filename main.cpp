#include <iostream>
#include <iomanip>
#include <vector>
#include <cuda_runtime.h>
#include "LUT.h"
#include "ManagePhotonAbsorption.h"
#include "MyCudaToolkit.h"
#include "RandomLUT.h"
#include <helper_cuda.h>
#include <helper_timer.h>

using namespace std;

/// <summary>
/// Print the everage probability of the given depth from LUT
/// </summary>
/// <param name="lut"></param>	LookUpTable address
/// <param name="depthid"></param>	Depth number
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
	cout << "	 " << setw(20) << "--sets-num=<N>" << "Specify number of sets in each simulation" << endl;
	cout << "    " << setw(20) << "--photon-num=<N>" << "Specify photon number in each sets" << endl;
	cout << "    " << setw(20) << "--verification" << "Enable verfication mode" << endl;
	cout << "    " << setw(20) << "--block-size=<N>" << "Specify number of threads per block" << endl;
	cout << "    " << setw(20) << "--seed=<N>" << "Specify the seed to use for the random number generator" << endl;
	cout << endl;

}

int main(int argc, char **argv) {

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

		if (getCmdLineArgumentString(argc, argv, "lut", &value))
		{
			table = new LUT(2, dim_size, value);
		}
		else
		{
			table = new LUT();
		}

		if (getCmdLineArgumentString(argc, argv, "sim-num", &value))
		{
			sim_num = (int)atoi(value);
		}
		else
		{
			sim_num = 1;
		}

		if (getCmdLineArgumentString(argc, argv, "sets-num", &value))
		{
			sets_num = (int)atoi(value);
		}
		else
		{
			sets_num = 1;
		}

		if (getCmdLineArgumentString(argc, argv, "photon-num", &value))
		{
			photon_num = (int)atoi(value);
		}
		else
		{
			photon_num = 1e7;
		}

		if (checkCmdLineFlag(argc, (const char**)argv, "verification"))
		{
			printf("Verification mode enabled\n");
			verification_mode = true;
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
	//incident_depth.push_back(0.5);	incident_depth.push_back(1.7);	incident_depth.push_back(20.7);
	//incident_photon_num.push_back(1e7);	incident_photon_num.push_back(1e7);	incident_photon_num.push_back(1e7);
	for (int i = 0; i < 100; i++) {
		/*incident_depth.push_back(0.5);*/	incident_depth.push_back(20.5);
		/*incident_photon_num.push_back(1e8);	*/incident_photon_num.push_back(1e7);
	}

	double MAXDPT = 100., MINDPT = 0.;

	ManagePhotonAbsorption* mpa = new ManagePhotonAbsorption(table, MAXDPT, MINDPT);


	PrintMeanProbOfDepth(table, 20);
	//PrintMeanProbOfDepth(table, 0);
	//PrintMeanProbOfDepth(table, 1);
	//PrintMeanProbOfDepth(table, 2);
	//PrintMeanProbOfDepth(table, 3);
	//PrintMeanProbOfDepth(table, 4);


	//simulation start
	double iStart, iElaps;
	iStart = cpuSecond();
	StopWatchInterface* timer = NULL;
	sdkCreateTimer(&timer);
	sdkStartTimer(&timer);

	auto result = mpa->getAbsorbedPhotonNum(incident_depth, incident_photon_num, 512);
	
	iElaps = cpuSecond() - iStart;
	sdkStopTimer(&timer);
	double elapsedTime = sdkGetAverageTimerValue(&timer) / 1000.0f;
	printf("Time elapsed (by my toolkit): %f s\n", iElaps);
	printf("Time (by helper_time.h) = % .2f(ms),\n", elapsedTime * 1000.0f);
	return 1;
}