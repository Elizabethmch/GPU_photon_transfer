#ifndef MANAGEPHOTONABSORPTION
#define MANAGEPHOTONABSORPTION
#include <iostream>
#include <vector>
#include "LUT.h"


using std::vector;

/// <summary>
/// Simulate photon absorption
/// </summary>
class ManagePhotonAbsorption
{
public:
	ManagePhotonAbsorption(LUT* lut, double maxdepth, double mindepth, int blocksize = 512) : look_up_table(lut),
		max_depth(maxdepth), min_depth(mindepth), threadBlockSize(blocksize) {
	}
	~ManagePhotonAbsorption() {
	}

	vector<long> getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num, unsigned int rndmseed = 0);


private:
	LUT* look_up_table;
	double max_depth;
	double min_depth;
	int threadBlockSize;

};
#endif

