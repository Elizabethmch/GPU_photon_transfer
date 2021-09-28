#ifndef MANAGEPHOTONABSORPTION
#define MANAGEPHoTONABSORPTION
#include <iostream>
#include <vector>
#include "LUT.h"

using std::vector;

class ManagePhotonAbsorption
{
public:
	ManagePhotonAbsorption(LUT* lut, double maxdepth, double mindepth): look_up_table(lut),
					max_depth(maxdepth), min_depth(mindepth) {
	}
	~ManagePhotonAbsorption() {
	}

	long getAbsorbedPhotonNum(vector<double> depth, vector<long> incident_photon_num);


private:
	LUT* look_up_table;
	double max_depth;
	double min_depth;
};

#endif

