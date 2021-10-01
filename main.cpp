#include <iostream>
#include <vector>
#include "LUT.h"
#include "ManagePhotonAbsorption.h"
#include "RandomLUT.h"


using namespace std;

int main(int argc, char **argv) {
	//PrintRndmLUT();
	LUT* table = new LUT();
	//const vector<double>* test = table->getLUTAddress();
	//table->PrintLUT();
	vector<double> incident_depth;	incident_depth.push_back(0.5);
	vector<long> incident_photon_num;	incident_photon_num.push_back(1e7);
	double MAXDPT = 99., MINDPT = 0.;
	ManagePhotonAbsorption* mpa = new ManagePhotonAbsorption(table, MAXDPT, MINDPT);
	mpa->getAbsorbedPhotonNum(incident_depth, incident_photon_num);
	return 1;
}