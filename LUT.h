#ifndef LUT_H
#define LUT_H

#include <vector>
#include <string>

using std::string;
using std::vector;
/// <summary>
/// class LUT: store look up table data
/// 
/// variables:
/// int dim: dimension number of the LUT (e.g. depth & angle means dim=2)
/// vector<int> size: size of each dimension (default 100 depth * 50 angle)
/// vector<double> prob_table: LUT data
/// 
/// member function:
/// getLUTAddress(): return the address for the vector of LUT data
/// void LoadLUTData(): Load the LUT from datafile
/// void PrintLUT(): Print detailed LUT
/// </summary>
class LUT
{
public:
	LUT();
	LUT(int dim,  vector<int> dim_size, string LutDataFile);
	~LUT();
	void LoadLUTData(string DataFileName = "lut.dat");
	void PrintLUT();
	int getLUTDimension() { return dim; }
	vector<int> getLUTSize() { return size; }
	string getLUTName() { return lut_name; }
	const vector<double>* getLUTAddress() { return &prob_table; }

private:
	int dim;			//dimension number of LUT
	vector<int> size;	//size of each dimension
	vector<double> prob_table;	//look up table for probability
	string lut_name;

};


#endif
