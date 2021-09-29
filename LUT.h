#ifndef LUT_H
#define LUT_H

#include <vector>
#include <string>

using std::string;
using std::vector;

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
