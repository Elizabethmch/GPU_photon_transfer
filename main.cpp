#include <iostream>
#include "LUT.h"
#include "RandomLUT.h"


using namespace std;

int main(int argc, char **argv) {
	//PrintRndmLUT();
	LUT* table = new LUT();
	const vector<double>* test = table->getLUTAddress();
	//tabel->PrintLUT();
	return 1;
}