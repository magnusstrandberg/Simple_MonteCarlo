#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
#include "timer.h"

using namespace std;
void fileReader(string);


int main() {
	timer test;
	test.startTime();
	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data(filelocation);
	if (input_data.checkInputCompletness() == 1) { cout << "Check input data file! There are errors"; }
	
	double time = test.calcStop();
	//cout << time;
	cin.ignore(); //So we can see shit
	return 0;
}
