#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"

using namespace std;
void fileReader(string);


int main() {
	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data(filelocation);

	cin.ignore(); //So we can see shit
	return 0;
}
