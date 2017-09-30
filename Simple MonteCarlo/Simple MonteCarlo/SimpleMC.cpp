#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
#include "timer.h"
#include "coin.h"

using namespace std;
void fileReader(string);


int main() {

	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data(filelocation);
	if (input_data.checkInputCompletness() == 1) { cout << "Check input data file! There are errors"; }
	coin coin_flip(input_data);
	coin_flip.flipflip();
	cout << coin_flip.Timed();
	cin.ignore(); //So we can see shit
	return 0;
}
