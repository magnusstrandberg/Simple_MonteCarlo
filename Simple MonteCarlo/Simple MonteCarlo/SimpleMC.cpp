#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
#include "timer.h"
#include "coin.h"
#include "needle.h"
#include "stat.h"

using namespace std;
void fileReader(string);


int main() {

	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data(filelocation);
	ofstream log;
	log.open("read.txt", ios::out | ios::app);
	log << "*************************\n";
	if (input_data.checkInputCompletness() == 1) { log << "Check input data file! There are errors! \n"; }
	coin coin_flip(input_data);
	coin_flip.flipflip();
	log << "N: " << input_data.N << "\n";
	log << "M: " << input_data.M << "\n";
	log << "Calculating Pi with the coin toss method\n";
	log << "Mean Pi value: " << coin_flip.Mean() << "\n";
	log << "Standard diviation: "<< coin_flip.StdDiviation() << "\n";
	log << "Figure of merit: " << coin_flip.FOM() << "\n";
	log << "Time to calculate: "<< coin_flip.Timed() << "s \n";
	log << "A value:" << coin_flip.A() << "\n";

	needle needle_toss(input_data);
	needle_toss.throwthrow();
	log << "\n";
	log << "Calculating Pi with the Buffon's needle method\n";
	log << "l: " << input_data.l << "\n";
	log << "L: " << input_data.L << "\n";
	log << "Mean Pi value: " << needle_toss.Mean() << "\n";
	log << "Standard diviation: " << needle_toss.StdDiviation() << "\n";
	log << "Figure of merit: " << needle_toss.FOM() << "\n";
	log << "Time to calculate: " << needle_toss.Timed() << "s \n";
	log << "A value:" << needle_toss.A() << "\n";
	log.close();
	ofstream mat;
	mat.open("log.txt", ios::out | ios::app);
	mat << input_data.N << ";";
	mat << coin_flip.Mean() << ";";
	mat << coin_flip.StdDiviation() << ";";
	mat << coin_flip.FOM() << ";";
	mat << coin_flip.Timed() << ";";
	mat << needle_toss.Mean() << ";";
	mat << needle_toss.StdDiviation() << ";";
	mat << needle_toss.FOM() << ";";
	mat << needle_toss.Timed() << ";\n";
	mat.close();

	return 0;
}
