#include "Input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;


Input::Input()
{
}


Input::~Input()
{
}

void Input::defaultValues() {
	N = -999;
	M = -999;
	Needle = -999;
	Coin = -999;
	l = -999;
	L = -999;
}


void Input::fileReader(const std::string filelocation)
{
	Input::defaultValues();
	string line;
	ifstream input(filelocation);
	if (input.is_open())
	{
		while (getline(input, line))
		{

			if (!line.empty() && line.at(0) != '%')
			{
				//cout << line << "\n";
				dataParser(line);
			}

		}
		input.close();
	}

	else cout << "Unable to open file";

	return;
}

void Input::dataParser(const std::string line)
{
	bool parameter = true;
	string para, value;
	
	string::size_type split = line.find('=');

	para = line.substr(0, split);
	para.erase(std::remove(para.begin(), para.end(), ' '), para.end());

	value = line.substr(split + 1);
	value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
	value.erase(std::remove(value.begin(), value.end(), ';'), value.end());


	if (!para.empty() && !value.empty()) { createDataSet(para, value); }
	return;
}

void Input::createDataSet(const std::string para, const std::string value)
{
	if (para == "N") { 
		N = stoi(value);
	}
	
	if (para == "M") { 
		M = stoi(value);
	}
	if (para == "Needle") {
		Needle = stoi(value);
	}
	if (para == "Coin") {
		Coin = stoi(value);
	}
	if (para == "L") {
		L = stod(value);
	}
	if (para == "l") {
		l = stod(value);
	}
	return;

}

void Input::printData()
{
	cout << N << "\n";
	cout << M << "\n";
}

int Input::checkInputCompletness()
{
	if (N == -999) {return 1;}
	if (M == -999) {return 1;}
	if (Needle == -999) { return 1; }
	if (Coin == -999) { return 1; }
	if (L == -999) { return 1; }
	if (l == -999) { return 1; }
	return 0;
}

void Input::createspheredata(double * centre, double R)
{
	double parameters[10];
	parameters[0] = 1;
	parameters[1] = 1;
	parameters[2] = -1;
	parameters[3] = 0;
	parameters[4] = 0;
	parameters[5] = 0;
	parameters[6] = -2 * centre[0];
	parameters[7] = -2 * centre[1];
	parameters[8] = +2 * centre[2];
	parameters[9] = (-pow(R, 2) - pow(centre[2], 2)
		+ pow(centre[0], 2) + pow(centre[0], 2));
}
