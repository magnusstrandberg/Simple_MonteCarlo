#include "Input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;


Input::Input(const std::string input)
{
	N = 1;
	M = 1;
	input_location = input;
	fileReader(input_location);
}


Input::~Input()
{
}

void Input::fileReader(const std::string filelocation)
{
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
	if (para == "N") { N = stoi(value);
	}
	
	if (para == "M") { M = stoi(value);
	}
	
	return;

}

void Input::printData()
{
	cout << N << "\n";
	cout << M << "\n";
}
