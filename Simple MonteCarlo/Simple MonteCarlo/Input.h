#pragma once
#include <iostream>
#include <string>

class Input
{
public:
	int N, M, Needle, Coin;
	double l, L;
	std::string input_location;
	Input (const std::string);
	~Input ();
	void defaultValues();
	void fileReader(const std::string);
	void dataParser(const std::string);
	void createDataSet(const std::string, const std::string);
	void printData();
	int checkInputCompletness();
};


