#pragma once
#include <iostream>
#include <string>

class Input
{
public:
	int N, M;
	std::string input_location;
	Input (const std::string);
	~Input ();
	void fileReader(const std::string);
	void dataParser(const std::string);
	void createDataSet(const std::string, const std::string);
	void printData();
};


