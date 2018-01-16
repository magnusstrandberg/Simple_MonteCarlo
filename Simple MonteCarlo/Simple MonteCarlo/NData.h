#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <map>

struct MT
{
	int MT_number;
	double Q_value;
	std::map<double, double> Ec_pair;
};

struct materialdata
{
	int Z, A;
	double T, Aw;
	std::string Symbol;
	int Fissionable;
	std::map<double, double> Fission_nubar;
	std::vector <MT> MT_data;

};



class NData
{
public:
	NData();
	~NData();
	void readdatafile(std::string filename, int it);
	void readlist(std::string filename);
	std::vector <materialdata> materials;
};

