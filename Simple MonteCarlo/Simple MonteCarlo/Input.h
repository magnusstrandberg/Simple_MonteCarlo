#pragma once
#include <iostream>
#include <string>
#include <vector>

enum Surf_type
{
	plane,
	cubid,
	sphere,
	cylinder,
	prism,
};

struct Cell_input
{
	int Cell_id;
	std::string cell_name;
	std::vector <int> members_id;
	std::vector <bool> complement;

};
struct Surf_input
{
	int complex_id;
	std::vector <std::vector <double> > Surfparams;
	std::vector <int> surf_ids;
	std::vector <bool> complement;
};
class Input
{
public:
	int N, M, Needle, Coin;
	double l, L;
	std::string input_location;
	std::vector <Surf_input> Complex_surf_input;
	std::vector <Cell_input> Cell_input_data;

	Input ();
	~Input ();
	void defaultValues();
	void fileReader(const std::string);
	void dataParser(const std::string);
	void createDataSet(const std::string, const std::string);
	void printData();
	int checkInputCompletness();
};


