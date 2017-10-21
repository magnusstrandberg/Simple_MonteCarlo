#pragma once
#include <iostream>
#include <string>
#include <vector>

enum Surf_type
{
	plane,
	cubid,
	prism_square_inf,
	sphere,
	cylinder,
	cylinder_inf,
	prism_hex,
	prism_hex_inf,
};

//Make a subspace struct that contains ID, and boundry box indications, preferebly a Cell ID comp.

struct cell_comp
{
	int comp_surface_id;
	bool cell_complement;
};

struct subspace_input
{
	int subspace_id;
	int boundery_cell_id; //This cell is the outside.
};

struct Cell_input
{
	int Cell_id;
	int subspace_rank;
	std::string cell_name;
	std::vector <cell_comp> cell_complements;
};

struct Surf_input
{
	int complex_id;
	Surf_type type;
	int subspace_rank;
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
	std::vector <subspace_input> subspace_input_data;
	Input ();
	~Input ();
	void defaultValues();
	void fileReader(const std::string);
	void dataParser(const std::string);
	void createDataSet(const std::string, const std::string);
	void printData();
	int checkInputCompletness();
	void createspheredata(double *, double);
};


