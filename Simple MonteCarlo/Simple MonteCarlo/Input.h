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

struct Latice
{
	int type;
	int subspace_inside;
	int hexa_pitch; //placeholder

};

struct ranges
{
	std::vector <double> x;
	std::vector <double> y;
	std::vector <double> z;
};

struct cell_comp
{
	int comp_surface_id;
	bool cell_complement;
};

struct subspace_input
{
	int subspace_id;
	int boundery_cell_id; //This cell is the outside.
	ranges subspacerange;
};

struct Cell_input
{
	int Cell_id;
	int subspace_rank;
	std::string cell_name;
	std::vector <cell_comp> cell_complements;
	Latice latice_info;
	bool moved = false;
	std::vector <double> place;
	bool rotated = false;
	std::vector <double> angles;

};

struct Surf_input
{
	int complex_id;
	Surf_type type;
	int subspace_rank;
	std::vector <std::vector <double> > Surfparams;
	std::vector <int> surf_ids;
	std::vector <bool> complement;
	bool moved = false;
	std::vector <double> place;
	bool rotated = false;
	std::vector <double> angles;
};
class Input
{
public:
	int N, M;
	std::string input_location;
	std::vector <Surf_input> Complex_surf_input;
	std::vector <Cell_input> Cell_input_data;
	std::vector <subspace_input> subspace_input_data;
	Input ();
	~Input ();
	void defaultValues();
	void fileReader(const std::string);
	void printData();
	int checkInputCompletness();
private:
	void dataParser(const std::string);
	void subspaceCreator(const std::string);
	void cellCreator(const std::string);
	void surfCreator(const std::string);
	void surfSphere(const std::string, Surf_input);
	void surfCylinder(const std::string, Surf_input);
	void surfParPlane(const std::string, Surf_input);
	void surfPlane(const std::string, Surf_input);
	void surfHexPrism(const std::string, Surf_input);
	std::vector <double> generalPlanParams(double[], double[], double[]);
	void surfCubid(const std::string, Surf_input);
	void surfSquarePrismInf(const std::string, Surf_input);
	void createDataSet(const std::string, const std::string);
	void transformCreator(const std::string data);

};


