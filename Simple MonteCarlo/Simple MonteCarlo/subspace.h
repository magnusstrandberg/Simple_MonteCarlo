#pragma once
#include "Surface.h"
#include "Input.h"
#include <vector>
#include <string>
#include <algorithm>

class cell;
class complex_surf;

class subspace {
protected:
	Input input_data;
	int Id,boundry_id;
	std::vector <cell> cells;
	std::vector <complex_surf> complex_surfs;
	

public:
	subspace();
	~subspace();
	void makeSubspace(subspace_input, std::vector <Surf_input>,	std::vector <Cell_input>);
};



class complex_surf
{
	friend class cell;
public:
	complex_surf();
	~complex_surf();
	void createComplexSurface(Surf_input);
	int insideComplexSurface(double *, bool);
	double distanceComplexSurface(double *, double *);
private:
	int complex_id;
	Surf_type type;
	std::vector <Surface> surfs;
	std::vector <bool> complements;
};

class cell
{
public:
	cell();
	~cell();
	void createCellfromCompSurf(Cell_input, std::vector <complex_surf>);
	int insideCell(double *, bool);
	double distanceCell(double *, double *);
private:
	int cell_id;
	std::string cell_name;
	std::vector <int> members_id;
	std::vector <bool> comp_cell_lvl;
	std::vector <complex_surf> components;

};
