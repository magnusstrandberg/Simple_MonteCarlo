#pragma once
#include "Surface.h"
#include "Input.h"
#include <vector>
#include <string>
#include <algorithm>
#include "Rotations.h"


class cell;
class complex_surf;

class subspace {
	friend class Universe;
public:
	Input input_data;
	int Id, boundry_id, boundry_index;
	double x_r, y_r, z_r, totalV;
	std::vector <cell> cells;
	std::vector <complex_surf> complex_surfs;
	ranges subspaceranges;

	subspace();
	~subspace();
	void makeSubspace(subspace_input, std::vector <Surf_input>,	std::vector <Cell_input>);
	int findCellatpoint(double *);
	int printCellID(int);
	std::string printCellName(int);
	int printCellMaterial(int);
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
	int insideComplexSurfaceNRM(double *, bool);
private:
	int complex_id;
	Surf_type type;
	std::vector <Surface> surfs;
	std::vector <bool> complements;
	bool moved;
	std::vector <double> place;
	bool rotated;
	std::vector <double> angles;
	Rotations Rota;
};

class cell
{
	friend class subspace;
public:
	cell();
	~cell();
	void createCellfromCompSurf(Cell_input, std::vector <complex_surf>);
	int insideCell(double *, bool);
	int insideCellNRM(double *, bool);
	double distanceCell(double * position, double * direction);
	Latice internalsubspace;
protected:
	std::string cell_name;
	int cell_id, material;
	bool rotated;
	bool moved;
private:
	std::vector <int> members_id;
	std::vector <bool> comp_cell_lvl;
	std::vector <complex_surf> components;
	std::vector <double> place;
	std::vector <double> angles;
	Rotations Rota;

};
