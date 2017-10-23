#include "subspace.h"
#include "Surface.h"
#include "Input.h"
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>


subspace::subspace()
{
}


subspace::~subspace()
{
}



void subspace::makeSubspace(subspace_input subspace_data, 
							std::vector <Surf_input> Complex_surf_input,
							std::vector <Cell_input> Cell_input_data)
{
	Id = subspace_data.subspace_id;
	boundry_id = subspace_data.boundery_cell_id;
	subspaceranges = subspace_data.subspacerange;
	for (int i = 0; i < Complex_surf_input.size(); i++)
	{
		//Okey, as all different sub-spaces have their own cells and surfaces.
		if (Complex_surf_input[i].subspace_rank == Id)
		{
			complex_surfs.push_back(complex_surf());
			complex_surfs.back().createComplexSurface(Complex_surf_input[i]);
		}
	}


	for (int i = 0; i < Cell_input_data.size(); i++)
	{
		if (Cell_input_data[i].subspace_rank == Id)
		{
			cells.push_back(cell());
			cells.back().createCellfromCompSurf(Cell_input_data[i], complex_surfs);
		}
	}
}


int subspace::findCellatpoint(double * point)
{
	//returns cell vector index to the cell
	std::vector <int> inside_values;
	std::vector<int>::iterator it;
	for (size_t i = 0; i < cells.size(); i++)
	{
		inside_values.push_back(cells[i].insideCell(point,0));
	}

	int insides = std::count(inside_values.begin(), inside_values.end(), -1);

	if (insides == 1) 
	{
		for (int i = 0; i < inside_values.size(); i++)
		{
			if (inside_values[i]== -1)
			{
				return i;
			}
		}
	}

	return -1;
}

int subspace::printCellID(int order) 
{
	return cells[order].cell_id;
}

std::string subspace::printCellName(int order)
{
	return cells[order].cell_name;
}

complex_surf::complex_surf()
{
}

complex_surf::~complex_surf()
{
}

void complex_surf::createComplexSurface(Surf_input input)
{
	complex_id = input.complex_id;
	type = input.type;
	for (int i = 0; i < input.surf_ids.size(); i++)
		{
			surfs.push_back(Surface());
			surfs.back().CreateSurface(input.Surfparams[i], input.surf_ids[i], 0);
			complements.push_back(input.complement[i]);
		}

}

int complex_surf::insideComplexSurface(double * position, bool complemnt)
{
		bool inside = true;
		int value = 0;
		for (int i = 0; i < surfs.size(); i++)
		{
			if (surfs[i].insideSurf(position, complements[i]) != -1)
			{
				inside = false;
				break;
			}
		}
		
		if (inside) value = -1;
		else value = 1;
		
		if (complemnt)
		{
			return (-1 * value);
		}
		else
		{
			return value;
		}
}


double complex_surf::distanceComplexSurface(double * position, double * direction)
{
	double tmp[2];
	switch (type)
	{
	case plane:
		{
		if (surfs[0].distToSurf(position, direction, tmp) == 1) 
		{
			return tmp[0];
			
		}
		return -1.0;
		break;
		}
	case cylinder_inf:
		{
		if (surfs[0].distToSurf(position, direction, tmp) == 1)
		{
			return tmp[0];
		}
		return -1.0;
		break;
		}
	case sphere:
		{
		if (surfs[0].distToSurf(position, direction, tmp) == 1)
		{
			return tmp[0];
		}
		return -1.0;
		break;
		}
	case prism_square_inf:
	{
		std::vector <double> tries;
		for (int i = 0; i < surfs.size(); i++)
		{
			if(surfs[i].distToSurf(position, direction, tmp) ==1 )
				tries.push_back(tmp[0]);
		}
		std::sort(tries.begin(),tries.end());
		double epsilon = 10e-10;
		double newpos[3] = { 0,0,0 };
		for (int i = 0; i < tries.size(); i++)
		{
			if (tries[i] < 0)
			{
				continue;
			}
			newpos[0] = position[0] + tries[i] * direction[0]+epsilon;
			newpos[1] = position[1] + tries[i] * direction[1]+epsilon;
			newpos[2] = position[2] + tries[i] * direction[2]+epsilon;

			if (insideComplexSurface(position,0) != insideComplexSurface(newpos,0))
			{
				return tries[i];
			}
		}
		return -1;
		break;
	}

	case cubid:
	{
		std::vector <double> tries;
		for (int i = 0; i < surfs.size(); i++)
		{
			if (surfs[i].distToSurf(position, direction, tmp) == 1)
				tries.push_back(tmp[0]);
		}
		std::sort(tries.begin(), tries.end());
		double epsilon = 10e-10;
		double newpos[3] = { 0,0,0 };
		for (int i = 0; i < tries.size(); i++)
		{
			if (tries[i] < 0)
			{
				continue;
			}
			newpos[0] = position[0] + tries[i] * direction[0] + epsilon;
			newpos[1] = position[1] + tries[i] * direction[1] + epsilon;
			newpos[2] = position[2] + tries[i] * direction[2] + epsilon;

			if (insideComplexSurface(position, 0) != insideComplexSurface(newpos, 0))
			{
				return tries[i];
			}
		}
		return -1;
		break;
	}
	case cylinder:
	{
		std::vector <double> tries;
		for (int i = 0; i < surfs.size(); i++)
		{
			if (surfs[i].distToSurf(position, direction, tmp) == 1)
				tries.push_back(tmp[0]);
		}
		std::sort(tries.begin(), tries.end());
		double epsilon = 10e-10;
		double newpos[3] = { 0,0,0 };
		for (int i = 0; i < tries.size(); i++)
		{
			if (tries[i] < 0)
			{
				continue;
			}
			newpos[0] = position[0] + tries[i] * direction[0] + epsilon;
			newpos[1] = position[1] + tries[i] * direction[1] + epsilon;
			newpos[2] = position[2] + tries[i] * direction[2] + epsilon;

			if (insideComplexSurface(position, 0) != insideComplexSurface(newpos, 0))
			{
				return tries[i];
			}
		}
		return -1;
		break;
	}
	
	case prism_hex_inf:
	{
		std::vector <double> tries;
		for (int i = 0; i < surfs.size(); i++)
		{
			if (surfs[i].distToSurf(position, direction, tmp) == 1)
				tries.push_back(tmp[0]);
		}
		std::sort(tries.begin(), tries.end());
		double epsilon = 10e-10;
		double newpos[3] = { 0,0,0 };
		for (int i = 0; i < tries.size(); i++)
		{
			if (tries[i] < 0)
			{
				continue;
			}
			newpos[0] = position[0] + tries[i] * direction[0] + epsilon;
			newpos[1] = position[1] + tries[i] * direction[1] + epsilon;
			newpos[2] = position[2] + tries[i] * direction[2] + epsilon;

			if (insideComplexSurface(position, 0) != insideComplexSurface(newpos, 0))
			{
				return tries[i];
			}
		}
		return -1;
		break;
	}
	case prism_hex:
	{
		std::vector <double> tries;
		for (int i = 0; i < surfs.size(); i++)
		{
			if (surfs[i].distToSurf(position, direction, tmp) == 1)
				tries.push_back(tmp[0]);
		}
		std::sort(tries.begin(), tries.end());
		double epsilon = 10e-10;
		double newpos[3] = { 0,0,0 };
		for (int i = 0; i < tries.size(); i++)
		{
			if (tries[i] < 0)
			{
				continue;
			}
			newpos[0] = position[0] + tries[i] * direction[0] + epsilon;
			newpos[1] = position[1] + tries[i] * direction[1] + epsilon;
			newpos[2] = position[2] + tries[i] * direction[2] + epsilon;

			if (insideComplexSurface(position, 0) != insideComplexSurface(newpos, 0))
			{
				return tries[i];
			}
		}
		return -1;
		break;
	}

	default:
		return -1.0;
	}
}





cell::cell()
{
}

cell::~cell()
{
}

void cell::createCellfromCompSurf(Cell_input input, std::vector <complex_surf> comp)
{
	cell_id = input.Cell_id;
	cell_name = input.cell_name;

	for (int i = 0; i < comp.size(); i++)
	{
		//Need to be like this so compsurfuces can be part of more than one cell
		for (int j = 0; j < input.cell_complements.size() ; j++)
		{
			if (comp[i].complex_id == input.cell_complements[j].comp_surface_id )
			{
				components.push_back(comp[i]);
				members_id.push_back(input.cell_complements[j].comp_surface_id);
				comp_cell_lvl.push_back(input.cell_complements[j].cell_complement);
			}
		}
	}

	return;
}


int cell::insideCell(double * position, bool complement)
{
	bool inside = true;
	int value = 0;
	for (int i = 0; i < components.size(); i++)
	{
		if (components[i].insideComplexSurface(position, comp_cell_lvl[i]) != -1)
		{
			inside = false;
			break;
		}
	}

	if (inside) value = -1;
	else value = 1;

	if (complement)
	{
		return (-1 * value);
	}
	else
	{
		return value;
	}

}

double cell::distanceCell(double * position, double * direction)
{
	std::vector <double> tries;
	for (int i = 0; i < components.size(); i++)
	{
			tries.push_back(components[i].distanceComplexSurface(position, direction));
	}
	std::sort(tries.begin(), tries.end());
	double epsilon = 10e-10;
	double newpos[3] = { 0,0,0 };
	for (int i = 0; i < tries.size(); i++)
	{
		if (tries[i] < 0)
		{
			continue;
		}
		newpos[0] = position[0] + tries[i] * direction[0] + epsilon;
		newpos[1] = position[1] + tries[i] * direction[1] + epsilon;
		newpos[2] = position[2] + tries[i] * direction[2] + epsilon;

		if (insideCell(position, 0) != insideCell(newpos, 0))
		{
			return tries[i];
		}
	}
	return -1;
	return 0.0;
}