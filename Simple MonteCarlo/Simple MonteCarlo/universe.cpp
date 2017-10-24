#include "universe.h"
#include "Surface.h"
#include "subspace.h"
#include "stat.h"
#include "timer.h"
#include <random>
#include <iostream>
#include <fstream>
#include <string>




Universe::Universe()
{
}

Universe::~Universe()
{
}

void Universe::buildSubspaces(Input data)
{
	N = data.N;
	M = data.M;


	for (int i = 0; i < data.subspace_input_data.size(); i++)
	{
		subspace a_subspace;
		a_subspace.makeSubspace(data.subspace_input_data[i],
			data.Complex_surf_input, data.Cell_input_data);
		subspaces.push_back(a_subspace);
	}

	return;
}

void Universe::calculateVolumes(int subspacerank)
{
	double totalV = (subspaces[subspacerank].subspaceranges.x[1]
					- subspaces[subspacerank].subspaceranges.x[0])
				* (subspaces[subspacerank].subspaceranges.y[1]
					- subspaces[subspacerank].subspaceranges.y[0])
				* (subspaces[subspacerank].subspaceranges.z[1]
					- subspaces[subspacerank].subspaceranges.z[0]);
	
	timer t;
	t.startTime();

	std::vector <std::vector <int> > hits
			(M, std::vector<int>(subspaces[subspacerank].cells.size()));

	

	for (int i = 0; i < M; i++)
	{
		hits[i] = pointVolume(subspacerank);
	}

	std::vector<double> means(subspaces[subspacerank].cells.size());
	std::vector<double> stdiv(subspaces[subspacerank].cells.size());
	std::vector<double> FOM(subspaces[subspacerank].cells.size());

	double T = t.calcStop();

	for (int i = 0; i < subspaces[subspacerank].cells.size(); i++)
	{
		std::vector <double> tmp_ratios(M);
		for (int j = 0; j < M; j++)
		{
			tmp_ratios[j] = ((double)hits[j][i]/(double)N);
		}
		means[i] = getMean(tmp_ratios);
		stdiv[i] = getStddiv(tmp_ratios, means[i]);
		FOM[i] = getFOM(stdiv[i], T);
		
	}

	std::string logfile = "Pointsvol.txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::app);
	log << "*************************\n";
	log << "Volumes in subspace: " << subspacerank + 1 << "\n";


		for (int i = 0; i < means.size(); i++)
	{	
		log << "Component " << subspaces[subspacerank].printCellName(i)
			<< " has a mean volume procentage of " << means[i]* totalV;
		log << " with a std_div of: " << stdiv[i]*totalV << " ";
		log << "FOM: " << FOM[i] << " \n";
	}
	
	log.close();



}

void Universe::plotSlice(double z0, int ID)
{
	double point[3];

	point[2] = z0;

	point[0] = subspaces[ID].subspaceranges.x[0];
	point[1] = subspaces[ID].subspaceranges.y[0];

	double save = point[1];

	double x_r = subspaces[ID].subspaceranges.x[1]
		- subspaces[ID].subspaceranges.x[0];
	double y_r = subspaces[ID].subspaceranges.y[1]
		- subspaces[ID].subspaceranges.y[0];

	int cell_id;
	
	std::string logfile;

	logfile = "plot" + std::to_string(ID) + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::app);
	log << "*************************\n";
	log << "Ranges: " << x_r << " " << y_r << "\n";
	log << "x;y;cellID; \n";

	int Nx, Ny;
	Nx = 100;
	Ny = 100;

	double x_add = x_r / Nx;
	double y_add = y_r / Ny;

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			cell_id = subspaces[ID].findCellatpoint(point);
			log	<< cell_id << ";";

			point[1] = point[1] + y_add;
		}
		log << "\n";
		point[0] = point[0] + x_add;
		point[1] = save;
	}

	log.close();

	return;
}


//Returns a vector of points that are inside that cell in the order the cells are listed.
std::vector<int> Universe::pointVolume(int subspaceRank)
{
	std::vector<int> cell_indexes(N);
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	double x_r = subspaces[subspaceRank].subspaceranges.x[1]
		- subspaces[subspaceRank].subspaceranges.x[0];
	double y_r = subspaces[subspaceRank].subspaceranges.y[1]
		- subspaces[subspaceRank].subspaceranges.y[0];
	double z_r = subspaces[subspaceRank].subspaceranges.z[1]
		- subspaces[subspaceRank].subspaceranges.z[0];

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < N; i++)
	{
		double point[3];
		point[0] = (dis(gen) * x_r) + subspaces[subspaceRank].subspaceranges.x[0];
		point[1] = (dis(gen) * y_r) + subspaces[subspaceRank].subspaceranges.y[0];
		point[2] = (dis(gen) * z_r) + subspaces[subspaceRank].subspaceranges.z[0];

		cell_indexes[i] = subspaces[subspaceRank].findCellatpoint(point);

	}
	
	std::vector<int> hit(subspaces[subspaceRank].cells.size());

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < subspaces[subspaceRank].cells.size(); i++)
	{
		hit[i] = std::count(cell_indexes.begin(), cell_indexes.end(), i);
	}
	
	return hit;
}

