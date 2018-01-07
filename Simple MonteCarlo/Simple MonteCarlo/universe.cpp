#include "universe.h"
#include "Surface.h"
#include "subspace.h"
#include "stat.h"
#include "timer.h"
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>



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
	double totalV = subspaces[subspacerank].totalV;
	
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
	log.open(logfile, std::ios::out | std::ios::trunc);
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

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < N; i++)
	{
		double point[3];
		point[0] = (dis(gen) * subspaces[subspaceRank].x_r) + subspaces[subspaceRank].subspaceranges.x[0];
		point[1] = (dis(gen) * subspaces[subspaceRank].y_r) + subspaces[subspaceRank].subspaceranges.y[0];
		point[2] = (dis(gen) * subspaces[subspaceRank].z_r) + subspaces[subspaceRank].subspaceranges.z[0];

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

void Universe::CalculateLineVolume(int subspacerank)
{

	double Volume = subspaces[subspacerank].totalV;

	timer t;
	t.startTime();

	std::vector <std::vector <double> > volpercentages
	(M, std::vector<double>(subspaces[subspacerank].cells.size()));

	for (int i = 0; i < M; i++)
	{
		volpercentages[i];
	}

	double T = t.calcStop();

	std::vector<double> means (subspaces[subspacerank].cells.size());
	std::vector<double> stdiv (subspaces[subspacerank].cells.size());
	std::vector<double> FOM (subspaces[subspacerank].cells.size());

	for (size_t i = 0; i < subspaces[subspacerank].cells.size(); i++)
	{
		std::vector <double> tmp_percents (M,0);
		for (size_t j = 0; j < M; j++)
		{
			tmp_percents[j] = volpercentages[j][i];
		}
		means[i] = getMean(tmp_percents);
		stdiv[i] = getStddiv(tmp_percents, means[i]);
		FOM[i] = getFOM(stdiv[i], T);
	}

	
	std::string logfile = "linevol.txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::app);
	log << "*************************\n";
	log << "The Volume calculated by lines in subspace " << subspacerank + 1 << "\n \n";

	for (int i = 0; i < means.size(); i++)
	{
		log << "Mean volume of " << subspaces[subspacerank].printCellName(i);
		log << " is " << means[i] << " % \n";
		log << "With a STD Div of " << stdiv[i] <<"and a FOM of: " << FOM[i] <<"\n";
	}
	log.close(); 
	return;
}

std::vector<double> Universe::lineCalc(int subspacerank)
{
	// Vector containg the total length, and individual lengths of each cell
	std::vector <std::vector <double> > Nlengths
	(N, std::vector <double>(subspaces[subspacerank].cells.size()+1));
	//Container for per centages
	std::vector <double> volpercentages(subspaces[subspacerank].cells.size(), 0);
	
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < N; i++)
	{

		int startside = floor(dis(gen) * 3);

		double direction[3];
		Randomdirfrombound(startside, direction);
		double epsilon = 10e-10;
		std::vector <std::pair<double, int>> lengths;
		std::pair <double, int> tmp_data;
		std::vector <std::pair<double, int>> tmp_list;

		double point[] = { 0,0,0 };

		if (startside == 0)
		{
			point[0] = subspaces[subspacerank].subspaceranges.x[0] + epsilon;
			point[1] = dis(gen) * subspaces[subspacerank].y_r - subspaces[subspacerank].subspaceranges.y[0];
			point[2] = dis(gen) * subspaces[subspacerank].z_r - subspaces[subspacerank].subspaceranges.z[0];
		}
		else if (startside == 1)
		{
			point[0] = dis(gen) * subspaces[subspacerank].x_r - subspaces[subspacerank].subspaceranges.x[0];
			point[1] = subspaces[subspacerank].subspaceranges.y[0] + epsilon;
			point[2] = dis(gen) * subspaces[subspacerank].z_r - subspaces[subspacerank].subspaceranges.z[0];
		}
		else
		{
			point[0] = dis(gen) * subspaces[subspacerank].x_r - subspaces[subspacerank].subspaceranges.x[0];
			point[1] = dis(gen) * subspaces[subspacerank].y_r - subspaces[subspacerank].subspaceranges.y[0];
			point[2] = subspaces[subspacerank].subspaceranges.z[0] + epsilon;
		}


		do
		{
			for (size_t i = 0; i < subspaces[subspacerank].cells.size(); i++)
			{
				double tmp = subspaces[subspacerank].cells[i].distanceCell(point, direction);
				if (0.0 < tmp)
				{
					tmp_data.first = tmp;
					tmp_data.second = -1;
					tmp_list.push_back(tmp_data);
				}
			}

			std::sort(tmp_list.begin(), tmp_list.end());
			tmp_list[0].second = subspaces[subspacerank].findCellatpoint(point);
			lengths.push_back(tmp_list[0]);

			point[0] = point[0] + tmp_list[0].first * direction[0] + epsilon;
			point[1] = point[1] + tmp_list[0].first * direction[1] + epsilon;
			point[2] = point[2] + tmp_list[0].first * direction[2] + epsilon;

			tmp_list.clear();

		} while (subspaces[subspacerank].cells[subspaces[subspacerank].boundry_index].insideCell(point, 0) == 1);

		std::vector <double> tmp_lengths(subspaces[subspacerank].cells.size() + 1, 0);

		for (size_t j = 0; j < lengths.size(); j++)
		{
			tmp_lengths[0] += lengths[j].first;
			tmp_lengths[lengths[j].second+1] += lengths[j].first;
		}

		Nlengths[i] = tmp_lengths;

	}

	 //sum of all line segments in N runs 
	std::vector <double> sumlengths (subspaces[subspacerank].cells.size() + 1,0);
	
	for (int i = 0; i < N; i++)
	{
		#pragma omp parallel for schedule(dynamic,1)
		for (int j = 0; j < sumlengths.size(); j++)
		{
			sumlengths[j] += Nlengths[i][j];
		}
	}

	for (size_t i = 0; i < volpercentages.size(); i++)
	{
		volpercentages[i] = sumlengths[i + 1] / sumlengths[0];
	}
	
	return volpercentages;
}

void Universe::Randomdirfrombound(int startside,double * dir)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	
	// startside 0= -x boundery surf, 1 = Y, 2 = Z

	if (startside == 0)
	{
		dir[0] = dis(gen);
		dir[1] = (dis(gen) * 2) -1;
		dir[2] = (dis(gen) * 2)- 1;
	}
	else if (startside == 1)
	{
		dir[0] = (dis(gen) * 2) - 1;
		dir[1] = dis(gen); 
		dir[2] = (dis(gen) * 2) - 1;
	}
	else
	{
		dir[0] = (dis(gen) * 2) - 1;
		dir[1] = (dis(gen) * 2) - 1;
		dir[2] = dis(gen);
	}

	double norm = sqrt(pow(dir[0], 2) + pow(dir[1], 2) + pow(dir[2], 2));

	dir[0] = dir[0] / norm;
	dir[1] = dir[1] / norm;
	dir[2] = dir[2] / norm;

	return;
}