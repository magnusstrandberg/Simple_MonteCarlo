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

	for (size_t i = 0; i < subspaces.size(); i++)
	{
		for (size_t j = 0; j < subspaces[i].cells.size(); j++)
		{
			if (subspaces[i].printCellMaterial(j) == 0)
			{
				int ID = subspaces[i].cells[j].internalsubspace.subspace_inside;
				if (subspaces[i].cells[j].internalsubspace.type == 1)
				{
					double x = subspaces[ID].x_r;
					double y = subspaces[ID].y_r;
					double z = subspaces[ID].z_r;
					double x_tot = subspaces[i].x_r;
					double y_tot = subspaces[i].y_r;
					double z_tot = subspaces[i].z_r;
					subspaces[i].cells[j].internalsubspace.pitch.push_back(x);
					subspaces[i].cells[j].internalsubspace.pitch.push_back(y);
					subspaces[i].cells[j].internalsubspace.pitch.push_back(z);
					// amount of cells. 
					subspaces[i].cells[j].internalsubspace.nr_x = (int) floor((x + 1e-10) / x_tot);
					subspaces[i].cells[j].internalsubspace.nr_y = (int) floor((y + 1e-10) / y_tot);
					subspaces[i].cells[j].internalsubspace.nr_z = (int) floor((z + 1e-10) / z_tot);
				}
				if (subspaces[i].cells[j].internalsubspace.type == 2)
				{
					
				}
			}
		}
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
	int target[2];

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			CellInUniverse(point,ID,target);
			log	<< subspaces[target[0]].printCellID(target[1]) << ";";

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
		volpercentages[i] = lineCalc(subspacerank);
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
		log << " is " << means[i]*100 << " %, and " << means[i]*Volume << " cubics \n";
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
			point[1] = ( dis(gen) * subspaces[subspacerank].y_r ) + subspaces[subspacerank].subspaceranges.y[0];
			point[2] = ( dis(gen) * subspaces[subspacerank].z_r ) + subspaces[subspacerank].subspaceranges.z[0];
		}
		else if (startside == 1)
		{
			point[0] = (dis(gen) * subspaces[subspacerank].x_r) + subspaces[subspacerank].subspaceranges.x[0];
			point[1] = subspaces[subspacerank].subspaceranges.y[0] + epsilon;
			point[2] = (dis(gen) * subspaces[subspacerank].z_r ) + subspaces[subspacerank].subspaceranges.z[0];
		}
		else
		{
			point[0] = (dis(gen) * subspaces[subspacerank].x_r ) + subspaces[subspacerank].subspaceranges.x[0];
			point[1] = (dis(gen) * subspaces[subspacerank].y_r ) + subspaces[subspacerank].subspaceranges.y[0];
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


void Universe::CellInUniverse(double * point, int top_subspace, int * target)
//Target should be a int [2], giving subspace and ID of location.
{
	int target_subspace, target_index, tmp_material;

	double tmp_point[3];

	tmp_point[0] = point[0];
	tmp_point[1] = point[1];
	tmp_point[2] = point[2];

	target_subspace = top_subspace;
	
	do
	{
		target_index = subspaces[target_subspace].findCellatpoint(tmp_point);
		tmp_material = subspaces[target_subspace].printCellMaterial(target_index);

		if (tmp_material == 0)
		{
			//Calculate point in relation to bottom corner of cell containing subspace.
			tmp_point[0] = tmp_point[0] - subspaces[target_subspace].cells[target_index].internalsubspace.bottom[0];
			tmp_point[1] = tmp_point[1] - subspaces[target_subspace].cells[target_index].internalsubspace.bottom[1];
			tmp_point[2] = tmp_point[2] - subspaces[target_subspace].cells[target_index].internalsubspace.bottom[2];
			//
			tmp_point[0] = tmp_point[0] - 
				(floor(tmp_point[0] / subspaces[target_subspace].cells[target_index].internalsubspace.pitch[0])
				*subspaces[target_subspace].cells[target_index].internalsubspace.pitch[0]);
			tmp_point[1] = tmp_point[1] -
				(floor(tmp_point[1] / subspaces[target_subspace].cells[target_index].internalsubspace.pitch[1])
					*subspaces[target_subspace].cells[target_index].internalsubspace.pitch[1]);
			tmp_point[2] = tmp_point[2] -
				(floor(tmp_point[2] / subspaces[target_subspace].cells[target_index].internalsubspace.pitch[2])
					*subspaces[target_subspace].cells[target_index].internalsubspace.pitch[2]);
			//Put in relation to bottom of target subspace. OBS! rounding errors!
			int internal_subspace = subspaces[target_subspace].cells[target_index].internalsubspace.subspace_inside;
			tmp_point[0] = tmp_point[0] + subspaces[internal_subspace].subspaceranges.x[0] + 1e-8;
			tmp_point[1] = tmp_point[1] + subspaces[internal_subspace].subspaceranges.y[0] + 1e-8;
			tmp_point[2] = tmp_point[2] + subspaces[internal_subspace].subspaceranges.z[0] + 1e-8;
			//move point to next subspace at that index.
			target_subspace = internal_subspace;
		}
	} while (tmp_material == 0);

	target[0] = target_subspace;
	target[1] = target_index;
	return;
}