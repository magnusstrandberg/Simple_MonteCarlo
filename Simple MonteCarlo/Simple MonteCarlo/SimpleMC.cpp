#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
#include "timer.h"
#include "coin.h"
#include "stat.h"
#include "Surface.h"
#include "subspace.h"
#include "universe.h"
#include "Rotations.h"
#include <vector>
#include "NData.h"

using namespace std;
void fileReader(string);
void testInputCompSurf(Input);
void testInputCell(Input);
void testInputSubspace(Input);



int main()
{
	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data;
	input_data.fileReader(filelocation);
	
	string log = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\listofmaterial.txt";
	NData Nucleardata;
	Nucleardata.readlist(log);

	//testInputCompSurf(input_data);
	//testInputCell(input_data); 
	//testInputSubspace(input_data);
	
	//Universe uni;
	//uni.buildSubspaces(input_data);
	//uni.calculateVolumes(0);
	//uni.calculateVolumes(1);
	//uni.calculateVolumes(2);
	//uni.plotSlice(0.0, 0);
	
	//testRotations();
	
	//uni.calculateVolumes(1);

	/*
	double direction[3];
	direction[0] = 1;
	direction[1] = 0;
	direction[2] = 0; 
	*/

	//double point[] = { 1, 0, 0.0 };
	//double test = uni.subspaces[0].complex_surfs[1].distanceComplexSurface(point, direction);
	//double tmp = uni.subspaces[0].y_r;
	//std:vector <double> vittu = uni.lineCalc(0);
	//int target[2];
	//uni.CellInUniverse(point, 1, target);
	
	//int test;
	//uni.CalculateLineVolume(0);
	//uni.plotSlice(0,1);
	return 0;
}


void testInputCompSurf(Input input_data) 
{
	ofstream log;
	log.open("read.txt", ios::out | ios::app);
	log << "*************************\n";

	log << "Input reader:\n";
	log << "How many Complex surfaces: " << input_data.Complex_surf_input.size() << "\n";

	for (int i = 0; i < input_data.Complex_surf_input.size(); i++)
	{
		log << "\n Surface nr: " << i + 1 << "\n";
		log << "Complex Id: " << input_data.Complex_surf_input[i].complex_id << "\n";
		log << "Subspace rank: " << input_data.Complex_surf_input[i].subspace_rank << "\n";
		log << "Type: " << input_data.Complex_surf_input[i].type << "\n";
		log << "How many surfaces components: " << input_data.Complex_surf_input[i].Surfparams.size() << "\n";
		for (int j = 0; j < input_data.Complex_surf_input[i].Surfparams.size(); j++)
		{
			log << "Surface nr" << j + 1 << " components \n";
			log << "Surface is complement: " << input_data.Complex_surf_input[i].complement[j] << "\n";
			for (int k = 0; k < 10; k++)
			{
				log << " " << input_data.Complex_surf_input[i].Surfparams[j].at(k) << "\n";
			}
		}
	}
	log.close();
	return;
}

void testInputCell(Input input_data)
{
	ofstream log;
	log.open("read.txt", ios::out | ios::app);
	log << "*************************\n";

	log << "Input reader:\n";
	log << "How many Cells: " << input_data.Cell_input_data.size() << "\n";

	for (int i = 0; i < input_data.Cell_input_data.size(); i++)
	{
		log << "\n Cell nr: " << i + 1 << "\n";
		log << "Cell Id: " << input_data.Cell_input_data[i].Cell_id << "\n";
		log << "Subspace rank: " << input_data.Cell_input_data[i].subspace_rank << "\n";
		log << "Name: " << input_data.Cell_input_data[i].cell_name << "\n";
		log << "How many surfaces components: " << input_data.Cell_input_data[i].cell_complements.size() << "\n";
		for (int j = 0; j < input_data.Cell_input_data[i].cell_complements.size(); j++)
		{
			log << "Surface ID: " << 
				input_data.Cell_input_data[i].cell_complements[j].comp_surface_id << " \n";
			log << "Surface is complement: " << 
					input_data.Cell_input_data[i].cell_complements[j].cell_complement << "\n";
		}
		log << "Contains a Latice: ";
		if (input_data.Cell_input_data[i].latice_info.type == 0)
		{
			log << "No\n";
		}
		else
		{
			log << "Yes, of type " << input_data.Cell_input_data[i].latice_info.type << "\n";
			log << "Inside is subspace: " << input_data.Cell_input_data[i].latice_info.subspace_inside << "\n";
		}
	}
	log.close();
	return;
}

void testInputSubspace(Input input_data)
{
	ofstream log;
	log.open("read.txt", ios::out | ios::app);
	log << "*************************\n";

	log << "Input reader:\n";
	log << "How many Subspaces: " << input_data.subspace_input_data.size() << "\n";

	for (int i = 0; i < input_data.subspace_input_data.size(); i++)
	{
		log << "\n Cell nr: " << i + 1 << "\n";
		log << "Subspace Id/rank: " << input_data.subspace_input_data[i].subspace_id << "\n";
		log << "Boundry cell ID: " << input_data.subspace_input_data[i].boundery_cell_id << "\n";
		log << "Boundry x range from: "
			<< input_data.subspace_input_data[i].subspacerange.x[0]
			<< " to "
			<< input_data.subspace_input_data[i].subspacerange.x[1]
			<< "\n";
		log << "Boundry y range from: "
			<< input_data.subspace_input_data[i].subspacerange.y[0]
			<< " to "
			<< input_data.subspace_input_data[i].subspacerange.y[1]
			<< "\n";
		log << "Boundry z range from: "
			<< input_data.subspace_input_data[i].subspacerange.z[0]
			<< " to "
			<< input_data.subspace_input_data[i].subspacerange.z[1]
			<< "\n";
	}
	log.close();
	return;
}


