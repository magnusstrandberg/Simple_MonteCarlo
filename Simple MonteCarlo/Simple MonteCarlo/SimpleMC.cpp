#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
#include "timer.h"
#include "coin.h"
#include "needle.h"
#include "stat.h"
#include "Surface.h"
#include "subspace.h"
#include "universe.h"

using namespace std;
void fileReader(string);


int main()
{
	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	Input input_data;
	input_data.fileReader(filelocation);
	ofstream log;
	log.open("read.txt", ios::out | ios::app);
	log << "*************************\n";
	log << "Testing Cells";
	
	double data[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 2};
	double data2[] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, -2 };
	double data3[] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 2 };
	double data4[] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, -2 };
	double data5[] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };
	vector <double> datav1(data, data + sizeof(data) / sizeof(double));
	vector <double> datav2(data2, data2 + sizeof(data2) / sizeof(double));
	vector <double> datav3(data3, data3 + sizeof(data3) / sizeof(double));
	vector <double> datav4(data4, data4 + sizeof(data4) / sizeof(double));
	vector <double> datav5(data5, data5 + sizeof(data5) / sizeof(double));
	double point[3] = {0, 0, -1};
	double dir[] = {1,0,1};
	double dist[] = {0,0};


	Surf_input test2;
	test2.complex_id = 100;
	test2.surf_ids.push_back(1);
	test2.surf_ids.push_back(2);
	test2.surf_ids.push_back(3);
	test2.surf_ids.push_back(4);
	test2.type = prism_square_inf;
	test2.complement.push_back(1);
	test2.complement.push_back(0);
	test2.complement.push_back(1);
	test2.complement.push_back(0);
	test2.Surfparams.push_back(datav1);
	test2.Surfparams.push_back(datav2);
	test2.Surfparams.push_back(datav3);
	test2.Surfparams.push_back(datav4);
	
	Surf_input test3;
	test3.complex_id = 200;
	test3.surf_ids.push_back(1);
	test3.type = plane;
	test3.complement.push_back(0);
	test3.Surfparams.push_back(datav5);

	complex_surf comptest1;
	comptest1.createComplexSurface(test2);
	complex_surf comptest2;
	comptest2.createComplexSurface(test3);

	vector <complex_surf> celltestsurf;
	celltestsurf.push_back(comptest1);
	celltestsurf.push_back(comptest2);

	cell_comp comp1;
	comp1.cell_complement = 0;
	comp1.comp_surface_id = 100;

	cell_comp comp2;
	comp2.cell_complement = 1;
	comp2.comp_surface_id = 200;

	

	Cell_input celltest;
	celltest.Cell_id = 10000;
	celltest.subspace_rank = 1;
	celltest.cell_name = "test";
	celltest.cell_complements.push_back(comp1);
	celltest.cell_complements.push_back(comp2);

	cell cell1;
	cell1.createCellfromCompSurf(celltest, celltestsurf);


	log << "Testing the complex surface of type: " << test2.type << "\n";
	log << "If negative inside: " << comptest1.insideComplexSurface(point, 0) << "\n";

	double dist2 = comptest1.distanceComplexSurface(point, dir);
	log << "Distance from complex: " << dist2 << "\n";

	log << "Testing the cell \n";
	log << "Negative if inside cell: " << cell1.insideCell(point, 0) << "\n";
	dist2 = cell1.distanceCell(point, dir);
	log << "Distance from cell: " << dist2 << "\n";
	

	log.close();

	return 0;
}
