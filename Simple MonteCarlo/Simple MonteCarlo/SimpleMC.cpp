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
	log << "Testing simple surface \n";
	Surface test;
	double data[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 2};
	double data2[] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, -2 };
	double data3[] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 2 };
	double data4[] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, -2 };
	vector <double> datav1(data, data + sizeof(data) / sizeof(double));
	vector <double> datav2(data2, data2 + sizeof(data2) / sizeof(double));
	vector <double> datav3(data3, data3 + sizeof(data3) / sizeof(double));
	vector <double> datav4(data4, data4 + sizeof(data4) / sizeof(double));

	test.CreateSurface(datav1,0,0);
	double point[3] = {-1, 0, 0};
	double dir[] = {1,0,0};
	double dist[] = {0,0};
	log << "If negative inside: " << test.insideSurf(point, 0) << "\n";
	int tmp = test.distToSurf(point, dir, dist);
	if(tmp == 2)
	{
		log << "D1: " << dist[0] << " D2: " << dist[1] << "\n";
	}
	else if (tmp == 1)
	{
		log << "D1: " << dist[0] << " only one! \n";
	}
	else
	{
		log << "No hit with direction \n";
	}

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

	complex_surf comptest;
	comptest.createComplexSurface(test2);

	log << "Testing the complex surface of type: " << test2.type << "\n";
	log << "If negative inside: " << comptest.insideComplexSurface(point, 0) << "\n";

	double dist2 = comptest.distanceComplexSurface(point, dir);
	log << "Distance from complex: " << dist2 << "\n";


	log.close();

	return 0;
}
