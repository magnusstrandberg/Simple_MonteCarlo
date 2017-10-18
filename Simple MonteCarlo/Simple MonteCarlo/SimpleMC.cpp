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
	Surface test;
	double data[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, -2};
	test.CreateSurface(data,0,0,0);
	double point[3] = {0.0, 0, 0};
	double dir[] = {1,0,0};
	double dist[] = {0,0};
	log << test.insideSurf(point) << "\n";
	int tmp = test.distToSurf(point, dir, dist);
	if(tmp == 2)
	{
		log << "D1:" << dist[0] << " D2:" << dist[1] << "\n";
	}
	else if (tmp == 1)
	{
		log << "D1:" << dist[0] << " only one! \n";
	}
	else
	{
		log << "No hit with direction \n";
	}
	log.close();

	return 0;
}
