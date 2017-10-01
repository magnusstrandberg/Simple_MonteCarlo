#include "coin.h"
#include "Input.h"
#include "stat.h"
#include "timer.h"
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>


using namespace std;
coin::coin(const Input data):input_data(data)
{
}



coin::~coin()
{
}

void coin::flipflip()
{
	ratio = new double[input_data.M];
	timer t;
	t.startTime();
	for (int i = 0; i < input_data.M; i++) {
		ratio[i] = flipNtimes();
	}
	T = t.calcStop();
	mean = getMean(ratio, input_data.M);
	stddiv = getStddiv(ratio, mean, input_data.M);
	fom = getFOM(stddiv, T);
	
}

double coin::Timed()
{
	return T;
}

double coin::Mean()
{
	return mean;
}

double coin::StdDiviation()
{
	return stddiv;
}

double coin::flipNtimes(){
	int incircle = 0;
	double * tmp;
	tmp = new double[input_data.N];
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < input_data.N; i++)
	{
		double x = dis(gen);
		double y = dis(gen);
		tmp[i] = sqrt(pow(x,2) + pow(y,2));
	}
	for (int i = 0; i < input_data.N; i++)
	{
		if (tmp[i] < 1)
		{
			incircle += 1;
		}
	}
	delete[] tmp;
	double ratio = (4 * (double)incircle )/(double)input_data.N;

	return ratio;
}

double coin::FOM()
{
	return fom;
}
