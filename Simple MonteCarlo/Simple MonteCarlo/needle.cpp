#include "needle.h"
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

needle::needle(const Input data) :input_data(data)
{
}

needle::~needle()
{
}

//Main function for the Buffon Needle

void needle::throwthrow()
{
	ratio = new double[input_data.M];

	//Timer
	timer t;
	t.startTime();

	// For M sets of N Needles.
	for (int i = 0; i < input_data.M; i++) {
		ratio[i] = throwNtimes();
	}
	//Stop Timer
	T = t.calcStop();
	//Get statistical data
	mean = getMean(ratio, input_data.M);
	stddiv = getStddiv(ratio, mean, input_data.M);
	fom = getFOM(stddiv, T);
}



double needle::throwNtimes()
{
	/*Idea is that if needles are thrown on a stripped floor,
	the amount of needles that cross a stripline will be proportional to Pi.

	From input data the length of the needle "l" and width of stripes "L"

	Even though problem is 2D in real life in simulation only the x dimension matter.
	*/

	//Needed variables.
	int overLine = 0;
	bool * tmp;
	tmp = new bool[input_data.N];
	double pi_std = 4 * atan(1);
	double l = input_data.l;
	double L = input_data.L;
	//So to achive good random values.
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, L);

	//For loop for "needle tossing" and seeing end point.
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < input_data.N; i++)
	{
		double x = dis(gen)* L/2; //distance from needle centre to wall
		double phi = dis(gen) * pi_std; //random radian

		double needlelength_orth = (l / 2) * sin(phi);

		if (x < needlelength_orth)
		{
			tmp[i] = true;
		}
		else
		{
			tmp[i] = false;
		}
		
	}
	//Checking endpoints.
	for (int i = 0; i < input_data.N; i++)
	{
		if (tmp[i])
		{
			overLine++;
		}
	}
	delete[] tmp;

	//Calculating the Pi aproximation.
	double P = (double)overLine / (double)input_data.N;
	double pi_aprox=0;
	if (l <= L) {
		pi_aprox = (2* l)/(L*P);
	}
	else if( l == L)
	{
		pi_aprox = 2/P;
	}
	else
	{
		double C = (2 * acos(L / l));
		C += (l * 2 / L)*(1 - sqrt(1 - pow(L / l, 2)));
		pi_aprox = (1/P)*C;
	}
	return pi_aprox;
}


//Functions to return values

double needle::Timed()
{
	return T;
}

double needle::Mean()
{
	return mean;
}

double needle::StdDiviation()
{
	return stddiv;
}

double needle::FOM()
{
	return fom;
}
