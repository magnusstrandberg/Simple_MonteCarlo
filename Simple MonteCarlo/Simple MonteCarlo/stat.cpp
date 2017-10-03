#include "stat.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

double getMean(const double ratios[], const int M)
{
	double mean = 0.0;
	for (int i = 0; i < M; i++) {
		mean += ratios[i];
	}
	return (mean/(double)M);
}


double getStddiv(const double ratios[], const double mean, const int M)
{
	double stddiv=0;
	for (int i = 0; i < M; i++)
	{
		stddiv += pow((ratios[i]-mean),2);
	}

	stddiv = stddiv / M;
	stddiv = sqrt(stddiv);
	return stddiv;
}

double getFOM(const double stddiv, const double T)
{
	return (1 / (pow(stddiv, 2)*T));
}

/*using the fancy relations between the complementary error function,
and CDF of a nice and ordered data set. */

double giveCFD(double value) 
{
	
	return 0.5*erfc(-value/sqrt(2));
}

double AndersonDarlingTest(const double sample[], const int M, const double mean, const double stddiv) 
{
	//Need to make a copy of the sample as it will be heavily modified
	double * X;
	X = new double[M];
	memcpy(X, sample, M * sizeof(double));
	
	std::sort(X, X+M);
	
	//Standardize X
	for (int i = 0; i < M; i++)
	{
		X[i] = (X[i] - mean) / stddiv;
	}


	double S =  0; //Readability
	for (int i = 0; i < M; i++)
	{
		/*Need to be careful with loss of resolution, moving N outside of for loop*/
		S += ((2*(i+1) -1)*(log(giveCFD(X[i]))+log(1-giveCFD(X[M-1-i]))));
	}
	S = S/M;
	double A = (-M - S);




	delete[] X;
	return A;
}