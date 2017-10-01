#include "stat.h"
#include <cmath>


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
