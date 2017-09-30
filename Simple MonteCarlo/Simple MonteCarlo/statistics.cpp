#include "statistics.h"
#include <cmath>

double Mean(const double ratios[], const int M)
{
	double mean = 0.0;
	for (int i = 0; i < M; i++) {
		mean += ratios[i];
	}
	return (mean/(double)M);
}

double stddiv(const double ratios[], const double mean, const int M)
{
	return 0.0;
}
