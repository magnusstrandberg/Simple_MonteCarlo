#include "subspace.h"
#include "Surface.h"
#include "Input.h"
#include <cmath>
#include <math.h>



subspace::subspace()
{
}


subspace::~subspace()
{
}

void subspace::makeSubspace(Input data)
{

}




void subspace::makeSphere(double * centre, double R)
{
	double parameters[10];
	parameters[0] = 1;
	parameters[1] = 1;
	parameters[2] = -1;
	parameters[3] = 0;
	parameters[4] = 0;
	parameters[5] = 0;
	parameters[6] = -2 * centre[0];
	parameters[7] = -2 * centre[1];
	parameters[8] = +2 * centre[2];
	parameters[9] = (-pow(R,2) -pow(centre[2],2)
					+ pow(centre[0], 2) + pow(centre[0], 2));
	
}
