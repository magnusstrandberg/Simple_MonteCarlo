#pragma once
#include <iostream>
#include <string>
#include "Input.h"
#include <cmath>
#include <math.h>


class Surface {

	
	double surf_param[10];
	int type;

public:
	Surface();
	~Surface();
	void CreateSurface(double *, int);
	int insideSurf(double *);
	int distToSurf(double *, double *, double *);
	
};