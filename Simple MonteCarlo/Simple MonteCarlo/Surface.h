#pragma once
#include <iostream>
#include <string>
#include "Input.h"
#include <cmath>
#include <math.h>
#include <vector>


class Surface {

	
	double surf_param[10];
	int ID,type;

public:
	Surface();
	~Surface();
	void CreateSurface(std::vector <double>, int, int);
	int insideSurf(double *, bool);
	int distToSurf(double *, double *, double *);
	int showID();
	
};