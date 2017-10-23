#pragma once
#include "Input.h"
#include "Surface.h"
#include "subspace.h"
#include "stat.h"
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <fstream>



class Universe
{
public:
	Universe();
	~Universe();
	void buildSubspaces(Input);
	void calculateVolumes(int);
	std::vector <int> pointVolume(int);

	std::vector <subspace> subspaces;
private:
	int N, M;
};


