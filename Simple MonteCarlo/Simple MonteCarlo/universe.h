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
	void plotSlice(double, int);
	std::vector<int> pointVolume(int);

	void CalculateLineVolume(int subspacerank);

	std::vector<double> lineCalc(int subspacerank);

	void Randomdirfrombound(int startside, double * dir);

	void CellInUniverse(double * point, int top_subspace, int * target);

	std::vector <subspace> subspaces;
private:
	int N, M;
};


