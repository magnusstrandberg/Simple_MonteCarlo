#pragma once
#include "Surface.h"
#include "Input.h"
#include <vector>


class subspace
{
	Input input_data;
	int Id;
	std::vector <Surface> surfs;



public:
	subspace();
	~subspace();
	void makeSubspace(Input);
	void makeSphere(double *, double);
};

