#pragma once
#include <vector>
class Rotations
{
public:
	Rotations();
	~Rotations();
	void make(std::vector<double>);
	void Armultiplier(std::vector <double>, double*, double *);
	void Rotation(double*);
	std::vector <double> Rx;
	std::vector <double> Ry;
	std::vector <double> Rz;
};

