#include "Rotations.h"
#include <cmath>
#include <vector>



Rotations::Rotations()
{

}


Rotations::~Rotations()
{
}

void Rotations::Armultiplier(std::vector<double> A, double *r, double *r_prime)
{
	for (int i = 0; i < 3; i++)
	{
		r_prime[i] = (A[i] * r[0]) + (A[i + 3] * r[1]) + (A[i + 6] * r[2]);
	}

	return;
}

void Rotations::make(std::vector<double> ang) 
{
	if (ang.size()!= 3)
	{
		return;
	}
	double Rx_tmp[] = { 1, 0, 0,
		0, cos(ang[0]), sin(ang[0]),
		0, -sin(ang[0]),cos(ang[0]) };
	Rx.assign(Rx_tmp, Rx_tmp + 9);

	double Ry_tmp[] = { cos(ang[1]),0, -sin(ang[1]),
		0,			1,		0,
		sin(ang[1]),0, cos(ang[1]) };
	Ry.assign(Ry_tmp, Ry_tmp + 9);
	double Rz_tmp[] = { cos(ang[2]), sin(ang[2]), 0,
		-sin(ang[2]), cos(ang[2]), 0,
		0, 0, 1 };
	Rz.assign(Rz_tmp, Rz_tmp + 9);
	
	return;
}


void Rotations::Rotation(double * r)
{
	//Rotates a vector by the amount specified in ang for the different rotational axis.

	double r_prime[3];

	Armultiplier(Rx, r, r_prime);

	double tmp[3];
	for (int i = 0; i < 3; i++)
	{
		tmp[i] = r_prime[i];
	}

	Armultiplier(Ry, tmp, r_prime);

	for (int i = 0; i < 3; i++)
	{
		tmp[i] = r_prime[i];
	}
	Armultiplier(Rz, tmp, r_prime);

	r[0] = r_prime[0];
	r[1] = r_prime[1];
	r[2] = r_prime[2];
	return;

}


