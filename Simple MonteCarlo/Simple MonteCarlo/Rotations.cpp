#include "Rotations.h"
#include <cmath>



Rotations::Rotations()
{
}


Rotations::~Rotations()
{
}

void Rotations::Armultiplier(double * A, double * r, double * r_prime)
{
	for (int i = 0; i < 3; i++)
	{
		r_prime[i] = (A[i] * r[0]) + (A[i + 3] * r[1]) + (A[i + 6] * r[2]);
	}

	return;
}

void Rotations::Rotation(double * r, double * r_prime, double * ang)
{
	//Rotates a vector by the amount specified in ang for the different rotational axis.

	double Rx[] = {1, 0, 0, 
				0, cos(ang[0]), sin(ang[0]),
				0, -sin(ang[0]),cos(ang[0])};

	double Ry[] = { cos(ang[1]),0, -sin(ang[1]),
					0,			1,		0,
					sin(ang[1]),0, cos(ang[1]) };
	
	double Rz[] = { cos(ang[2]), sin(ang[2]), 0,
					-sin(ang[2]), cos(ang[2]), 0,
					0, 0, 1 };

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

	return;

}


