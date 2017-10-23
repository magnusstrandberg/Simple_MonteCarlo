#include "Surface.h"
#include "Input.h"
#include <cmath>
#include <math.h>
#include <vector>

Surface::Surface()
{
}

Surface::~Surface()
{
}



void Surface::CreateSurface(std::vector <double> data, int identifier, int placeholder)
{
	ID = identifier;
	int tmp = 6;
	for (int i = 0; i < 10; i++)
	{
		surf_param[i] = data[i];
		if (surf_param[i] == 0 && i < 6)
		{
			tmp--;
		}
	}
	if (tmp == 0) 
	{
		type = 1;
	}
	else
	{
		type = 2;
	}

	return;
}

int Surface::insideSurf(double * position, bool complement)
{
				//Square parameters
	double S = surf_param[0]*pow(position[0],2) 
		+ (surf_param[1]*pow(position[1],2)) 
		+ (surf_param[2] * pow(position[2],2))
		//Cross parameters
		+ (surf_param[3] * position[0] * position[1]) 
		+ (surf_param[4] * position[1] * position[2])
		+ (surf_param[5] * position[2] * position[0])
		//Linear parameters
		+ (surf_param[6] * position[0])
		+ (surf_param[7] * position[1]) 
		+ (surf_param[8] * position[2])
		//Constant
		+ surf_param[9];
	
	int inside;
	
	if (S < 0)
	{
		inside = -1;

	}
	else if(S > 10e-14)
	{
		inside = 1;
	}
	else
	{
		return 0;
	}

	if (complement) {
		return -1 * inside;
	}
	else
	{
		return inside;
	}
}
/*
A=0	B=1
C=2	D=3
E=4 F=5
G=6 H=7
I=8 J=9
*/

int Surface::distToSurf(double * position, double * direction, double * distance) 
{
	if (type == 1)
	{
		double dominator = (surf_param[6] * position[0])
							+ (surf_param[7] * position[1])
							+ (surf_param[8] * position[2])
							+ surf_param[9];
		double nominator = +(surf_param[6] * direction[0])
							+ (surf_param[7] * direction[1])
							+ (surf_param[8] * direction[2]);
		if (nominator != 0)
		{
			distance[0] = -1*(dominator / nominator);
			return 1;
		}
		return 4;
	}
	else
	{
		double K, L, M;
		K = (surf_param[0] * pow(position[0], 2))
			+ (surf_param[1] * pow(position[1], 2))
			+ (surf_param[2] * pow(position[2], 2))
			+ (surf_param[3] * position[0] * position[1])
			+ (surf_param[4] * position[1] * position[2])
			+ (surf_param[5] * position[2] * position[0])
			+ (surf_param[6] * position[0])
			+ (surf_param[7] * position[1])
			+ (surf_param[8] * position[2])
			+ surf_param[9];

		L = (2 * ((surf_param[0] * position[0] * direction[0]))
			+ (surf_param[1] * position[1] * direction[1])
			+ (surf_param[2] * position[2] * direction[2]))
			+ (surf_param[3] * ((position[0] * direction[1]) + (position[1] * direction[0])))
			+ (surf_param[4] * ((position[1] * direction[2]) + (position[2] * direction[1])))
			+ (surf_param[5] * ((position[0] * direction[2]) + (position[2] * direction[0])))
			+ (surf_param[6] * position[0])
			+ (surf_param[7] * position[1])
			+ (surf_param[8] * position[2]);

		M = (surf_param[0] * pow(direction[0], 2))
			+ (surf_param[1] * pow(direction[1], 2))
			+ (surf_param[2] * pow(direction[2], 2))
			+ (surf_param[3] * direction[0] * direction[1])
			+ (surf_param[4] * direction[1] * direction[2])
			+ (surf_param[5] * direction[2] * direction[0]);
		//sqrt(L^2-4MK)
		if (M == 0)
		{
			return 3;
		}
		double square_value = (pow(L, 2) - (4 * M*K));

		if (square_value > 0)
		{
			double tmp[] = { ((-M + sqrt(square_value)) / (2 * M)), ((-M - sqrt(square_value)) / (2 * M)) };
			if (abs(tmp[0]) < abs(tmp[1])) {
				distance[0] = tmp[0];
				distance[1] = tmp[1];
			}
			else {
				distance[0] = tmp[1];
				distance[1] = tmp[0];
			}
			return 2;
		}
		else if (square_value < 10e-14 && square_value > 0)
		{
			distance[0] = (-M) / (2 * M);
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

void Surface::surfaceNorm(double * position,double * norm, bool side)
{
	double nabS_x;
	double nabS_y;
	double nabS_z;
	double abs_xyz;

	nabS_x = ((2 * surf_param[0] * position[0])
			+ (surf_param[3] * position[1])
			+ (surf_param[5] * position[2])
			+ (surf_param[6]));
	
	nabS_y = ((2 * surf_param[1] * position[1])
			+ (surf_param[3] * position[0])
			+ (surf_param[4] * position[2])
			+ (surf_param[7]));
	nabS_z = ((2 * surf_param[2] * position[2])
			+ (surf_param[4] * position[1])
			+ (surf_param[5] * position[0])
			+ (surf_param[8]));

	abs_xyz = sqrt(pow(nabS_x, 2) + pow(nabS_y, 2) + pow(nabS_z, 2));

	//normal on the other side.
	if (!side)
	{
		nabS_x = -1 * nabS_x;
		nabS_y = -1 * nabS_y;
		nabS_z = -1 * nabS_z;
	}

	norm[0] = nabS_x / abs_xyz;
	norm[1] = nabS_y / abs_xyz;
	norm[2] = nabS_z / abs_xyz;

	return;
}

int Surface::showID()
{
	return ID;
}


