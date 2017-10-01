#pragma once
#include <iostream>
#include <string>
#include "Input.h"
#include "stat.h"

class coin
{
	Input input_data;
	double  T;
	double * ratio;
	double  mean;
	double  stddiv;
	double fom;
public:
	coin(Input);
	~coin();

	void flipflip();
	double Timed();
	double Mean();
	double StdDiviation();
	double flipNtimes();
	double FOM();

};

