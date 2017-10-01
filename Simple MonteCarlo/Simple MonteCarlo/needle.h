#pragma once
#include <iostream>
#include <string>
#include "Input.h"
#include "stat.h"

class needle
{
	Input input_data;
	double  T;
	double * ratio;
	double  mean;
	double  stddiv;
	double fom;
public:
	needle(Input);
	~needle();

	void throwthrow();
	double Timed();
	double Mean();
	double StdDiviation();
	double throwNtimes();
	double FOM();

};
