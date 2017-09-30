#pragma once
#include <iostream>
#include <string>
#include "Input.h"
#include "statistics.h"

class coin
{
	Input input_data;
	double  T;
	double * ratio;
	double  mean;
	double  stddiv;
public:
	coin(Input);
	~coin();

	void flipflip();
	double Timed();
	double flipNtimes();

};

