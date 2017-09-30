#include "coin.h"
#include "Input.h"
#include "statistics.h"
#include "timer.h"
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>


using namespace std;
coin::coin(const Input data):input_data(data)
{
}



coin::~coin()
{
}

void coin::flipflip()
{
	ratio = new double[input_data.M];
	timer t;
	t.startTime();
	double sum=0.0;
	for (int i = 0; i < input_data.M; i++) {
		ratio[i] = flipNtimes();
		//cout << ratio[i] << "\n";
		sum += ratio[i];
	}
	mean = sum/input_data.M;
	T = t.calcStop();
	cout << mean<<"\n";
}

double coin::Timed()
{
	return T;
}

double coin::flipNtimes(){
	int incircle = 0;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	#pragma omp parallel for
	for (int i = 0; i < input_data.N; i++)
	{
		double x = dis(gen);
		double y = dis(gen);
		if(sqrt(pow(x,2)+pow(y,2)) <= 1){
			incircle += 1;
			}
	}
	double ratio = 4*(double)incircle /(double)input_data.N;

	return ratio;
}
