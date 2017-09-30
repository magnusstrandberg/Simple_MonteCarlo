#include "timer.h"
#include <time.h>

timer::timer()
{
}

timer::~timer()
{
}

void timer::startTime() {
	t_start = clock();
}

double timer::calcStop()
{
	t_stop = clock();
	time = ((double)t_stop - (double)t_start)/CLOCKS_PER_SEC;
	return time;
}
