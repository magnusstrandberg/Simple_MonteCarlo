#include <time.h>


class timer
{
public:
	timer();
	~timer();
	clock_t t_start,t_stop;
	void startTime();
	double time;
	double calcStop();

private:

};


