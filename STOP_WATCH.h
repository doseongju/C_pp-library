#pragma once

#include <iostream>
#include <time.h>
#include "ARRAY.h"
#include <assert.h>

class STOP_WATCH {
private: 
	clock_t start_time; // time the stop watch was started
	clock_t  stop_time; // stop the stop watch
	int num_toc;

public:
	STOP_WATCH()
	{
		num_toc = 0;
	}
	~STOP_WATCH(){}
	void Tic()
	{
		start_time = clock();
		num_toc = 0;
	}
	double Toc()
	{
		stop_time  = clock();
		num_toc ++;
		return (double)(stop_time - start_time) / CLOCKS_PER_SEC;
	}
	void TocAndReadDuration(const char* name_work)
	{
		double duration = Toc();
		printf( "%2.3f seconds in %s", duration, name_work);
	}
};