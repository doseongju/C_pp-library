#pragma once

#include <iostream>
#include <time.h>

class STOPWATCH {
private: 
	clock_t startTime; // time the stop watch was started
	clock_t  stopTime; // stop the stop watch
	int total_ticks;
public:
	STOPWATCH() {total_ticks=0;}
	void stop (){ stopTime =clock(); } // stop the stopwatch
	void start(){ startTime=clock(); } // start the stopwatch
	void read_duration( const char* name_Of_work) 	{ 
		double duration = (double)(stopTime - startTime) / CLOCKS_PER_SEC;
		printf( "%2.3f seconds in %s\n", duration, name_Of_work );	}
	void read_ticks( const char* name_Of_work) 	{ 
		int number_of_ticks = (int)(stopTime - startTime);
		printf( "%d ticks in %s\n", number_of_ticks, name_Of_work );	}
	void add_into_total_ticks()	{total_ticks += (int)(stopTime - startTime);	}
	void read_total_ticks( const char* name_Of_work) { printf( "%d ticks in %s\n", total_ticks, name_Of_work );}
};