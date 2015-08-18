#include <omp.h>
#include "MACROS.h"

class OMP_MANAGER
{
public:
	int num_thread;
	double* arr_double;
	int* arr_int;

public:
	OMP_MANAGER():num_thread(0), arr_double(0), arr_int(0){}
	OMP_MANAGER(const int& num_thread_):num_thread(0), arr_double(0), arr_int(0)
	{
		Initialize(num_thread_);
	}
	~OMP_MANAGER()
	{
		SAFE_DELETE_ARRAY(arr_double);
		SAFE_DELETE_ARRAY(arr_int);
	}

public:
	void Initialize(const int& num_thread_)
	{
		num_thread = num_thread_;

		SAFE_DELETE_ARRAY(arr_double);
		SAFE_DELETE_ARRAY(arr_int);

		arr_double = new double[num_thread];
		arr_int = new int[num_thread];

		omp_set_num_threads(num_thread);
	}

	void Input(const double& val)
	{
		arr_double[omp_get_thread_num()] = val;
	}
	void Input(const int& val)
	{
		arr_int[omp_get_thread_num()] = val;
	}

	double Max_double()
	{
		double val = arr_double[0];
		for(int id = 1; id < num_thread; id ++) val = MAX(val ,arr_double[id]);
		return val;
	}
	int Max_int()
	{
		int val = arr_int[0];
		for(int id = 1; id < num_thread; id ++) val = MAX(val ,arr_int[id]);
		return val;
	}

	double Min_double()
	{
		double val = arr_double[0];
		for(int id = 1; id < num_thread; id ++) val = MIN(val ,arr_double[id]);
		return val;
	}
	int Min_int()
	{
		int val = arr_int[0];
		for(int id = 1; id < num_thread; id ++) val = MIN(val ,arr_int[id]);
		return val;
	}

	double Sum_double()
	{
		double val = arr_double[0];
		for(int id = 1; id < num_thread; id ++) val += arr_double[id];
		return val;
	}
	int Sum_int()
	{
		int val = arr_int[0];
		for(int id = 1; id < num_thread; id ++) val += arr_int[id];
		return val;
	}
};