#pragma once

#include <iostream>
#include <cassert>
#include <fstream>
#include "MACROS.h"

template<class T>
class ARRAY_2D
{
public:
	int i_start, i_end,
		j_start, j_end,
		i_res, j_res, ij_res;
	T* val;

public:
	ARRAY_2D():val(0), i_start(0), i_end(0), j_start(0), j_end(0), i_res(0), j_res(0), ij_res(0){}
	ARRAY_2D(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_, const bool& initialize = false):val(0)
	{
		Initialize(i_start_, j_start_, i_res_, j_res_, initialize);
	}
	ARRAY_2D(const ARRAY_2D<T>& arr):val(0)
	{
		Initialize(arr.i_start, arr.j_start, arr.i_res, arr.j_res, true);
		for(int i = 0; i < arr.ij_res; i++)
			val[i] = arr.val[i];
	}

	~ARRAY_2D(){if(val != 0) delete[] val; val = 0;}

public:
	void Initialize(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_, const bool& initialize = false)
	{
		i_start = i_start_;
		i_res = i_res_;
		i_end = i_start + i_res - 1;
		j_start = j_start_;
		j_res = j_res_;
		j_end = j_start + j_res - 1;
		ij_res = i_res * j_res;

		if(val != 0) delete [] val;
		val = new T[ij_res];

		if(initialize == true)
		{
			for(int i = 0; i < ij_res; i++) val[i] = T();
		}
	}

	void Initialize(const ARRAY_2D<T>& arr, const bool& initialize = false)
	{
		Initialize(arr.i_start, arr.j_start, arr.i_res, arr.j_res, initialize);
	}

	void AssignAllValues(const T& value)
	{
		for(int i = 0; i < ij_res; i++) val[i] = value;
	}

	T FindMax() const
	{
		assert(val != 0);
		T result = val[0];
		for(int i = 1; i < ij_res; i++)
		{
			if(val[i] > result) result = val[i];
		}
		return result;
	}

	T FindMin() const
	{
		assert(val != 0);
		T result = val[0];
		for(int i = 1; i < ij_res; i++)
		{
			if(val[i] < result) result = val[i];
		}
		return result;
	}

	T FindMaxMagnitude() const
	{
		assert(val != 0);
		T result = ABS(val[0]);
		for(int i = 1; i < ij_res; i++)
		{
			if(ABS(val[i]) > result) result = ABS(val[i]);
		}
		return result;
	}

	T FindMinMagnitude() const
	{
		assert(val != 0);
		T result = ABS(val[0]);
		for(int i = 1; i < ij_res; i++)
		{
			if(ABS(val[i]) < result) result = ABS(val[i]);
		}
		return result;
	}

public:
	T& operator()(const int& i, const int& j)
	{
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);

		return val[i - i_start + i_res * (j - j_start)];
	}

	const T& operator()(const int& i, const int& j) const 
	{
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);

		return val[i - i_start + i_res * (j - j_start)];
	}

	void operator=(const ARRAY_2D<T>& arr_)
	{
		Initialize(arr_.i_start, arr_.j_start, arr_.i_res, arr_.j_res);
		for(int i = 0; i < ij_res; i++) val[i] = arr_.val[i];
	}


public:
	void Print(const char* filename)
	{
		std::ofstream fout;
		fout.open(filename);
		for(int j = j_start; j <= j_end; j++)
		{
			for(int i = i_start; i <= i_end; i++)
			{
				fout<<(*this)(i,j)<<" ";
			}
			fout<<std::endl;
		}
		fout.close();
	}
};

template<class T>
std::ostream& operator<<(std::ostream& output, const ARRAY_2D<T>& arr)
{
	for(int i = 0; i < arr.ij_res; i++)
		output<<arr.val[i]<<" ";
	return output;
}