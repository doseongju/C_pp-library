#pragma once
#include <iostream>
#include <cassert>
#include <fstream>
#include "MACROS.h"

template<class T>
class ARRAY_1D
{
public:
	int i_start, i_end, i_res;
	T* val;

public:
	ARRAY_1D():val(0), i_start(0), i_end(0), i_res(0){}
	ARRAY_1D(const int& i_start_, const int& i_res_): val(0)
	{
		Initialize(i_start_, i_res_, true);
	}
	ARRAY_1D(const ARRAY_1D<T>& arr_): val(0)
	{
		Initialize(arr_.i_start, arr_.i_res);
	}
	~ARRAY_1D()
	{
		SAFE_DELETE_ARRAY(val);
	}

public:
	void Initialize(const int& i_start_, const int& i_res_, const bool& initialize = false)
	{
		i_start = i_start_;
		i_res   = i_res_;
		i_end   = i_start + i_res - 1;

		SAFE_DELETE_ARRAY(val);
		val = new T [i_res];

		if(initialize)
		{
			for(int i = 0; i < i_res; i++) val[i] = T();
		}
	}
	
	T& operator()(const int& i)
	{
		assert(i >= i_start && i <= i_end);
		return val[i - i_start];
	}

	const T& operator()(const int& i) const
	{
		assert(i >= i_start && i <= i_end);
		return val[i - i_start];
	}

	void operator=(const ARRAY_1D<T>& arr_)
	{
		Initialize(arr_.i_start, arr_.i_res);
		for(int i = 0; i < i_res; i++) val[i] = arr_.val[i];
	}

};