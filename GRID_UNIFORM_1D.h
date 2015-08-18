#pragma once
#include <iostream>
#include <fstream>

class GRID_UNIFORM_1D
{
public:
	double x_min, x_max;
	double dx, one_over_dx;
	int i_start, i_end, i_res;

public:
	GRID_UNIFORM_1D(){}
	GRID_UNIFORM_1D(const int& i_start_, const int& i_res_, const double& x_min_, const double& x_max_)
	{
		Initialize(i_start_, i_res_, x_min_, x_max_);
	}
	GRID_UNIFORM_1D(const GRID_UNIFORM_1D& grid)
	{
		i_start = grid.i_start;
		i_end   = grid.i_end;
		i_res   = grid.i_res;
		x_min   = grid.x_min;
		x_max   = grid.x_max;
		dx      = grid.dx;
		one_over_dx = grid.one_over_dx;
	}
	~GRID_UNIFORM_1D(){}

public:
	void Initialize(const int& i_start_, const int& i_res_, const double& x_min_, const double& x_max_)
	{
				i_start = i_start_;
		i_res   = i_res_;
		i_end   = i_start + i_res - 1;
		x_min   = x_min_;
		x_max   = x_max_;
		dx      = (x_max - x_min) / (i_res - 1);
		one_over_dx = 1/dx;
	}
	void Initialize(const GRID_UNIFORM_1D& grid)
	{
		Initialize(grid.i_start, grid.i_res, grid.x_min, grid.x_max);
	}
	double operator()(const int& i) const
	{
		return x_min + dx * (i - i_start);
	}

public:
	void PrintGrid(const char* filename)
	{
		std::ofstream fout(filename);
		for(int i = i_start; i <= i_end; i++)
		{
			fout<<x_min + dx * (i - i_start)<<" ";
		}
		fout.close();
	}

public:
	int ClampIndexI(const int& i) const
	{
		if(i < i_start) return i_start;
		else if(i > i_end) return i_end;
		else return i;
	}

	bool IsInside(const int& i) const
	{
		if(i_start <= i && i <= i_end) return true;
		else return false;
	}
};

std::ostream& operator<<(std::ostream& input, const GRID_UNIFORM_1D& grid)
{
	return input<<"i_start : "<<grid.i_start<<"  i_res : "<<grid.i_res<<"  x_min : "<<grid.x_min<<"  x_max : "<<grid.x_max<<
		"  dx : "<<grid.dx<<std::endl;
}