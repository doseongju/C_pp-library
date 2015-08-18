#pragma once
#include <fstream>
#include <cassert>
#include "VECTOR_2D.h"

class GRID_UNIFORM_2D
{
public:
	int i_start, i_res, i_end,
		j_start, j_res, j_end;
	double x_min, x_max,
		   y_min, y_max;
	double dx, one_over_dx,
		   dy, one_over_dy;

public:
	GRID_UNIFORM_2D():i_start(0), i_res(0), i_end(0), j_start(0), j_res(0), j_end(0),
					  x_min(0), x_max(0), y_min(0), y_max(0), dx(0), dy(0), one_over_dx(0), one_over_dy(0){}
	GRID_UNIFORM_2D(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_,
				    const double& x_min_, const double& y_min_, const double& x_max_, const double& y_max_)
	{
		Initialize( i_start_,  j_start_,  i_res_,  j_res_, x_min_,  y_min_,  x_max_,  y_max_);
	}

	~GRID_UNIFORM_2D(){}

public:
	void Initialize(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_,
				    const double& x_min_, const double& y_min_, const double& x_max_, const double& y_max_)
	{
		i_start = i_start_;
		i_res = i_res_;
		i_end = i_start + i_res - 1;
		j_start = j_start_;
		j_res = j_res_;
		j_end = j_start + j_res - 1;

		x_min = x_min_; x_max = x_max_;
		y_min = y_min_; y_max = y_max_;
		dx = (x_max - x_min) / (i_res - 1); one_over_dx = 1 / dx;
		dy = (y_max - y_min) / (j_res - 1); one_over_dy = 1 / dy;
	}

	void Initialize(const GRID_UNIFORM_2D& grid)
	{
		Initialize(grid.i_start, grid.j_start, grid.i_res, grid.j_res, grid.x_min, grid.y_min, grid.x_max, grid.y_max);
	}

public:
	VECTOR_2D<double> operator()(const int& i, const int& j) const
	{
		return VECTOR_2D<double>(x_min + (i - i_start) * dx, y_min + (j - j_start) * dy);
	}

public:
	const void PrintGridX(const char* filename) const
	{
		std::ofstream fout(filename);
		for(int i = i_start; i <= i_end; i++)
		{
			fout<<x_min + dx * (i - i_start)<<" ";
		}
		fout.close();
	}
	const void PrintGridY(const char* filename) const
	{
		std::ofstream fout(filename);
		for(int j = j_start; j <= j_end; j++)
		{
			fout<<y_min + dy * (j - j_start)<<" ";
		}
		fout.close();
	}

	GRID_UNIFORM_2D MakeMACGrid(const bool& x_MAC = false, const bool& y_MAC = false) const
	{
		int i_start_(i_start), j_start_(j_start), i_res_(i_res), j_res_(j_res);
		double x_min_(x_min), y_min_(y_min), x_max_(x_max), y_max_(y_max);
		
		if(x_MAC == true)
		{
			i_start_ --;
			i_res_ ++;
			x_min_ -= 0.5 * dx;
			x_max_ += 0.5 * dx;
		}
		if(y_MAC == true)
		{
			j_start_ --;
			j_res_ ++;
			y_min_ -= 0.5 * dy;
			y_max_ += 0.5 * dy;
		}

		return GRID_UNIFORM_2D(i_start_, j_start_, i_res_, j_res_, x_min_, y_min_, x_max_, y_max_);
	}

	int ClampIndexI(const int& i) const
	{
		if(i < i_start) return i_start;
		else if(i > i_end) return i_end;
		else return i;
	}

	int ClampIndexJ(const int& j) const
	{
		if(j < j_start) return j_start;
		else if(j > j_end) return j_end;
		else return j;
	}

	VECTOR_2D<int> LeftBottomCell(const double& x, const double& y) const
	{
		int i, j;
		i = int((x - x_min) / dx + i_start);
		j = int((y - y_min) / dy + j_start);

		return VECTOR_2D<int>(i,j);
	}

	void IsInside(const int& i, const int& j) const
	{
		assert(i_start <= i && i <= i_end);
		assert(j_start <= j && j <= j_end);
	}

	bool IsInside(const int& i, const int& j, const int& width) const
	{
		if(i_start - width <= i && i <= i_end + width )
		{
			if(j_start - width <= j && j <= j_end + width)
			{
				return true;
			}
		}
		return false;
	}
};

std::ostream& operator<< (std::ostream& output, const GRID_UNIFORM_2D grid)
{
	return output<< "i_start : " << grid.i_start<< " i_res : "<<grid.i_res<< " x_min : "<<grid.x_min<<" x_max : "<<grid.x_max<<" dx : "<<grid.dx<<std::endl<<
		"j_start : " << grid.j_start << " j_res : "<<grid.j_res<<" y_min : "<<grid.y_min<<" y_max : "<<grid.y_max<<" dy : "<<grid.dy<<std::endl;

}