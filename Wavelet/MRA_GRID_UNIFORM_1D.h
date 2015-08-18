#pragma once

#include "ARRAY.h"
#include "ARRAY_1D.h"
#include <iostream>
#include <math.h>

#define ITERATION_GRID_INT_1D(grid, width) for(int i = grid.i_start - width * grid.dx; i <= grid.i_end + width * grid.dx; i += grid.dx)
class GRID_INT_1D
{
public:
	int i_start, i_end, i_res;
	int dx, dx2;

public:
	GRID_INT_1D()
	{
		i_start = 0;
		i_end = 0;
		i_res = 0;
		dx = 0;
		dx2 = 0;
	}
	GRID_INT_1D(const int& i_start_, const int& i_end_, const int& i_res_)
	{
		Initialize(i_start_, i_end_, i_res_);
	}
	~GRID_INT_1D(){}

public:
	void Initialize(const int& i_start_, const int& i_end_, const int& i_res_)
	{
		i_start = i_start_;
		i_end = i_end_;
		i_res = i_res_;
		dx = (i_end - i_start) / (i_res-1);
		dx2 = dx * 2;
	}

	bool IsInside(const int& i, const int& width = 0) const
	{
		if(i < i_start - width || i_end + width < i) return false;
		else return true;
	}
};

std::ostream& operator<< (std::ostream& input, const GRID_INT_1D& grid)
{
	return input<<grid.i_start<<" "<<grid.i_end<<" "<<grid.i_res<<" "<<grid.dx;
}

//NOTE : no need to know physical coordinate.
class MRA_GRID_INT_1D
{
public:
	int num_levels;
	int Nx_coarse;
	int Nx;

	ARRAY<int> dx;
	ARRAY<GRID_INT_1D> grid, grid_x;

	GRID_INT_1D grid_finest;

public:
	MRA_GRID_INT_1D(){}
	~MRA_GRID_INT_1D(){}

public:
	void Initialize(const int& num_levels_, const int& Nx_coarse_)
	{
		num_levels = num_levels_;

		//the number of grid points.
		// x = 0, 1, 2, 3, ..., Nx-1.
		Nx_coarse = Nx_coarse_;
		Nx = (Nx_coarse-1) * int(pow(2.0, num_levels-1)) + 1;

		dx.Initialize(num_levels);
		grid.Initialize(num_levels);
		grid_x.Initialize(num_levels-1);
		for(int lev = num_levels-1; lev >= 0; lev --)
		{
			dx[lev] = int(pow(2.0, num_levels - lev - 1));
			//NOTE : start indices of the grid is not '0' but '1'!!!! In FWT, I want to not manipulate boundary layer with width 1.
			//       But, grid_finest starts from 0.
			if(lev == num_levels-1)
			{
				grid[lev].Initialize(0,(Nx-1), Nx);
			}
			else
			{
				grid[lev].Initialize(0, Nx-1, (Nx_coarse-1) * int(pow(2.0, lev)) + 1);
				grid_x[lev].Initialize(dx[lev+1], (Nx-1) - dx[lev+1], int((Nx-1) / pow(2.0, num_levels-1-lev)));
			}
		}
		grid_finest = grid[num_levels-1];
	}
};