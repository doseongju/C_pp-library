#pragma once

#include "ARRAY.h"
#include "ARRAY_2D.h"
#include <iostream>
#include <math.h>

#define ITERATION_GRID_INT_2D(subgrid, width) for(int j = subgrid.j_start - width*subgrid.dy; j <= subgrid.j_end + width*subgrid.dy; j += subgrid.dy)\
												for(int i = subgrid.i_start - width*subgrid.dx; i <= subgrid.i_end + width*subgrid.dx; i += subgrid.dx)
class GRID_INT_2D
{
public:
	int i_start, i_end, i_res;
	int j_start, j_end, j_res;
	int dx, dx2;
	int dy, dy2;

public:
	GRID_INT_2D()
	{
		i_start = 0;
		j_start = 0;
		i_end = 0;
		j_end = 0;
		i_res = 0;
		j_res = 0;
		dx = 0;
		dy = 0;
		dx2 = 0;
		dy2 = 0;
	}
	GRID_INT_2D(const int& i_start_, const int& j_start_, const int& i_end_, const int& j_end_, const int& i_res_, const int& j_res_)
	{
		Initialize(i_start_, j_start_, i_end_, j_end_, i_res_, j_res_);
	}
	~GRID_INT_2D(){}

public:
	void Initialize(const int& i_start_, const int& j_start_, const int& i_end_, const int& j_end_, const int& i_res_, const int& j_res_)
	{
		i_start = i_start_;
		j_start = j_start_;
		i_end = i_end_;
		j_end = j_end_;
		i_res = i_res_;
		j_res = j_res_;
		dx = (i_end - i_start) / (i_res-1);
		dy = (j_end - j_start) / (j_res-1);
		dx2 = dx * 2;
		dy2 = dy * 2;

	}

	bool IsInside(const int& i, const int& j, const int& width = 0) const
	{
		if(i < i_start - width || i_end + width < i) return false;
		else if(j < j_start - width || j_end + width < j) return false;
		else return true;
	}
};

std::ostream& operator<< (std::ostream& input, const GRID_INT_2D& grid)
{
	return input<<grid.i_start<<" "<<grid.j_start<<" "<<grid.i_end<<" "<<grid.j_end<<" "<<grid.i_res<<" "<<grid.j_res<<" "<<grid.dx<<" "<<grid.dy;
}

//NOTE : no need to know physical coordinate.
class MRA_GRID_INT_2D
{
public:
	int num_levels;
	int Nx_coarse;
	int Ny_coarse;
	int Nx, Ny;

	ARRAY<int> mra;
	ARRAY<GRID_INT_2D> grid;
	ARRAY_2D<GRID_INT_2D> subgrid;

	GRID_INT_2D grid_finest;

public:
	MRA_GRID_INT_2D(){}
	~MRA_GRID_INT_2D(){}

public:
	void Initialize(const int& num_levels_, const int& Nx_coarse_, const int& Ny_coarse_)
	{
		num_levels = num_levels_;
		//the number of grid points.
		// x = 0, 1, 2, 3, ..., Nx-1.
		// y = 0, 1, 2, 3, ..., Ny-1.
		Nx_coarse = Nx_coarse_;
		Ny_coarse = Ny_coarse_;
		Nx = (Nx_coarse-1) * int(pow(2.0, num_levels-1)) + 1;
		Ny = (Ny_coarse-1) * int(pow(2.0, num_levels-1)) + 1;

		mra.Initialize(num_levels);		
		grid.Initialize(num_levels);
		subgrid.Initialize(0,1, num_levels-1, 3);
		for(int lev = num_levels-1; lev >= 0; lev --)
		{
			mra[lev] = int(pow(2.0, num_levels - lev - 1));
			//NOTE : start indices of the grid is not '0' but '1'!!!! In FWT, I want to not manipulate boundary layer with width 1.
			//       But, grid_finest starts from 0.
			if(lev == num_levels-1)
			{
				grid[lev].Initialize(0, 0, (Nx-1), (Ny-1), int((Nx_coarse-1) * pow(2.0, lev) + 1), int((Ny_coarse-1) * pow(2.0, lev) + 1));
			}
			else
			{
				grid[lev].Initialize(mra[lev], mra[lev], (Nx-1) - mra[lev] , (Ny-1) - mra[lev], int((Nx_coarse-1) * pow(2.0, lev) + 1) - 2, int((Ny_coarse-1) * pow(2.0, lev) + 1) - 2);
			}
			for(int mu = 1; mu <= 3; mu ++)
			{
				if(lev == num_levels-1) continue;
				if(mu == 1)
				{
					subgrid(lev, mu).Initialize(int(0.5*mra[lev]), mra[lev], int(Nx-0.5*mra[lev]-1), Ny-mra[lev]-1, int(Nx/(mra[lev])), int(Ny/(mra[lev])-1));

				}
				else if(mu == 2)
				{
					subgrid(lev, mu).Initialize(mra[lev], int(0.5*mra[lev]), Nx-mra[lev]-1, int(Ny-0.5*mra[lev]-1), int(Nx/(mra[lev])-1), int(Ny/(mra[lev])));
				}
				else
				{
					subgrid(lev, mu).Initialize(int(0.5*mra[lev]), int(0.5*mra[lev]), int(Nx-0.5*mra[lev]-1), int(Ny-0.5*mra[lev]-1), int(Nx/(mra[lev])), int(Ny/(mra[lev])));
				}

			}
		}
		grid_finest = grid[num_levels-1];
	}
};