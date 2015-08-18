#pragma once

#include "FIELD_UNIFORM_1D.h"
#include "MRA_GRID_UNIFORM_1D.h"
#include "VECTOR_3D.h"

class AWCM_MANAGER_1D
{
public:
	AWCM_MANAGER_1D(){}
	~AWCM_MANAGER_1D(){}

public:
	static void AddAdjacentZone(FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra)
	{
		const int max_lev = mra.num_levels-1;
		FIELD_UNIFORM_1D<bool> mask_temp(mask);
		
		ITERATION_GRID_INT_1D(mra.grid_x[max_lev-1], 0)
		{
			if(mask_temp(i) == true)
			{
				mask(i-1) = true;
				mask(i+1) = true;
			}
		}


		for (int lev = max_lev - 2; lev >= 0; lev--)
		{
			const int dx = mra.grid[lev].dx;
			const int dx_half = dx/2;
			const int dx_quarter = dx/4;
			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)//TODO : 0 is correct. but -1 is not sufficient? Must check this.
			{
				if (mask_temp(i) == true)
				{
					mask(i - dx_quarter) = true;
					mask(i + dx_quarter) = true;
				}
			}
		}
	}

	static void CheckReconstruction(FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra)
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end;

		for(int lev = max_lev-1; lev >= 0; lev --)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask(i) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				if(i_start < 0) i_start = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					mask(i_start + cnt*dx) = true;
				}
			}
		}
	}

	static void AddPointForDifferentiation(FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<int>& level)
	{
		FIELD_UNIFORM_1D<bool> mask_temp = mask;
		const int Nx = mra.grid_finest.i_end;
		const int max_lev = mra.num_levels-1;
		ITERATION_GRID_1D(mask, 0)
		{
			if(mask_temp(i) == false) continue;
			const int lev = level(i);
			const int dx = mra.grid[lev].dx;
			int i_start = i - 2*dx, i_end = i + 2*dx;
			
			if(lev == max_lev)
			{
				for(int cnt = -3; cnt < 4; cnt ++) mask(i + cnt) = true;
			}

			else
			{
				if(i_start <= 0)
				{
					for(int cnt = 0; cnt < 6; cnt ++) mask(cnt*dx) = true;
				}
				else if(i_end >= Nx)
				{
					for(int cnt = 0; cnt < 6; cnt ++) mask(Nx - cnt*dx) = true;
				}
				else
				{

					for(int cnt = 0; cnt < 5; cnt ++) mask(i_start + cnt*dx) = true;
				}
			}
		}
	}

	static void ComputeLevel(FIELD_UNIFORM_1D<int>& level, const FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra)
	{
		const int max_lev = mra.num_levels-1;
		ITERATION_GRID_1D(mask, -1)
		{
			if(mask(i) == false) 
			{
				level(i) = 0;
				continue;
			}

			int dx = 1;
			for(int lev = max_lev; lev >= 0; lev --)
			{
				if(mask(i-dx) == true || mask(i+dx) == true)
				{
					level(i) = lev;
					break;
				}
				dx *= 2;
			}
		}
		level(0) = max_lev;
		level(mra.grid_finest.i_end) = max_lev;
	}

	static void RemoveResidue(FIELD_UNIFORM_1D<double>& q, const FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra)
	{
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask(i) == false)
			{
				q(i) = 0;
			}
		}
	}

	static void RemoveResidue(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const FIELD_UNIFORM_1D<bool>& mask, const MRA_GRID_INT_1D& mra)
	{
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask(i) == false)
			{
				q(i)[0] = 0;
				q(i)[1] = 0;
				q(i)[2] = 0;
			}
		}
	}

	template<class T>
	static void Print(const FIELD_UNIFORM_1D<T>& q, const FIELD_UNIFORM_1D<bool>& mask, const char* filename)
	{
		std::ofstream fout(filename);
		const GRID_UNIFORM_1D& grid(q.base_grid);
		ITERATION_GRID_1D(q, 0)
		{
			if(mask(i) == true)
			{
				fout<<grid(i)<<" "<<q(i)<<std::endl;
			}
		}
		fout.close();
	}

	static void Print(const FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const FIELD_UNIFORM_1D<bool>& mask, const char* filename)
	{
		std::ofstream fout(filename);
		const GRID_UNIFORM_1D& grid(q.base_grid);
		ITERATION_GRID_1D(q, 0)
		{
			if(mask(i) == true)
			{
				fout<<grid(i)<<" "<<q(i)<<std::endl;
			}
		}
		fout.close();
	}

};