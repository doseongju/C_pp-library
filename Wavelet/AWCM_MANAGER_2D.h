#pragma once

#include "MRA_GRID_UNIFORM_2D.h"
#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"

class AWCM_MANAGER_2D
{
public:
	AWCM_MANAGER_2D(){}
	~AWCM_MANAGER_2D(){}

public:
	static void AddAdjacentZone(FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid)
	{
		const int max_lev = mra_grid.num_levels - 1;
		FIELD_UNIFORM_2D<bool> mask_temp(mask);

		for (int mu = 1; mu <= 3; mu++)
		{
#pragma omp parallel for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(max_lev - 1, mu), 0)
			{
				if (mask_temp(i, j) == true)
				{
					mask(i - 1, j    ) = true;
					mask(i + 1, j    ) = true;
					mask(i    , j - 1) = true;
					mask(i    , j + 1) = true;
				}
			}
		}

		for (int lev = max_lev - 2; lev >= 0; lev--)
		{
			const int dx = mra_grid.grid[lev].dx;
			const int dy = mra_grid.grid[lev].dy;
			const int dx_half = dx / 2, dx_quarter = dx / 4;
			const int dy_half = dy / 2, dy_quarter = dy / 4;

			for (int mu = 1; mu <= 3; mu++)
			{
#pragma omp parallel for
				ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, mu), 0)
				{
					if (mask_temp(i, j) == true)
					{
						mask(i - dx_quarter, j             ) = true;
						mask(i + dx_quarter, j             ) = true;
						mask(i             , j - dy_quarter) = true;
						mask(i             , j + dy_quarter) = true;

						mask(i - dx_half, j          ) = true;
						mask(i + dx_half, j          ) = true;
						mask(i          , j - dy_half) = true;
						mask(i          , j + dy_half) = true;
					}
				}
			}
		}
	}

	static void CheckReconstruction(FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid)
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;
		int i_start, j_start, i_end, j_end;

		for (int lev = max_lev - 1; lev >= 0; lev--)
		{
			const int dx = mra[lev];
			const int dx_half = dx/2;
			const int dy = dx;
			const int dy_half = dy/2;

			//d3
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask(i, j) == false) continue;

				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				if (i_start < 0) i_start = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx;
				if (j_start < 0) j_start = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy;


				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					mask(m, j) = true;
					m += dx;
				}
				for (int cnt = 0; cnt < 4; cnt++)
				{
					mask(i, n) = true;
					n += dy;
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						mask(m, n) = true;
						m += dx;
					}
					n += dy;
				}
			}


			//d1
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask(i, j) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;

				if (i_start < 0) i_start = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx;

				int m = i_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					mask(m, j) = true;
					m += dx;
				}
			}

			//d2
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask(i, j) == false) continue;

				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;

				if (j_start < 0) j_start = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy;

				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					mask(i, n) = true;
					n += dy;
				}
			}

		}
	}

	static void AddPointForDifferentiation(FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<int>& level)
	{
		FIELD_UNIFORM_2D<bool> mask_temp(mask);
		const int Nx = mra_grid.grid_finest.i_end;
		const int Ny = mra_grid.grid_finest.j_end;
		const int max_lev = mra_grid.num_levels - 1;
#pragma omp parallel for
		ITERATION_GRID_2D(mask, -1)
		{
			if (mask_temp(i, j) == false) continue;
			const int lev = level(i, j);
			const int dx = mra_grid.grid[lev].dx;
			const int i_start = i - 2*dx, i_end = i + 2*dx;
			const int j_start = j - 2*dx, j_end = j + 2*dx;

			if (lev == max_lev)
			{
				for (int cnt = -3; cnt <= 3; cnt++)
				{
					mask(i + cnt, j) = true;
					mask(i, j + cnt) = true;
				}
			}
			else
			{
				if (i_start < 0)
				{
					for (int cnt = -1; cnt <= 2; cnt++) mask(i + dx*cnt, j) = true;
				}
				else if (i_end > Nx)
				{
					for (int cnt = -2; cnt <= 1; cnt++) mask(i + dx*cnt, j) = true;
				}
				else
				{
					for (int cnt = -2; cnt <= 2; cnt++) mask(i + dx*cnt, j) = true;
				}

				if (j_start < 0)
				{
					for (int cnt = -1; cnt <= 2; cnt++) mask(i, j + dx*cnt) = true;
				}
				else if (j_end > Ny)
				{
					for (int cnt = -2; cnt <= 1; cnt++) mask(i, j + dx*cnt) = true;
				}
				else
				{
					for (int cnt = -2; cnt <= 2; cnt++) mask(i, j + dx*cnt) = true;
				}
				
			}
		}
	}

	static void ComputeLevel(FIELD_UNIFORM_2D<int>& level, const FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid)
	{
		const int max_lev = mra_grid.num_levels - 1;
#pragma omp parallel for
		ITERATION_GRID_2D(mask, 0)
		{
			if (mask(i, j) == false)
			{
				level(i, j) = 0;
				continue;
			}

			int dx = 1;
			int lev;
			for (lev = max_lev; lev >= 0; lev--)
			{
				if (mask(i - dx, j) == true || mask(i + dx, j) == true ||
					mask(i, j - dx) == true || mask(i, j + dx) == true)
				{
					level(i, j) = lev;
					break;
				}
				dx *= 2;
			}
		}

	}

	template<class T>
	static void Print(const FIELD_UNIFORM_2D<T>& q, const FIELD_UNIFORM_2D<bool>& mask, char* filename)
	{
		std::ofstream fout(filename);
		fout.precision(15);
		const GRID_UNIFORM_2D& grid(q.grid);
		ITERATION_GRID_2D(q, 0)
		{
			if (mask(i, j) == false) continue;
			fout << grid(i, j).i << " " << grid(i, j).j << " " << q(i, j) << std::endl;
		}
		fout.close();
	}

	static void RemoveResidue(FIELD_UNIFORM_2D<double>& q, const FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid)
	{
		ITERATION_GRID_INT_2D(mra_grid.grid_finest, -1)
		{
			if (mask(i, j) == false)
			{
				q(i, j) = 0.0;
			}
		}
	}
	static void RemoveResidueEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const FIELD_UNIFORM_2D<bool>& mask, const MRA_GRID_INT_2D& mra_grid)
	{
#pragma omp parallel for
		ITERATION_GRID_INT_2D(mra_grid.grid_finest, -1)
		{
			if (mask(i, j) == false)
			{
				q(i, j)[0] = 0.0;
				q(i, j)[1] = 0.0;
				q(i, j)[2] = 0.0;
				q(i, j)[3] = 0.0;
			}
		}
	}
};