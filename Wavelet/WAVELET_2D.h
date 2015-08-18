#pragma once

#include "MRA_GRID_UNIFORM_2D.h"
#include "AWCM_MANAGER_2D.h"
#include "VECTOR_4D.h"
#include "ARRAY_1D.h"
#include <omp.h>
using std::cout;
using std::endl;
class WAVELET_2D
{
public:
	const int N;
	ARRAY_2D<double> weight;

public:
	WAVELET_2D(const int& N_ = 2) :N(N_)
	{
		//evaluate Lagrangian basis function l_i(j+1/2).
		weight.Initialize(0, -1, 2 * N, 2 * N + 1);
		weight.AssignAllValues(1);
		for (int j = -1; j < 2 * N; j++)
		{
			for (int i = 0; i < 2 * N; i++)
			{
				for (int k = 0; k < 2 * N; k++)
				{
					if (k == i) continue;
					weight(i, j) *= double(j + 0.5 - k) / double(i - k);
				}
			}
		}
	}
	~WAVELET_2D(){}

public:
	void ForwardWaveletTransformation(FIELD_UNIFORM_2D<double>& q, const MRA_GRID_INT_2D& mra_grid, FIELD_UNIFORM_2D<bool>& mask, const double& tol) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;

		int i_start, j_start, i_end, j_end, k, k1, k2;

		for(int lev = max_lev-1; lev >= 0; lev--)
		{
			const int dx = mra[lev];
			const int dx_half = dx/2;
			const int dy = dx;
			const int dy_half = dy/2;

			//d3
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev,3), 0)
			{
				if (mask(i,j) == false) continue;

				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end   = i + dx + dx_half;
				j_end   = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;
				

				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) -= q(m, j) * weight(cnt, k1);
					m += dx;
				}
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) -= q(i, n) * weight(cnt, k2);
					n += dy;
				}


				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						q(i, j) += q(m, n) * weight(cnt1, k1) * weight(cnt2, k2);
						m += dx;
					}
					n += dy;
				}
				q(i, j) *= 0.25;
				if (ABS(q(i, j)) < tol) mask(i, j) = false;
			}


			//d1
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask(i,j) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				int m = i_start;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i, j) -= q(m, j) * weight(cnt, k);
					m += dx;
				}
				q(i,j) *= 0.5;

				if (ABS(q(i, j)) < tol) mask(i, j) = false;
			}


			//d2
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev,2), 0)
			{
				if (mask(i,j) == false) continue;

				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				int n = j_start;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i, j) -= q(i, n) * weight(cnt, k);
					n += dy;
				}
				q(i,j) *= 0.5;

				if (ABS(q(i, j)) < tol) mask(i, j) = false;
			}
			


			//c
			ITERATION_GRID_INT_2D(mra_grid.grid[lev], 0)
			{
				if (mask(i,j) == false) continue;

				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end   = i + dx + dx_half;
				j_end   = j + dy + dy_half;
				k1 = 1;
				k2 = 1;

				if (i_start < 0)                    i_start = dx_half, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx - dx_half, k1 = 2;
				if (j_start < 0)                    j_start = dy_half, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy - dy_half, k2 = 2;


				int m = i_start;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i,j) += q(m,j) * weight(cnt, k1);
					m += dx;
				}

				int n = j_start;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i,j) += q(i,n) * weight(cnt, k2);
					n += dy;
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						q(i, j) += q(m, n) * weight(cnt1, k1) * weight(cnt2, k2);
						m += dx;
					}
					n += dy;
				}
			}

		}
	}


	void InverseWaveletTransformation(FIELD_UNIFORM_2D<double>& q, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<bool>& mask) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;

		int i_start, j_start, i_end, j_end, k, k1, k2;

		for (int lev = 0; lev < max_lev; lev++)
		{
			const int dx = mra[lev];
			const int dx_half = dx/2;
			const int dy = dx;
			const int dy_half = dy/2;

			//c
			ITERATION_GRID_INT_2D(mra_grid.grid[lev], 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j) = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = dx_half, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx - dx_half, k1 = 2;
				if (j_start < 0) j_start = dy_half, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy - dy_half, k2 = 2;


				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) -= q(m, j) * weight(cnt, k1);
					m += dx;
				}
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) -= q(i, n) * weight(cnt, k2);
					n += dx;
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						q(i, j) -= q(m, n) * weight(cnt1, k1) * weight(cnt2, k2);
						m += dx;
					}
					n += dx;
				}
			}




			//d1
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j) = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				q(i, j) *= 2.0;
				int m = i_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(m, j) * weight(cnt, k);
					m += dx;
				}

			}

			//d2
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j) = 0;
					continue;
				}
				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				q(i, j) *= 2.0;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(i, n) * weight(cnt, k);
					n += dy;
				}
			}


			//d3
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j) = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

				q(i, j) *= 4.0;
				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(m, j) * weight(cnt, k1);
					m += dx;
				}
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(i, n) * weight(cnt, k2);
					n += dy;
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						q(i, j) -= q(m, n) * weight(cnt1, k1) * weight(cnt2, k2);
						m += dx;
					}
					n += dy;
				}
			}
			
		}
	}

	void Interpolate(FIELD_UNIFORM_2D<double>& q, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<bool>& mask_from, const FIELD_UNIFORM_2D<bool>& mask_to) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;
		FIELD_UNIFORM_2D<bool> mask_temp = mask_to;
		ITERATION_GRID_INT_2D(grid_finest, 0)
		{
			if (mask_from(i, j) == true) mask_temp(i, j) = false;
		}
		AWCM_MANAGER_2D::CheckReconstruction(mask_temp, mra_grid);

		int i_start, j_start, i_end, j_end, k, k1, k2;
		for (int lev = 0; lev < max_lev; lev++)
		{
			const int dx = mra[lev];
			const int dx_half = dx / 2;
			const int dy = dx;
			const int dy_half = dy / 2;

			//d1
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				q(i, j) = 0.0;
				int m = i_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(m, j) * weight(cnt, k);
					m += dx;
				}

			}

			//d2
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				q(i, j) = 0.0;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(i, n) * weight(cnt, k);
					n += dy;
				}
			}

			//d3
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

				q(i, j) = 0.0;
				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(m, j) * weight(cnt, k1);
					m += dx;
				}
				for (int cnt = 0; cnt < 4; cnt++)
				{
					q(i, j) += q(i, n) * weight(cnt, k2);
					n += dy;
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++)
					{
						q(i, j) -= q(m, n) * weight(cnt1, k1) * weight(cnt2, k2);
						m += dx;
					}
					n += dy;
				}
			}

		}
	}



	///////////////////////////////////////Euler equation////////////////////////////////////////
	void ForwardWaveletTransformationEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const MRA_GRID_INT_2D& mra_grid, FIELD_UNIFORM_2D<bool>& mask, const double& tol) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;

		for (int lev = max_lev - 1; lev >= 0; lev--)
		{
			const int dx = mra[lev];
			const int dx_half = dx / 2;
			const int dy = dx;
			const int dy_half = dy / 2;

#pragma omp parallel
{
			int i_start, j_start, i_end, j_end, k, k1, k2;

			//d3
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask(i, j) == false) continue;

				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] -= q(m, j)[0] * weight(cnt, k1);

				}

				for (int cnt = 0; cnt < 4; cnt++, n += dy)
				{
					q(i, j)[0] -= q(i, n)[0] * weight(cnt, k2);
				}
				
				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++, n += dy)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++, m += dx)
					{
						q(i, j)[0] += q(m, n)[0] * weight(cnt1, k1) * weight(cnt2, k2);
					}
				}
				q(i, j)[0] *= 0.25;
				if (ABS(q(i, j)[0]) < tol) mask(i, j) = false;
			}

			//d1
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask(i, j) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				int m = i_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] -= q(m, j)[0] * weight(cnt, k);
				}
				q(i, j)[0] *= 0.5;

				if (ABS(q(i, j)[0]) < tol) mask(i, j) = false;

			}
			//d2
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask(i, j) == false) continue;

				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, n += dy)
				{
					q(i, j)[0] -= q(i, n)[0] * weight(cnt, k);
				}
				q(i, j)[0] *= 0.5;

				if (ABS(q(i, j)[0]) < tol) mask(i, j) = false;
			}
			//c
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.grid[lev], 0)
			{
				if (mask(i, j) == false) continue;

				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = dx_half, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx - dx_half, k1 = 2;
				if (j_start < 0) j_start = dy_half, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy - dy_half, k2 = 2;


				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] += q(m, j)[0] * weight(cnt, k1);
				}
				for (int cnt = 0; cnt < 4; cnt++, n += dx)
				{
					q(i, j)[0] += q(i, n)[0] * weight(cnt, k2);
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++, n += dx)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++, m += dx)
					{
						q(i, j)[0] += q(m, n)[0] * weight(cnt1, k1) * weight(cnt2, k2);
					}
				}
			}
}
		}
	}

	void InverseWaveletTransformationEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<bool>& mask) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;

		for (int lev = 0; lev < max_lev; lev++)
		{
			const int dx = mra[lev];
			const int dx_half = dx / 2;
			const int dy = dx;
			const int dy_half = dy / 2;

#pragma omp parallel
{
			int i_start, j_start, i_end, j_end, k, k1, k2;

			//c
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.grid[lev], 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j)[0] = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = dx_half, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx - dx_half, k1 = 2;
				if (j_start < 0) j_start = dy_half, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy - dy_half, k2 = 2;


				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] -= q(m, j)[0] * weight(cnt, k1);
				}
				for (int cnt = 0; cnt < 4; cnt++, n += dx)
				{
					q(i, j)[0] -= q(i, n)[0] * weight(cnt, k2);
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++, n += dx)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++, m += dx)
					{
						q(i, j)[0] -= q(m, n)[0] * weight(cnt1, k1) * weight(cnt2, k2);
					}
				}
			}


			//d2
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j)[0] = 0;
					continue;
				}
				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				q(i, j)[0] *= 2.0;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, n += dy)
				{
					q(i, j)[0] += q(i, n)[0] * weight(cnt, k);
				}
			}


			//d1
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j)[0] = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				q(i, j)[0] *= 2.0;
				int m = i_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] += q(m, j)[0] * weight(cnt, k);
				}
			}
			

			//d3
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask(i, j) == false)
				{
					q(i, j)[0] = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

				q(i, j)[0] *= 4.0;
				int m = i_start;
				int n = j_start;
				for (int cnt = 0; cnt < 4; cnt++, m += dx)
				{
					q(i, j)[0] += q(m, j)[0] * weight(cnt, k1);
				}
				for (int cnt = 0; cnt < 4; cnt++, n += dy)
				{
					q(i, j)[0] += q(i, n)[0] * weight(cnt, k2);
				}

				n = j_start;
				for (int cnt2 = 0; cnt2 < 4; cnt2++, n += dy)
				{
					m = i_start;
					for (int cnt1 = 0; cnt1 < 4; cnt1++, m += dx)
					{
						q(i, j)[0] -= q(m, n)[0] * weight(cnt1, k1) * weight(cnt2, k2);
					}
				}
			}
}

		}
	}

	void InterpolateEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<bool>& mask_from, const FIELD_UNIFORM_2D<bool>& mask_to) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;
		FIELD_UNIFORM_2D<bool> mask_temp = mask_to;
		ITERATION_GRID_INT_2D(grid_finest, 0)
		{
			if (mask_from(i, j) == true) mask_temp(i, j) = false;
		}
		AWCM_MANAGER_2D::CheckReconstruction(mask_temp, mra_grid);

		for (int lev = 0; lev < max_lev; lev++)
		{
			const int dx = mra[lev];
			const int dx_half = dx / 2;
			const int dy = dx;
			const int dy_half = dy / 2;

#pragma omp parallel
{
			int i_start, j_start, i_end, j_end, k, k1, k2;

			//d1
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;

				if (i_start < 0) i_start = 0, k = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

				for (int p = 1; p < 4; p++)
				{
					q(i, j)[p] = 0.0;
					int m = i_start;
					for (int cnt = 0; cnt < 4; cnt++, m += dx)
					{
						q(i, j)[p] += q(m, j)[p] * weight(cnt, k);
					}
				}

			}

			//d2
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				j_start = j - dy - dy_half;
				j_end = j + dy + dy_half;
				k = 1;

				if (j_start < 0) j_start = 0, k = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

				for (int p = 1; p < 4; p++)
				{
					q(i, j)[p] = 0.0;
					int n = j_start;
					for (int cnt = 0; cnt < 4; cnt++, n += dy)
					{
						q(i, j)[p] += q(i, n)[p] * weight(cnt, k);
					}
				}
			}

			//d3
#pragma omp for
			ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
			{
				if (mask_temp(i, j) == false || mask_from(i, j) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				j_start = j - dy - dy_half;
				i_end = i + dx + dx_half;
				j_end = j + dy + dy_half;
				k1 = 1;
				k2 = 1;
				if (i_start < 0) i_start = 0, k1 = 0;
				else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
				if (j_start < 0) j_start = 0, k2 = 0;
				else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

				for (int p = 1; p < 4; p++)
				{
					q(i, j)[p] = 0.0;
					int m = i_start;
					int n = j_start;
					for (int cnt = 0; cnt < 4; cnt++, m += dx)
					{
						q(i, j)[p] += q(m, j)[p] * weight(cnt, k1);
					}
					for (int cnt = 0; cnt < 4; cnt++, n += dy)
					{
						q(i, j)[p] += q(i, n)[p] * weight(cnt, k2);
					}

					n = j_start;
					for (int cnt2 = 0; cnt2 < 4; cnt2++, n += dy)
					{
						m = i_start;
						for (int cnt1 = 0; cnt1 < 4; cnt1++, m += dx)
						{
							q(i, j)[p] -= q(m, n)[p] * weight(cnt1, k1) * weight(cnt2, k2);
						}
					}
				}
			}
}

		}
	}




	void InterpolateAllGridEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const MRA_GRID_INT_2D& mra_grid, const FIELD_UNIFORM_2D<bool>& mask_from) const
	{
		const int max_lev = mra_grid.num_levels - 1;
		const GRID_INT_2D& grid_finest = mra_grid.grid_finest;
		const ARRAY<int>& mra = mra_grid.mra;
		FIELD_UNIFORM_2D<bool> mask_temp(mask_from);
		mask_temp.AssignAllValues(true);
		ITERATION_GRID_INT_2D(grid_finest, 0)
		{
			if (mask_from(i, j) == true) mask_temp(i, j) = false;
		}
		AWCM_MANAGER_2D::CheckReconstruction(mask_temp, mra_grid);

		for (int lev = 0; lev < max_lev; lev++)
		{
			const int dx = mra[lev];
			const int dx_half = dx / 2;
			const int dy = dx;
			const int dy_half = dy / 2;

#pragma omp parallel
			{
				int i_start, j_start, i_end, j_end, k, k1, k2;

				//d1
#pragma omp for
				ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 1), 0)
				{
					if (mask_temp(i, j) == false || mask_from(i, j) == true)
					{
						continue;
					}
					i_start = i - dx - dx_half;
					i_end = i + dx + dx_half;
					k = 1;

					if (i_start < 0) i_start = 0, k = 0;
					else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k = 2;

					for (int p = 0; p < 4; p++)
					{
						q(i, j)[p] = 0.0;
						int m = i_start;
						for (int cnt = 0; cnt < 4; cnt++)
						{
							q(i, j)[p] += q(m, j)[p] * weight(cnt, k);
							m += dx;
						}
					}

				}

				//d2
#pragma omp for
				ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 2), 0)
				{
					if (mask_temp(i, j) == false || mask_from(i, j) == true)
					{
						//q(i, j) = 0;
						continue;
					}
					j_start = j - dy - dy_half;
					j_end = j + dy + dy_half;
					k = 1;

					if (j_start < 0) j_start = 0, k = 0;
					else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k = 2;

					for (int p = 0; p < 4; p++)
					{
						q(i, j)[p] = 0.0;
						int n = j_start;
						for (int cnt = 0; cnt < 4; cnt++)
						{
							q(i, j)[p] += q(i, n)[p] * weight(cnt, k);
							n += dy;
						}
					}
				}

				//d3
#pragma omp for
				ITERATION_GRID_INT_2D(mra_grid.subgrid(lev, 3), 0)
				{
					if (mask_temp(i, j) == false || mask_from(i, j) == true)
					{
						//q(i, j) = 0;
						continue;
					}
					i_start = i - dx - dx_half;
					j_start = j - dy - dy_half;
					i_end = i + dx + dx_half;
					j_end = j + dy + dy_half;
					k1 = 1;
					k2 = 1;
					if (i_start < 0) i_start = 0, k1 = 0;
					else if (i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3 * dx, k1 = 2;
					if (j_start < 0) j_start = 0, k2 = 0;
					else if (j_end > grid_finest.j_end) j_start = grid_finest.j_end - 3 * dy, k2 = 2;

					for (int p = 0; p < 4; p++)
					{
						q(i, j)[p] = 0.0;
						int m = i_start;
						int n = j_start;
						for (int cnt = 0; cnt < 4; cnt++)
						{
							q(i, j)[p] += q(m, j)[p] * weight(cnt, k1);
							m += dx;
						}
						for (int cnt = 0; cnt < 4; cnt++)
						{
							q(i, j)[p] += q(i, n)[p] * weight(cnt, k2);
							n += dy;
						}

						n = j_start;
						for (int cnt2 = 0; cnt2 < 4; cnt2++)
						{
							m = i_start;
							for (int cnt1 = 0; cnt1 < 4; cnt1++)
							{
								q(i, j)[p] -= q(m, n)[p] * weight(cnt1, k1) * weight(cnt2, k2);
								m += dx;
							}
							n += dy;
						}
					}
				}
			}

		}
	}
};