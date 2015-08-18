#pragma once

#include "AWCM_MANAGER_1D.h"
#include "ARRAY_2D.h"
#include "VECTOR_3D.h"

using std::cout;
using std::endl;

class AWCM_1D;

class WAVELET_1D
{
public:
	const int N;
	ARRAY_2D<double> weight;

public:
	WAVELET_1D(const int& N_ = 2):N(N_)
	{
		Initialize();
	}
	~WAVELET_1D(){}

public:
	void Initialize()
	{
		//evaluate Lagrangian basis function l_i(j+1/2).
		weight.Initialize(0,-1, 2*N, 2*N+1);
		weight.AssignAllValues(1);
		for(int j = -1; j < 2*N; j ++)
		{
			for(int i = 0; i < 2*N; i ++)
			{
				for(int k = 0; k < 2*N; k ++)
				{
					if(k == i) continue;
					weight(i, j) *= (j + 0.5 - k) / (i - k);
				}
			}
		}
	}

	void ForwardWaveletTransform(FIELD_UNIFORM_1D<double>& q, const MRA_GRID_INT_1D& mra) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;

		for(int lev = max_lev-1; lev >= 0; lev --)
		{

			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) -= weight(cnt, k) * q(i_start + cnt*dx);
				}
				q(i) *= 0.5;
			}

			ITERATION_GRID_INT_1D(mra.grid[lev], -1)
			{
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0)
				{
					i_start = dx_half, k = 0;
				}
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx - dx_half, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) += weight(cnt, k) * q(i_start + cnt * dx);
				}
			}
		}
	}

	void ForwardWaveletTransform(FIELD_UNIFORM_1D<double>& q, const MRA_GRID_INT_1D& mra, FIELD_UNIFORM_1D<bool>& mask, const double& tol) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;

		for(int lev = max_lev-1; lev >= 0; lev --)
		{
			
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask(i) == false) continue;
				
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) -= weight(cnt, k) * q(i_start + cnt*dx);
				}
				q(i) *= 0.5;

				if(ABS(q(i)) < tol)
				{
					mask(i) = false;
					//q(i) = 0;
				}
			}

			ITERATION_GRID_INT_1D(mra.grid[lev], -1)
			{
				if(mask(i) == false) continue;
				
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0)
				{
					i_start = dx_half, k = 0;
				}
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx - dx_half, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) += weight(cnt, k) * q(i_start + cnt * dx);
				}
			}
		}
	}

	void InverseWaveletTransform(FIELD_UNIFORM_1D<double>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;

		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid[lev], -1)
			{
				if(mask(i) == false)
				{
					q(i) = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = dx_half, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx - dx_half, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) -= weight(cnt, k) * q(i_start + cnt * dx);
				}
			}

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask(i) == false)
				{
					q(i) = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				q(i) *= 2;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) += weight(cnt, k) * q(i_start + cnt*dx);
				}
			}
		}
	}

	void Interpolate(FIELD_UNIFORM_1D<double>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask_from, const FIELD_UNIFORM_1D<bool>& mask_to) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;
		FIELD_UNIFORM_1D<bool> mask_temp = mask_to;
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask_from(i) == true) mask_temp(i) = false;
		}
		AWCM_MANAGER_1D::CheckReconstruction(mask_temp, mra);
		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask_temp(i) == false || mask_from(i) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				q(i) = 0;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) += weight(cnt, k) * q(i_start + cnt*dx);
				}
			}
		}
	}

	void Interpolate(FIELD_UNIFORM_1D<double>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask_from) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;
		FIELD_UNIFORM_1D<bool> mask_temp = mask_from;
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask_from(i) == true) mask_temp(i) = false;
			else
			{
				mask_temp(i) = true;
				//q(i) = 0;
			}
		}
		AWCM_MANAGER_1D::CheckReconstruction(mask_temp, mra);
		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask_temp(i) == false || mask_from(i) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				q(i) = 0;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i) += weight(cnt, k) * q(i_start + cnt*dx);
				}
			}
		}
	}


	////////////////////////////////////////////////Euler equation////////////////////////////////////////



	void ForwardWaveletTransformEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const MRA_GRID_INT_1D& mra, FIELD_UNIFORM_1D<bool>& mask, const double& tol) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;

		for(int lev = max_lev-1; lev >= 0; lev --)
		{

			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask(i) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i)[0] -= weight(cnt, k) * q(i_start + cnt*dx)[0];
				}
				q(i)[0] *= 0.5;

				if(ABS(q(i)[0]) < tol)
				{
					mask(i) = false;
				}
			}

			ITERATION_GRID_INT_1D(mra.grid[lev], -1)
			{
				if(mask(i) == false) continue;

				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0)
				{
					i_start = dx_half, k = 0;
				}
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx - dx_half, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i)[0] += weight(cnt, k) * q(i_start + cnt * dx)[0];
				}
			}
		}
	}

	void InverseWaveletTransformEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;

		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid[lev], -1)
			{
				if(mask(i) == false)
				{
					q(i)[0] = 0;
					//q(i)[1] = 0;
					//q(i)[2] = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = dx_half, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx - dx_half, k = 2;

				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i)[0] -= weight(cnt, k) * q(i_start + cnt * dx)[0];
				}
			}

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask(i) == false)
				{
					q(i)[0] = 0;
					//q(i)[1] = 0;
					//q(i)[2] = 0;
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				q(i)[0] *= 2;
				for(int cnt = 0; cnt < 4; cnt ++)
				{
					q(i)[0] += weight(cnt, k) * q(i_start + cnt*dx)[0];
				}
			}
		}
	}

	void InterpolateEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask_from, const FIELD_UNIFORM_1D<bool>& mask_to) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;
		FIELD_UNIFORM_1D<bool> mask_temp = mask_to;
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask_from(i) == true) mask_temp(i) = false;
		}
		AWCM_MANAGER_1D::CheckReconstruction(mask_temp, mra);
		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask_temp(i) == false || mask_from(i) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				for(int p = 1; p < 3; p ++)
				{
					q(i)[p] = 0;
					for(int cnt = 0; cnt < 4; cnt ++)
					{
						q(i)[p] += weight(cnt, k) * q(i_start + cnt*dx)[p];
					}
				}
			}
		}
	}

	void InterpolateEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const MRA_GRID_INT_1D& mra, const FIELD_UNIFORM_1D<bool>& mask_from) const
	{
		const int max_lev = mra.num_levels-1;
		const GRID_INT_1D& grid_finest = mra.grid_finest;
		int i_start, i_end, k;
		FIELD_UNIFORM_1D<bool> mask_temp = mask_from;
		ITERATION_GRID_INT_1D(mra.grid_finest, 0)
		{
			if(mask_from(i) == true) mask_temp(i) = false;
			else
			{
				mask_temp(i) = true;
				//q(i) = 0;
			}
		}
		AWCM_MANAGER_1D::CheckReconstruction(mask_temp, mra);
		for(int lev = 0; lev < max_lev; lev ++)
		{
			const int dx = mra.grid_x[lev].dx;
			const int dx_half = dx/2;

			ITERATION_GRID_INT_1D(mra.grid_x[lev], 0)
			{
				if(mask_temp(i) == false || mask_from(i) == true)
				{
					continue;
				}
				i_start = i - dx - dx_half;
				i_end = i + dx + dx_half;
				k = 1;
				if(i_start < 0) i_start = 0, k = 0;
				else if(i_end > grid_finest.i_end) i_start = grid_finest.i_end - 3*dx, k = 2;

				// maybe p = 1 will be OK.
				for(int p = 0; p < 3; p ++)
				{
					q(i)[p] = 0;
					for(int cnt = 0; cnt < 4; cnt ++)
					{
						q(i)[p] += weight(cnt, k) * q(i_start + cnt*dx)[p];
					}
				}
			}
		}
	}

};