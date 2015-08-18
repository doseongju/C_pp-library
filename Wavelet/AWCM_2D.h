#pragma once

#include "WAVELET_2D.h"
#include "AWCM_MANAGER_2D.h"
#include "EULER_EQUATION.h"
#include "Problems/2D/PROBLEM_WINDTUNNEL.h"
#include "Problems/2D/PROBLEM_DOUBLEMACH.h"
#include "Problems/2D/PROBLEM_SHOCK_VORTEX_INTERACTION.h"

class AWCM_2D
{
public:
	WAVELET_2D wavelet;
	MRA_GRID_INT_2D& mra_grid;
	FIELD_UNIFORM_2D<bool> mask_temp, mask0, mask1, mask2, mask3;
	FIELD_UNIFORM_2D<int> level0, level1, level2, level3;
	const double d1, d2;
	double d00, d01, d02;
	double d10, d11, d12, d13;
	double d20, d21, d22, d23, d24;
	double d30, d31, d32, d33, d34, d35;

public:
	AWCM_2D(MRA_GRID_INT_2D& mra_grid_, const double& gamma = 1.4) :mra_grid(mra_grid_), d1(2.0 / 3.0), d2(-1.0 / 12.0)
	{
		mask_temp.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny, 3);
		mask0.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny, 3);
		mask1.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny, 3);
		mask2.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny, 3);
		mask3.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny, 3);

		level0.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny);
		level1.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny);
		level2.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny);
		level3.Initialize(0, 0, mra_grid.Nx, mra_grid.Ny, 0, 0, mra_grid.Nx, mra_grid.Ny);

		mask0.AssignAllValues(true);
		level0.AssignAllValues(mra_grid.num_levels - 1);
		level1.AssignAllValues(mra_grid.num_levels - 1);
		level2.AssignAllValues(mra_grid.num_levels - 1);
		level3.AssignAllValues(mra_grid.num_levels - 1);

		d00 = -11.0 / 6.0; d01 = -1.0 / 3.0; d02 =  1.0 / 12.0;
		d10 =  3.0       ; d11 = -1.0 / 2.0; d12 = -2.0 / 3.0 ; d13 =  1.0 / 12.0;
		d20 = -3.0 / 2.0 ; d21 =  1.0      ; d22 =  0         ; d23 = -2.0 / 3.0 ; d24 =  1.0 / 12.0;
		d30 =  1.0 / 3.0 ; d31 = -1.0 / 6.0; d32 =  2.0 / 3.0 ; d33 =  0         ; d34 = -2.0 / 3.0 ; d35 = 1.0 / 12.0;

	}
	~AWCM_2D(){}

public:
	void Prepare(FIELD_UNIFORM_2D<double>& q, const double& tol = 1e-6)
	{
		wavelet.ForwardWaveletTransformation(q, mra_grid, mask0, tol);

		//!!NOTE : Banding boundary process should be handled properly according to the simulations.
		BandingBoundary(mask0);

		AWCM_MANAGER_2D::AddAdjacentZone(mask0, mra_grid);
		AWCM_MANAGER_2D::CheckReconstruction(mask0, mra_grid);

		mask1 = mask0;
		AWCM_MANAGER_2D::ComputeLevel(level0, mask0, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask1, mra_grid, level0);
		AWCM_MANAGER_2D::CheckReconstruction(mask1, mra_grid);


		mask2 = mask1;
		AWCM_MANAGER_2D::ComputeLevel(level1, mask1, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask2, mra_grid, level1);
		AWCM_MANAGER_2D::CheckReconstruction(mask2, mra_grid);

		mask3 = mask2;
		AWCM_MANAGER_2D::ComputeLevel(level2, mask2, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask3, mra_grid, level2);
		AWCM_MANAGER_2D::CheckReconstruction(mask3, mra_grid);

		wavelet.InverseWaveletTransformation(q, mra_grid, mask3);
	}

	void PostProcess(FIELD_UNIFORM_2D<double>& q)
	{
		AWCM_MANAGER_2D::RemoveResidue(q, mask0, mra_grid);
	}
	void BandingBoundary(FIELD_UNIFORM_2D<bool>& mask)
	{
		ITERATION_GRID_2D(mask, 0)
		{
			// 2 is a parameter can be tuned.
			if (2 <= i && i <= mask.i_end - 2 && j >= 2 && j <= mask.j_end - 2) continue;
			mask(i, j) = true;
		}
	}

	void Differentiate(const double& ll, const double& l, const double& c, const double& r, const double& rr, double& value)
	{
		value = d1 * (r - l) + d2 * (rr - ll);
	}
	void DifferentiateAtLeftBoundary(const double& l, const double& c, const double& r, const double& rr, double& value)
	{
		//value = d01 * l + d11 * c + d21 * r + d31 * rr;
		value = d1 * (r - l) + d2 * (rr - l);

	}
	void DifferentiateAtRightBoundary(const double& ll, const double& l, const double& c, const double& r, double& value)
	{
		//value = d01 * r + d11 * c + d21 * l + d31 * ll;
		value = d1 * (r - l) + d2 * (r - ll);

	}

	/////////////////////////////////////////////Euler equation/////////////////////////////////////////////

	void PrepareEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q, const double& tol = 1e-6)
	{
		mask_temp = mask0;
		wavelet.ForwardWaveletTransformationEuler(q, mra_grid, mask0, tol);

		//!!Boundary must have the max_level!!
		//!!NOTE : Banding boundary process should be handled properly according to the simulations.
		//BOUNDARY_MANAGER_WT boundary_manager(q);
		BOUNDARY_MANAGER_SHOCK_VORTEX boundary_manager(q);
		//BOUNDARY_MANAGER_DMACH boundary_manager(q);
		boundary_manager.BandingMask(mask0);

		AWCM_MANAGER_2D::AddAdjacentZone(mask0, mra_grid);
		AWCM_MANAGER_2D::CheckReconstruction(mask0, mra_grid);

		mask1 = mask0;
		AWCM_MANAGER_2D::ComputeLevel(level0, mask0, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask1, mra_grid, level0);
		AWCM_MANAGER_2D::CheckReconstruction(mask1, mra_grid);


		mask2 = mask1;
		AWCM_MANAGER_2D::ComputeLevel(level1, mask1, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask2, mra_grid, level1);
		AWCM_MANAGER_2D::CheckReconstruction(mask2, mra_grid);

		mask3 = mask2;
		AWCM_MANAGER_2D::ComputeLevel(level2, mask2, mra_grid);
		AWCM_MANAGER_2D::AddPointForDifferentiation(mask3, mra_grid, level2);
		AWCM_MANAGER_2D::CheckReconstruction(mask3, mra_grid);

		wavelet.InverseWaveletTransformationEuler(q, mra_grid, mask3);
		wavelet.InterpolateEuler(q, mra_grid, mask_temp, mask3);
	}

	void PostProcessEuler(FIELD_UNIFORM_2D<VECTOR_4D<double>>& q)
	{
		AWCM_MANAGER_2D::RemoveResidueEuler(q, mask0, mra_grid);
	}

	void Differentiate(const VECTOR_4D<double>& ll, const VECTOR_4D<double>& l, const VECTOR_4D<double>& c, const VECTOR_4D<double>& r, const VECTOR_4D<double>& rr, VECTOR_4D<double>& value)
	{
		value = d1 * (r - l) + d2 * (rr - ll);

		/*euler_eqn.ComputeFu(ll, F_ll);
		euler_eqn.ComputeFu( l, F_l );
		euler_eqn.ComputeFu( c, F_c );
		euler_eqn.ComputeFu( r, F_r );
		euler_eqn.ComputeFu(rr, F_rr);

		value = d1 * (F_r - F_l) + d2 * (F_rr - F_ll);*/
	}
	void DifferentiateAtLeftBoundary(const VECTOR_4D<double>& l, const VECTOR_4D<double>& c, const VECTOR_4D<double>& r, const VECTOR_4D<double>& rr, VECTOR_4D<double>& value)
	{
		value = d1 * (r - l) + d2 * (rr - l);

	/*	euler_eqn.ComputeFu(l , F_l);
		euler_eqn.ComputeFu(c , F_c);
		euler_eqn.ComputeFu(r , F_r);
		euler_eqn.ComputeFu(rr, F_rr);
		
		value = d1 * (F_r - F_l) + d2 * (F_rr - F_l);*/
	}
	void DifferentiateAtRightBoundary(const VECTOR_4D<double>& ll, const VECTOR_4D<double>& l, const VECTOR_4D<double>& c, const VECTOR_4D<double>& r, VECTOR_4D<double>& value)
	{
		value = d1 * (r - l) + d2 * (r - ll);

		/*euler_eqn.ComputeFu(ll, F_ll);
		euler_eqn.ComputeFu( l, F_l );
		euler_eqn.ComputeFu( c, F_c );
		euler_eqn.ComputeFu( r, F_r );

		value = d1 * (F_r - F_l) + d2 * (F_r - F_ll);*/
	}
};