#pragma once

#include "WAVELET_1D.h"
#include "AWCM_MANAGER_1D.h"
#include "EULER_EQUATION.h"

class AWCM_1D
{
public:
	WAVELET_1D wavelet;
	MRA_GRID_INT_1D& mra;
	FIELD_UNIFORM_1D<bool> mask_temp, mask0, mask1, mask2, mask3;
	FIELD_UNIFORM_1D<int> level0, level1, level2, level3;
	const double d1, d2;
	double d00, d01, d02;
	double d10, d11, d12, d13;
	double d20, d21, d22, d23, d24;
	double d30, d31, d32, d33, d34, d35;

public:
	AWCM_1D(MRA_GRID_INT_1D& mra_):mra(mra_), d1(2.0/3.0), d2(-1.0/12.0)
	{
		mask_temp.Initialize(0, mra.Nx, 0, mra.Nx, 3);
		mask0.Initialize(0, mra.Nx, 0, mra.Nx, 3);
		mask1.Initialize(0, mra.Nx, 0, mra.Nx, 3);
		mask2.Initialize(0, mra.Nx, 0, mra.Nx, 3);
		mask3.Initialize(0, mra.Nx, 0, mra.Nx, 3);

		level0.Initialize(0, mra.Nx, 0, mra.Nx);
		level1.Initialize(0, mra.Nx, 0, mra.Nx);
		level2.Initialize(0, mra.Nx, 0, mra.Nx);
		level3.Initialize(0, mra.Nx, 0, mra.Nx);

		mask0.AssignAllValues(true);

		d00 = -11.0 / 6.0; d01 = -1.0 / 3.0; d02 =  1.0 / 12.0;
		d10 =   3.0      ; d11 = -1.0 / 2.0; d12 = -2.0 /  3.0; d13 =  1.0 / 12.0;
		d20 =  -3.0 / 2.0; d21 =  1.0      ; d22 = 0          ; d23 = -2.0 / 3.0 ; d24 =  1.0 / 12.0;
		d30 =   1.0 / 3.0; d31 = -1.0 / 6.0; d32 =  2.0 /  3.0; d33 = 0          ; d34 = -2.0 /  3.0; d35 = 1.0 / 12.0;

	}
	~AWCM_1D(){}

public:
	void Prepare(FIELD_UNIFORM_1D<double>& q, const double& tol = 1e-6)
	{

		wavelet.ForwardWaveletTransform(q, mra, mask0, tol);
		BadingBoundary(mask0);

		AWCM_MANAGER_1D::AddAdjacentZone(mask0, mra);
		AWCM_MANAGER_1D::CheckReconstruction(mask0, mra);

		mask1 = mask0;
		AWCM_MANAGER_1D::ComputeLevel(level0, mask0, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask1, mra, level0);
		AWCM_MANAGER_1D::CheckReconstruction(mask1, mra);

		mask2 = mask1;
		AWCM_MANAGER_1D::ComputeLevel(level1, mask1, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask2, mra, level1);
		AWCM_MANAGER_1D::CheckReconstruction(mask2, mra);

		mask3 = mask2;
		AWCM_MANAGER_1D::ComputeLevel(level2, mask2, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask3, mra, level2);
		AWCM_MANAGER_1D::CheckReconstruction(mask3, mra);

		AWCM_MANAGER_1D::ComputeLevel(level3, mask3, mra);

		wavelet.InverseWaveletTransform(q, mra, mask3);
	}
	void BadingBoundary(FIELD_UNIFORM_1D<bool>& mask)
	{
		ITERATION_GRID_1D(mask, 0)
		{
			if (3 <= i && i <= mask.i_end-3) continue;
			mask(i) = true;
		}
	}
	void PostProcess(FIELD_UNIFORM_1D<double>& q)
	{
		AWCM_MANAGER_1D::RemoveResidue(q, mask0, mra);
	}

	void PrepareEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q, const double& tol = 1e-6)
	{

		mask_temp = mask0;

		wavelet.ForwardWaveletTransformEuler(q, mra, mask0, tol);
		BadingBoundary(mask0);

		AWCM_MANAGER_1D::AddAdjacentZone(mask0, mra);
		AWCM_MANAGER_1D::CheckReconstruction(mask0, mra);

		mask1 = mask0;
		AWCM_MANAGER_1D::ComputeLevel(level0, mask0, mra);
		AWCM_MANAGER_1D::AddAdjacentZone(mask1, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask1, mra, level0);
		AWCM_MANAGER_1D::CheckReconstruction(mask1, mra);


		mask2 = mask1;
		AWCM_MANAGER_1D::ComputeLevel(level1, mask1, mra);
		AWCM_MANAGER_1D::AddAdjacentZone(mask2, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask2, mra, level1);
		AWCM_MANAGER_1D::CheckReconstruction(mask2, mra);

		mask3 = mask2;
		AWCM_MANAGER_1D::ComputeLevel(level2, mask2, mra);
		AWCM_MANAGER_1D::AddAdjacentZone(mask3, mra);
		AWCM_MANAGER_1D::AddPointForDifferentiation(mask3, mra, level2);
		AWCM_MANAGER_1D::CheckReconstruction(mask3, mra);

		AWCM_MANAGER_1D::ComputeLevel(level3, mask3, mra);

		wavelet.InverseWaveletTransformEuler(q, mra, mask3);
		wavelet.InterpolateEuler(q, mra, mask_temp, mask3);
	}

	void PostProcessEuler(FIELD_UNIFORM_1D<VECTOR_3D<double>>& q)
	{
		AWCM_MANAGER_1D::RemoveResidue(q, mask0, mra);
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
		//value = d02 * ll + d12 * l + d22 * c + d32 * r + d2 * rr;
		value = d1 * (r - l) + d2 * (r - ll);
	}

	void DifferentiateEuler(const VECTOR_3D<double>& ll, const VECTOR_3D<double>& l, const VECTOR_3D<double>& c, const VECTOR_3D<double>& r, const VECTOR_3D<double>& rr, VECTOR_3D<double>& value)
	{
		EULER_EQUATION_1D euler_eqn(1.4);
		VECTOR_3D<double> f_ll, f_l, f_c, f_r, f_rr;
		euler_eqn.ComputeFu(ll, f_ll);
		euler_eqn.ComputeFu(l , f_l );
		euler_eqn.ComputeFu(c , f_c );
		euler_eqn.ComputeFu(r , f_r );
		euler_eqn.ComputeFu(rr, f_rr);

		for(int p = 0; p < 3; p ++)
		{
			value[p] = d1 * (f_r[p] - f_l[p]) + d2 * (f_rr[p] - f_ll[p]);
		}
	}

	inline void DifferentiateEulerAtLeftBoundary(const VECTOR_3D<double>& l, const VECTOR_3D<double>& c, const VECTOR_3D<double>& r, const VECTOR_3D<double>& rr, VECTOR_3D<double>& value)
	{
		EULER_EQUATION_1D euler_eqn(1.4);
		VECTOR_3D<double> f_l, f_c, f_r, f_rr;
		euler_eqn.ComputeFu(l , f_l );
		euler_eqn.ComputeFu(c , f_c );
		euler_eqn.ComputeFu(r , f_r );
		euler_eqn.ComputeFu(rr, f_rr);

		value = d01 * f_l + d11 * f_c + d21 * f_r + d31 * f_rr;
	}
	inline void DifferentiateEulerAtRightBoundary(const VECTOR_3D<double>& ll, const VECTOR_3D<double>& l, const VECTOR_3D<double>& c, const VECTOR_3D<double>& r, VECTOR_3D<double>& value)
	{
		DifferentiateEulerAtLeftBoundary(r, c, l, ll, value);
	}

	inline void DifferentiateEulerAtLeftBoundary(const VECTOR_3D<double>& ll, const VECTOR_3D<double>& l, const VECTOR_3D<double>& c, const VECTOR_3D<double>& r, const VECTOR_3D<double>& rr, VECTOR_3D<double>& value)
	{
		EULER_EQUATION_1D euler_eqn(1.4);
		VECTOR_3D<double> f_ll, f_l, f_c, f_r, f_rr;
		euler_eqn.ComputeFu(ll, f_ll);
		euler_eqn.ComputeFu(l , f_l );
		euler_eqn.ComputeFu(c , f_c );
		euler_eqn.ComputeFu(r , f_r );
		euler_eqn.ComputeFu(rr, f_rr);

		value = d02 * f_ll + d12 * f_l + d22 * f_c + d32 * f_r + d2 * f_rr;
	}
	inline void DifferentiateEulerAtRightBoundary(const VECTOR_3D<double>& ll, const VECTOR_3D<double>& l, const VECTOR_3D<double>& c, const VECTOR_3D<double>& r, const VECTOR_3D<double>& rr, VECTOR_3D<double>& value)
	{
		DifferentiateEulerAtLeftBoundary(rr, r, c, l, ll, value);
	}
};