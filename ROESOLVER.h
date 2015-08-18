#pragma once
#include "VECTOR_3D.h"
#include "VECTOR_4D.h"
#include "IDEAL_GAS_EQUATIONS.h"

class ROE_SOLVER_1D
{
public:
	static void RoeAverage(const VECTOR_3D<double>& Q_l, const VECTOR_3D<double>& Q_r, double& u, double& H, double& c, const double& gamma = 1.4)
	{
		IDEAL_GAS_1D EOS(gamma);
		double rho_l, rho_r, u_l, u_r, H_l, H_r;
		rho_l = Q_l[0]         , rho_r = Q_r[0];
		u_l   = Q_l[1] / Q_l[0], u_r   = Q_r[1] / Q_r[0];
		H_l = EOS.H(Q_l), H_r = EOS.H(Q_r);
		double denominator = sqrt(rho_l) + sqrt(rho_r);

		u  = sqrt(rho_l) * u_l + sqrt(rho_r) * u_r;
		u /= denominator;

		H  = sqrt(rho_l) * H_l + sqrt(rho_r) * H_r;
		H /= denominator;

		c = sqrt((gamma-1) * (H-0.5*u*u));
	}

	static void Inversion(const double& u, const double& H, const double& c,
					      const double& b1, const double& b2, const double& b3,
						  double& x1, double& x2, double& x3)
	{
		x2 = (H-u*u)*b1 + u*b2 - b3;
		x2 *= 0.4 / (c*c);
		x3 = b2 + (c-u)*b1 - c*x2;
		x3 /= 2*c;
		x1 = b1 - x2 - x3;
	}

	static void Inversion(const double& u, const double& H, const double& c, const VECTOR_3D<double>& b, VECTOR_3D<double>& x)
	{
		Inversion(u, H, c, b[0], b[1], b[2], x[0], x[1], x[2]);
	}
};

class ROE_SOLVER_2D
{
public:
	//left state, right state -> Roe-average value
	//bottom state, top state -> Roe-average value
	static void RoeAverage(const VECTOR_4D<double>& Q_l, const VECTOR_4D<double>& Q_r, double& u, double& v, double& H, double& c, const double& gamma = 1.4)
	{
		IDEAL_GAS_2D EOS(gamma);
		//Left & Right state
		double rho_l, rho_r, H_l, H_r, u_l, u_r, v_l, v_r;
		rho_l = Q_l[0], rho_r = Q_r[0];
		H_l = EOS.H(Q_l), H_r = EOS.H(Q_r);
		u_l = Q_l[1] / Q_l[0], u_r = Q_r[1] / Q_r[0];
		v_l = Q_l[2] / Q_l[0], v_r = Q_r[2] / Q_r[0];
		double denominator = sqrt(rho_l) + sqrt(rho_r);

		u  = sqrt(rho_l) * u_l + sqrt(rho_r) * u_r;
		u /= denominator;

		v  = sqrt(rho_l) * v_l + sqrt(rho_r) * v_r;
		v /= denominator;

		H  = sqrt(rho_l) * H_l + sqrt(rho_r) * H_r;
		H /= denominator;

		c = sqrt(0.4 * (H - 0.5 * (u*u+v*v)));
	}

	//x = R_inv * b
	static void Inversion_X(const double& u, const double & v, const double & H, const double & c,
		             const double& b1, const double& b2, const double& b3, const double& b4,
				     double      & x1, double      & x2, double      & x3, double      & x4)
	{
		x3 = (H-u*u-v*v)*b1 + u*b2 + v*b3 - b4;
		x3 *= 0.4 / c / c;
		x2 = b3 - v*b1;
		x4 = b2 + (c-u)*b1 - c* x3;
		x4 /= 2*c;
		x1 = b1 - x3 - x4;
	}

	//x = R_inv * b
	static void Inversion_Y(const double& u, const double & v, const double & H, const double & c,
					 const double& b1, const double& b2, const double& b3, const double& b4,
					 double      & x1, double      & x2, double      & x3, double      & x4)
	{
		x3 = (H-u*u-v*v)*b1 + v*b3 + u*b2 - b4;
		x3 *= 0.4 / c / c;
		x2 = b2 - u*b1;
		x4 = b3 + (c-v)*b1 - c*x3;
		x4 /= 2*c;
		x1 = b1 - x3 - x4;
	}

	//x = R_inv * b
	static void Inversion_X(const double& u, const double & v, const double & H, const double & c,
						    const VECTOR_4D<double>& b, VECTOR_4D<double>& x)
	{
		Inversion_X(u, v, H, c, b[0], b[1], b[2], b[3], x[0], x[1], x[2], x[3]);		
	}

	//x = R_inv * b
	static void Inversion_Y(const double& u, const double & v, const double & H, const double & c,
		const VECTOR_4D<double>& b, VECTOR_4D<double>& x)
	{
		Inversion_Y(u, v, H, c, b[0], b[1], b[2], b[3], x[0], x[1], x[2], x[3]);	
	}

};