#pragma once

#include "VECTOR_8D.h"
#include "IDEAL_GAS_EQUATIONS.h"

//Actually 1.5D.
class MHD_EQUATION_1D
{
public:
	const double gamma;
	const IDEAL_GAS_3D EOS;

public:
	MHD_EQUATION_1D(const double& gamma_):gamma(gamma_), EOS(gamma){}
	~MHD_EQUATION_1D(){}

public:
	double p_star(const VECTOR_8D<double>& q) const
	{
		return p(q) + 0.5 * (SQR(q[5])+SQR(q[6])+SQR(q[7]));
	}
	double E(const VECTOR_8D<double>& q) const
	{
		return q[4];
	}
	double E(const double& rho, const double& u, const double& v, const double& w, const double& p, const double& B1, const double& B2, const double& B3) const
	{
		return EOS.E(rho, u, v, w, p) + 0.5 * (SQR(B1)+SQR(B2)+SQR(B3));
	}
	double p(const VECTOR_8D<double>& q) const
	{
		return (gamma-1) * (E(q) - 0.5 * (SQR(q[1])+SQR(q[2])+SQR(q[3])) / q[0] - 0.5 * (SQR(q[5])+SQR(q[6])+SQR(q[7])));
	}

	void PrimitiveToConservative(const double& rho, const double& u, const double& v, const double& w,
		const double& p, const double& B1, const double& B2, const double& B3, VECTOR_8D<double>& q) const
	{
		q[0] = rho;
		q[1] = rho * u;
		q[2] = rho * v;
		q[3] = rho * w;
		q[4] = E(rho, u, v, w, p, B1, B2, B3);
		q[5] = B1;
		q[6] = B2;
		q[7] = B3;
	}

	void ConservativeToPrimitive(const VECTOR_8D<double>& q, double& rho, double& u, double& v, double& w,
		 double& p, double& B1, double& B2, double& B3) const
	{
		rho = q[0];
		u = q[1] / q[0];
		v = q[2] / q[0];
		w = q[3] / q[0];
		p = EOS.p(VECTOR_5D<double>(q[0],q[1],q[2],q[3],q[4]));
		B1 = q[5];
		B2 = q[6];
		B3 = q[7];
	}

	void ComputeFu(const VECTOR_8D<double>& q, VECTOR_8D<double>& F)
	{
		double rho, u, v, w, p, p_star_, B1, B2, B3, E_;
		ConservativeToPrimitive(q, rho, u, v, w, p, B1, B2, B3);
		p_star_ = p_star(q);
		E_ = E(q);

		F[0] = rho*u ;
		F[1] = rho*u*u + p_star_ ;
		F[2] = rho*u*v - B1*B2 ;
		F[3] = rho*u*w - B1*B3 ;
		F[4] = (E_+p_star_)*u - (u*B1 + v*B2 + w*B3)*B1 ;
		F[5] = 0 ;
		F[6] = u*B2 - B1*v ;
		F[7] = u*B3 - B1*w ;
	}

	double SpectralRadius(const VECTOR_8D<double>& q) const
	{
		const double rho = q[0];
		const double u = q[1] / q[0];
		const double p_ = p(q);
		const double B_sqrnorm = SQR(q[5]) + SQR(q[6]) + SQR(q[7]);
		double c_f;

		c_f = SQR((gamma*p_ + B_sqrnorm) / rho) - 4*gamma*p_*SQR(q[5])/SQR(rho);
		c_f = sqrt(c_f);
		c_f += (gamma*p_ + B_sqrnorm) / rho;
		c_f = sqrt(0.5 * c_f);

		return ABS(u)+c_f;
	}

	double MaxEigenvalue(const VECTOR_8D<double>& q) const
	{
		const double rho = q[0];
		const double u = q[1] / q[0];
		const double p_ = p(q);
		const double B_sqrnorm = SQR(q[5]) + SQR(q[6]) + SQR(q[7]);
		double c_f;

		c_f = SQR((gamma*p_ + B_sqrnorm) / rho) - 4*gamma*p_*SQR(q[5])/SQR(rho);
		c_f = sqrt(c_f);
		c_f += (gamma*p_ + B_sqrnorm) / rho;
		c_f = sqrt(0.5 * c_f);

		return u + c_f;
	}

	double MinEigenvalue(const VECTOR_8D<double>& q) const
	{
		const double rho = q[0];
		const double u = q[1] / q[0];
		const double p_ = p(q);
		const double B_sqrnorm = SQR(q[5]) + SQR(q[6]) + SQR(q[7]);
		double c_f;

		c_f = SQR((gamma*p_ + B_sqrnorm) / rho) - 4*gamma*p_*SQR(q[5])/SQR(rho);
		c_f = sqrt(c_f);
		c_f += (gamma*p_ + B_sqrnorm) / rho;
		c_f = sqrt(0.5 * c_f);

		return u - c_f;
	}

};
