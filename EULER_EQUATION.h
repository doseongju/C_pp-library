#pragma once
#include "FIELD_UNIFORM_1D.h"
#include "VECTOR_3D.h"
#include "IDEAL_GAS_EQUATIONS.h"

class EULER_EQUATION_1D
{
public:
	const double gamma;
	const IDEAL_GAS_1D EOS;
public:
	EULER_EQUATION_1D(const double& gamma_ = 1.4):gamma(gamma_), EOS(gamma_){}
	~EULER_EQUATION_1D(){}

public:
	void ComputeFu(const VECTOR_3D<double>& Q, VECTOR_3D<double>& fu) const
	{
		double rho = Q[0];
		double u = Q[1] / Q[0];
		double E = Q[2];
		double p = EOS.p(Q);

		fu[0] = rho * u;
		fu[1] = rho * u * u + p;
		fu[2] = (E + p) * u;
	}

	void ConservedToPrimitive(const VECTOR_3D<double>& Q, double& rho, double& u, double& p) const
	{
		rho = Q[0];
		u = Q[1] / Q[0];
		p = EOS.p(Q);
	}

	void PrimitiveToConserved(const double& rho, const double& u, const double& p, VECTOR_3D<double>& Q) const
	{
		Q[0] = rho;
		Q[1] = rho * u;
		Q[2] = p / (gamma-1) + 0.5 * rho * (u*u);
	}

	double Spectral_Radius(const VECTOR_3D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return ABS(u) + c;
	}

	double MaxEigenvalue(const VECTOR_3D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return u + c;
	}

	double MinEigenvalue(const VECTOR_3D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return u - c;
	}
};

class EULER_EQUATION_2D
{
public:
	const double gamma;
	const IDEAL_GAS_2D EOS;
public:
	EULER_EQUATION_2D(const double& gamma_ = 1.4):gamma(gamma_), EOS(gamma_){}
	~EULER_EQUATION_2D(){}

public:
	void ComputeFu(const VECTOR_4D<double>& Q, VECTOR_4D<double>& fu) const
	{
		double rho = Q[0];
		double u = Q[1] / Q[0];
		double v = Q[2] / Q[0];
		double E = Q[3];
		double p = EOS.p(Q);

		fu[0] = rho     * u;
		fu[1] = rho * u * u + p;
		fu[2] = rho * v * u;
		fu[3] = (E + p) * u;
	}
	void ComputeGu(const VECTOR_4D<double>& Q, VECTOR_4D<double>& gu) const
	{
		double rho = Q[0];
		double u = Q[1] / Q[0];
		double v = Q[2] / Q[0];
		double E = Q[3];
		double p = EOS.p(Q);

		gu[0] = rho     * v;
		gu[1] = rho * u * v;
		gu[2] = rho * v * v + p;
		gu[3] = (E + p) * v;
	}


	void ConservedToPrimitive(const VECTOR_4D<double>& Q, double& rho, double& u, double& v, double& p) const
	{
		rho = Q[0];
		u = Q[1] / Q[0];
		v = Q[2] / Q[0];
		p = EOS.p(Q);
	}
	void PrimitiveToConserved(const double& rho, const double& u, const double& v, const double& p, VECTOR_4D<double>& Q) const
	{
		Q[0] = rho;
		Q[1] = rho * u;
		Q[2] = rho * v;
		Q[3] = p / (gamma-1) + 0.5 * rho * (u*u + v*v);
	}


	double Spectral_RadiusX(const VECTOR_4D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return ABS(u) + c;
	}
	double MaxEigenvalueX(const VECTOR_4D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return u + c;
	}
	double MinEigenvalueX(const VECTOR_4D<double>& q) const
	{
		double u = q[1] / q[0];
		double c = EOS.c(q);
		return u - c;
	}


	double Spectral_RadiusY(const VECTOR_4D<double>& q) const
	{
		double v = q[2] / q[0];
		double c = EOS.c(q);
		return ABS(v) + c;
	}
	double MaxEigenvalueY(const VECTOR_4D<double>& q) const
	{
		double v = q[2] / q[0];
		double c = EOS.c(q);
		return v + c;
	}
	double MinEigenvalueY(const VECTOR_4D<double>& q) const
	{
		double v = q[2] / q[0];
		double c = EOS.c(q);
		return v - c;
	}

};