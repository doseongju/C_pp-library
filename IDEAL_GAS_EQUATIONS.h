#pragma once
#include "MACROS.h"
#include "VECTOR_3D.h"
#include "VECTOR_4D.h"
#include "VECTOR_5D.h"
#include <iostream>

class IDEAL_GAS_1D
{
public:
	const double gamma;
public:
	IDEAL_GAS_1D(const double& gamma_ = 1.4):gamma(gamma_){}
	~IDEAL_GAS_1D(){}
public:
	//From the primitives
	double E(const double& rho, const double& u, const double& p) const
	{
		return p / (gamma-1) + 0.5 * rho * SQR(u);
	}

	double e(const double& rho, const double& u, const double& p) const
	{
		return E(rho, u, p) / rho - 0.5 * SQR(u);
	}

	double H(const double& rho, const double& u, const double& p) const
	{
		return (E(rho, u, p) + p) / rho;
	}

	double c(const double& rho, const double& p) const
	{
		return sqrt(gamma * p / rho);
	}

	//From the conservatives
	double E(const VECTOR_3D<double>& q) const
	{
		return q[2];
	}
	double e(const VECTOR_3D<double>& q) const
	{
		return q[2] / q[0] - 0.5 * SQR(q[1]/q[0]);
	}
	double p(const VECTOR_3D<double>& q) const
	{
		return (gamma-1) * (q[2] - 0.5 * SQR(q[1]) / q[0]);
	}
	double H(const VECTOR_3D<double>& q) const
	{
		return (q[2] + p(q)) / q[0];
	}
	double c(const VECTOR_3D<double>& q) const
	{
		return sqrt(gamma * p(q) / q[0]);
	}
};

class IDEAL_GAS_2D
{
public:
	const double gamma;
public:
	IDEAL_GAS_2D(const double& gamma_ = 1.4):gamma(gamma_){}
	~IDEAL_GAS_2D(){}
public:
	//From the primitives
	double E(const double& rho, const double& u, const double& v, const double& p) const
	{
		return p / (gamma-1) + 0.5 * rho * (u*u + v*v);
	}

	double e(const double& rho, const double& u, const double& v, const double& p) const
	{
		return E(rho, u, v, p) / rho - 0.5 * (u*u + v*v);
	}

	double H(const double& rho, const double& u, const double& v, const double& p) const
	{
		return (E(rho, u, v,p) + p) / rho;
	}

	double c(const double& rho, const double& p) const
	{
		return sqrt(gamma * p / rho);
	}

	//From the conservatives
	double E(const VECTOR_4D<double>& q) const
	{
		return q[3];
	}
	double e(const VECTOR_4D<double>& q) const
	{
		return q[3] / q[0] - 0.5 * (SQR(q[1]) + SQR(q[2])) / SQR(q[0]);
	}
	double p(const VECTOR_4D<double>& q) const
	{
		return (gamma-1) * (q[3] - 0.5 * (SQR(q[1]) + SQR(q[2])) / q[0]);
	}
	double H(const VECTOR_4D<double>& q) const
	{
		return (q[3] + p(q)) / q[0];
	}
	double c(const VECTOR_4D<double>& q) const
	{
		return sqrt(gamma * p(q) / q[0]);
	}

};

class IDEAL_GAS_3D
{
public:
	const double gamma;

public:
	IDEAL_GAS_3D(const double& gamma_ = 1.4):gamma(gamma_){}
	~IDEAL_GAS_3D(){}

public:
	//From the primitives
	double E(const double& rho, const double& u, const double& v, const double& w, const double& p) const
	{
		return p / (gamma-1) + 0.5 * rho * (u*u + v*v + w*w);
	}
	double e(const double& rho, const double& u, const double& v, const double& w, const double& p) const
	{
		return E(rho, u, v, w, p) / rho - 0.5 * (u*u + v*v + w*w);
	}
	double H(const double& rho, const double& u, const double& v, const double& w, const double& p) const
	{
		return (E(rho, u, v, w, p) + p) / rho;
	}
	double c(const double& rho, const double& p) const 
	{
		return sqrt(gamma * p / rho);
	}

	//From the conservatives
	double E(const VECTOR_5D<double>& q) const
	{
		return q[4];
	}
	double e(const VECTOR_5D<double>& q) const
	{
		return q[4] / q[0] - 0.5 * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]) / SQR(q[0]);
	}
	double p(const VECTOR_5D<double>& q) const
	{
		return q[0] * e(q) * (gamma-1);
	}
	double H(const VECTOR_5D<double>& q) const
	{
		return (q[4] + p(q)) / q[0];
	}
	double c(const VECTOR_5D<double>& q) const
	{
		return sqrt(p(q)*gamma/q[0]);
	}
};