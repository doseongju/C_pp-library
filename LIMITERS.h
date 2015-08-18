#pragma once
#include "MACROS.h"

enum LIMITER_TYPE
{
	NO_LIMITER,
	MINMOD_LIMITER,
	SUPERBEE_LIMITER,
	MC_LIMITER,
	VANLEER_LIMITER,
	THIRD_LIMITER,
	FIFTH_LIMITER, 
	MINMOD3_LIMITER,
	FIFTH_TEST_LIMITER
};

class LIMITER
{
public:
	static double Dummy(const double& theta)
	{
		return 0;
	}
	static double Minmod(const double& theta)
	{
		return MINMOD(1, theta);
	}

	static double Superbee(const double& theta)
	{
		return MAX3(0, MIN(1, 2*theta), MIN(2, theta));
	}

	static double MC(const double& theta)
	{
		return MAX(0, MIN3(0.5 * (1+theta), 2, 2*theta));
	}

	static double VanLeer(const double& theta)
	{
		return (theta + ABS(theta)) / (1 + ABS(theta));
	}

	
	static double ThirdOrder(const double& theta, const double& alpha)
	{
		const double beta = (1 + 2*theta) / 3.0;

		return MAX(0, MIN3(alpha, alpha*theta, beta));
	}

	static double FifthOrder(const double& ll, const double& l, const double& c, const double& r, const double& rr, const double& alpha)
	{
		const double beta = (-2 * (l-ll) + 11 * (c-l) + 24 *(r-c) -3 * (rr-r)) / (30.0 * (c-l));

		return MAX(0, MIN3(alpha, alpha*(r-c)/(c-l), beta));
	}

	static double FifthOrder(const double& ll, const double& l, const double& c, const double& r, const double& rr)
	{
		const double beta = (-2 * (l-ll) + 11 * (c-l) + 24 *(r-c) -3 * (rr-r)) / (30.0 * (c-l));

		return beta;
	}

	static double FifthOrderTest(const double& ll, const double& l, const double& c, const double& r, const double& rr, const double& alpha)
	{
		double ratio = (r - c) / (c - l);
		if (ratio <= 0) return 0;

		const double beta = (-2 * (l - ll) + 11 * (c - l) + 24 * (r - c) - 3 * (rr - r)) / (30.0 * (c - l));
		double A, B;

		A = MIN(1, ratio);
		B = MIN(alpha, alpha * ratio);
		return MAX(A, MIN(B, beta));

	}

	//MINMOD3
	static double MM3(const double& l, const double& c, const double& r, const double& alpha)
	{
		return MINMOD3(alpha * (c - l), alpha * (r - c), 0.5*(r - l));
	}
};