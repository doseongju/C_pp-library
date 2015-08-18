#pragma once
#include "LIMITERS.h"

template<class T>
class MUSCL
{
public:
	inline static void Interapolation(const T& lll, const T& ll, const T& l, T& left, T& right, const T& r, const T& rr, const T& rrr, LIMITER_TYPE limiter,
		const double& alpha_l = 1.7, const double& alpha_r = 1.7)
	{
		if(ABS(l-ll) < 1e-16)
		{
			left = l;
		}
		else
		{
			switch(limiter)
			{
			case NO_LIMITER:
				left = l;
				break;
			case MINMOD_LIMITER:
				left  = l + 0.5 * LIMITER::Minmod((r-l)/(l-ll)) * (l-ll);
				break;
			case SUPERBEE_LIMITER:
				left  = l + 0.5 * LIMITER::Superbee((r-l)/(l-ll)) * (l-ll);
				break;
			case VANLEER_LIMITER:
				left  = l + 0.5 * LIMITER::VanLeer((r-l)/(l-ll)) * (l-ll);
				break;
			case MC_LIMITER:
				left  = l + 0.5 * LIMITER::MC((r-l)/(l-ll)) * (l-ll);
				break;
			case THIRD_LIMITER:
				left  = l + 0.5 * LIMITER::ThirdOrder((r-l)/(l-ll), alpha_l) * (l-ll);
				break;
			case FIFTH_LIMITER:
				left  = l + 0.5 * LIMITER::FifthOrder(lll, ll, l, r, rr, alpha_l) * (l-ll);
				break;
			case MINMOD3_LIMITER:
				left = l + 0.5 * LIMITER::MM3(ll, l, r, alpha_l);
				break;
			case FIFTH_TEST_LIMITER:
				left = l + 0.5 * LIMITER::FifthOrderTest(lll, ll, l, r, rr, alpha_l) * (l - ll);
				break;
			}
		}
		
		if(ABS(rr-r) < 1e-16)
		{
			right = r;	
		}
		else
		{
			switch(limiter)
			{
			case NO_LIMITER:
				right = r;
				break;
			case MINMOD_LIMITER:
				right  = r - 0.5 * LIMITER::Minmod((r-l)/(rr-r)) * (rr-r);
				break;
			case SUPERBEE_LIMITER:
				right  = r - 0.5 * LIMITER::Superbee((r-l)/(rr-r)) * (rr-r);
				break;
			case VANLEER_LIMITER:
				right  = r - 0.5 * LIMITER::VanLeer((r-l)/(rr-r)) * (rr-r);
				break;
			case MC_LIMITER:
				right  = r - 0.5 * LIMITER::MC((r-l)/(rr-r)) * (rr-r);
				break;
			case THIRD_LIMITER:
				right  = r - 0.5 * LIMITER::ThirdOrder((r-l)/(rr-r), alpha_r) * (rr-r);
				break;
			case FIFTH_LIMITER:
				right  = r - 0.5 * LIMITER::FifthOrder(rrr, rr, r, l, ll, alpha_r) * (rr-r);
				break;
			case MINMOD3_LIMITER:
				right = r - 0.5 * LIMITER::MM3(l, r, rr, alpha_r);
				break;
			case FIFTH_TEST_LIMITER:
				right = r - 0.5 * LIMITER::FifthOrderTest(rrr, rr, r, l, ll, alpha_r) * (rr - r);
				break;

			}
		}
	}
};