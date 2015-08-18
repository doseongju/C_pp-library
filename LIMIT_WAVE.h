#pragma once
#include "LIMITERS.h"
#include "VECTOR_ND.h"


template<class T>
class LIMIT_WAVE
{
public:
	//NOTE : It also works under the variables, wave_bottom & wave_center & wave_top
	static void LimitWave(const T& wave_left, const T& wave_center, const T& wave_right,
		                  const double& s, const LIMITER_TYPE& limiter_type, T& wave_limited)
	{
		double dot_left(0), dot_right(0);
		double norm;
		double theta = 0;
		double coef = 0;

		norm      = DotProduct(wave_center, wave_center);
		dot_left  = DotProduct(wave_center, wave_left  );
		dot_right = DotProduct(wave_center, wave_right );

		if(norm == 0)
		{
			wave_limited = wave_center;
			return;
		}

		if(s > 0) theta = dot_left  / norm;
		else      theta = dot_right / norm;

		switch (limiter_type)
		{
		case MINMOD_LIMITER:
			coef = LIMITER::Minmod(theta);
			break;
		case SUPERBEE_LIMITER:
			coef = LIMITER::Superbee(theta);
			break;
		case VANLEER_LIMITER:
			coef = LIMITER::VanLeer(theta);
			break;
		case MC_LIMITER:
			coef = LIMITER::MC(theta);
			break;
		case THIRD_LIMITER:
			coef = LIMITER::ThirdOrder(theta, 1.8);
		}

		wave_limited = coef * wave_center;

	}
};