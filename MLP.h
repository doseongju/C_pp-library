#pragma once
#include "MACROS.h"
#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "LIMITERS.h"

class MLP_2D
{
public:
	const double eps;
public:
	MLP_2D(const double& eps_ = 1e-16):eps(eps_){}
	~MLP_2D(){}

public:
	//MLP-2005
	void ComputeAlpha(const double& ll, const double& l,
					  const double& r, const double& rr,
					  const double& tl, const double& tr, const double& bl, const double& br,
					  double& alpha_l, double& alpha_r)

	{
		double tantheta_l, tantheta_r;
		double rL, rR;

		tantheta_l = ABS((tl - bl) / (r - ll + eps));
		tantheta_r = ABS((tr - br) / (rr - l + eps));
		rL = (r - l) / (l - ll + eps);
		rR = (r - l) / (rr - r + eps);
		alpha_l  = 1.0 + MAX(0.0, tantheta_r / (rR + eps));
		alpha_l *= 2.0 * MAX(1.0, rL);
		alpha_l /= 1.0 + tantheta_l;
		alpha_l = MAX(1.0, MIN(2.0, alpha_l));

		alpha_r  = 1.0 + MAX(0.0, tantheta_l / (rL + eps));
		alpha_r *= 2.0 * MAX(1.0, rR);
		alpha_r /= 1.0 + tantheta_r;
		alpha_r = MAX(1.0, MIN(2.0, alpha_r));
		
	}

	//MLP-2005 modified
	void ComputeAlphaTest(const double& ll, const double& l,
					  const double& r, const double& rr,
					  const double& tl, const double& tr, const double& bl, const double& br,
					  double& alpha_l, double& alpha_r)
	{
		double tantheta_l, tantheta_r;
		double rL, rR;

		tantheta_l = ABS((tl - bl) / (r - ll + eps));
		tantheta_r = ABS((tr - br) / (rr - l + eps));
		rL = (r - l) / (l - ll + eps);
		rR = (r - l) / (rr - r + eps);

		alpha_l  = 1.0 + tantheta_r;
		alpha_l *= 2.0 * MAX(1.0, rL);
		alpha_l /= 1.0 + tantheta_l;
		alpha_l = MAX(1.0, MIN(2.0, alpha_l));

		alpha_r  = 1.0 + tantheta_l;
		alpha_r *= 2.0 * MAX(1.0, rR);
		alpha_r /= 1.0 + tantheta_r;
		alpha_r = MAX(1.0, MIN(2.0, alpha_r));
	}

	//Proposed by DSJ
	void ComputeAlphaYoonTest(const double& ll, const double& l,
		const double& r, const double& tl, const double& tr, const double& bl, const double& br,
		double& alpha_l)
	{
		double r_x, r_xy, min_max;
		if ((tl - l) * (l - bl) <= 0)
		{
			alpha_l = 2.0;
			return;
		}
		if (r - l >= 0)
		{
			if (tl - l >= 0) 
			{
				min_max = MAX4(l, r, tr, tl);
				r_xy = tl - l;
			}
			else 
			{
				min_max = MAX4(l, r, br, bl);
				r_xy = l - bl;
			}
		}
		else
		{
			if (tl - l >= 0) 
			{
				min_max = MIN4(l, r, br, bl);
				r_xy = l - bl;
			}
			else 
			{
				min_max = MIN4(l, r, tr, tl);
				r_xy = tl - l;
			}
		}

		r_x = (r - l) / (l - ll + eps);
		r_xy /= r - l + eps;

		alpha_l = (2 * MAX(1, r_x) / (1 + ABS(r_xy)));
		alpha_l *= ABS((min_max - l) / (r - l + eps));
		alpha_l = MAX(1, MIN(alpha_l, 2));
	}

	//MLP-2008
	void ComputeAlphaYoon(const double& ll, const double& l, const double& c, const double& r, const double& rr,
		                  const double& bb, const double& b,                  const double& t, const double& tt,
						  const double& tl, const double& tc, const double& tr, const double& bl, const double& bc, const double& br,
						  double& alpha_x, double& alpha_y)
	{
		int k1, k2, k1_, k2_;
		double r_x, r_xy, r_y, r_yx;
		double min, max;
		double deviation_x, deviation_y;
		ARRAY_2D<double> q(-1,-1,3,3);
		q(-1, 1) = tl; q(0, 1) = t; q(1, 1) = tr;
		q(-1, 0) = l ; q(0, 0) = c; q(1, 0) = r ;
		q(-1,-1) = bl; q(0,-1) = b; q(1,-1) = br;
		
		r_x = (r-c) / (c-l + eps);
		r_y = (t-c) / (c-b + eps);

		//compute deviation along x, y
		if(ABS(c-l) > eps) deviation_x = 0.5 * LIMITER::FifthOrder(ll, l, c, r, rr) * (c-l);
		else deviation_x = 0;
		if(ABS(c-b) > eps) deviation_y = 0.5 * LIMITER::FifthOrder(bb, b, c, t, tt) * (c-b);
		else deviation_y = 0;

		//find max/min direction
		if(deviation_x >= 0) k1 = 1;
		else k1 = -1;
		if(deviation_y >= 0) k2 = 1;
		else k2 = -1;
		k1_ = -1 * k1;
		k2_ = -1 * k2;

/*
		//for symmetricity
		double deviation_x_, deviation_y_;
		if(ABS(c-r) > eps) deviation_x_ = 0.5 * LIMITER::FifthOrder(rr, r, c, l, ll) * (c-r);
		else deviation_x_ = 0;
		if(ABS(c-t) > eps) deviation_y_ = 0.5 * LIMITER::FifthOrder(tt, t, c, b, bb) * (c-t);
		else deviation_y_ = 0;

		if(ABS(deviation_x_) > ABS(deviation_x))
		{
			deviation_x = SIGN(deviation_x) * ABS(deviation_x_);
			r_x = (c-l) / (r-c + eps);
		}
		if(ABS(deviation_y_) > ABS(deviation_y))
		{
			deviation_y = SIGN(deviation_y) * ABS(deviation_y_);
			r_y = (c-b) / (t-c + eps);
		}*/

		max = MAX4(q(k1,k2), q(0,0), q(k1,0), q(0,k2));
		min = MIN4(q(k1_,k2_), q(0,0), q(k1_,0), q(0,k2_));

		r_xy = ABS(deviation_y / (deviation_x + eps));
		r_yx = ABS(deviation_x / (deviation_y + eps));

		if(ABS(r-c) > eps)alpha_x = ABS(2 * MAX(1, r_x) / ((1+r_xy)*(r-c))) * (MIN(ABS(max-c), ABS(min-c)));
		else alpha_x = ABS(2 * MAX(1, r_x) / ((1+r_xy)*eps)) * (MIN(ABS(max-c), ABS(min-c)));
		//alpha_x = 0.5 * alpha_x + 1;
		alpha_x = MIN(alpha_x, 2);

		if(ABS(t-c) > eps) alpha_y = ABS(2*MAX(1, r_y) / ((1+r_yx) * (t-c))) * (MIN(ABS(max-c), ABS(min-c)));
		else alpha_y = ABS(2*MAX(1, r_y) / ((1+r_yx) * eps)) * (MIN(ABS(max-c), ABS(min-c)));
		//alpha_y = 0.5 * alpha_y + 1;
		alpha_y = MIN(alpha_y, 2);
	}

	//5th order interpolation(from ENO coefficient)
	static void ReconstructRight(const double& ll, const double& l, const double& c, const double& r, const double& rr, const double& alpha, double& right)
	{
		if(ABS(c-l) < 1e-16)
		{
			right = c;
		}
		else
		{
			const double beta = (-2 * (l-ll) + 11 * (c-l) + 24 *(r-c) -3 * (rr-r)) / (30.0 * (c-l));
			double phi = MAX(0, MIN3(alpha, alpha*(r-c)/(c-l), beta));

			right = c + 0.5 * phi * (c-l);
		}
	}

	static void ReconstructLeft(const double& ll, const double& l, const double& c, const double& r, const double& rr, const double& alpha, double& left)
	{
		ReconstructRight(rr, r, c, l, ll, alpha, left);
	}
};