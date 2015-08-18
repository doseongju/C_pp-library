#pragma once
#include "MACROS.h"
#include "ARRAY_1D.h"
#include "ARRAY_2D.h"
#include "ROESOLVER.h"
#include "EULER_EQUATION.h"
#include "MLP.h"
#include "MUSCL.h"

class WENOZ5
{
public:
	double eps;
	ARRAY_2D<double> coef;
public:
	WENOZ5()
	{
		eps = 1e-40;
		coef.Initialize(-1,0,4,3, true);
		coef(-1,0) =  11.0/6.0, coef(-1,1) = -7.0/6.0, coef(-1,2) =  1.0 /3.0;
		coef( 0,0) =  1.0 /3.0, coef( 0,1) =  5.0/6.0, coef( 0,2) = -1.0 /6.0;
		coef( 1,0) = -1.0 /6.0, coef( 1,1) =  5.0/6.0, coef( 1,2) =  1.0 /3.0;
		coef( 2,0) =  1.0 /3.0, coef( 2,1) = -7.0/6.0, coef( 2,2) =  11.0/6.0;
	}
	~WENOZ5(){}
public:

	//TODO : Optimization.
	inline void Reconstruct(const double& ll, const double& l, const double& c, const double& r, const double& rr, double& left, double& right)
	{
		//right
		double d0 = 0.3, d1 = 0.6, d2 = 0.1;
		double b0, b1, b2;
		double tau;

		b0 = 13.0/12.0 * SQR(c - 2*r + rr) + 0.25 * SQR(3*c - 4*r + rr);
		b1 = 13.0/12.0 * SQR(l - 2*c + r ) + 0.25 * SQR(l - r);
		b2 = 13.0/12.0 * SQR(ll - 2*l + c) + 0.25 * SQR(ll - 4*l +3*c);
		tau = ABS(b0 - b2);

		double a0, a1, a2;
		double a_sum(0);
		double w0, w1, w2;

		a0 = d0 * (1 + tau / (b0 + eps));
		a1 = d1 * (1 + tau / (b1 + eps));
		a2 = d2 * (1 + tau / (b2 + eps));
		a_sum = a0 + a1 + a2;

		w0 = a0 / a_sum;
		w1 = a1 / a_sum;
		w2 = a2 / a_sum;

		right  = w0 * ( coef(0,0)*c + coef(0,1)*r + coef(0,2)*rr );
		right += w1 * ( coef(1,0)*l + coef(1,1)*c + coef(1,2)*r  );
		right += w2 * ( coef(2,0)*ll+ coef(2,1)*l + coef(2,2)*c  );

		//left
		d0 = 0.1, d2 = 0.3;

		a0 = d0 / SQR(eps + b0);
		a1 = d1 / SQR(eps + b1);
		a2 = d2 / SQR(eps + b2);
		a_sum = a0 + a1 + a2;

		w0 = a0 / a_sum;
		w1 = a1 / a_sum;
		w2 = a2 / a_sum;

		left  = w0 * ( coef(-1,0)*c  + coef(-1,1)*r + coef(-1,2)*rr );
		left += w1 * ( coef( 0,0)*l  + coef( 0,1)*c + coef( 0,2)*r  );
		left += w2 * ( coef( 1,0)*ll + coef( 1,1)*l + coef( 1,2)*c  );
	}

	inline void ReconstructLeft(const double& ll, const double& l, const double& c, const double& r, const double& rr, double& left)
	{
		double d0 = 0.1, d1 = 0.6, d2 = 0.3;
		double b0, b1, b2;
		double tau;

		b0 = 13.0/12.0 * SQR(c - 2*r + rr) + 0.25 * SQR(3*c - 4*r + rr);
		b1 = 13.0/12.0 * SQR(l - 2*c + r ) + 0.25 * SQR(l - r);
		b2 = 13.0/12.0 * SQR(ll - 2*l + c) + 0.25 * SQR(ll - 4*l +3*c);
		tau = ABS(b0 - b2);

		double a0, a1, a2;
		double a_sum(0);
		double w0, w1, w2;

		d0 = 0.1, d2 = 0.3;

		a0 = d0 * (1 + tau / (b0 + eps));
		a1 = d1 * (1 + tau / (b1 + eps));
		a2 = d2 * (1 + tau / (b2 + eps));
		a_sum = a0 + a1 + a2;

		w0 = a0 / a_sum;
		w1 = a1 / a_sum;
		w2 = a2 / a_sum;

		left  = w0 * ( coef(-1,0)*c  + coef(-1,1)*r + coef(-1,2)*rr );
		left += w1 * ( coef( 0,0)*l  + coef( 0,1)*c + coef( 0,2)*r  );
		left += w2 * ( coef( 1,0)*ll + coef( 1,1)*l + coef( 1,2)*c  );
	}

	inline void ReconstructRight(const double& ll, const double& l, const double& c, const double& r, const double& rr, double& right)
	{
		ReconstructLeft(rr, r, c, l, ll, right);
		/*//right
		double d0 = 0.3, d1 = 0.6, d2 = 0.1;
		double b0, b1, b2;
		double tau;

		b0 = 13.0/12.0 * SQR(c - 2*r + rr) + 0.25 * SQR(3*c - 4*r + rr);
		b1 = 13.0/12.0 * SQR(l - 2*c + r ) + 0.25 * SQR(l - r);
		b2 = 13.0/12.0 * SQR(ll - 2*l + c) + 0.25 * SQR(ll - 4*l +3*c);
		tau = ABS(b0 - b2);

		double a0, a1, a2;
		double a_sum(0);
		double w0, w1, w2;

		a0 = d0 * (1 + tau / (b0 + eps));
		a1 = d1 * (1 + tau / (b1 + eps));
		a2 = d2 * (1 + tau / (b2 + eps));
		a_sum = a0 + a1 + a2;

		w0 = a0 / a_sum;
		w1 = a1 / a_sum;
		w2 = a2 / a_sum;

		right  = w0 * ( coef(0,0)*c + coef(0,1)*r + coef(0,2)*rr );
		right += w1 * ( coef(1,0)*l + coef(1,1)*c + coef(1,2)*r  );
		right += w2 * ( coef(2,0)*ll+ coef(2,1)*l + coef(2,2)*c  );*/
	}
};

class WENOZ_EULER_SOLVER_1D
{
public:
	WENOZ5 wenoz;
	const double gamma;
	//local array for characteristic decompositions.
	ARRAY_1D<VECTOR_3D<double>> R_Q, R_fu;
	//local array for vector-splitting.
	ARRAY_1D<VECTOR_3D<double>> fu, fu_plus, fu_minus;

public:
	WENOZ_EULER_SOLVER_1D(const double& gamma_):gamma(gamma_)
	{
		R_Q.Initialize(-3, 6);
		R_fu.Initialize(-3, 6);
		fu.Initialize(-3, 6);
		fu_plus.Initialize(-3, 6);
		fu_minus.Initialize(-3, 6);
	}
	~WENOZ_EULER_SOLVER_1D(){}

public:
	void ComputeFlux(const VECTOR_3D<double>& lll, const VECTOR_3D<double>& ll, const VECTOR_3D<double>& l,
		const VECTOR_3D<double>& r, const VECTOR_3D<double>& rr, const VECTOR_3D<double>& rrr,
		const double& alpha, VECTOR_3D<double>& flux)
	{
		double u, H, c;
		EULER_EQUATION_1D euler_eqn(gamma);
		euler_eqn.ComputeFu(lll, fu(-3));
		euler_eqn.ComputeFu(ll , fu(-2));
		euler_eqn.ComputeFu(l  , fu(-1));
		euler_eqn.ComputeFu(r  , fu( 0));
		euler_eqn.ComputeFu(rr , fu( 1));
		euler_eqn.ComputeFu(rrr, fu( 2));

		ROE_SOLVER_1D::RoeAverage(l, r, u, H, c, gamma);

		VECTOR_3D<double> r1, r2, r3;
		r1[0] = 1,     r2[0] = 1            , r3[0] = 1    ;
		r1[1] = u-c,   r2[1] = u            , r3[1] = u+c  ;
		r1[2] = H-u*c, r2[2] = 0.5*(u*u)    , r3[2] = H+u*c;


		ROE_SOLVER_1D::Inversion(u, H, c, lll, R_Q(-3));
		ROE_SOLVER_1D::Inversion(u, H, c, ll , R_Q(-2));
		ROE_SOLVER_1D::Inversion(u, H, c, l  , R_Q(-1));
		ROE_SOLVER_1D::Inversion(u, H, c, r  , R_Q( 0));
		ROE_SOLVER_1D::Inversion(u, H, c, rr , R_Q( 1));
		ROE_SOLVER_1D::Inversion(u, H, c, rrr, R_Q( 2));

		for(int p = -3; p < 3; p ++) ROE_SOLVER_1D::Inversion(u, H, c, fu(p), R_fu(p));
		for (int p = -3; p < 3; p ++)
		{
			fu_plus (p) = 0.5 * (R_fu(p) + alpha*R_Q(p));
			fu_minus(p) = 0.5 * (R_fu(p) - alpha*R_Q(p));
		}

		VECTOR_4D<double> left_l, left_r, right_l, right_r;
		for(int p = 0; p < 3; p ++)
		{
			wenoz.ReconstructRight(fu_plus(-3)[p], fu_plus(-2)[p], fu_plus(-1)[p], fu_plus(0)[p], fu_plus(1)[p], left_r[p]);
			wenoz.ReconstructLeft(fu_minus(-2)[p], fu_minus(-1)[p], fu_minus(0)[p], fu_minus(1)[p], fu_minus(2)[p], right_l[p]);

			//MLP_2D::ReconstructRight(fu_plus(-3)[p], fu_plus(-2)[p], fu_plus(-1)[p], fu_plus(0)[p], fu_plus(1)[p], 1.7, left_r[p]);
			//MLP_2D::ReconstructLeft(fu_minus(-2)[p], fu_minus(-1)[p], fu_minus(0)[p], fu_minus(1)[p], fu_minus(2)[p], 1.7, right_l[p]);
		}
		flux  = (left_r[0]+right_l[0]) * r1 + (left_r[1]+right_l[1]) * r2 + (left_r[2]+right_l[2]) * r3;
	}
};

class WENOZ_EULER_SOLVER_2D
{
public:
	WENOZ5 wenoz;
	ARRAY_1D<VECTOR_4D<double>> fu, gu;
	ARRAY_1D<VECTOR_4D<double>> R_Q, R_fu, R_gu;
	ARRAY_1D<VECTOR_4D<double>> R_fu_plus, R_fu_minus, R_gu_plus, R_gu_minus;
	const double gamma;

	MLP_2D mlp;
public:
	WENOZ_EULER_SOLVER_2D(const double& gamma_ = 1.4):gamma(gamma_)
	{
		fu.Initialize(-3,6);
		gu.Initialize(-3,6);

		R_Q.Initialize (-3,6);
		R_fu.Initialize(-3,6);
		R_gu.Initialize(-3,6);

		R_fu_plus.Initialize (-3,6);
		R_fu_minus.Initialize(-3,6);
		R_gu_plus.Initialize (-3,6);
		R_gu_minus.Initialize(-3,6);
	}
	~WENOZ_EULER_SOLVER_2D(){}

public:
	void ComputeFluxX(const VECTOR_4D<double>& lll, const VECTOR_4D<double>& ll, const VECTOR_4D<double>& l,
		const VECTOR_4D<double>& r, const VECTOR_4D<double>& rr, const VECTOR_4D<double>& rrr,
		const double& alpha, VECTOR_4D<double>& flux)
	{
		double u, v, c, H;
		EULER_EQUATION_2D euler_eqn(gamma);
		euler_eqn.ComputeFu(lll, fu(-3));
		euler_eqn.ComputeFu(ll , fu(-2));
		euler_eqn.ComputeFu(l  , fu(-1));
		euler_eqn.ComputeFu(r  , fu( 0));
		euler_eqn.ComputeFu(rr , fu( 1));
		euler_eqn.ComputeFu(rrr, fu( 2));

		ROE_SOLVER_2D::RoeAverage(l, r, u, v, H, c, gamma);
		VECTOR_4D<double> r1, r2, r3, r4;
		r1[0] = 1,     r2[0] = 0, r3[0] = 1            , r4[0] = 1    ;
		r1[1] = u-c,   r2[1] = 0, r3[1] = u            , r4[1] = u+c  ;
		r1[2] = v,     r2[2] = 1, r3[2] = v            , r4[2] = v    ;
		r1[3] = H-u*c, r2[3] = v, r3[3] = 0.5*(u*u+v*v), r4[3] = H+u*c;

		ROE_SOLVER_2D::Inversion_X(u, v, H, c, lll, R_Q(-3));
		ROE_SOLVER_2D::Inversion_X(u, v, H, c, ll , R_Q(-2));
		ROE_SOLVER_2D::Inversion_X(u, v, H, c, l  , R_Q(-1));
		ROE_SOLVER_2D::Inversion_X(u, v, H, c, r  , R_Q( 0));
		ROE_SOLVER_2D::Inversion_X(u, v, H, c, rr , R_Q( 1));
		ROE_SOLVER_2D::Inversion_X(u, v, H, c, rrr, R_Q( 2));

		for(int p = -3; p < 3; p ++) ROE_SOLVER_2D::Inversion_X(u, v, H, c, fu(p), R_fu(p));
		for (int p = -3; p < 3; p ++)
		{
			R_fu_plus (p) = 0.5 * (R_fu(p) + alpha*R_Q(p));
			R_fu_minus(p) = 0.5 * (R_fu(p) - alpha*R_Q(p));
		}

		VECTOR_4D<double> left_r, right_l;
		for(int p = 0; p < 4; p ++)
		{
			wenoz.ReconstructRight(R_fu_plus(-3)[p], R_fu_plus(-2)[p], R_fu_plus(-1)[p], R_fu_plus(0)[p], R_fu_plus(1)[p], left_r[p]);
			wenoz.ReconstructLeft(R_fu_minus(-2)[p], R_fu_minus(-1)[p], R_fu_minus(0)[p], R_fu_minus(1)[p], R_fu_minus(2)[p], right_l[p]);
		}
		flux  = (left_r[0]+right_l[0]) * r1 + (left_r[1]+right_l[1]) * r2 + (left_r[2]+right_l[2]) * r3 + (left_r[3]+right_l[3]) * r4;
	}

	void ComputeFluxY(const VECTOR_4D<double>& bbb, const VECTOR_4D<double>& bb, const VECTOR_4D<double>& b,
		const VECTOR_4D<double>& t, const VECTOR_4D<double>& tt, const VECTOR_4D<double>& ttt,
		const double& alpha, VECTOR_4D<double>& flux)
	{
		double u, v, c, H;
		EULER_EQUATION_2D euler_eqn(gamma);
		euler_eqn.ComputeGu(bbb, gu(-3));
		euler_eqn.ComputeGu(bb , gu(-2));
		euler_eqn.ComputeGu(b  , gu(-1));
		euler_eqn.ComputeGu(t  , gu( 0));
		euler_eqn.ComputeGu(tt , gu( 1));
		euler_eqn.ComputeGu(ttt, gu( 2));

		ROE_SOLVER_2D::RoeAverage(b, t, u, v, H, c, gamma);
		VECTOR_4D<double> r1, r2, r3, r4;
		r1[0] = 1  ,   r2[0] = 0, r3[0] = 1            , r4[0] = 1    ;
		r1[1] = u  ,   r2[1] = 1, r3[1] = u            , r4[1] = u    ;
		r1[2] = v-c,   r2[2] = 0, r3[2] = v            , r4[2] = v+c  ;
		r1[3] = H-v*c, r2[3] = u, r3[3] = 0.5*(u*u+v*v), r4[3] = H+v*c;

		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, bbb, R_Q(-3));
		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, bb , R_Q(-2));
		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, b  , R_Q(-1));
		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, t  , R_Q( 0));
		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, tt , R_Q( 1));
		ROE_SOLVER_2D::Inversion_Y(u, v, H, c, ttt, R_Q( 2));

		for(int p = -3; p < 3; p ++) ROE_SOLVER_2D::Inversion_Y(u, v, H, c, gu(p), R_gu(p));
		for (int p = -3; p < 3; p ++)
		{
			R_gu_plus (p) = 0.5 * (R_gu(p) + alpha*R_Q(p));
			R_gu_minus(p) = 0.5 * (R_gu(p) - alpha*R_Q(p));
		}

		VECTOR_4D<double> bottom_t, top_b;
		for(int p = 0; p < 4; p ++)
		{
			wenoz.ReconstructRight(R_gu_plus (-3)[p], R_gu_plus (-2)[p], R_gu_plus(-1)[p], R_gu_plus (0)[p], R_gu_plus (1)[p], bottom_t[p]);
			wenoz.ReconstructLeft(R_gu_minus(-2)[p], R_gu_minus(-1)[p], R_gu_minus(0)[p], R_gu_minus(1)[p], R_gu_minus(2)[p], top_b[p]);
		}
		flux = (bottom_t[0] + top_b[0]) * r1 + (bottom_t[1] + top_b[1]) * r2 + (bottom_t[2] + top_b[2]) * r3 + (bottom_t[3] + top_b[3]) * r4;
	}
};