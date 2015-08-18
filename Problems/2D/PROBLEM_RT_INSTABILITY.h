#pragma once

#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "EULER_EQUATION.h"

class RT_INSTABILITY
{
public:
	RT_INSTABILITY(){}
	~RT_INSTABILITY(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = 0, x_upper = 0.25;
		const double y_lower = 0, y_upper = 1;
		gamma   = 5.0/3.0;
		t_final = 1.95;

		IDEAL_GAS_2D EOS(gamma);

		grid_xy.Initialize(0,0, resolution+1, 4*resolution+1, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid_xy.dx;
		const double dy = grid_xy.dy;
		grid_x.Initialize(0,0, resolution+1, 4*resolution  , x_lower, y_lower+0.5*dx, x_upper, y_upper-0.5*dx);
		grid_y.Initialize(0,0, resolution  , 4*resolution+1, x_lower+0.5*dy, y_lower, x_upper-0.5*dy, y_upper);
		grid.Initialize(0,0, resolution, 4*resolution, x_lower+0.5*dx, y_lower+0.5*dy, x_upper-0.5*dx, y_upper-0.5*dy);

		Q.Initialize(grid,3);
		ITERATION_GRID_2D(grid, 3)
		{

			const double x = grid(i,j).x;
			const double y = grid(i,j).y;
			double rho, u, v, p;

			if(y <= 0.5)
			{
				rho = 2;
				p = 2*y + 1;
			}
			else
			{
				rho = 1;
				p = y + 1.5;
			}
			u = 0;
			v = -0.025 * EOS.c(rho, p) * cos(8*PI*x);

			Q(i,j)[0] = rho;
			Q(i,j)[1] = rho*u;
			Q(i,j)[2] = rho*v;
			Q(i,j)[3] = EOS.E(rho, u, v, p);
		}
	}
};

//Zero-order extrapolation
class BOUNDARY_MANAGER_RTI
{
public:
	int i_start;
	int i_start_g;
	int i_end;
	int i_end_g;
	int j_start;
	int j_start_g;
	int j_end;
	int j_end_g;
	int ghost_width;

	FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q;

public:
	BOUNDARY_MANAGER_RTI(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_):Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER_RTI(){}

public:
	void Initialize(const int i_start_, const int i_start_g_, const int i_end_, const int i_end_g_,
		const int j_start_, const int j_start_g_, const int j_end_, const int j_end_g_,
		const int ghost_width_)
	{
		i_start   = i_start_;
		i_start_g = i_start_g_;
		i_end     = i_end_;
		i_end_g   = i_end_g_;
		j_start   = j_start_;
		j_start_g = j_start_g_;
		j_end     = j_end_;
		j_end_g   = j_end_g_;
		ghost_width = ghost_width_;
	}

	void SetLeftBoundaryCondition()
	{
		for(int j = j_start; j <= j_end; j ++)
		{
			int cnt = 0;
			for(int i = i_start-1; i >= i_start_g; i --)
			{
				Q(i,j)[0] =  Q(i_start/*+cnt*/, j)[0];
				Q(i,j)[1] = -Q(i_start/*+cnt*/, j)[1];
				Q(i,j)[2] =  Q(i_start/*+cnt*/, j)[2];
				Q(i,j)[3] =  Q(i_start/*+cnt*/, j)[3];
				cnt ++;
			}
		}
	}

	void SetRightBoundaryCondition()
	{
		for(int j = j_start; j <= j_end; j ++)
		{
			int cnt = 0;
			for(int i = i_end+1; i <= i_end_g; i ++)
			{
				Q(i,j)[0] =  Q(i_end/*-cnt*/, j)[0];
				Q(i,j)[1] = -Q(i_end/*-cnt*/, j)[1];
				Q(i,j)[2] =  Q(i_end/*-cnt*/, j)[2];
				Q(i,j)[3] =  Q(i_end/*-cnt*/, j)[3];
				cnt ++;
			}
		}
	}

	void SetBottomBoundaryCondition()
	{
		EULER_EQUATION_2D euler_eqn(5.0/3.0);
		double rho, u, v, p;
		rho = 1; u = 0; v = 0; p = 1;

		for(int i = i_start; i <= i_end; i ++)
		{
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				euler_eqn.PrimitiveToConserved(rho, u, v, p, Q(i,j));
			}
		}
	}

	void SetTopBoundaryCondition()
	{
		EULER_EQUATION_2D euler_eqn(5.0/3.0);
		double rho, u, v, p;
		rho = 2; u = 0; v = 0; p = 2.5;

		for(int i = i_start; i <= i_end; i ++)
		{
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				euler_eqn.PrimitiveToConserved(rho, u, v, p, Q(i,j));
			}
		}
	}

	void SetDiagonalBoundaryCondition()
	{
		for(int i = i_start_g; i < i_start; i ++)
		{
			for(int j = j_start_g; j < j_start; j ++)
			{
				Q(i,j) = Q(i_start,j_start);
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_start,j_end);
			}
		}
		for(int i = i_end+1; i <= i_end_g; i ++)
		{
			for(int j = j_start_g; j < j_start; j ++)
			{
				Q(i,j) = Q(i_end,j_start);
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_end,j_end);
			}
		}
	}

	void SetBoundaryCondition()
	{
		SetLeftBoundaryCondition();
		SetRightBoundaryCondition();
		SetBottomBoundaryCondition();
		SetTopBoundaryCondition();
		//SetDiagonalBoundaryCondition();
	}


	//For AWCM!!
	void BaningMask(FIELD_UNIFORM_2D<bool>& mask)
	{
		ITERATION_GRID_2D(mask, 0)
		{
			// 2 is a parameter can be tuned.
			if (2 <= i && i <= mask.i_end - 2 && j >= 2 && j <= mask.j_end - 2) continue;
			mask(i, j) = true;
		}
	}
};