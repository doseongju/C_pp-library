#pragma once

#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "IDEAL_GAS_EQUATIONS.h"
#include <omp.h>

class SHOCK_VORTEX
{
public:
	SHOCK_VORTEX(){}
	~SHOCK_VORTEX(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = -16, x_upper = 8;
		const double y_lower = -12, y_upper = 12;
		gamma   = 1.4;
		t_final = 6;

		IDEAL_GAS_2D EOS(gamma);
		grid.Initialize(0, 0, resolution, resolution, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid.dx;
		const double dy = grid.dy;
		grid_x.Initialize (0, 0, resolution + 1, resolution    , x_lower - 0.5*dx, y_lower         , x_upper + 0.5*dx, y_upper         );
		grid_y.Initialize (0, 0, resolution    , resolution + 1, x_lower         , y_lower - 0.5*dy, x_upper         , y_upper + 0.5*dy);
		grid_xy.Initialize(0, 0, resolution + 1, resolution + 1, x_lower - 0.5*dx, y_lower - 0.5*dy, x_upper + 0.5*dx, y_upper + 0.5*dy);

		Q.Initialize(grid,3);
#pragma omp parallel for
		ITERATION_GRID_2D(grid, 3)
		{

			const double x = grid(i, j).x;
			const double y = grid(i, j).y;
			const double M_v = 0.5;
			const double M_s = 1.2;
			const double x_v = 2, y_v = 0;
			const double r = sqrt(SQR(x - x_v) + SQR(y - y_v));
			double rho, u, v, p;

			if (x <= 0)
			{
				rho = 1.394708041144614;
				u = -0.860403729339060;
				v = 0;
				p = 1.121834000000000;
			}
			else
			{
				rho = 1 - (gamma - 1) / 2 * SQR(M_v) * exp(1 - r*r);
				rho = pow(rho, 1 / (gamma - 1));

				p = 1 - (gamma - 1) / 2 * SQR(M_v) * exp(1 - r*r);
				p = pow(p, gamma / (gamma - 1));
				p = 1 / gamma * p;

				u = M_v * exp(0.5*(1 - r*r)) * (-y) - 1.200011999940001;

				v = M_v * exp(0.5*(1 - r*r)) * (x - x_v);
			}

			
			Q(i,j)[0] = rho;
			Q(i,j)[1] = rho*u;
			Q(i,j)[2] = rho*v;
			Q(i,j)[3] = EOS.E(rho, u, v, p);

		}
	}
};

//Zero-order extrapolation
class BOUNDARY_MANAGER_SHOCK_VORTEX
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
	BOUNDARY_MANAGER_SHOCK_VORTEX(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_):Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER_SHOCK_VORTEX(){}

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
#pragma omp parallel for
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_start-1; i >= i_start_g; i --)
			{
				Q(i,j) = Q(i_start, j);
			}
		}
	}

	void SetRightBoundaryCondition()
	{
#pragma omp parallel for
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_end+1; i <= i_end_g; i ++)
			{
				Q(i,j) = Q(i_end, j);
			}
		}
	}

	void SetBottomBoundaryCondition()
	{
#pragma omp parallel for
		for(int i = i_start; i <= i_end; i ++)
		{
			for (int cnt = 0; cnt < 3; cnt ++)
			{
				Q(i, j_start - cnt - 1)[0] =  Q(i, j_start)[0];
				Q(i, j_start - cnt - 1)[1] =  Q(i, j_start)[1];
				Q(i, j_start - cnt - 1)[2] = -Q(i, j_start)[2];
				Q(i, j_start - cnt - 1)[3] =  Q(i, j_start)[3];
			}
		}
	}

	void SetTopBoundaryCondition()
	{
#pragma omp parallel for
		for(int i = i_start; i <= i_end; i ++)
		{
			for (int cnt = 0; cnt < 3; cnt ++)
			{
				Q(i, j_end + cnt + 1)[0] =  Q(i, j_end)[0];
				Q(i, j_end + cnt + 1)[1] =  Q(i, j_end)[1];
				Q(i, j_end + cnt + 1)[2] = -Q(i, j_end)[2];
				Q(i, j_end + cnt + 1)[3] =  Q(i, j_end)[3];
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
	}

	//For AWCM!!
	void BandingMask(FIELD_UNIFORM_2D<bool>& mask)
	{
		ITERATION_GRID_2D(mask, 0)
		{
			// 2 is a parameter can be tuned.
			if (2 <= i && i <= mask.i_end - 2 && j >= 2 && j <= mask.j_end - 2) continue;
			mask(i, j) = true;
		}
	}
};