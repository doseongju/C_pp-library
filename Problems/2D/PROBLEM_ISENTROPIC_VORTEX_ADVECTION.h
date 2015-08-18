#pragma once

#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "IDEAL_GAS_EQUATIONS.h"
#include <omp.h>

class ISENTROPIC_VORTEX_ADVECTION
{
public:
	ISENTROPIC_VORTEX_ADVECTION(){}
	~ISENTROPIC_VORTEX_ADVECTION(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = 0 , x_upper = 10;
		const double y_lower = -5, y_upper =  5;
		gamma = 1.4;

		IDEAL_GAS_2D EOS(gamma);

		grid.Initialize(0, 0, resolution, resolution, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid.dx;
		const double dy = grid.dy;
		grid_x.Initialize(0, 0, resolution + 1, resolution, x_lower - 0.5*dx, y_lower, x_upper + 0.5*dx, y_upper);
		grid_y.Initialize(0, 0, resolution, resolution + 1, x_lower, y_lower - 0.5*dy, x_upper, y_upper + 0.5*dy);
		grid_xy.Initialize(0, 0, resolution + 1, resolution + 1, x_lower - 0.5*dx, y_lower - 0.5*dy, x_upper + 0.5*dx, y_upper + 0.5*dy);

		//t_final = 10 * (10 + dx);	
		t_final = (10 + dx);

		Q.Initialize(grid, 3);
#pragma omp parallel for
		ITERATION_GRID_2D(grid, 3)
		{

			const double x = grid(i, j).x;
			const double y = grid(i, j).y;
			const double rr = SQR(x-5) + SQR(y);
			const double beta = 5;	//vortex strength
			double rho, rhou, rhov, p;

			rho = 1 - (gamma - 1) *SQR(beta) / (8 * gamma*PI*PI) * exp(1 - rr);
			rho = pow(rho, 1.0 / (gamma - 1.0));
			

			rhou = 0 - beta / (2 * PI)*exp((1 - rr) * 0.5) * y;
			rhou *= rho;

			rhov = 1 + beta / (2 * PI) * exp((1 - rr) * 0.5) * (x - 5);
			rhov *= rho;

			p = pow(rho, gamma);

			Q(i, j)[0] = rho;
			Q(i, j)[1] = rhou;
			Q(i, j)[2] = rhov;
			Q(i, j)[3] = EOS.E(rho, rhou/rho, rhov/rho, p);

		}
	}
};

//periodic boundary conditions
class BOUNDARY_MANAGER_ISENVORTEX
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
	BOUNDARY_MANAGER_ISENVORTEX(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_) :Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER_ISENVORTEX(){}

public:
	void Initialize(const int i_start_, const int i_start_g_, const int i_end_, const int i_end_g_,
		const int j_start_, const int j_start_g_, const int j_end_, const int j_end_g_,
		const int ghost_width_)
	{
		i_start = i_start_;
		i_start_g = i_start_g_;
		i_end = i_end_;
		i_end_g = i_end_g_;
		j_start = j_start_;
		j_start_g = j_start_g_;
		j_end = j_end_;
		j_end_g = j_end_g_;
		ghost_width = ghost_width_;
	}

	void SetLeftBoundaryCondition()
	{
		for (int j = j_start; j <= j_end; j++)
		{
			int cnt = 0;
			for (int i = i_start - 1; i >= i_start_g; i--)
			{
				Q(i, j) = Q(i_end - cnt, j);
				cnt++;
				//Q(i, j) = Q(i_start, j);
			}
		}
	}

	void SetRightBoundaryCondition()
	{
		for (int j = j_start; j <= j_end; j++)
		{
			int cnt = 0;
			for (int i = i_end + 1; i <= i_end_g; i++)
			{
				Q(i, j) = Q(i_start + cnt, j);
				cnt++;
				//Q(i, j) = Q(i_end, j);
			}
		}
	}

	void SetBottomBoundaryCondition()
	{

		for (int i = i_start; i <= i_end; i++)
		{
			int cnt = 0;
			for (int j = j_start - 1; j >= j_start_g; j--)
			{
				Q(i, j) = Q(i, j_end - cnt);
				cnt++;
				//Q(i, j) = Q(i, j_start);
			}
		}
	}

	void SetTopBoundaryCondition()
	{

		for (int i = i_start; i <= i_end; i++)
		{
			int cnt = 0;
			for (int j = j_end + 1; j <= j_end_g; j++)
			{
				Q(i, j) = Q(i, j_start + cnt);
				cnt++;
				//Q(i, j) = Q(i, j_end);
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