#pragma once

#include "IDEAL_GAS_EQUATIONS.h"
#include "FIELD_UNIFORM_2D.h"

class DOUBLEMACH
{
public:
	DOUBLEMACH(){}
	~DOUBLEMACH(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = 0, x_upper = 4;
		const double y_lower = 0, y_upper = 1;
		gamma   = 1.4;
		t_final = 0.2;

		IDEAL_GAS_2D EOS(gamma);

		grid.Initialize(0, 0, (resolution-1)*4+1, resolution, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid.dx;
		const double dy = grid.dy;
		grid_x.Initialize(0, 0, (resolution-1) * 4 + 2, resolution, x_lower - 0.5*dx, y_lower, x_upper + 0.5*dx, y_upper);
		grid_y.Initialize(0, 0, (resolution - 1) * 4 + 1, resolution + 1, x_lower, y_lower - 0.5*dy, x_upper, y_upper + 0.5*dy);
		grid_xy.Initialize(0, 0, (resolution - 1) * 4 + 2, resolution + 1, x_lower - 0.5*dx, y_lower - 0.5*dy, x_upper + 0.5*dx, y_upper + 0.5*dy);
		

		Q.Initialize(grid,3);
		ITERATION_GRID_2D(grid, 3)
		{
			const double x = grid(i,j).x;
			const double y = grid(i,j).y;
			double rho;
			double u;
			double v;
			double p;

			if(y > sqrt(3.0) * (x - 1.0/6.0))
			{
				rho = 8;
				u = 8.25 * sqrt(3.0) * 0.5;
				v = -8.25 * 0.5;
				p = 116.5;
			}
			else
			{
				rho = 1.5;
				u = 0;
				v = 0;
				p = 1;
			}

			Q(i,j)[0] = rho;
			Q(i,j)[1] = rho * u;
			Q(i,j)[2] = rho * v;
			Q(i,j)[3] = EOS.E(rho, u, v, p);
		}
	}
};

class BOUNDARY_MANAGER_DMACH
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
	IDEAL_GAS_2D EOS;

	FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q;

public:
	BOUNDARY_MANAGER_DMACH(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_):Q(Q_), EOS(1.4)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER_DMACH(){}

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

	//Do nothing(inflow boundary condition)
	void SetLeftBoundaryCondition()
	{
		
		double rho;
		double u;
		double v;
		double p;

		rho = 8;
		u = 8.25 * sqrt(3.0) * 0.5;
		v = -8.25 * 0.5;
		p = 116.5;

		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_start-1; i >= i_start_g; i--)
			{
				Q(i, j)[0] = rho;
				Q(i, j)[1] = rho * u;
				Q(i, j)[2] = rho * v;
				Q(i, j)[3] = EOS.E(rho, u, v, p);
			}
		}
	}

	void SetRightBoundaryCondition()
	{
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
		const GRID_UNIFORM_2D& grid(Q.grid);
		double sqrt3 = sqrt(3.0);

		for(int i = i_start; i <= i_end; i ++)
		{
			int cnt = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				const double x = grid(i,j).x;
				if(x > 1.0/6.0)
				{
					Q(i,j)[0] =  Q(i, j_start)[0];
					Q(i,j)[1] =  Q(i, j_start)[1];
					Q(i,j)[2] = -Q(i, j_start)[2];
					Q(i,j)[3] =  Q(i, j_start)[3];
				}
				else
				{
					Q(i,j)[0] =  8;
					Q(i,j)[1] =  8 * 8.25 * sqrt3 * 0.5; 
					Q(i,j)[2] = -8 * 8.25 * 0.5; 
					Q(i,j)[3] = EOS.E(8, 8.25 * sqrt3 * 0.5, -8.25 * 0.5, 116.5);
				}
				cnt ++;
			}
		}
	}

	void SetTopBoundaryCondition(const double& t)
	{
		const GRID_UNIFORM_2D& grid(Q.grid);
		double sqrt3 = sqrt(3.0);

		for(int i = i_start; i <= i_end; i ++)
		{
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				const double x = grid(i,j).x;
				if(x < 1.0/6.0 + (1 + 20*t) / sqrt3)
				{
					Q(i,j)[0] =  8;
					Q(i,j)[1] =  8 * 8.25 * sqrt3 * 0.5; 
					Q(i,j)[2] = -8 * 8.25 * 0.5; 
					Q(i,j)[3] = EOS.E(8, 8.25 * sqrt3 * 0.5, -8.25 * 0.5, 116.5);
				}
				else
				{
					Q(i,j)[0] = 1.4;
					Q(i,j)[1] = 1.4 * 0;
					Q(i,j)[2] = 1.4 * 0;
					Q(i,j)[3] = EOS.E(1.4, 0, 0, 1);
				}
			}
		}
	}

	void SetDiagonalBoundaryCondition()
	{
		int cnt1 = 0;
		for(int i = i_start-1; i >= i_start_g; i --)
		{
			int cnt2 = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				Q(i,j) = Q(i_end-cnt1,j_end-cnt2);
				cnt2 ++;
			}

			cnt2 = 0;
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_end-cnt1,j_start+cnt2);
				cnt2 ++;
			}
			cnt1 ++;
		}

		cnt1 = 0;
		for(int i = i_end+1; i <= i_end_g; i ++)
		{
			int cnt2 = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				Q(i,j) = Q(i_start+cnt1, j_end-cnt2);
				cnt2 ++;
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_start+cnt1, j_start+cnt2);
				cnt2 ++;
			}
			cnt1 ++;
		}
	}

	void SetBoundaryCondition(const double& t)
	{
		SetLeftBoundaryCondition();
		SetRightBoundaryCondition();
		SetBottomBoundaryCondition();
		SetTopBoundaryCondition(t);
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