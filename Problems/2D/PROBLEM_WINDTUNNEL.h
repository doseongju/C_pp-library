#pragma once

#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "IDEAL_GAS_EQUATIONS.h"

class WINDTUNNEL
{
public:
	WINDTUNNEL(){}
	~WINDTUNNEL(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = 0, x_upper = 3;
		const double y_lower = 0, y_upper = 1;
		gamma   = 1.4;
		t_final = 4;

		IDEAL_GAS_2D EOS(gamma);

		grid.Initialize(0, 0, 3 * resolution - 2, resolution, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid.dx;
		const double dy = grid.dy;
		grid_x.Initialize(0, 0, 3 * resolution - 1, resolution, x_lower - 0.5*dx, y_lower, x_upper + 0.5*dx, y_upper);
		grid_y.Initialize(0, 0, 3 * resolution - 2, resolution + 1, x_lower, y_lower - 0.5*dy, x_upper, y_upper + 0.5*dx);
		grid_xy.Initialize(0, 0, 3 * resolution - 1, resolution + 1, x_lower - 0.5*dx, y_lower - 0.5*dy, x_upper + 0.5*dx, y_upper + 0.5*dy);

		Q.Initialize(grid,4);
		ITERATION_GRID_2D(grid, 4)
		{
			const double x = grid(i,j).x;
			const double y = grid(i,j).y;
			double rho = 1.4;
			double u = 3;
			double v = 0;
			double p = 1;

			Q(i,j)[0] = rho;
			Q(i,j)[1] = rho * u;
			Q(i,j)[2] = rho * v;
			Q(i,j)[3] = EOS.E(rho, u, v, p);

		}
	}
};

class BOUNDARY_MANAGER_WT
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

	int i_bdry, j_bdry;

	FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q;

public:
	BOUNDARY_MANAGER_WT(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_):
	  Q(Q_)
	  {
		  Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
		  const double dx = Q.grid.dx;
		  const double dy = Q.grid.dy;
		  i_bdry = int(0.6/dx);
		  j_bdry = int(0.2/dy); 
	  }
	  ~BOUNDARY_MANAGER_WT(){}

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
		double rho = 1.4;
		double u = 3;
		double v = 0;
		double p = 1;

		IDEAL_GAS_2D EOS(1.4);
		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_start - 1; i >= i_start_g; i--)
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
			int cnt = 0;
			for(int i = i_end+1; i <= i_end_g; i ++)
			{
				Q(i,j) = Q(i_end, j);
			}
		}
	}


	void SetBottomBoundaryCondition(const bool x_direction, const bool y_direction)
	{
		if (x_direction == true && y_direction == true) assert(false);
		if (x_direction == false && y_direction == false) assert(false);

		if (x_direction == true)
		{
			for (int j = j_start_g; j <= j_bdry; j++)
			{
				for (int cnt = 0; cnt < 3; cnt++)
				{
					Q(i_bdry + cnt, j)[0] =  Q(i_bdry - 1, j)[0];
					Q(i_bdry + cnt, j)[1] = -Q(i_bdry - 1, j)[1];
					Q(i_bdry + cnt, j)[2] =  Q(i_bdry - 1, j)[2];
					Q(i_bdry + cnt, j)[3] =  Q(i_bdry - 1, j)[3];
				}
			}
		}

		else
		{
			//for (int i = i_start_g; i < i_bdry; i++)
			for (int i = i_start_g; i < i_end_g; i++)
			{
				for (int cnt = 0; cnt < 3; cnt++)
				{
					Q(i, j_start - 1 - cnt)[0] =  Q(i, j_start)[0];
					Q(i, j_start - 1 - cnt)[1] =  Q(i, j_start)[1];
					Q(i, j_start - 1 - cnt)[2] = -Q(i, j_start)[2];
					Q(i, j_start - 1 - cnt)[3] =  Q(i, j_start)[3];
				}
			}
			for (int i = i_bdry; i <= i_end_g; i++)
			{
				for (int cnt = 0; cnt < 3; cnt++)
				{
					Q(i, j_bdry - cnt)[0] =  Q(i, j_bdry + 1)[0];
					Q(i, j_bdry - cnt)[1] =  Q(i, j_bdry + 1)[1];
					Q(i, j_bdry - cnt)[2] = -Q(i, j_bdry + 1)[2];
					Q(i, j_bdry - cnt)[3] =  Q(i, j_bdry + 1)[3];
				}
			}
		}
	}
	void SetTopBoundaryCondition()
	{

		for(int i = i_start; i <= i_end; i ++)
		{
			int cnt = 0;
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j)[0] =  Q(i, j_end)[0];
				Q(i,j)[1] =  Q(i, j_end)[1];
				Q(i,j)[2] = -Q(i, j_end)[2];
				Q(i,j)[3] =  Q(i, j_end)[3];
				cnt ++;
			}
		}
	}

	//Except bottom boundary!!
	void SetBoundaryCondition()
	{
		SetLeftBoundaryCondition();
		SetRightBoundaryCondition();
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
		
		ITERATION_GRID_2D(mask, 0)
		{
			if (ABS(j - j_bdry) <= 2 && i >= i_bdry - 2)
			{
				mask(i, j) = true;
			}
			if (ABS(i - i_bdry) <= 2 && j <= j_bdry + 2)
			{
				mask(i, j) = true;
			}
		}
	}
};