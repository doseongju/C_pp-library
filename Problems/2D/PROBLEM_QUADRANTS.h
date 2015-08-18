#pragma once

#include "FIELD_UNIFORM_2D.h"
#include "VECTOR_4D.h"
#include "IDEAL_GAS_EQUATIONS.h"
#include <omp.h>

class QUADRANTS
{
public:
	QUADRANTS(){}
	~QUADRANTS(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q, double& t_final, double& gamma)
	{

		const double x_lower = 0, x_upper = 1;
		const double y_lower = 0, y_upper = 1;
		gamma   = 1.4;
		t_final = 0.8;

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

			const double x = grid(i,j).x;
			const double y = grid(i,j).y;
			double rho, u, v, p;

			//1-quadrant
			if(x >= 0.8 && y >= 0.8)
			{
				rho = 1.5;
				u   = 0;
				v   = 0;
				p   = 1.5;

			}
			//2-quadrant
			else if(x < 0.8 && y >= 0.8)
			{
				rho = 0.532258064516129;
				u   = 1.206045378311055;
				v   = 0;
				p   = 0.3;
			}
			//3-quadrant
			else if(x < 0.8 && y < 0.8)
			{
				rho = 0.137992831541219;
				u   = 1.206045378311055;
				v   = 1.206045378311055;
				p   = 0.029032258064516;
			}
			//4-quadrant
			else
			{
				rho = 0.532258064516129;
				u   = 0;
				v   = 1.206045378311055;
				p   = 0.3;
			}
			Q(i,j)[0] = rho;
			Q(i,j)[1] = rho*u;
			Q(i,j)[2] = rho*v;
			Q(i,j)[3] = EOS.E(rho, u, v, p);

		}
	}
};

//Zero-order extrapolation
class BOUNDARY_MANAGER_QUAD
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
	BOUNDARY_MANAGER_QUAD(FIELD_UNIFORM_2D<VECTOR_4D<double>>& Q_):Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER_QUAD(){}

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
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				Q(i,j) = Q(i, j_start);
			}
		}
	}

	void SetTopBoundaryCondition()
	{
#pragma omp parallel for
		for(int i = i_start; i <= i_end; i ++)
		{
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i, j_end);
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