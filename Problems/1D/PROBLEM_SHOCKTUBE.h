#pragma once

#include "FIELD_UNIFORM_1D.h"
#include "VECTOR_3D.h"
#include "IDEAL_GAS_EQUATIONS.h"

class SHOCKTUBE
{
public:
	SHOCKTUBE(){}
	~SHOCKTUBE(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_1D& grid, GRID_UNIFORM_1D& grid_x,
		FIELD_UNIFORM_1D<VECTOR_3D<double>>& Q, double& t_final, double& gamma)
	{
		//t_final = 0.1677;
		t_final = 2;
		gamma = 1.4;
		const double x_lower = 0, x_upper = 10;

		grid_x.Initialize(0, resolution+1, x_lower, x_upper);
		const double dx = grid_x.dx;
		grid.Initialize(0, resolution, x_lower+0.5*dx, x_upper-0.5*dx);

		/*grid.Initialize(0, resolution, x_lower, x_upper);
		const double dx = grid.dx;
		grid_x.Initialize(0, resolution + 1, x_lower - 0.5*dx, x_upper + 0.5*dx);*/

		const IDEAL_GAS_1D EOS(gamma);

		Q.Initialize(grid, 3);
		ITERATION_GRID_1D(grid, 3)
		{
			const double x = grid(i);
			double rho, p, u;
			if(x <= 5)
			{
				rho = 1;
				u = 0;
				p = 1;
			}
			else
			{
				rho = 0.125;
				u = 0;
				p = 0.1;
			}
			Q(i)[0] = rho;
			Q(i)[1] = rho*u;
			Q(i)[2] = EOS.E(rho, u, p);
		}
	}
};

class BOUNDARY_MANAGER
{
public:
	int i_start;
	int i_start_g;
	int i_end;
	int i_end_g;
	int ghost_width;

	FIELD_UNIFORM_1D<VECTOR_3D<double>>& Q;

public:
	BOUNDARY_MANAGER(FIELD_UNIFORM_1D<VECTOR_3D<double>>& Q_):Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER(){}

public:
	void Initialize(const int i_start_, const int i_start_g_, const int i_end_, const int i_end_g_, const int ghost_width_)
	{
		i_start       = i_start_;
		i_start_g     = i_start_g_;
		i_end         = i_end_;
		i_end_g       = i_end_g_;
		ghost_width   = ghost_width_;
	}

public:
	void SetLeftBoundaryCondition()
	{
		for(int i = i_start-1; i >= i_start_g; i --) 
		{
			Q(i) =  Q(i_start);

		}
	}

	void SetRightBoundaryCondition()
	{
		for(int i = i_end+1; i <= i_end_g; i ++) 
		{
			Q(i) =  Q(i_end);
		}
	}

	void SetBoundaryCondition()
	{
		SetLeftBoundaryCondition ();
		SetRightBoundaryCondition();
	}
};