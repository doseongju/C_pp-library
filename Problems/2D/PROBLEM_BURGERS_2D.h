#pragma once

#include "FIELD_UNIFORM_2D.h"

class BURGERS_2D
{
public:
	BURGERS_2D(){}
	~BURGERS_2D(){}

public:
	static void Initialize(const int& resolution,
		GRID_UNIFORM_2D& grid, GRID_UNIFORM_2D& grid_x, GRID_UNIFORM_2D& grid_y, GRID_UNIFORM_2D& grid_xy,
		FIELD_UNIFORM_2D<double>& Q, double& t_final)
	{
		const double x_lower = 0, x_upper = 1;
		const double y_lower = 0, y_upper = 1;
		//t_final = 2;

		t_final = 0.65;

		grid_xy.Initialize(0,0, resolution+1, resolution+1, x_lower, y_lower, x_upper, y_upper);
		const double dx = grid_xy.dx;
		const double dy = grid_xy.dy;
		grid_x.Initialize(0,0, resolution+1, resolution  , x_lower, y_lower+0.5*dx, x_upper, y_upper-0.5*dx);
		grid_y.Initialize(0,0, resolution  , resolution+1, x_lower+0.5*dy, y_lower, x_upper-0.5*dy, y_upper);
		grid.Initialize(0,0, resolution, resolution, x_lower+0.5*dx, y_lower+0.5*dy, x_upper-0.5*dx, y_upper-0.5*dy);

		Q.Initialize(grid,3);
		ITERATION_GRID_2D(grid, 3)
		{

			const double x = grid(i,j).x;
			const double y = grid(i,j).y;
			const double r = sqrt(x*x + y*y);

			Q(i,j) = sin(2 * PI*r) + sin(PI*r) / 3.0;

			/*if (ABS(x - 0.35) <= 0.25 && ABS(y - 0.35) <= 0.25)
			{
				Q(i, j) = 1;
			}
			else
			{
				Q(i, j) = 0;
			}*/
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
	int j_start;
	int j_start_g;
	int j_end;
	int j_end_g;
	int ghost_width;

	FIELD_UNIFORM_2D<double>& Q;

public:
	BOUNDARY_MANAGER(FIELD_UNIFORM_2D<double>& Q_):Q(Q_)
	{
		Initialize(Q.i_start, Q.i_start_g, Q.i_end, Q.i_end_g, Q.j_start, Q.j_start_g, Q.j_end, Q.j_end_g, Q.ghost_width);
	}
	~BOUNDARY_MANAGER(){}

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
				//Q(i,j) = Q(i_end-cnt, j);
				//cnt ++;
				Q(i, j) = Q(i_start, j);
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
				//Q(i,j) = Q(i_start+cnt, j);
				//cnt ++;
				Q(i, j) = Q(i_end, j);
			}
		}
	}

	void SetBottomBoundaryCondition()
	{

		for(int i = i_start; i <= i_end; i ++)
		{
			int cnt = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				//Q(i,j) = Q(i, j_end-cnt);
				//cnt ++;
				Q(i, j) = Q(i, j_start);
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
				//Q(i,j) = Q(i, j_start+cnt);
				//cnt ++;
				Q(i, j) = Q(i, j_end);
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
				//cnt2 ++;
			}

			cnt2 = 0;
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_end-cnt1,j_start+cnt2);
				//cnt2 ++;
			}
			//cnt1 ++;
		}

		cnt1 = 0;
		for(int i = i_end+1; i <= i_end_g; i ++)
		{
			int cnt2 = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				Q(i,j) = Q(i_start+cnt1, j_end-cnt2);
				//cnt2 ++;
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				Q(i,j) = Q(i_start+cnt1, j_start+cnt2);
				//cnt2 ++;
			}
			//cnt1 ++;
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