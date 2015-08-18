#include "Problems/2D/PROBLEM_BOX_ADVECTION.h"
#include "AWCM_2D.h"
#include "STOP_WATCH.h"
#include "WENO.h"


void main()
{
	
	const int num_lev = 6;
	const int Nx_coarse = 15;
	const int Ny_coarse = 15;
	MRA_GRID_INT_2D mra_grid;
	mra_grid.Initialize(num_lev, Nx_coarse, Ny_coarse);
	const double tol = 1e-4;
	const int max_lev = num_lev - 1;

	const int resolution = mra_grid.grid_finest.i_res;

	GRID_UNIFORM_2D grid, grid_x, grid_y, grid_xy;
	FIELD_UNIFORM_2D<double> Q;
	double t_final;

	BOX_ADVECTION::Initialize(resolution, grid, grid_x, grid_y, grid_xy, Q, t_final);

	FIELD_UNIFORM_2D<double> Q_0(grid, 3), Q_1(grid, 3), Q_2(grid, 3);
	grid.PrintGridX("x");
	grid.PrintGridY("y");

	Q.Print("Q");

	AWCM_2D awcm(mra_grid);
	WENO5 weno;
	const double dt = 0.4 * grid.dx;
	BOUNDARY_MANAGER boundary_manager(Q);
	BOUNDARY_MANAGER boundary_manager_0(Q_0);
	BOUNDARY_MANAGER boundary_manager_1(Q_1);
	BOUNDARY_MANAGER boundary_manager_2(Q_2);


	STOP_WATCH stop_watch;
	stop_watch.Tic();
	for (int iter = 1; iter <= 100; iter++)
	{
		std::cout << iter << std::endl;
		awcm.Prepare(Q, tol);

		boundary_manager.SetBoundaryCondition();

		ITERATION_GRID_2D(grid, 0)
		{
			double dF = 0, dG = 0;
			if (awcm.mask2(i,j) == false) continue;
			const int lev = awcm.level2(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				double a, b;
				weno.ReconstructRight(Q(i - 3, j), Q(i - 2, j), Q(i - 1, j), Q(i    , j), Q(i + 1, j), a);
				weno.ReconstructRight(Q(i - 2, j), Q(i - 1, j), Q(i    , j), Q(i + 1, j), Q(i + 2, j), b);
				dF = b - a;
				weno.ReconstructRight(Q(i, j - 3), Q(i, j - 2), Q(i, j - 1), Q(i, j    ), Q(i, j + 1), a);
				weno.ReconstructRight(Q(i, j - 2), Q(i, j - 1), Q(i, j    ), Q(i, j + 1), Q(i, j + 2), b);
				dG = b - a;
			}
			else
			{
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q(i - d, j), Q(i, j), Q(i + d, j), Q(i + 2*d, j), dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q(i - 2*d, j), Q(i - d, j), Q(i, j), Q(i + d, j), dF);
				}
				else
				{
					awcm.Differentiate(Q(i - 2*d, j), Q(i - d, j), Q(i, j), Q(i + d, j), Q(i + 2*d, j), dF);
				}

				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q(i, j-d), Q(i, j), Q(i, j+d), Q(i, j+2*d), dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q(i, j - 2*d), Q(i, j - d), Q(i, j), Q(i, j + d), dG);
				}
				else
				{
					awcm.Differentiate(Q(i, j - 2*d), Q(i, j - d), Q(i, j), Q(i, j + d), Q(i, j + 2*d), dG);
				}
			}
			Q_2(i, j) = Q(i, j) - dt / dx * (dF+dG);
		}


		boundary_manager_2.SetBoundaryCondition();
		ITERATION_GRID_2D(grid, 0)
		{
			double dF = 0, dG = 0;
			if (awcm.mask1(i, j) == false) continue;
			const int lev = awcm.level1(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				double a, b;
			
				weno.ReconstructRight(Q_2(i - 3, j), Q_2(i - 2, j), Q_2(i - 1, j), Q_2(i, j), Q_2(i + 1, j), a);
				weno.ReconstructRight(Q_2(i - 2, j), Q_2(i - 1, j), Q_2(i, j), Q_2(i + 1, j), Q_2(i + 2, j), b);
				dF = b - a;
				weno.ReconstructRight(Q_2(i, j - 3), Q_2(i, j - 2), Q_2(i, j - 1), Q_2(i, j), Q_2(i, j + 1), a);
				weno.ReconstructRight(Q_2(i, j - 2), Q_2(i, j - 1), Q_2(i, j), Q_2(i, j + 1), Q_2(i, j + 2), b);
				dG = b - a;
			}
			else
			{
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q_2(i - d, j), Q_2(i, j), Q_2(i + d, j), Q_2(i + 2 * d, j), dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q_2(i - 2 * d, j), Q_2(i - d, j), Q_2(i, j), Q_2(i + d, j), dF);
				}
				else
				{
					awcm.Differentiate(Q_2(i - 2 * d, j), Q_2(i - d, j), Q_2(i, j), Q_2(i + d, j), Q_2(i + 2 * d, j), dF);
				}

				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q_2(i, j - d), Q_2(i, j), Q_2(i, j + d), Q_2(i, j + 2 * d), dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q_2(i, j - 2 * d), Q_2(i, j - d), Q_2(i, j), Q_2(i, j + d), dG);
				}
				else
				{
					awcm.Differentiate(Q_2(i, j - 2 * d), Q_2(i, j - d), Q_2(i, j), Q_2(i, j + d), Q_2(i, j + 2 * d), dG);
				}
			}
			Q_1(i, j) = 0.75 * Q(i,j) + 0.25 * (Q_2(i, j) - dt / dx * (dF + dG));
		}

		boundary_manager_1.SetBoundaryCondition();
		ITERATION_GRID_2D(grid, 0)
		{
			double dF = 0, dG = 0;
			if (awcm.mask0(i, j) == false) continue;
			const int lev = awcm.level0(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				double a, b;

				weno.ReconstructRight(Q_1(i - 3, j), Q_1(i - 2, j), Q_1(i - 1, j), Q_1(i, j), Q_1(i + 1, j), a);
				weno.ReconstructRight(Q_1(i - 2, j), Q_1(i - 1, j), Q_1(i, j), Q_1(i + 1, j), Q_1(i + 2, j), b);
				dF = b - a;
				weno.ReconstructRight(Q_1(i, j - 3), Q_1(i, j - 2), Q_1(i, j - 1), Q_1(i, j), Q_1(i, j + 1), a);
				weno.ReconstructRight(Q_1(i, j - 2), Q_1(i, j - 1), Q_1(i, j), Q_1(i, j + 1), Q_1(i, j + 2), b);
				dG = b - a;
			}
			else
			{
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q_1(i - d, j), Q_1(i, j), Q_1(i + d, j), Q_1(i + 2 * d, j), dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q_1(i - 2 * d, j), Q_1(i - d, j), Q_1(i, j), Q_1(i + d, j), dF);
				}
				else
				{
					awcm.Differentiate(Q_1(i - 2 * d, j), Q_1(i - d, j), Q_1(i, j), Q_1(i + d, j), Q_1(i + 2 * d, j), dF);
				}

				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(Q_1(i, j - d), Q_1(i, j), Q_1(i, j + d), Q_1(i, j + 2 * d), dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(Q_1(i, j - 2 * d), Q_1(i, j - d), Q_1(i, j), Q_1(i, j + d), dG);
				}
				else
				{
					awcm.Differentiate(Q_1(i, j - 2 * d), Q_1(i, j - d), Q_1(i, j), Q_1(i, j + d), Q_1(i, j + 2 * d), dG);
				}
			}
			Q_0(i, j) = Q_1(i, j) - dt / dx * (dF + dG);
		}


		ITERATION_GRID_2D(grid, 0)
		{
			if (awcm.mask0(i,j) == true) 
			{
				Q(i, j) = (1.0/3.0*Q(i, j) + 2.0/3.0*Q_0(i, j));
				//Q(i, j) = Q_2(i, j);
			}
		}

		awcm.PostProcess(Q);
	}
	std::cout << stop_watch.Toc();

	AWCM_MANAGER_2D::Print(Q, awcm.mask0, "data");
	FIELD_UNIFORM_2D<bool> mask_init(grid);
	mask_init.AssignAllValues(true);
	awcm.wavelet.Interpolate(Q, mra_grid, awcm.mask0, mask_init);
	Q.Print("Q_");
	awcm.level0.Print("level");

}