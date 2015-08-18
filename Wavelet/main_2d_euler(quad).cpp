#include "OMP_MANAGER.h"
#include "Problems/2D/PROBLEM_QUADRANTS.h"
#include "AWCM_2D.h"
#include "STOP_WATCH.h"
#include "WENO.h"
#include "WENOZ.h"
#include "EULER_EQUATION.h"

void main()
{
	OMP_MANAGER omp_manager(12);
	const int num_lev = 6;
	const int Nx_coarse = 15;
	const int Ny_coarse = 15;
	MRA_GRID_INT_2D mra_grid;
	mra_grid.Initialize(num_lev, Nx_coarse, Ny_coarse);
	const double tol = 1e-4;
	const int max_lev = num_lev - 1;
	double gamma;
	const double cfl = 0.5;
	const int num_print = 9;

	const int resolution = mra_grid.grid_finest.j_res;

	GRID_UNIFORM_2D grid, grid_x, grid_y, grid_xy;
	FIELD_UNIFORM_2D<VECTOR_4D<double>> Q;
	double t_final;

	QUADRANTS::Initialize(resolution, grid, grid_x, grid_y, grid_xy, Q, t_final, gamma);

	EULER_EQUATION_2D euler_eqn(gamma);
	FIELD_UNIFORM_2D<VECTOR_4D<double>> Q_0(grid, 3), Q_1(grid, 3), Q_2(grid, 3);
	ARRAY_1D<double> alpha_x(0, resolution), alpha_y(0,resolution);
	grid.PrintGridX("x");
	grid.PrintGridY("y");
	Q.TagName("Q");

	AWCM_2D awcm(mra_grid);
	double t(0), dt;
	static int cnt = 0;
	BOUNDARY_MANAGER_QUAD boundary_manager(Q);
	BOUNDARY_MANAGER_QUAD boundary_manager_0(Q_0);
	BOUNDARY_MANAGER_QUAD boundary_manager_1(Q_1);
	BOUNDARY_MANAGER_QUAD boundary_manager_2(Q_2);
	const double dx = grid.dx;
	WENO_EULER_SOLVER_2D weno_solver;

	ARRAY<double> arr_printing_times(num_print);
	for (int p = 0; p < num_print; p++)
	{
		double increment = t_final / (num_print - 1);
		arr_printing_times[p] = increment * p;
	}
	Q.Print();
	cnt++;
	while (true)
	{
		awcm.PrepareEuler(Q, tol);

		//Determine the time step.
		double max_velocity = 0;

#pragma omp parallel
		{
			double max_velocity = 0;

#pragma omp for
			ITERATION_GRID_2D(grid, 0)
			{
				if (awcm.mask3(i, j) == false) continue;
				max_velocity = MAX(max_velocity, euler_eqn.Spectral_RadiusX(Q(i, j)) + euler_eqn.Spectral_RadiusY(Q(i, j)));
			}
			omp_manager.Input(max_velocity);
		}
		max_velocity = omp_manager.Max_double();
		dt = cfl * dx / max_velocity;

		static bool do_print = false;
		if (t + dt > arr_printing_times[cnt])
		{
			dt = arr_printing_times[cnt] - t;
			do_print = true;
		}

		boundary_manager.SetBoundaryCondition();
		alpha_x.AssignAllValue(0), alpha_y.AssignAllValue(0);
		
#pragma omp parallel
{

#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			if (awcm.mask3(i, j) == 0) continue;

			alpha_x(j) = MAX(alpha_x(j), euler_eqn.Spectral_RadiusX(Q(i, j)));
		}
#pragma omp for
		ITERATION_GRID_2D_IJ(grid, 0)
		{
			if (awcm.mask3(i, j) == 0) continue;

			alpha_y(i) = MAX(alpha_y(i), euler_eqn.Spectral_RadiusY(Q(i, j)));
		}

		WENO_EULER_SOLVER_2D weno_solver;
		WENOZ_EULER_SOLVER_2D wenoz_solver;
#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			VECTOR_4D<double> dF, dG;
			VECTOR_4D<double> F_ll, F_l, F_c, F_r, F_rr;
			VECTOR_4D<double> G_bb, G_b	, G_c, G_t, G_tt;
			if (awcm.mask2(i, j) == false) continue;
			const int lev = awcm.level2(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				VECTOR_4D<double> a, b;

				wenoz_solver.ComputeFluxX(Q(i - 3, j), Q(i - 2, j), Q(i - 1, j), Q(i, j), Q(i + 1, j), Q(i + 2, j), alpha_x(j), a);
				wenoz_solver.ComputeFluxX(Q(i - 2, j), Q(i - 1, j), Q(i, j), Q(i + 1, j), Q(i + 2, j), Q(i + 3, j), alpha_x(j), b);
				dF = b - a;
				wenoz_solver.ComputeFluxY(Q(i, j - 3), Q(i, j - 2), Q(i, j - 1), Q(i, j), Q(i, j + 1), Q(i, j + 2), alpha_y(i), a);
				wenoz_solver.ComputeFluxY(Q(i, j - 2), Q(i, j - 1), Q(i, j), Q(i, j + 1), Q(i, j + 2), Q(i, j + 3), alpha_y(i), b);
				dG = b - a;
			}
			else
			{
				euler_eqn.ComputeFu(Q(i - 2 * d, j), F_ll), euler_eqn.ComputeFu(Q(i - d, j), F_l), euler_eqn.ComputeFu(Q(i, j), F_c);
				euler_eqn.ComputeFu(Q(i + d, j), F_r), euler_eqn.ComputeFu(Q(i + 2 * d, j), F_rr);
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(F_l, F_c, F_r, F_rr, dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(F_ll, F_l, F_c, F_r, dF);
				}
				else
				{
					awcm.Differentiate(F_ll, F_l, F_c, F_r, F_rr, dF);
				}

				euler_eqn.ComputeGu(Q(i, j - 2 * d), G_bb), euler_eqn.ComputeGu(Q(i, j - d), G_b), euler_eqn.ComputeGu(Q(i, j), G_c);
				euler_eqn.ComputeGu(Q(i, j + d), G_t), euler_eqn.ComputeGu(Q(i, j + 2 * d), G_tt);
				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(G_b, G_c, G_t,G_tt, dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(G_bb, G_b, G_c, G_t, dG);
				}
				else
				{
					awcm.Differentiate(G_bb, G_b, G_c, G_t, G_tt, dG);
				}
			}
			Q_2(i, j) = Q(i, j) - dt / dx * (dF + dG);
		}
}

		boundary_manager_2.SetBoundaryCondition();
		alpha_x.AssignAllValue(0), alpha_y.AssignAllValue(0);
		
#pragma omp parallel
{

#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			if (awcm.mask2(i, j) == 0) continue;
			alpha_x(j) = MAX(alpha_x(j), euler_eqn.Spectral_RadiusX(Q_2(i, j)));
		}
#pragma omp for
		ITERATION_GRID_2D_IJ(grid, 0)
		{
			if (awcm.mask2(i, j) == 0) continue;
			alpha_y(i) = MAX(alpha_y(i), euler_eqn.Spectral_RadiusY(Q_2(i, j)));
		}


		WENO_EULER_SOLVER_2D weno_solver;
		WENOZ_EULER_SOLVER_2D wenoz_solver;
#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			VECTOR_4D<double> dF, dG;
			VECTOR_4D<double> F_ll, F_l, F_c, F_r, F_rr;
			VECTOR_4D<double> G_bb, G_b, G_c, G_t, G_tt;
			if (awcm.mask1(i, j) == false) continue;
			const int lev = awcm.level1(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				VECTOR_4D<double> a, b;
				wenoz_solver.ComputeFluxX(Q_2(i - 3, j), Q_2(i - 2, j), Q_2(i - 1, j), Q_2(i, j), Q_2(i + 1, j), Q_2(i + 2, j), alpha_x(j), a);
				wenoz_solver.ComputeFluxX(Q_2(i - 2, j), Q_2(i - 1, j), Q_2(i, j), Q_2(i + 1, j), Q_2(i + 2, j), Q_2(i + 3, j), alpha_x(j), b);
				dF = b - a;
				wenoz_solver.ComputeFluxY(Q_2(i, j - 3), Q_2(i, j - 2), Q_2(i, j - 1), Q_2(i, j), Q_2(i, j + 1), Q_2(i, j + 2), alpha_y(i), a);
				wenoz_solver.ComputeFluxY(Q_2(i, j - 2), Q_2(i, j - 1), Q_2(i, j), Q_2(i, j + 1), Q_2(i, j + 2), Q_2(i, j + 3), alpha_y(i), b);
				dG = b - a;

			}
			else
			{
				euler_eqn.ComputeFu(Q_2(i - 2 * d, j), F_ll), euler_eqn.ComputeFu(Q_2(i - d, j), F_l), euler_eqn.ComputeFu(Q_2(i, j), F_c);
				euler_eqn.ComputeFu(Q_2(i + d, j), F_r), euler_eqn.ComputeFu(Q_2(i + 2 * d, j), F_rr);
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(F_l, F_c, F_r, F_rr, dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(F_ll, F_l, F_c, F_r, dF);
				}
				else
				{
					awcm.Differentiate(F_ll, F_l, F_c, F_r, F_rr, dF);
				}

				euler_eqn.ComputeGu(Q_2(i, j - 2 * d), G_bb), euler_eqn.ComputeGu(Q_2(i, j - d), G_b), euler_eqn.ComputeGu(Q_2(i, j), G_c);
				euler_eqn.ComputeGu(Q_2(i, j + d), G_t), euler_eqn.ComputeGu(Q_2(i, j + 2 * d), G_tt);
				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(G_b, G_c, G_t, G_tt, dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(G_bb, G_b, G_c, G_t, dG);
				}
				else
				{
					awcm.Differentiate(G_bb, G_b, G_c, G_t, G_tt, dG);
				}
			}
			Q_1(i, j) = 0.75 * Q(i, j) + 0.25 * (Q_2(i, j) - dt / dx * (dF + dG));
		}
}

		boundary_manager_1.SetBoundaryCondition();
		alpha_x.AssignAllValue(0), alpha_y.AssignAllValue(0);
		
#pragma omp parallel
{

#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			if (awcm.mask1(i, j) == 0) continue;
			alpha_x(j) = MAX(alpha_x(j), euler_eqn.Spectral_RadiusX(Q_1(i, j)));
		}
#pragma omp for
		ITERATION_GRID_2D_IJ(grid, 0)
		{
			if (awcm.mask1(i, j) == 0) continue;
			alpha_y(i) = MAX(alpha_y(i), euler_eqn.Spectral_RadiusY(Q_1(i, j)));
		}

		WENO_EULER_SOLVER_2D weno_solver;
		WENOZ_EULER_SOLVER_2D wenoz_solver;
#pragma omp for
		ITERATION_GRID_2D(grid, 0)
		{
			VECTOR_4D<double> dF, dG;
			VECTOR_4D<double> F_ll, F_l, F_c, F_r, F_rr;
			VECTOR_4D<double> G_bb, G_b, G_c, G_t, G_tt;
			if (awcm.mask0(i, j) == false) continue;
			const int lev = awcm.level0(i, j);
			const int d = mra_grid.grid[lev].dx;
			const double dx = grid.dx * d;

			if (lev == max_lev)
			{
				VECTOR_4D<double> a, b;

				wenoz_solver.ComputeFluxX(Q_1(i - 3, j), Q_1(i - 2, j), Q_1(i - 1, j), Q_1(i, j), Q_1(i + 1, j), Q_1(i + 2, j), alpha_x(j), a);
				wenoz_solver.ComputeFluxX(Q_1(i - 2, j), Q_1(i - 1, j), Q_1(i, j), Q_1(i + 1, j), Q_1(i + 2, j), Q_1(i + 3, j), alpha_x(j), b);
				dF = b - a;
				wenoz_solver.ComputeFluxY(Q_1(i, j - 3), Q_1(i, j - 2), Q_1(i, j - 1), Q_1(i, j), Q_1(i, j + 1), Q_1(i, j + 2), alpha_y(i), a);
				wenoz_solver.ComputeFluxY(Q_1(i, j - 2), Q_1(i, j - 1), Q_1(i, j), Q_1(i, j + 1), Q_1(i, j + 2), Q_1(i, j + 3), alpha_y(i), b);
				dG = b - a;
			}
			else
			{
				euler_eqn.ComputeFu(Q_1(i - 2 * d, j), F_ll), euler_eqn.ComputeFu(Q_1(i - d, j), F_l), euler_eqn.ComputeFu(Q_1(i, j), F_c);
				euler_eqn.ComputeFu(Q_1(i + d, j), F_r), euler_eqn.ComputeFu(Q_1(i + 2 * d, j), F_rr);
				if (i <= d)
				{
					awcm.DifferentiateAtLeftBoundary(F_l, F_c, F_r, F_rr, dF);
				}
				else if (i >= grid.i_end - d)
				{
					awcm.DifferentiateAtRightBoundary(F_ll, F_l, F_c, F_r, dF);
				}
				else
				{
					awcm.Differentiate(F_ll, F_l, F_c, F_r, F_rr, dF);
				}

				euler_eqn.ComputeGu(Q_1(i, j - 2 * d), G_bb), euler_eqn.ComputeGu(Q_1(i, j - d), G_b), euler_eqn.ComputeGu(Q_1(i, j), G_c);
				euler_eqn.ComputeGu(Q_1(i, j + d), G_t), euler_eqn.ComputeGu(Q_1(i, j + 2 * d), G_tt);
				if (j <= d)
				{
					awcm.DifferentiateAtLeftBoundary(G_b, G_c, G_t, G_tt, dG);
				}
				else if (j >= grid.j_end - d)
				{
					awcm.DifferentiateAtRightBoundary(G_bb, G_b, G_c, G_t, dG);
				}
				else
				{
					awcm.Differentiate(G_bb, G_b, G_c, G_t, G_tt, dG);
				}
			}
			Q_0(i, j) = Q_1(i, j) - dt / dx * (dF + dG);
		}
}

#pragma omp parallel for
		ITERATION_GRID_2D(grid, 0)
		{
			if (awcm.mask0(i, j) == true)
			{
				Q(i, j) = (Q(i, j) + 2.0*Q_0(i, j)) / 3.0;
			}
		}

		awcm.PostProcessEuler(Q);
		
		t += dt;
		std::cout.precision(7);
		std::cout << t << std::endl;

		if (do_print == true)
		{
			Q.Print();
			char file_name[10];
			sprintf_s(file_name, "data_%4.4d",cnt);
			AWCM_MANAGER_2D::Print(Q, awcm.mask0, file_name);
			cnt++;
			do_print = false;
		}
		if (t >= t_final) break;
	}
}