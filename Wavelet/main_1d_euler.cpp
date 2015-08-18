#include "WENO.h"
#include "WENOZ.h"
#include "ARRAY.h"
#include "ARRAY_1D.h"
#include "Problems/1D/PROBLEM_SHOCKTUBE.h"
#include "Wavelet/AWCM_1D.h"
#include "OMP_MANAGER.h"

using std::cout;
using std::endl;

void main()
{
	const int Nx_coarse = 16;
	const int num_level = 11;
	const int max_level = num_level-1;
	const double tol = 1e-6;

	MRA_GRID_INT_1D mra;
	mra.Initialize(num_level, Nx_coarse);
	AWCM_1D awcm(mra);

	const int resolution = mra.grid_finest.i_res;
	const double cfl = 0.5;
	const int num_print = 41;

	GRID_UNIFORM_1D grid, grid_x;
	FIELD_UNIFORM_1D<VECTOR_3D<double>> Q;
	double t_final, gamma;

	SHOCKTUBE::Initialize(resolution, grid, grid_x, Q, t_final, gamma);

	const double dx = grid.dx;
	const double one_over_dx = 1/dx;

	const IDEAL_GAS_1D EOS(gamma);

	double alpha_x;	//for LF-splitting.		
	grid.PrintGrid("x");

	FIELD_UNIFORM_1D<VECTOR_3D<double>> F(grid_x);

	ARRAY<double> arr_printing_times(num_print);
	for(int p = 0; p < num_print; p++)
	{
		double increment = t_final / (num_print-1);
		arr_printing_times[p] = increment * p;
	}

	BOUNDARY_MANAGER bc_manager(Q);

	double t = 0;
	double dt = 0;
	double max_velocity = 0;
	WENO_EULER_SOLVER_1D weno(gamma);
	WENOZ_EULER_SOLVER_1D wenoz(gamma);

	Q.Print();

	FIELD_UNIFORM_1D<VECTOR_3D<double>> Q_temp, Q_copy;
	while(true)
	{
		awcm.PrepareEuler(Q, tol);
		static int cnt = 1;
		static bool do_print = false;

		max_velocity = 0;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask3(i) == false) continue;
			const double c = EOS.c(Q(i));
			const double u = Q(i)[1] / Q(i)[0];
			max_velocity = MAX(max_velocity, (ABS(u)+c));
		}
		dt = cfl * dx / max_velocity;

		if(t+dt > arr_printing_times[cnt])
		{
			dt = arr_printing_times[cnt] - t;
			do_print = true;
		}

		bc_manager.SetBoundaryCondition();

		//preparation for LF splitting.
		alpha_x = 0;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask3(i) == false || awcm.level3(i) != max_level) continue;
			double u, v, c;
			u = Q(i)[1] / Q(i)[0];
			c = EOS.c(Q(i));

			alpha_x = MAX(alpha_x, ABS(u)+c);
		}

		Q_copy = Q;
		Q_temp = Q;
		ITERATION_GRID_1D(grid, 0)
		{
			
			if(awcm.mask2(i) == false) continue;
			const int lev = awcm.level2(i);
			const int dx_ = mra.dx[lev];

			VECTOR_3D<double> dF, F_l, F_r;
			if(lev == max_level)
			{
				weno.ComputedF(Q_temp(i-3), Q_temp(i-2), Q_temp(i-1), Q_temp(i), Q_temp(i+1), Q_temp(i+2), Q_temp(i+3), alpha_x, dF);
			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.DifferentiateEuler(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if (i == 2*dx_)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i > grid.i_end - 2*dx_)
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), dF);
				}
				else
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
			}
			
			Q(i) = Q(i) - dt / (dx*dx_) * dF;
		}

		bc_manager.SetBoundaryCondition();


		//preparation for LF splitting.
		alpha_x = 0;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask2(i) == false || awcm.level2(i) != max_level) continue;
			double u, v, c;
			u = Q(i)[1] / Q(i)[0];
			c = EOS.c(Q(i));

			alpha_x = MAX(alpha_x, ABS(u)+c);
		}

		Q_temp = Q;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask1(i) == false) continue;
			const int lev = awcm.level1(i);
			const int dx_ = mra.dx[lev];

			VECTOR_3D<double> dF, F_l, F_r;
			if(lev == max_level)
			{
				weno.ComputedF(Q_temp(i-3), Q_temp(i-2), Q_temp(i-1), Q_temp(i), Q_temp(i+1), Q_temp(i+2), Q_temp(i+3), alpha_x, dF);

			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.DifferentiateEuler(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if (i == 2*dx_)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i > grid.i_end - 2*dx_)
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), dF);
				}
				else
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
			}
			Q(i) = 0.75 * Q_copy(i) + 0.25 * (Q(i) - dt / (dx*dx_) * dF);

		}

		bc_manager.SetBoundaryCondition();


		//preparation for LF splitting.
		alpha_x = 0;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask1(i) == false || awcm.level1(i) != max_level) continue;
			double u, v, c;
			u = Q(i)[1] / Q(i)[0];
			c = EOS.c(Q(i));
		
			alpha_x = MAX(alpha_x, ABS(u)+c);
		}

		Q_temp = Q;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask0(i) == false) continue;
			const int lev = awcm.level0(i);
			const int dx_ = mra.dx[lev];

			VECTOR_3D<double> dF, F_l, F_r;
			if(lev == max_level)
			{
				weno.ComputedF(Q_temp(i-3), Q_temp(i-2), Q_temp(i-1), Q_temp(i), Q_temp(i+1), Q_temp(i+2), Q_temp(i+3), alpha_x, dF);
			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.DifferentiateEuler(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if (i == 2*dx_)
				{
					awcm.DifferentiateEulerAtLeftBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
				else if(i > grid.i_end - 2*dx_)
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), dF);
				}
				else
				{
					awcm.DifferentiateEulerAtRightBoundary(Q_temp(i-2*dx_), Q_temp(i-dx_), Q_temp(i), Q_temp(i+dx_), Q_temp(i+2*dx_), dF);
				}
			}
			Q(i) = (Q_copy(i) + 2.0 * (Q(i) - dt / (dx*dx_) * dF))/3.0;
			
		}

		awcm.PostProcessEuler(Q);

		t += dt;
		cout<<t<<endl;
		if(do_print == true)
		{
			char filename[10];
			sprintf_s(filename, "%4.4d", cnt);
			AWCM_MANAGER_1D::Print(Q, awcm.mask0, filename);
			cnt ++;
			do_print = false;
		}
		if(t >= t_final) break;
	}

	Q.Print("q");
	awcm.mask0.Print("mask");
	awcm.level0.Print("level");

}