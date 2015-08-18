#include "FIELD_UNIFORM_2D.h"
#include "ARRAY.h"
#include "ARRAY_1D.h"
#include "MRA_GRID_UNIFORM_2D.h"
#include "WAVELET_2D.h"
#include "AWCM_2D.h"
#include "BoundaryConditions/2D/BOUNDARY_QUADRANT.h"
#include "WENO.h"
#include "ROESOLVER.h"
#include "EULER_EQUATION.h"


void main()
{
	//wavelet parameters.
	const int num_levels = 7;		//lev = 0, ... ,num_levels-1
	const int x_res_coarse = 13;
	const int y_res_coarse = 13;
	const int N = 2;	//TODO : make N tunable. Only the case N = 2 is working!
	const double tol = 1e-4;

	//problem parameters
	const double x_lower = 0, x_upper = 1;
	const double y_lower = 0, y_upper = 1;
	const int resolution = (x_res_coarse-1) * pow(2.0, num_levels-1) + 1;
	const double cfl = 0.5;
	const int num_print = 81;
	const double t_final = 0.8;
	const double gamma = 1.4;
	const IDEAL_GAS_2D EOS(gamma);
	const EULER_EQUATION_2D euler_eqn(gamma);

	//Preparing MRA grid.
	MRA_GRID_INT_2D mra_grid;
	mra_grid.Initialize(num_levels, x_res_coarse, y_res_coarse);

	GRID_UNIFORM_2D grid_macxy(0,0, resolution+1, resolution+1, x_lower, y_lower, x_upper, y_upper);
	const double dx = grid_macxy.dx;
	const double dy = grid_macxy.dy;
	const double one_over_dx = 1/dx;
	const double one_over_dy = 1/dy;
	GRID_UNIFORM_2D grid_macx (0,0, resolution+1, resolution  , x_lower, y_lower+0.5*dx, x_upper, y_upper-0.5*dx);
	GRID_UNIFORM_2D grid_macy (0,0, resolution  , resolution+1, x_lower+0.5*dy, y_lower, x_upper-0.5*dy, y_upper);
	GRID_UNIFORM_2D grid(0,0, resolution, resolution, x_lower+0.5*dx, y_lower+0.5*dy, x_upper-0.5*dx, y_upper-0.5*dy);
	ARRAY_1D<double> alpha_x(0, resolution), alpha_y(0, resolution);	//for LF-splitting.	
	FIELD_UNIFORM_2D<bool> mask(grid,300);
	FIELD_UNIFORM_2D<int> level(grid,300);
	grid.PrintGridX("x");
	grid.PrintGridY("y");


	FIELD_UNIFORM_2D<VECTOR_4D<double>> dF(grid), dG(grid), Q(grid, 3);
	FIELD_UNIFORM_2D<VECTOR_4D<double>> fu(grid, 3), gu(grid, 3);

	ARRAY<double> arr_printing_times(num_print);
	for(int p = 0; p < num_print; p++)
	{
		double increment = t_final / (num_print-1);
		arr_printing_times[p] = increment * p;
	}

	ITERATION_GRID_2D(grid, 3)
	{

		const double x = grid(i,j).x;
		const double y = grid(i,j).y;

		//1-quadrant
		if(x >= 0.8 && y >= 0.8)
		{
			Q(i,j)[0] = 1.5;
			Q(i,j)[1] = 1.5 * 0;
			Q(i,j)[2] = 1.5 * 0;
			Q(i,j)[3] = EOS.E(1.5, 0, 0, 1.5);

		}
		//2-quadrant
		else if(x < 0.8 && y >= 0.8)
		{
			Q(i,j)[0] = 0.532258064516129;
			Q(i,j)[1] = 0.532258064516129 * 1.206045378311055;
			Q(i,j)[2] = 0.532258064516129 * 0;
			Q(i,j)[3] = EOS.E(0.532258064516129, 1.206045378311055, 0, 0.3);
		}
		//3-quadrant
		else if(x < 0.8 && y < 0.8)
		{
			Q(i,j)[0] = 0.137992831541219;
			Q(i,j)[1] = 0.137992831541219 * 1.206045378311055;
			Q(i,j)[2] = 0.137992831541219 * 1.206045378311055;
			Q(i,j)[3] = EOS.E(0.137992831541219, 1.206045378311055, 1.206045378311055, 0.029032258064516);
		}
		//4-quadrant
		else
		{
			Q(i,j)[0] = 0.532258064516129;
			Q(i,j)[1] = 0.532258064516129 * 0;
			Q(i,j)[2] = 0.532258064516129 * 1.206045378311055;
			Q(i,j)[3] = EOS.E(0.532258064516129, 0, 1.206045378311055, 0.3);
		}
	}

	double t = 0;
	double dt = 0;
	double max_velocity = 0;
	WENO_EULER_SOLVER_2D weno(gamma);
	BOUNDARY_MANAGER_2D bc_manager(Q);
	mask.AssignAllValues(true);

	AWCM_2D awcm(mra_grid, N);

	Q.TagName("rho");
	awcm.level1.TagName("level");
	mask.TagName("mask");

	Q.Print();
	awcm.mask1.Print();
	awcm.level1.Print();

	FIELD_UNIFORM_2D<VECTOR_4D<double>> Q_copy;

	while(true)
	{
		static int cnt = 1;
		static bool do_print = false;

		double max_velocity = 0;
		ITERATION_GRID_2D(grid, 0)
		{
			max_velocity = MAX(max_velocity, euler_eqn.Spectral_RadiusX(Q(i,j)) + euler_eqn.Spectral_RadiusY(Q(i,j)));
		}
		dt = cfl * dx / max_velocity;

		if(t+dt > arr_printing_times[cnt])
		{
			dt = arr_printing_times[cnt] - t;
			do_print = true;
		}

		awcm.PrepareAWCM(Q, mask, tol);

		Q_copy = Q;

		for(int iter = 1; iter <= 3; iter ++)
		{
			bc_manager.SetBoundaryCondition();

			if(iter == 1)
			{
				mask = awcm.mask3;
				level = awcm.level3;
			}
			else if(iter == 2)
			{
				mask = awcm.mask2;
				level = awcm.level2;
			}
			else
			{
				mask = awcm.mask1;
				level = awcm.level1;
			}
			ITERATION_GRID_2D(grid, 3)
			{
				if(mask(i,j) == true)
				{
					euler_eqn.ComputeFu(Q(i,j), fu(i,j));
					euler_eqn.ComputeGu(Q(i,j), gu(i,j));
				}
				else
				{
					fu(i,j).AssignAllValue(0);
					gu(i,j).AssignAllValue(0);
				}
			}

			ITERATION_GRID_2D(grid, 0)
			{
				alpha_x(j) = MAX(alpha_x(j), euler_eqn.Spectral_RadiusX(Q(i,j)));
				alpha_y(i) = MAX(alpha_y(i), euler_eqn.Spectral_RadiusY(Q(i,j)));
			}
			ITERATION_GRID_2D(grid, 0)
			{
				VECTOR_4D<double> flux_l, flux_r;
				if(mask(i,j) == true)
				{
					const int lev = level(i,j);
					if(lev < num_levels-1)
					{
						int res = mra_grid.grid[lev].dx;

						weno.ComputeFluxX(Q(grid.ClampIndexI(i-3*res),j), Q(grid.ClampIndexI(i-2*res),j), Q(grid.ClampIndexI(i-res),j), Q(i,j), 
							Q(grid.ClampIndexI(i+res),j), Q(grid.ClampIndexI(i+2*res),j), alpha_x(j), flux_l);
						weno.ComputeFluxX(Q(grid.ClampIndexI(i-2*res),j), Q(grid.ClampIndexI(i-res),j), Q(i,j), Q(grid.ClampIndexI(i+res),j),
							Q(grid.ClampIndexI(i+2*res),j), Q(grid.ClampIndexI(i+3*res),j), alpha_x(j), flux_r);
						dF(i,j) = (flux_r - flux_l) / (res*dx);

						weno.ComputeFluxY(Q(i,grid.ClampIndexJ(j-3*res)), Q(i,grid.ClampIndexJ(j-2*res)), Q(i,grid.ClampIndexJ(j-res)), Q(i,j),
							Q(i,grid.ClampIndexJ(j+res)), Q(i,grid.ClampIndexJ(j+2*res)), alpha_y(i), flux_l);
						weno.ComputeFluxY(Q(i,grid.ClampIndexJ(j-2*res)), Q(i,grid.ClampIndexJ(j-res)), Q(i,j), Q(i,grid.ClampIndexJ(j+res)),
							Q(i,grid.ClampIndexJ(j+2*res)), Q(i,grid.ClampIndexJ(j+3*res)), alpha_y(i), flux_r);
						dG(i,j) = (flux_r - flux_l) / (res*dy);

						// 					dF(i,j).AssignAllValue(0);
						// 					dG(i,j).AssignAllValue(0);
					}
					else
					{
						weno.ComputeFluxX(Q(i-3,j), Q(i-2,j), Q(i-1,j), Q(i,j), Q(i+1,j), Q(i+2,j), alpha_x(j), flux_l);
						weno.ComputeFluxX(Q(i-2,j), Q(i-1,j), Q(i,j), Q(i+1,j), Q(i+2,j), Q(i+3,j), alpha_x(j), flux_r);
						dF(i,j) = (flux_r - flux_l) / dx;

						weno.ComputeFluxY(Q(i,j-3), Q(i,j-2), Q(i,j-1), Q(i,j), Q(i,j+1), Q(i,j+2), alpha_y(i), flux_l);
						weno.ComputeFluxY(Q(i,j-2), Q(i,j-1), Q(i,j), Q(i,j+1), Q(i,j+2), Q(i,j+3), alpha_y(i), flux_r);
						dG(i,j) = (flux_r - flux_l) / dy;
					}
				}
				else
				{
					dF(i,j).AssignAllValue(0);
					dG(i,j).AssignAllValue(0);
				}

			}

			if(iter == 1)
			{
				ITERATION_GRID_2D(grid, 0)
				{
					if(mask(i,j) == true) Q(i,j) = Q(i,j) - dt * (dF(i,j) + dG(i,j));
					else Q(i,j).AssignAllValue(0);
				}
			}
			else if(iter == 2)
			{
				ITERATION_GRID_2D(grid, 0)
				{
					if(mask(i,j) == true)
					{
						Q(i,j) = 0.75*Q_copy(i,j) + 0.25*(Q(i,j) - dt * (dF(i,j) + dG(i,j)));
					}
					else
					{
						Q(i,j).AssignAllValue(0);
					}
				}
			}
			else
			{
				ITERATION_GRID_2D(grid, 0)
				{
					if(mask(i,j) == true)
					{
						Q(i,j) = 1.0/3.0 * Q_copy(i,j) + 2.0/3.0 * (Q(i,j) - dt * (dF(i,j) + dG(i,j)));
					}
					else
					{
						Q(i,j).AssignAllValue(0);
					}
				}
			}
		}


		t += dt;
		std::cout<<t<<std::endl;
		if(do_print == true)
		{
			Q.Print();
			cnt ++;
			do_print = false;
		}

//		Q.Print();
//		awcm.mask1.Print();
//		awcm.level1.Print();

		if(t >= 0.8) break;
	}
	awcm.level1.Print("level1");
	mask.Print("mask");

	std::ofstream fout;
	fout.open("Q");
	ITERATION_GRID_2D(grid, 0)
	{
		if(mask(i,j) == true)
		{
			double x = grid(i,j).x;
			double y = grid(i,j).y;
			fout<<x<<" "<<y<<" "<<Q(i,j)[0]<<std::endl;
		}
	}
	fout.close();

}