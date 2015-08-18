#include "Wavelet/MRA_GRID_UNIFORM_1D.h"
#include "Wavelet/WAVELET_1D.h"
#include "Wavelet/AWCM_1D.h"
#include "WENO.h"
#include "STOP_WATCH.h"

using std::cout;
using std::endl;

void main()
{
	const int num_lev = 13;
	const int Nx_coarse = 9;
	MRA_GRID_INT_1D mra;
	mra.Initialize(num_lev, Nx_coarse);
	const double tol = 1e-6;
	const int max_lev = num_lev - 1;

	const int i_res = mra.grid_finest.i_res; 

	GRID_UNIFORM_1D grid(0, i_res, -0.5, 0.5);
	FIELD_UNIFORM_1D<double> q(grid, 3), q_0, q_1, q_2;
	FIELD_UNIFORM_1D<bool> mask(grid);

	grid.PrintGrid("x");
	ITERATION_GRID_1D(grid, 3)
	{
		double x = grid(i);
		//q(i) = exp(-3200*SQR(x));
		//q(i) = exp(-100*x*x);
		if(x <= 0) q(i) = 1;
		else q(i) = 2;
	}

	q.Print("q");

	AWCM_1D awcm(mra);
	

	const double dx = grid.dx;
	const double dt = 0.5 * dx;
	WENO5 weno;
	
	
	STOP_WATCH stop_watch;
	stop_watch.Tic();
	for(int iter = 1; iter <= 1000; iter ++)
	{
		awcm.Prepare(q, tol);

		q_2 = q;
		double dF = 0;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask2(i) == false)
			{
				continue;
			}

			const int lev = awcm.level2(i);
			const int dx_ = mra.grid[lev].dx;

			if(lev == max_lev)
			{
				double F_r, F_l;

				weno.ReconstructRight(q(i-3), q(i-2), q(i-1), q(i), q(i+1), F_l);
				weno.ReconstructRight(q(i-2), q(i-1), q(i), q(i+1), q(i+2), F_r);

				q_2(i) = q_2(i) - dt * (F_r - F_l) / (dx*dx_);
				
			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.Differentiate(q(i-2*dx_), q(i-dx_), q(i), q(i+dx_), q(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateAtLeftBoundary(q(i-dx_), q(i), q(i+dx_), q(i+2*dx_), dF);
				}
				else if(i == 2*dx_)
				{
					awcm.DifferentiateAtLeftBoundary(q(i-2*dx_), q(i-dx_), q(i), q(i+dx_), q(i+2*dx_), dF);
				}
				else if(i +2 *dx_ > grid.i_end)
				{
					awcm.DifferentiateAtLeftBoundary(q(i+dx_), q(i), q(i-dx_), q(i-2*dx_), dF);
				}
				else
				{
					awcm.DifferentiateAtLeftBoundary(q(i+2*dx_), q(i+dx_), q(i), q(i-dx_), q(i-2*dx_), dF);
				}
				q_2(i) = q_2(i) - dt / (dx*dx_) * dF;
			}
		}

		q_2(-1) = q_2(0);
		q_2(-2) = q_2(0);
		q_2(-3) = q_2(0);
		q_2(i_res) = q_2(i_res-1);
		q_2(i_res+1) = q_2(i_res-1);
		q_2(i_res+2) = q_2(i_res-1);
		q_1 = q_2;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask1(i) == false)
			{
				continue;
			}
			const int lev = awcm.level1(i);
			const int dx_ = mra.grid[lev].dx;

			if(lev == max_lev)
			{
				double F_r, F_l;

				weno.ReconstructRight(q_2(i-3), q_2(i-2), q_2(i-1), q_2(i), q_2(i+1), F_l);
				weno.ReconstructRight(q_2(i-2), q_2(i-1), q_2(i), q_2(i+1), q_2(i+2), F_r);

				q_1(i) = q_1(i) - dt * (F_r - F_l) / (dx*dx_);
			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.Differentiate(q_2(i-2*dx_), q_2(i-dx_), q_2(i), q_2(i+dx_), q_2(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateAtLeftBoundary(q_2(i-dx_), q_2(i), q_2(i+dx_), q_2(i+2*dx_), dF);
				}
				else if(i == 2*dx_)
				{
					awcm.DifferentiateAtLeftBoundary(q_2(i-2*dx_), q_2(i-dx_), q_2(i), q_2(i+dx_), q_2(i+2*dx_), dF);
				}
				else if(i +2 *dx_ > grid.i_end)
				{
					awcm.DifferentiateAtLeftBoundary(q_2(i+dx_), q_2(i), q_2(i-dx_), q_2(i-2*dx_), dF);
				}
				else
				{
					awcm.DifferentiateAtLeftBoundary(q_2(i+2*dx_), q_2(i+dx_), q_2(i), q_2(i-dx_), q_2(i-2*dx_), dF);
				}
				q_1(i) = q_1(i) - dt / (dx*dx_) * dF;
			}
		}

		q_1(-1) = q_1(0);
		q_1(-2) = q_1(0);
		q_1(-3) = q_1(0);
		q_1(i_res) = q_1(i_res-1);
		q_1(i_res+1) = q_1(i_res-1);
		q_1(i_res+2) = q_1(i_res-1);

		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask1(i) == false) continue;
			q_1(i) = 0.75 * q(i) + 0.25 * q_1(i);
		}

		
		q_0 = q_1;
		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask0(i) == false) 
			{
				continue;
			}
			const int lev = awcm.level0(i);
			const int dx_ = mra.grid[lev].dx;

			if(lev == max_lev)
			{
				double F_r, F_l;

				weno.ReconstructRight(q_1(i-3), q_1(i-2), q_1(i-1), q_1(i), q_1(i+1), F_l);
				weno.ReconstructRight(q_1(i-2), q_1(i-1), q_1(i), q_1(i+1), q_1(i+2), F_r);

				q_0(i) = q_0(i) - dt * (F_r - F_l) / (dx*dx_);
			}
			else
			{
				if(i - 2*dx_ > 0 && i + 2*dx_ < grid.i_end)
				{
					awcm.Differentiate(q_1(i-2*dx_), q_1(i-dx_), q_1(i), q_1(i+dx_), q_1(i+2*dx_), dF);
				}
				else if(i - 2*dx_ < 0)
				{
					awcm.DifferentiateAtLeftBoundary(q_1(i-dx_), q_1(i), q_1(i+dx_), q_1(i+2*dx_), dF);
				}
				else if(i == 2*dx_)
				{
					awcm.DifferentiateAtLeftBoundary(q_1(i-2*dx_), q_1(i-dx_), q_1(i), q_1(i+dx_), q_1(i+2*dx_), dF);
				}
				else if(i +2 *dx_ > grid.i_end)
				{
					awcm.DifferentiateAtLeftBoundary(q_1(i+dx_), q_1(i), q_1(i-dx_), q_1(i-2*dx_), dF);
				}
				else
				{
					awcm.DifferentiateAtLeftBoundary(q_1(i+2*dx_), q_1(i+dx_), q_1(i), q_1(i-dx_), q_1(i-2*dx_), dF);
				}
				q_0(i) = q_0(i) - dt / (dx*dx_) * dF;
			}
		}
		q_0(-1) = q_0(0);
		q_0(-2) = q_0(0);
		q_0(-3) = q_0(0);
		q_0(i_res) = q_0(i_res-1);
		q_0(i_res+1) = q_0(i_res-1);
		q_0(i_res+2) = q_0(i_res-1);

		ITERATION_GRID_1D(grid, 0)
		{
			if(awcm.mask0(i) == false) 
			{
				continue;
			}

			const int lev = awcm.level0(i);
			const int dx_ = mra.grid[lev].dx;

			q(i) = 1.0/3.0*q(i) + 2.0/3.0*q_0(i);
		}
		awcm.PostProcess(q);

		if(iter % 100 == 0)
		{
			char filename[20];
			sprintf_s(filename, "%4.4d", iter/100);
			AWCM_MANAGER_1D::Print(q, awcm.mask0, filename);
		}

	}

	std::cout<<stop_watch.Toc()<<std::endl;
	//awcm.Prepare(q, tol);

	awcm.wavelet.Interpolate(q, mra, awcm.mask0);
	q.Print("q_new");
	awcm.mask0.Print("mask");
	awcm.level0.Print("level");


	AWCM_MANAGER_1D::Print(q, awcm.mask0, "q_adaptive");



}