#pragma once
#include "VECTOR_ND.h"
#include "VECTOR_4D.h"
#include "FIELD_UNIFORM_2D.h"

#ifndef BOUNDARY_CONDITION_DSJ
enum BOUNDARY_CONDITION
{
	SOLID_WALL = -4,
	PERIODIC,
	EXTRAPOLATION_ZERO_ORDER,
	SPECIFIED,
	FLUID
};
#endif
#define BOUNDARY_CONDITION_DSJ



//This class is not used now.
//TODO : Adjust this class to adapt more general situations.

template<class TT>
class BOUNDARY_MANAGER_2D
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

	const FIELD_UNIFORM_2D<BOUNDARY_CONDITION>& bc;

public:
	BOUNDARY_MANAGER_2D(const FIELD_UNIFORM_2D<BOUNDARY_CONDITION>& bc_):bc(bc_)
	{
		Initialize(bc.i_start, bc.i_start_g, bc.i_end, bc.i_end_g, bc.j_start, bc.j_start_g, bc.j_end, bc.j_end_g, bc.ghost_width);
	}
	~BOUNDARY_MANAGER_2D(){}

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

	void SetLeftBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		for(int j = j_start; j <= j_end; j ++)
		{
			int cnt = 0;
			for(int i = i_start-1; i >= i_start_g; i --)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_start,     j);
				if(bc(i,j) == PERIODIC                ) field(i,j) = field(i_end + cnt, j);
				cnt --;
				//TODO : Set other boundary conditions
			}
		}
	}

	void SetRightBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		for(int j = j_start; j <= j_end; j ++)
		{
			int cnt = 0;
			for(int i = i_end+1; i <= i_end_g; i ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_end        , j);
				if(bc(i,j) == PERIODIC                ) field(i,j) = field(i_start + cnt, j);
				cnt ++;
				//TODO : Set other boundary conditions
			}
		}
	}

	void SetBottomBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		
		for(int i = i_start; i <= i_end; i ++)
		{
			int cnt = 0;
			for(int j = j_start-1; j >= j_start_g; j --)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i,j_start    );
				if(bc(i,j) == PERIODIC                ) field(i,j) = field(i,j_end + cnt);
				cnt --;
			}
		}
		
	}

	void SetTopBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		
		for(int i = i_start; i <= i_end; i ++)
		{
			int cnt = 0;
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i,j_end        );
				if(bc(i,j) == PERIODIC                ) field(i,j) = field(i,j_start + cnt);
				cnt ++;
			}
		}
		
	}

	void SetDiagonalBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		
		for(int i = i_start_g; i < i_start; i ++)
		{
			for(int j = j_start_g; j < j_start; j ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_start,j_start);
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_start,j_end);
			}
		}
		for(int i = i_end+1; i <= i_end_g; i ++)
		{
			for(int j = j_start_g; j < j_start; j ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_end,j_start);
			}
			for(int j = j_end+1; j <= j_end_g; j ++)
			{
				if(bc(i,j) == EXTRAPOLATION_ZERO_ORDER) field(i,j) = field(i_end,j_end);
			}
		}
	}

	void SetBoundaryCondition(FIELD_UNIFORM_2D<TT>& field) const
	{
		SetLeftBoundaryCondition(field);
		SetRightBoundaryCondition(field);
		SetBottomBoundaryCondition(field);
		SetTopBoundaryCondition(field);
		SetDiagonalBoundaryCondition(field);
	}
};