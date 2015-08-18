#pragma once
#include "FIELD_UNIFORM_1D.h"

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


template<class T>
class BOUNDARY_MANAGER_1D
{
public:
	int i_start;
	int i_start_g;
	int i_end;
	int i_end_g;
	int ghost_width;

	const FIELD_UNIFORM_1D<BOUNDARY_CONDITION>& bc;

public:
	BOUNDARY_MANAGER_1D(const FIELD_UNIFORM_1D<BOUNDARY_CONDITION>& bc_):bc(bc_)
	{
		Initialize(bc.i_start, bc.i_start_g, bc.i_end, bc.i_end_g, bc.ghost_width);
	}
	~BOUNDARY_MANAGER_1D(){}

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
	void SetLeftBoundaryCondition(FIELD_UNIFORM_1D<T>& field) const
	{
		int cnt = 0;
		for(int i = i_start-1; i >= i_start_g; i --) 
		{
			if(bc(i) == EXTRAPOLATION_ZERO_ORDER) field(i) =  field(i_start);
			if(bc(i) == PERIODIC                ) field(i) =  field(i_end + cnt);
			cnt --;
			//TODO : Set other boundary conditions
		}
	}

	void SetRightBoundaryCondition(FIELD_UNIFORM_1D<T>& field) const
	{
		int cnt = 0;
		for(int i = i_end+1; i <= i_end_g; i ++) 
		{
			if(bc(i) == EXTRAPOLATION_ZERO_ORDER) field(i) =  field(i_end);
			if(bc(i) == PERIODIC                ) field(i) =  field(i_start + cnt);
			cnt ++;
			//TODO : Set other boundary conditions
		}
	}

	void SetBoundaryCondition(FIELD_UNIFORM_1D<T>& field)
	{
		SetLeftBoundaryCondition (field);
		SetRightBoundaryCondition(field);
	}
};