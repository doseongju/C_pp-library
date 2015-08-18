#pragma once
#include "MACROS.h"
#include "ARRAY_1D.h"
#include "GRID_UNIFORM_1D.h"
#include "EXTRAPOLATION.h"
#include "VECTOR_ND.h"

template <class T>
class FIELD_UNIFORM_1D
{
public:
	ARRAY_1D<T> arr;
	GRID_UNIFORM_1D base_grid;
	int i_start, i_res, i_end;
	int i_start_g, i_res_g, i_end_g;
	int ghost_width;

	char* name;

public:
	FIELD_UNIFORM_1D():name(0){}
	FIELD_UNIFORM_1D(const int& i_start_, const int& i_res_, const double& x_min, const double& x_max, const int& ghost_width_ = 0):name(0)
	{
		Initialize(i_start_, i_res_, x_min, x_max, ghost_width_);
	}
	FIELD_UNIFORM_1D(const GRID_UNIFORM_1D& grid_, const int& ghost_width_ = 0):name(0)
	{
		ghost_width = ghost_width_;
		Initialize(grid_, ghost_width);
	}
	FIELD_UNIFORM_1D(const FIELD_UNIFORM_1D<T>& field, const int& ghost_width_ = 0):name(0)
	{
		ghost_width = ghost_width_;
		Initialize(field.base_grid, ghost_width);
		ITERATION_GRID_1D(field.base_grid, 0)
		{
			arr(i) = field(i);
		}
	}
	~FIELD_UNIFORM_1D(){}

public:
	void Initialize(const int& i_start_, const int& i_res_, const double& x_min, const double& x_max, const int& ghost_width_ = 0)
	{
		ghost_width = ghost_width_;

		i_start = i_start_;
		i_res   = i_res_;
		i_end   = i_start + i_res - 1;

		i_start_g = i_start - ghost_width;
		i_res_g = i_res + 2 * ghost_width;
		i_end_g = i_end + ghost_width;

		arr.Initialize(i_start_g, i_res_g, true);

		base_grid.Initialize(i_start, i_res, x_min, x_max);
	}

	void Initialize(const GRID_UNIFORM_1D& grid_, const int& ghost_width_ = 0)
	{
		Initialize(grid_.i_start, grid_.i_res, grid_.x_min, grid_.x_max, ghost_width_);
	}

	void AssignAllValues(const T& values)
	{
		arr.AssignAllValue(values);
	}

	void Print(const char* filename, const bool ghost_cell = false) const
	{
		std::ofstream fout(filename);
		int start(0), end(0);
		if (ghost_cell == false)
		{
			start = i_start; end = i_end;
		}
		else
		{
			start = i_start_g; end = i_end_g;
		}
		for (int i = start; i <= end; i++)
		{
			fout << arr(i) << " ";
		}
		fout.close();
	}

	void Print(const bool& ghost_cell = false) const
	{
		static int cnt = 0;
		char filename[30];
		if(name != 0) sprintf(filename, "%s%4.4d", name, cnt);
		else sprintf(filename, "%4.4d", cnt);
		Print(filename, ghost_cell);
		cnt ++;
	}

	void TagName(char* name_)
	{
		name = name_;
	}

public:
	T& operator()(const int& i)
	{
		return arr(i);
	}
	const T& operator()(const int& i) const
	{
		return arr(i);
	}
};



template<class TT>
class FILL_GHOST_CELL_1D
{
public:
	FILL_GHOST_CELL_1D(){}
	~FILL_GHOST_CELL_1D(){}

public:
	static void Constant(FIELD_UNIFORM_1D<TT>& field)
	{
		const int i_start = field.i_start,
			      i_end = field.i_end,
				  i_start_g = field.i_start_g,
				  i_end_g = field.i_end_g;

		if(field.ghost_width == 0) return;
		for(int i =  i_start - 1; i >= i_start_g; i --)
		{
			EXTRAPOLATION<TT>::Constant(field(i_start), field(i));
		}
		
		for(int i = i_end + 1; i <= i_end_g; i ++)
		{
			EXTRAPOLATION<TT>::Constant(field(i_end), field(i));
		}
	}

	static void Linear(FIELD_UNIFORM_1D<TT>& field)
	{
		const int i_start = field.i_start,
			i_end = field.i_end,
			i_start_g = field.i_start_g,
			i_end_g = field.i_end_g;

		if(field.ghost_width == 0) return;
		for(int i =  i_start - 1; i >= i_start_g; i --)
		{
			EXTRAPOLATION<TT>::Linear(field(i+2), field(i+1), field(i));
		}

		for(int i = i_end + 1; i <= i_end_g; i ++)
		{
			EXTRAPOLATION<TT>::Linear(field(i-2), field(i-1), field(i));
		}
	}
};