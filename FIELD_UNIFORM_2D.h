#pragma once

#include "ARRAY_2D.h"
#include "GRID_UNIFORM_2D.h"
#include "EXTRAPOLATION.h"

template<class T>
class FIELD_UNIFORM_2D
{
public:
	ARRAY_2D<T> arr;
	GRID_UNIFORM_2D grid;
	int i_start, i_res, i_end,
		j_start, j_res, j_end,
		ij_res;

	int i_start_g, i_res_g, i_end_g,
		j_start_g, j_res_g, j_end_g,
		ij_res_g;

	int ghost_width;

	char* name;

public:
	FIELD_UNIFORM_2D():i_start(0), i_res(0), i_end(0), j_start(0), j_res(0), j_end(0), ij_res(0),
					   i_start_g(0), i_res_g(0), i_end_g(0), j_start_g(0), j_res_g(0), j_end_g(0), ij_res_g(0), name(0){}
	FIELD_UNIFORM_2D(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_,
		const double& x_min_, const double& y_min_, const double& x_max_, const double& y_max_, const int& ghost_width_ = 0) :name(0)
	{
		ghost_width = ghost_width_;
		Initialize(i_start_, j_start_, i_res_, j_res_, x_min_, y_min_, x_max_, y_max_, ghost_width);
	}
	FIELD_UNIFORM_2D(const GRID_UNIFORM_2D& grid_, const int& ghost_width_ = 0):name(0)
	{
		Initialize(grid_, ghost_width_);
	}

	~FIELD_UNIFORM_2D(){}

public:
	void Initialize(const int& i_start_, const int& j_start_, const int& i_res_, const int& j_res_,
					const double& x_min_, const double& y_min_, const double& x_max_, const double& y_max_, const int& ghost_width_ = 0)
	{
		ghost_width = ghost_width_;

		i_start = i_start_; i_res = i_res_; i_end = i_start + i_res - 1;
		j_start = j_start_; j_res = j_res_; j_end = j_start + j_res - 1;
		ij_res = i_res * j_res;

		i_start_g = i_start - ghost_width; i_res_g = i_res + 2*ghost_width; i_end_g = i_end + ghost_width;
		j_start_g = j_start - ghost_width; j_res_g = j_res + 2*ghost_width; j_end_g = j_end + ghost_width;
		ij_res_g = i_res_g * j_res_g;

		grid.Initialize(i_start, j_start, i_res, j_res, x_min_, y_min_, x_max_, y_max_);

		arr.Initialize(i_start_g, j_start_g, i_res_g, j_res_g, true);
	}

	void Initialize(const GRID_UNIFORM_2D& grid_, const int& ghost_width_ = 0)
	{
		ghost_width = ghost_width_;
		Initialize(grid_.i_start, grid_.j_start, grid_.i_res, grid_.j_res, grid_.x_min, grid_.y_min, grid_.x_max, grid_.y_max, ghost_width);
	}

	void AssignAllValues(const T& values)
	{
		arr.AssignAllValues(values);
	}

	void TagName(char* name_)
	{
		name = name_;
	}
public:
	T& operator()(const int& i, const int& j)
	{
		return arr(i,j);
	}
	const T& operator()(const int& i, const int& j) const
	{
		return arr(i,j);
	}


public:
	void Print(const char* filename, const bool print_ghost = false) const
	{
		std::ofstream fout;
		fout.open(filename);
		fout.precision(15);
		int i_start_(i_start), i_end_(i_end), j_start_(j_start), j_end_(j_end);
		if(print_ghost)
		{
			i_start_ = i_start_g, i_end_ = i_end_g, j_start_ = j_start_g, j_end_ = j_end_g;
		}

		for(int j = j_start_; j <= j_end_; j++)
		{
			for(int i = i_start_; i <= i_end_; i++)
			{
				fout<<(*this)(i,j)<<" ";
			}
			fout<<std::endl;
		}

		fout.close();
	}

	void Print(const bool& print_ghost = false) const
	{
		static int cnt = 0;
		char filename[30];
		if(name != 0) sprintf(filename, "%s_%4.4d" , name, cnt);
		else sprintf(filename, "%4.4d", cnt);
		Print(filename, print_ghost);
		cnt ++;
	}
};

template<class T, class TT>
void CheckFieldRange(const FIELD_UNIFORM_2D<T>& u, const FIELD_UNIFORM_2D<TT>& v)
{
	assert(u.i_start == v.i_start && u.i_end == v.i_end);
	assert(u.j_start == v.j_start && u.j_end == v.j_end);
}


template<class TT>
class FILL_GHOST_CELL_2D
{
	FILL_GHOST_CELL_2D(){}
	~FILL_GHOST_CELL_2D(){}

public:
	

	//Grid configuration
	//
	//   6  /      7        /  8
	//---------------------------
	//      /				/
	//		/				/
	//		/				/
	//	4	/				/  5
	//		/				/
	//		/				/
	//		/				/
	//---------------------------
	//  1   /		2		/  3
	//		/				/


	static void Constant(FIELD_UNIFORM_2D<TT>& field)
	{
		if(field.ghost_width == 0) return;

		int i_start, i_start_g,
			      j_start, j_start_g,
			      i_end  , i_end_g,
			      j_end  , j_end_g;
		i_start = field.i_start, i_start_g = field.i_start_g;
		j_start = field.j_start, j_start_g = field.j_start_g;
		i_end   = field.i_end  , i_end_g   = field.i_end_g  ;
		j_end   = field.j_end  , j_end_g   = field.j_end_g  ;

		//Region_1
		for(int j = j_start_g; j < j_start; j ++)
		{
			for(int i = i_start_g; i < i_start; i ++)
			{
				field(i,j) = field(i_start, j_start);
			}
		}

		//Region_2
		for(int j = j_start + 1; j >= j_start_g; j --)
		{
			for(int i = i_start; i <= i_end; i ++)
			{
				EXTRAPOLATION<TT>::Constant(field(i,j+1), field(i,j));
			}
		}

		//Region_3
		for(int j = j_start_g; j < j_start; j ++)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				field(i,j) = field(i_end, j_start);
			}
		}

		//Region_4
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_start-1; i >= i_start_g; i --)
			{
				EXTRAPOLATION<TT>::Constant(field(i+1,j), field(i,j));
			}
		}

		//Region_5
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				EXTRAPOLATION<TT>::Constant(field(i-1,j), field(i,j));
			}
		}

		//Region_6
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_start_g; i < i_start; i ++)
			{
				field(i,j) = field(i_start, j_end);
			}
		}

		//Region_7
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_start; i <= i_end	; i ++)
			{
				EXTRAPOLATION<TT>::Constant(field(i,j-1), field(i,j));
			}
		}

		//Region_8
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				field(i,j) = field(i_end, j_end);
			}
		}
	}


	static void Linear(FIELD_UNIFORM_2D<TT>& field)
	{
		if(field.ghost_width == 0) return;
		if(field.i_res < 3) return;

		int i_start, i_start_g,
			j_start, j_start_g,
			i_end  , i_end_g,
			j_end  , j_end_g;
		i_start = field.i_start, i_start_g = field.i_start_g;
		j_start = field.j_start, j_start_g = field.j_start_g;
		i_end   = field.i_end  , i_end_g   = field.i_end_g  ;
		j_end   = field.j_end  , j_end_g   = field.j_end_g  ;

		//Region_2
		for(int j = j_start - 1; j >= j_start_g; j --)
		{
			for(int i = i_start; i <= i_end; i ++)
			{
				EXTRAPOLATION<TT>::Linear(field(i,j+2), field(i,j+1), field(i,j));
			}
		}

		//Region_4
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_start - 1; i >= i_start_g; i --)
			{
				EXTRAPOLATION<TT>::Linear(field(i+2,j), field(i+1,j), field(i,j));
			}
		}

		//Region_5
		for(int j = j_start; j <= j_end; j ++)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				EXTRAPOLATION<TT>::Linear(field(i-2,j), field(i-1,j), field(i,j));
			}
		}

		//Region_7
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_start; i <= i_end	; i ++)
			{
				EXTRAPOLATION<TT>::Linear(field(i,j-2), field(i,j-1), field(i,j));
			}
		}



		//Region_1
		for(int j = j_start - 1; j >= j_start_g; j --)
		{ 
			for(int i = i_start - 1; i >= i_start_g; i --)
			{
				EXTRAPOLATION<TT>::Linear(field(i+2,j+2), field(i+1,j+1), field(i,j));
			}
		}

		//Region_3
		for(int j = j_start - 1; j >= j_start_g; j --)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				EXTRAPOLATION<TT>::Linear(field(i-2,j+2), field(i-1,j+1), field(i,j));
			}
		}

		//Region_6
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_start - 1; i >= i_start_g; i --)
			{
				EXTRAPOLATION<TT>::Linear(field(i+2,j-2), field(i+1,j-1), field(i,j));
			}
		}

		//Region_8
		for(int j = j_end + 1; j <= j_end_g; j ++)
		{
			for(int i = i_end + 1; i <= i_end_g; i ++)
			{
				EXTRAPOLATION<TT>::Linear(field(i-2,j-2), field(i-1,j-1), field(i,j));
			}
		}
	}
	
};