#pragma once
//#include "WENO_1D.h"
//#include "WENO_2D.h"

template<class TT>
class EXTRAPOLATION
{
	EXTRAPOLATION(){}
	~EXTRAPOLATION(){}

public:

	static void Constant(const TT& from, TT& to)
	{
		to = from;
	}



	//Extrapolate value of 'to' using 'from1' and 'from2'.
	//
	//grid configuration.
	// ----------------------------
	// from1       from2        to

	static void Linear(const TT& from_1, const TT& from_2, TT& to)
	{
		to = 2*from_2 - from_1;
	}

};