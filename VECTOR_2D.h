#pragma once

#include <iostream>

template<class T>
class VECTOR_2D
{
public:
	union
	{
		struct{T i, j;};
		struct{T x, y;};
		T values[2];
	};

public: 
	VECTOR_2D():i(T()), j(T()){}
	VECTOR_2D(const T& x_, const T& y_):x(x_), y(y_){}
	VECTOR_2D(const VECTOR_2D<T>& vector_)
	{
		i = vector_.i;
		j = vector_.j;
	}
	~VECTOR_2D(){}
public:
	void operator = (const VECTOR_2D& vector_)
	{
		i = vector_.i;
		j = vector_.j;
	}
	T& operator[](const int& i_)
	{
		return values[i_];
	}

	const T& operator[](const int& i_) const
	{
		return values[i_];
	}


	void operator +=(const T& a){i += a, j += a;}
	void operator -=(const T& a){i -= a, j -= a;}
	void operator *=(const T& a){i *= a, j *= a;}
	void operator /=(const T& a){i /= a, j /= a;}

	void operator +=(const VECTOR_2D<T>& vector_){i += vector_.i; j += vector_.j;}
	void operator -=(const VECTOR_2D<T>& vector_){i -= vector_.i; j -= vector_.j;}
	void operator *=(const VECTOR_2D<T>& vector_){i *= vector_.i; j *= vector_.j;}
	void operator /=(const VECTOR_2D<T>& vector_){i /= vector_.i; j /= vector_.j;}

	VECTOR_2D operator + (const VECTOR_2D<T>& vector_) const
	{ return VECTOR_2D(i + vector_.i, j + vector_.j); }
	VECTOR_2D operator - (const VECTOR_2D<T>& vector_) const
	{ return VECTOR_2D(i - vector_.i, j - vector_.j); }
	VECTOR_2D operator * (const VECTOR_2D<T>& vector_) const
	{ return VECTOR_2D(i * vector_.i, j * vector_.j); }
	VECTOR_2D operator / (const VECTOR_2D<T>& vector_) const
	{ return VECTOR_2D(i / vector_.i, j / vector_.j); }

	double SquareMagnitude()
	{
		return x*x + y*y;
	}

	double Magnitude()
	{
		return sqrt(SquareMagnitude());
	}

};

template<class T>
std::ostream& operator << (std::ostream& output_, const VECTOR_2D<T>& vector_)
{
	return output_<<vector_.i<<" "<<vector_.j;
}

template <class T>
T DotProduct(const VECTOR_2D<T>& vec1, const VECTOR_2D<T>& vec2)
{
	return vec1.i * vec2.i + vec1.j * vec2.j;
}

template<class T>
VECTOR_2D<T> operator *(const T& c, const VECTOR_2D<T>& vec)
{
	return VECTOR_2D<T>(c*vec.i, c*vec.j);
}