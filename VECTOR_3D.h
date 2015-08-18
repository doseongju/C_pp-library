#pragma once
#include <iostream>

template<class T>
class VECTOR_3D
{
public:
	union
	{
		struct{T i, j, k;};
		struct{T x, y, z;};
		T values[3];
	};

public:
	VECTOR_3D():i(0), j(0), k(0){}
	VECTOR_3D(const T& x_, const T& y_, const T& z_):x(x_), y(y_), z(z_){}
	VECTOR_3D(const VECTOR_3D<T>& vec)
	{
		i = vec.i;
		j = vec.j;
		k = vec.k;
	}

	~VECTOR_3D(){}

public:
	void operator=(const VECTOR_3D<T>& vec)
	{
		i = vec.i;
		j = vec.j;
		k = vec.k;
	}

	T& operator[](const int& i_)
	{
		return values[i_];
	}

	const T& operator[](const int& i_) const
	{
		return values[i_];
	}


	void operator +=(const T& a){i += a, j += a; k += a;}
	void operator -=(const T& a){i -= a, j -= a; k -= a;}
	void operator *=(const T& a){i *= a, j *= a; k *= a;}
	void operator /=(const T& a){i /= a, j /= a; k /= a;}

	void operator +=(const VECTOR_3D<T>& vector_){i += vector_.i; j += vector_.j; k += vector_.k;}
	void operator -=(const VECTOR_3D<T>& vector_){i -= vector_.i; j -= vector_.j; k -= vector_.k;}
	void operator *=(const VECTOR_3D<T>& vector_){i *= vector_.i; j *= vector_.j; k *= vector_.k;}
	void operator /=(const VECTOR_3D<T>& vector_){i /= vector_.i; j /= vector_.j; k /= vector_.k;}

	VECTOR_3D<T> operator + (const VECTOR_3D<T>& vector_) const
	{ return VECTOR_3D<T>(i + vector_.i, j + vector_.j, k + vector_.k); }
	VECTOR_3D<T> operator - (const VECTOR_3D<T>& vector_) const
	{ return VECTOR_3D<T>(i - vector_.i, j - vector_.j, k - vector_.k); }
	VECTOR_3D<T> operator * (const VECTOR_3D<T>& vector_) const
	{ return VECTOR_3D<T>(values[0]*vector_[0], values[1]*vector_[1], values[2]*vector_[2]); }

	VECTOR_3D<T> operator + (const T& c) const
	{ return VECTOR_3D<T>(i + c, j + c, k + c); }
	VECTOR_3D<T> operator - (const T& c) const
	{ return VECTOR_3D<T>(i - c, j - c, k - c); }
	VECTOR_3D<T> operator * (const T& c) const
	{ return VECTOR_3D<T>(i * c, j * c, k * c); }
	VECTOR_3D<T> operator / (const T& c) const
	{ return VECTOR_3D<T>(i / c, j / c, k / c); }

	double SquareMagnitude()
	{
		return x*x + y*y + z*z;
	}

	double Magnitude()
	{
		return sqrt(SquareMagnitude());
	}

	void Normalize()
	{
		T magnitude = Magnitude();
		i /= magnitude;
		j /= magnitude;
		k /= magnitude;
	}
};


template<class T>
std::ostream& operator << (std::ostream& output_, const VECTOR_3D<T>& vector_)
{
	return output_<<vector_.i<<" "<<vector_.j<<" "<<vector_.k;
}

template<class T>
T DotProduct(const VECTOR_3D<T>& v1, const VECTOR_3D<T>& v2)
{
	return(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template<class T>
VECTOR_3D<T> operator* (const T& c, const VECTOR_3D<T>& vec)
{
	return VECTOR_3D<T>(c*vec[0], c*vec[1], c*vec[2]);
}

template<class T>
VECTOR_3D<T> CrossProduct(const VECTOR_3D<T>& v1, const VECTOR_3D<T>& v2)
{
	VECTOR_3D<T> result;

	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y =-v1.x * v2.z + v1.z * v2.x;
	result.z = v1.x * v2.y - v1.y * v2.x;
	
	return result;

}