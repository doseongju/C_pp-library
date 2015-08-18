#pragma once
#include <iostream>

template<class T>
class VECTOR_5D
{
public:
	T values[5];

public:
	VECTOR_5D()
	{
		values[0] = T();
		values[1] = T();
		values[2] = T();
		values[3] = T();
		values[4] = T();
	}
	VECTOR_5D(const T& x_, const T& y_, const T& z_, const T& w_, const T& v_)
	{
		values[0] = x_;
		values[1] = y_;
		values[2] = z_;
		values[3] = w_;
		values[4] = v_;
	}
	VECTOR_5D(const VECTOR_5D<T>& vec)
	{
		values[0] = vec[0], values[1] = vec[1], values[2] = vec[2], values[3] = vec[3], values[4] = vec[4];
	}
	~VECTOR_5D(){}

public:

	T& operator[](const int& i_)
	{
		return values[i_];
	}
	const T& operator[](const int& i_) const
	{
		return values[i_];
	}

	void AssignAllValue(const T& c)
	{
		values[0] = c;
		values[1] = c;
		values[2] = c;
		values[3] = c;
		values[4] = c;
	}


	void operator = (const VECTOR_5D<T>& vec)
	{
		values[0] = vec[0], values[1] = vec[1], values[2] = vec[2], values[3] = vec[3], values[4] = vec[4];
	}

	void operator +=(const VECTOR_5D<T>& vec)
	{ values[0] += vec[0], values[1] += vec[1], values[2] += vec[2], values[3] += vec[3], values[4] += vec[4]; }
	void operator -=(const VECTOR_5D<T>& vec)
	{ values[0] -= vec[0], values[1] -= vec[1], values[2] -= vec[2], values[3] -= vec[3], values[4] -= vec[4]; }
	void operator *=(const VECTOR_5D<T>& vec)
	{ values[0] *= vec[0], values[1] *= vec[1], values[2] *= vec[2], values[3] *= vec[3], values[4] *= vec[4]; }
	void operator *=(const T& c)
	{ values[0] *= c, values[1] *= c, values[2] *= c, values[3] *= c, values[4] *= c; }
	void operator /=(const VECTOR_5D<T>& vec)
	{ values[0] /= vec[0], values[1] /= vec[1], values[2] /= vec[2], values[3] /= vec[3], values[4] /= vec[4]; }

	VECTOR_5D<T> operator + (const T& c) const
	{
		return VECTOR_5D<T>(values[0]+c, values[1]+c, values[2]+c, values[3]+c, values[4]+c);
	}
	VECTOR_5D<T> operator - (const T& c) const
	{
		return VECTOR_5D<T>(values[0]-c, values[1]-c, values[2]-c, values[3]-c, values[4]-c);
	}
	VECTOR_5D<T> operator * (const T& c) const
	{
		return VECTOR_5D<T>(values[0]*c, values[1]*c, values[2]*c, values[3]*c, values[4]*c);
	}
	VECTOR_5D<T> operator / (const T& c) const
	{
		return VECTOR_5D<T>(values[0]/c, values[1]/c, values[2]/c, values[3]/c, values[4]/c);
	}

	VECTOR_5D<T> operator + (const VECTOR_5D<T>& vec) const
	{
		return VECTOR_5D<T>(values[0]+vec[0], values[1]+vec[1], values[2]+vec[2], values[3]+vec[3], values[4]+vec[4]);
	}
	VECTOR_5D<T> operator - (const VECTOR_5D<T>& vec) const
	{
		return VECTOR_5D<T>(values[0]-vec[0], values[1]-vec[1], values[2]-vec[2], values[3]-vec[3], values[4]-vec[4]);
	}
	T operator * (const VECTOR_5D<T>& vec) const
	{
		return T(values[0] * vec[0] + values[1] * vec[1] + values[2] * vec[2] + values[3] * vec[3] + values[4] * vec[4]);
	}


	double SquareMagnitude() const
	{
		return values[0]*values[0] + values[1]*values[1] + values[2]*values[2] + values[3]*values[3] + values[4]*values[4];
	}
	double Magnitude() const
	{
		return sqrt(values[0]*values[0] + values[1]*values[1] + values[2]*values[2] + values[3]*values[3] + values[4]*values[4]);
	}
};


template<class T>
std::ostream& operator<<(std::ostream& output, const VECTOR_5D<T>& vec)
{
	return output<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" "<<vec[3]<<" "<<vec[4];
}

template<class T>
T DotProduct(const VECTOR_5D<T>& v1, const VECTOR_5D<T>& v2)
{
	return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3] + v1[4]*v2[4]);
}

template<class T>
VECTOR_5D<T> operator * (const T& c, const VECTOR_5D<T>& vec)
{
	return VECTOR_5D<T>(c*vec[0], c*vec[1], c*vec[2], c*vec[3], c*vec[4]);
}