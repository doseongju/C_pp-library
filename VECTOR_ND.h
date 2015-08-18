#pragma once
#include <iostream>

template <class T>
class VECTOR_ND
{
public:
	int size;		
	T* val;

public:
	VECTOR_ND():size(0), val(0){}
	VECTOR_ND(const int& size_):val(0)
	{
		Initialize(size_);
	}
	VECTOR_ND(const VECTOR_ND<T>& vec):val(0)
	{
		Initialize(vec.size);
		for(int i = 0; i < size; i++)
			val[i] = vec[i];
	}
	~VECTOR_ND()
	{
		SAFE_DELETE_ARRAY(val);
	}

public:
	void Initialize(const int& size_)
	{
		size = size_;
		SAFE_DELETE_ARRAY(val);
		val = new T[size];
		AssignAllValue(T());
	}

	T& operator[](const int& i)
	{
		return val[i];
	}

	const T& operator[](const int& i) const
	{
		return val[i];
	}

	void operator +=(const VECTOR_ND<T>& vec)
	{
		assert(size == vec.size);
		for(int i = 0; i < size; i++) val[i] += vec.val[i];
	}

	void operator -=(const VECTOR_ND<T>& vec)
	{
		assert(size == vec.size);
		for(int i = 0; i < size; i++) val[i] -= vec.val[i];
	}

	void operator *=(const T& c)
	{
		for(int i = 0; i < size; i++) val[i] *= c;
	}

	void operator =(const VECTOR_ND<T>& vec)
	{
		Initialize(vec.size);
		for(int i = 0; i < size; i++) val[i] = vec.val[i];
	}

	double operator * (const VECTOR_ND<T>& vec) const
	{
		assert(size == vec.size);
		double result(0);
		for(int p = 0; p < size; p ++) result += val[p] * vec[p];
		return result;
	}

	double Magnitude() const
	{
		return sqrt(SqrMagnitude());
	}

	double SqrMagnitude() const
	{
		double result(0);
		for(int i = 0; i < size; i++) result += val[i] * val[i];
		return result;
	}

	void AssignAllValue(const T& value)
	{
		for(int i = 0; i < size; i++) val[i] = value;
	}

	void Normalize()
	{
		T norm = Magnitude();
		if (norm < 1e-16)
		{
			val[0] = 1;
			for (int i = 1; i < size; i++) val[i] = 0;

		}
		else
		{
			for (int i = 0; i < size; i++) val[i] /= norm;
		}
	}

};

template<class T>
std::ostream& operator<<(std::ostream& output, const VECTOR_ND<T>& vec)
{
	for(int i = 0; i < vec.size; i++)
	{
		output<< vec[i]<<" ";
	}
	return output;
}

template<class TT>
double DotProduct(const VECTOR_ND<TT>& v1, const VECTOR_ND<TT>& v2)
{
	assert(v1.size == v2.size);
	TT result(0);
	for (int i = 0; i < v1.size; i++)
		result += v1[i] * v2[i];
	return result;
}

template<class T>
VECTOR_ND<T> operator*(const T& c, const VECTOR_ND<T>& vec)
{
	VECTOR_ND<T> vec_(vec.size);
	for(int i = 0; i < vec.size; i++)
	{
		vec_[i] = c * vec[i];
	}
	return vec_;
}

template<class T>
VECTOR_ND<T> operator+(const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2)
{
	assert(v1.size == v2.size);
	VECTOR_ND<T> vec(v1.size);
	for(int i = 0; i < vec.size; i++)
	{
		vec[i] = v1[i] + v2[i];
	}
	return vec;
}

template<class T>
VECTOR_ND<T> operator-(const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2)
{
	assert(v1.size == v2.size);
	VECTOR_ND<T> vec(v1.size);
	for(int i = 0; i < vec.size; i++)
	{
		vec[i] = v1[i] - v2[i];
	}
	return vec;
}