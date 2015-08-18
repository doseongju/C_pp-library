#pragma once

template<class T>
class VECTOR_8D
{
public:
	T values[8];

public:
	VECTOR_8D()
	{
		values[0] = T();
		values[1] = T();
		values[2] = T();
		values[3] = T();
		values[4] = T();
		values[5] = T();
		values[6] = T();
		values[7] = T();
	}
	VECTOR_8D(const T& y0, const T& y1, const T& y2, const T& y3, const T& y4, const T& y5, const T& y6, const T& y7)
	{
		values[0] = y0;
		values[1] = y1;
		values[2] = y2;
		values[3] = y3;
		values[4] = y4;
		values[5] = y5;
		values[6] = y6;
		values[7] = y7;
	}
	VECTOR_8D(const VECTOR_8D<T>& vec)
	{
		values[0] = vec.values[0];
		values[1] = vec.values[1];
		values[2] = vec.values[2];
		values[3] = vec.values[3];
		values[4] = vec.values[4];
		values[5] = vec.values[5];
		values[6] = vec.values[6];
		values[7] = vec.values[7];
	}
	~VECTOR_8D(){}

public:
	void operator = (const VECTOR_8D<T>& vec)
	{
		values[0] = vec.values[0];
		values[1] = vec.values[1];
		values[2] = vec.values[2];
		values[3] = vec.values[3];
		values[4] = vec.values[4];
		values[5] = vec.values[5];
		values[6] = vec.values[6];
		values[7] = vec.values[7];
	}
	void operator += (const VECTOR_8D<T>& vec)
	{
		values[0] += vec.values[0];
		values[1] += vec.values[1];
		values[2] += vec.values[2];
		values[3] += vec.values[3];
		values[4] += vec.values[4];
		values[5] += vec.values[5];
		values[6] += vec.values[6];
		values[7] += vec.values[7];
	}
	void operator -= (const VECTOR_8D<T>& vec)
	{
		values[0] -= vec.values[0];
		values[1] -= vec.values[1];
		values[2] -= vec.values[2];
		values[3] -= vec.values[3];
		values[4] -= vec.values[4];
		values[5] -= vec.values[5];
		values[6] -= vec.values[6];
		values[7] -= vec.values[7];
	}
	void operator *= (const T& c)
	{
		values[0] *= c;
		values[1] *= c;
		values[2] *= c;
		values[3] *= c;
		values[4] *= c;
		values[5] *= c;
		values[6] *= c;
		values[7] *= c;
	}
	VECTOR_8D<T> operator + (const VECTOR_8D<T>& vec) const
	{
		return VECTOR_8D<T> (values[0]+vec[0], values[1]+vec[1], values[2]+vec[2],
		values[3]+vec[3], values[4]+vec[4], values[5]+vec[5], values[6]+vec[6], values[7]+vec[7]);
	}
	VECTOR_8D<T> operator - (const VECTOR_8D<T>& vec) const
	{
		return VECTOR_8D<T> (values[0]-vec[0], values[1]-vec[1], values[2]-vec[2],
			values[3]-vec[3], values[4]-vec[4], values[5]-vec[5], values[6]-vec[6], values[7]-vec[7]);
	}
	VECTOR_8D<T> operator * (const T& c) const
	{
		return VECTOR_8D<T> (values[0]*c, values[1]*c, values[2]*c, values[3]*c,
		values[4]*c, values[5]*c, values[6]*c, values[7]*c);
	}

	T& operator [](const int& i)
	{
		return values[i];
	}
	const T& operator [](const int& i) const
	{
		return values[i];
	}
	
};

template<class T>
VECTOR_8D<T> operator * (const T& c, const VECTOR_8D<T>& vec)
{
	return VECTOR_8D<T> (c*vec[0], c*vec[1], c*vec[2], c*vec[3], c*vec[4], c*vec[5], c*vec[6], c*vec[7]);
}

template<class T>
std::ostream& operator << (std::ostream& output, const VECTOR_8D<T>& vec)
{
	return output<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" "<<vec[3]<<" "<<vec[4]<<" "
		<<vec[5]<<" "<<vec[6]<<" "<<vec[7];
}