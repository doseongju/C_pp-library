#pragma once

#include <assert.h>
#include <fstream>

template<class TT>
class ARRAY
{
public:
	int num_elements_; // number of all elements, i.e. last index is num_element - 1.
	TT *values_;

public:
	ARRAY(void)
		: num_elements_(0), values_(0)
	{};

	ARRAY(const int& num_elements_input)
		: num_elements_(0), values_(0)
	{
		Initialize(num_elements_input);
	}

	ARRAY(const int& num_elements_input, const TT& values_input)
		: num_elements_(0), values_(0)
	{
		Initialize(num_elements_input, values_input);
	}

	ARRAY(const ARRAY<TT>& arr) :num_elements_(0), values_(0)
	{
		Initialize(arr);
	}

	~ARRAY(void)
	{
		if(values_ != 0) delete [] values_;
		num_elements_ = 0;
	}

public:
	inline void Initialize(const int& num_elements_input)
	{
		num_elements_ = num_elements_input;

		if(values_ != 0) delete [] values_;
		values_ = new TT [num_elements_input];
	}

	void Initialize(const int& num_elements_input, const TT& values_input)
	{
		num_elements_ = num_elements_input;

		if(values_ != 0) delete [] values_;
		values_ = new TT [num_elements_];

		AssignAllValues(values_input);
	}

	void Initialize(const ARRAY<TT>& array_input)
	{
		Initialize(array_input.num_elements_);
		CopyFrom(array_input);
	}

	void AssignAllValues(const TT& constant)
	{
		for(int w=0; w < num_elements_; w++) values_[w] = constant;
	}

	void AssignValues(const int& start_ix, const int& end_ix, const TT& constant)
	{
		for(int w=start_ix; w <= end_ix; w++) values_[w] = constant;
	}

	void CopyFrom(const ARRAY<TT>& from)
	{
		assert(num_elements_ == from.num_elements_);

		TT *from_val = from.values_;

		for (int w = 0; w < num_elements_; w++) values_[w] = from_val[w];
	}

public:
	inline TT& operator[](const int& i)
	{
		return values_[i];
	}

	const inline TT& operator[](const int& i) const
	{
		return values_[i];
	}

	void operator = (const ARRAY<TT>& arr)
	{
		Initialize(arr);
	}

public:
	void Print(const char* filename) const
	{
		std::ofstream fout;
		fout.open(filename);
		for(int i = 0; i < num_elements_; i++)
		{
			fout<<values_[i]<<std::endl;
		}
		fout.close();
	}
};

