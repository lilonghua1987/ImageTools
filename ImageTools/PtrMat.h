#pragma once
#include <malloc.h>

class PtrMat
{
public:

	PtrMat(void)
		:row(0)
		,column(0)
		,channel(0)
		,data(nullptr)
	{
	}

	PtrMat(unsigned long row, unsigned long column, unsigned long channel = 1)
		:row(row)
		,column(column)
		,channel(channel)
	{
		creat();
	}

	PtrMat(const PtrMat& ptr)
	{
		this->row = ptr.row;
		this->column = ptr.column;
		this->channel = ptr.channel;
		destroy();
		creat();
		memcpy(this->data, ptr.data, sizeof(void *) * row * column * channel);
	}

	PtrMat& operator=(const PtrMat& ptr)
	{
		if (this == &ptr)
		{
			return *this;
		}

		this->row = ptr.row;
		this->column = ptr.column;
		this->channel = ptr.channel;
		destroy();
		creat();
		memcpy(this->data, ptr.data, sizeof(void *) * row * column * channel);
		return *this;
	}

	template<typename T>
	T* at(unsigned long row, unsigned long column , unsigned long channel = 0)
	{
		assert(row < this->row && column < this->column && channel < this->channel);
		return (T *)data+((this->column*row + column)*this->channel + channel);
	}

	~PtrMat(void)
	{
		destroy();
	}

private:
	void creat()
	{
		data = (void *)malloc(sizeof(void *)*row*column*channel);
	}
	void destroy()
	{
		if (data)
		{
			free(data);
			data = nullptr;
		}	
	}

public:
	unsigned long row;
	unsigned long column;
	unsigned long channel;
	void* data;
};

