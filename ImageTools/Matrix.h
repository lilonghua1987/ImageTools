#pragma once

#include <math.h>
#include <assert.h>
#include <exception>
#include <algorithm>
#include <limits>

#ifndef ulong
typedef unsigned long ulong;
#endif

template <typename T>
class Matrix
{
public:
	Matrix(void)
		:row(0)
		,column(0)
		,channel(channel)
		,data(nullptr)
	{

	}

	Matrix(unsigned long row, unsigned long column, unsigned long channel=1)
		:row(row)
		,column(column)
		,channel(channel)
	{
		data = create(row,column,channel);
	}

	Matrix(const Matrix& matrix)
	{
		row = matrix.row;
		column = matrix.column;
		channel = matrix.channel;
		data = create(row,column,channel);
		copyMenmory(&data, matrix.data);
    };

	~Matrix(void)
	{
		destory();
	}

private:
	T* create(unsigned long row, unsigned long column, unsigned long channel=1)
	{
		assert(channel >= 1);
		assert(typeid(void) != typeid(T));
		
		return new T[row*column*channel]();
	}
	void copyMenmory(T **dst, const T *src)
	{
		if (*dst) 
			destory();
		*dst = create(row,column,channel);
		memcpy(*dst, src, sizeof(T) * row * column * channel);
	}

public:
	unsigned long row;
	unsigned long column;
	unsigned long channel;
	T* data;

public:
	friend std::ofstream& operator << (std::ofstream& out, const Matrix& mat)
	{
		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (unsigned long c = 0; c < mat.channel; c++)
				{
					out << mat.at(i,j,c) << " ";
				}
				if (mat.channel > 1) out << ";";
			}
			out << "\n";
		}
		out << "\n";
		return out;
	}

	Matrix& operator = (const Matrix& matrix)
	{
		if (this == &matrix)
		{
			return *this;
		}
		row = matrix.row;
		column = matrix.column;
		channel = matrix.channel;
		copyMenmory(&data, matrix.data);
		return *this;
	}

	Matrix operator * (Matrix matrix)const
	{
		assert(this->channel == matrix.channel);
		assert(this->channel == 1);
		assert(this->column == matrix.row);

		Matrix result(row,matrix.column);

		for (ulong i = 0; i < row; i++)
		{
			ulong m = 0;
			while (m < matrix.column)
			{
				T t = NULL;
				for (ulong j = 0; j < column; j++)
				{					
					t+=get(i,j)*matrix.get(j,m);				
				}
				result.set(t,i,m);
				m++;
			}			
		}
		return result;
	}

	Matrix operator * (const T& t)const
	{
		assert(std::_Is_numeric<T>::value);

		Matrix mat(row, column, channel);
		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					mat.at(i, j, c) = at(i, j, c) * t;
				}
			}
		}

		return mat;
	}

	friend Matrix& operator/=(Matrix& mat, const T& t)
	{
		assert(std::_Is_numeric<T>::value);
		assert(t != 0);
		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					mat.at(i,j,c)/= t;
				}
			}
		}
		return mat;
	}

	Matrix operator-(const Matrix& mat)const
	{
		assert(row == mat.row && channel == mat.channel && column == mat.column);

		Matrix result(*this);

		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					result.at(i, j, c) = result.at(i, j, c) - mat.at(i, j, c);
				}
			}
		}

		return result;
	}

	Matrix operator+(const Matrix& mat)const
	{
		assert(row == mat.row && channel == mat.channel && column == mat.column);

		Matrix result(*this);

		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					result.at(i, j, c) = result.at(i, j, c) + mat.at(i, j, c);
				}
			}
		}

		return result;
	}

	Matrix dotMultip(const Matrix& mat)const
	{
		assert(row == mat.row && channel == mat.channel && column == mat.column);

		Matrix result(row, column, channel);

		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					result.at(i, j, c) = mat.at(i, j, c) * at(i, j, c);
				}
			}
		}

		return result;
	}

	Matrix dotDivsion(const Matrix& mat)const
	{
		assert(row == mat.row && channel == mat.channel && column == mat.column);

		Matrix result(row, column, channel);

		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					if (mat.at(i, j, c) == 0) throw(std::exception("dotDivsion: The matrix has zero !"));
					else result.at(i, j, c) = at(i, j, c) / mat.at(i, j, c);
				}
			}
		}

		return result;
	}

	Matrix equal(const T& t)const
	{
		assert(std::_Is_numeric<T>::value);

		Matrix mat(row, column, channel);
		for (ulong i = 0; i < mat.row; i++)
		{
			for (ulong j = 0; j < mat.column; j++)
			{
				for (ulong c = 0; c < mat.channel; c++)
				{
					if (at(i, j, c) == t) mat.at(i, j, c) = 1;
				}
			}
		}

		return mat;
	}

	T& get(unsigned long row, unsigned long column, unsigned long channel=0)const
	{
		assert(row>= 0);
		assert(column>= 0);
		assert(this->row > row);
		assert(this->column > column);
		assert(channel < this->channel && channel >= 0);

		return data[row*this->column*this->channel + column*this->channel + channel];
	}

	T* ptr(unsigned long row, unsigned long column, unsigned long channel=0)const
	{
		assert(row>= 0);
		assert(column>= 0);
		assert(this->row > row);
		assert(this->column > column);
		assert(channel < this->channel && channel >= 0);

		return &(data[row*this->column*this->channel + column*this->channel + channel]);
	}

	template<typename T2>
	void set(T2 t, unsigned long row, unsigned long column, unsigned long channel=0)
	{
		assert(row>= 0);
		assert(column>= 0);
		assert(this->row > row);
		assert(this->column > column);
		assert(channel < this->channel && channel >= 0);

		data[row*this->column*this->channel + column*this->channel + channel] = (T)t;
    };

	T& at(unsigned long row, unsigned long column, unsigned long channel=0)const
	{
		assert(channel < this->channel && channel >= 0);
		return this->data[row*this->column*this->channel + column*this->channel + channel];
	};

	static Matrix ones(ulong row, ulong column, ulong channel = 1)
	{
		assert(row > 0 && column > 0 && channel > 0);

		Matrix mat(row,column,channel);

		for (ulong r = 0; r < row; r++)
		{
			for (ulong c = 0; c < column; c++)
			{
				for (ulong i = 0; i < channel; i++)
				{
					mat.at(r, c, i) = 1;
				}
			}
		}

		return mat;
	}


	//create identity matrix
	static Matrix eye(const unsigned long m)
	{
		Matrix M(m,m);
		for (unsigned long i = 0; i < m; i++)
			M.at(i,i) = 1;
		return M;
	};
	void eye()
	{
		assert(row == column && channel == 1);
		for (unsigned long i = 0; i < row; i++)
			for (unsigned long j = 0; j < column; j++)
				at(i,j) = 0;
		for (unsigned long i = 0; i < std::min(row,column); i++)
			at(i,i) = 1;
	};

	Matrix Trans()
	{
		Matrix m(column,row,channel);

		for (unsigned long i = 0; i < row; i++)
		{
			for (unsigned long j = 0; j < column; j++)
			{
				for (unsigned long c = 0; c < channel; c++)
				{
					m.at(j,i,c) = get(i,j,c);
				}
			}
		}
		return m;
	}

	//forward diff and the first column elements are zereo (left to right)
	static void DiffHorizontal(const Matrix& src, Matrix& dest)
	{
		if (src.channel != dest.channel || src.row != dest.row || src.column != dest.column)
			throw("The Horizontal diff is unequal with original Matrix !");

		for (unsigned long i = 0; i < src.row; i++)
		{
			for (unsigned long j = 1; j < src.column; j++)
			{
				for (unsigned long c = 0; c < src.channel; c++)
				{
					dest.at(i,j,c) = src.at(i,j,c) - src.at(i,j-1,c);
				}
			}
		}
	}


	//The first row elements are zereo (top to bottom)
	static void DiffVertical(const Matrix& src, Matrix& dest)
	{
		if (src.channel != dest.channel || src.row != dest.row || src.column != dest.column)
			throw("The Vertical diff is unequal with original Matrix !");

		for (unsigned long i = 1; i < src.row; i++)
		{
			for (unsigned long j = 0; j < src.column; j++)
			{
				for (unsigned long c = 0; c < src.channel; c++)
				{
					dest.at(i,j,c) = src.at(i,j,c) - src.at(i-1,j,c);
				}
			}
		}
	}

	//forward sum (left to right)
	static void SumHorizontal(const Matrix& src, Matrix& dest)
	{
		if (src.channel != dest.channel || src.row != dest.row || src.column != dest.column)
			throw("The Horizontal sum is unequal with original Matrix !");

		for (unsigned long i = 0; i < src.row; i++)
		{
			for (unsigned long j = 0; j < src.column; j++)
			{
				for (unsigned long c = 0; c < src.channel; c++)
				{
					if (j == 0)
						dest.at(i,j,c) = src.at(i,j,c);
					else
					    dest.at(i,j,c) = src.at(i,j,c) + src.at(i,j-1,c);
				}
			}
		}
	}


	//The first row elements are zereo (top to bottom)
	static void SumVertical(const Matrix& src, Matrix& dest)
	{
		if (src.channel != dest.channel || src.row != dest.row || src.column != dest.column)
			throw("The Vertical sum is unequal with original Matrix !");

		for (unsigned long i = 0; i < src.row; i++)
		{
			for (unsigned long j = 0; j < src.column; j++)
			{
				for (unsigned long c = 0; c < src.channel; c++)
				{
					if (1 == 0)
						dest.at(i,j,c) = src.at(i,j,c);
					else
						dest.at(i,j,c) = src.at(i,j,c) + src.at(i-1,j,c);
				}
			}
		}
	}

	Matrix& norm(T minValue, T maxValue)
	{
		if(channel != 1) throw("The matrix channel must equal with 1 !");
		if(minValue >= maxValue) throw("The maxValue must bigger than minValue !");

		T minT = at(0,0),maxT = at(0,0);

		for (unsigned long i = 0; i < row; i++)
		{
			for (unsigned long j = 0; j < column; j++)
			{
				if(at(i,j) < minT) minT = at(i,j); 
				if(at(i,j) > maxT) maxT = at(i,j); 
			}
		}

		if(maxT <= minT) return *this;

		for (unsigned long i = 0; i < row; i++)
		{
			for (unsigned long j = 0; j < column; j++)
			{
				at(i,j) = (at(i,j) - minT) * (maxValue - minValue)/(maxT - minT) + minValue;
			}
		}

		return *this;
	}

	Matrix& abs()
	{
		for (unsigned long i = 0; i < row; i++)
		{
			for (unsigned long j = 0; j < column; j++)
			{
				for (unsigned long c = 0; c < channel; c++)
				{
					at(i,j,c) = static_cast<T>(std::abs(get(i,j,c)));
				}
			}
		}
		return *this;
	}

	T maxValue()
	{
		if(channel != 1) throw("The matrix channel must equal with 1 !");

		T maxV = (std::numeric_limits<T>::min)();

		for (ulong v = 0; v < row; v++)
		{
			for (ulong u = 0; u < column; u++)
			{
				if (at(v,u) > maxV) maxV = at(v,u);
			}
		}

		return maxV;
	}

	T minValue()
	{
		if(channel != 1) throw("The matrix channel must equal with 1 !");

		T minV = (std::numeric_limits<T>::max)();

		for (ulong v = 0; v < row; v++)
		{
			for (ulong u = 0; u < column; u++)
			{
				if (at(v,u) < minV) minV = at(v,u);
			}
		}

		return minV;
	}

	Matrix cumulativeRow()const
	{
		if (row < 1 || column < 1 || channel < 1) throw(std::exception("cumulativeRow: the matrix row or column value must be greater than 1 !"));

		Matrix out(*this);

		for (ulong v = 1; v < row; v++)
		{
			for (ulong u = 0; u < column; u++)
			{
				for (ulong c = 0; c < channel; c++)
				{
					out.at(v, u, c) += out.at(v - 1, u, c);
				}
			}
		}

		return out;
	}

	Matrix cumulativeCol()const
	{
		if (row < 1 || column < 1 || channel < 1) throw(std::exception("cumulativeCol: the matrix row or column value must be greater than 1 !"));

		Matrix out(*this);

		for (ulong u = 1; u < column; u++)
		{
			for (ulong v = 0; v < row; v++)
			{
				for (ulong c = 0; c < channel; c++)
				{
					out.at(v, u, c) += out.at(v, u - 1, c);
				}
			}
		}

		return out;
	}

	void destory()
	{
		delete[] data;
		data = nullptr;
	}

	static Matrix getChannel(const Matrix& mat, ulong index)
	{
		assert(index >= 0 && index < mat.channel);

		Matrix result(mat.row,mat.column);

		for (ulong v = 0; v < mat.row; v++)
		{
			for (ulong u = 0; u < mat.column; u++)
			{
				result.at(v, u) = mat.at(v,u,index);
			}
		}

		return result;
	}

	static void deleteRow(Matrix* matrix, long index)
	{
		for (long i = 0; i < matrix->row; i++)
		{
			if (i == index)
			{
				delete matrix->data[i];
				matrix->row--;
				break;
			}		
		}
	}

	static void deleteColumn(Matrix* matrix, long index)
	{
		Matrix m = matrix->Trans();
		for (long i = 0; i < m.row; i++)
		{
			if (i == index)
			{
				delete m.data[i];
				m.row--;
				break;
			}	
		}
		Matrix temp = m.Trans();
		matrix->row = temp.row;
		matrix->column = temp.column;
		matrix->data = temp.data;
	}

	Matrix& deleteRow (long index)
	{
		if (index < row)
		{
			delete data[index];
			row--;
		}

		return *this;
	}

	Matrix& deleteColumn (long index)
	{
		Matrix& m = this->Trans();
		if (index < m.row)
		{
			delete m.data[index];
			m.row--;
		}
		Matrix temp = m.Trans();
		row = temp.row;
		column = temp.column;
		data = temp.data;
		return *this;
	}

	// invert matrix M
	static Matrix inv(const Matrix &M)
	{
		assert(M.row == M.column && M.channel == 1);
		Matrix A(M);
		Matrix B = eye(M.row);
		B.solve(A);
		return B;
	};

	// invert this matrix
	bool inv()
	{
		assert(row == column && channel == 1);
		Matrix A(*this);
		eye();
		return solve(A);
	};

	static Matrix inverse(const Matrix& mat)
	{
		assert(std::_Is_numeric<T>::value);
		assert(mat.channel == 1 && mat.row == mat.column && mat.row >= 1);

		Matrix out(mat);
		Matrix<ulong> is(1,mat.row);
		Matrix<ulong> js(1, mat.row);

		T d, p;
		for (ulong k = 0; k < mat.row; k++)
		{
			d = 0.0;
			for (ulong i = k; i < mat.row; i++)
			{
				for (ulong j = k; j < mat.row; j++)
				{
					p = std::abs(out.at(i,j));
					if (p > d) 
					{ 
						d = p; 
						is.at(0, k) = i;
						js.at(0, k) = j;
					}
				}
			}

			if (0.0 == d) throw(std::exception("inverse: the matrix is singular matrix !"));

			if (is.at(0,k) != k)
				for (ulong j = 0; j < mat.row; j++)
				{
					std::swap(out.at(k, j), out.at(is.at(0, k), j));
				}

			if (js.at(0, k) != k)
				for (ulong i = 0; i < mat.row; i++)
				{
					std::swap(out.at(i, k), out.at(i, js.at(0, k)));
				}

			if (0.0 == out.at(k, k)) throw(std::exception("inverse: the matrix is singular matrix !"));
			else out.at(k, k) = 1.0 / out.at(k, k);

			for (ulong j = 0; j < mat.row; j++)
				if (j != k)
				{
					out.at(k, j) *= out.at(k, k);
				}

			for (ulong i = 0; i < mat.row; i++)
				if (i != k)
					for (ulong j = 0; j < mat.row; j++)
						if (j != k)
						{
							out.at(i, j) -= out.at(i, k) * out.at(k, j);
						}
			for (ulong i = 0; i < mat.row; i++)
				if (i != k)
				{
					out.at(i, k) = -out.at(i, k) * out.at(k, k);
				}
		}
		for (int k = (int)mat.row - 1; k >= 0; k--)
		{
			if (js.at(0, k) != k)
				for (ulong j = 0; j < mat.row; j++)
				{
					std::swap(out.at(k, j), out.at(js.at(0, k), j));
				}

			if (is.at(0, k) != k)
				for (ulong i = 0; i < mat.row; i++)
				{
					std::swap(out.at(i, k), out.at(i, is.at(0, k)));
				}
		}
		return out;
	}


	template <typename T2>
	static Matrix castType(const Matrix<T2>& srcMat)
	{
		assert(srcMat.column > 0 && srcMat.row > 0 && srcMat.channel > 0);
		assert(std::_Is_numeric<T>::value);
		assert(std::_Is_numeric<T2>::value);

		Matrix out(srcMat.row,srcMat.column,srcMat.channel);

		for (ulong i = 0; i < srcMat.row; i++)
		{
			for (ulong j = 0; j < srcMat.column; j++)
			{
				for (ulong c = 0; c < srcMat.channel; c++)
				{
					out.at(i, j, c) = static_cast<T>(srcMat.at(i,j,c));
				}
			}
		}
		return out;
	};

	// solve linear system M*x=B, replaces *this and M
	bool solve(const Matrix &mat,double eps=1e-20)
	{
		// substitutes
		const Matrix &A = mat;
		Matrix &B = *this;

		assert(A.row == A.column && A.column == B.row && A.row >= 1 && B.column >= 1);

		// index vectors for bookkeeping on the pivoting
		Matrix<unsigned long> indxc(row,1),indxr(row,1),ipiv(row,1);

		// loop variables
		unsigned long icol, irow;

		// main loop over the columns to be reduced
		for (unsigned long i = 0; i < row; i++)
		{
			T big=0;

			// search for a pivot element
			for (unsigned long j = 0; j < row; j++)
			{
				if (ipiv.at(j,0) != 1)
				{
					for (unsigned long k = 0; k < row; k++)
					{
						if (ipiv.at(k,0) == 0)
						{
							if (std::abs(A.at(j,k)) >= big) 
							{
								big = std::abs(A.at(j,k));
								irow=j;
								icol=k;
							}
							ipiv.at(icol,0)++;

							// We now have the pivot element, so we interchange rows, if needed, to put the pivot
							// element on the diagonal. The columns are not physically interchanged, only relabeled.
							if (irow != icol)
							{
								for (unsigned long l = 0; l < row; l++) std::swap(A.at(irow,l), A.at(icol,l));
								for (unsigned long l = 0; l < column; l++) std::swap(B.at(irow,l), B.at(icol,l));
							}

							indxr.at(i,0) = irow; // We are now ready to divide the pivot row by the
							indxc.at(i,0) = icol; // pivot element, located at irow and icol.

							// check for singularity
							if (std::abs(A.at(icol,icol)) < (T)eps) return false;

							T pivinv = (T)(1.f/A.at(icol,icol));
							A.at(icol,icol) = 1;
							for (unsigned long l = 0; l < row; l++) A.at(icol,l) *= pivinv;
							for (unsigned long l = 0; l < column; l++) B.at(icol,l) *= pivinv;

							// Next, we reduce the rows except for the pivot one
							for (unsigned long ll = 0; ll < row; ll++)
							{
								if (ll != icol) 
								{
									T dum = A.at(ll,icol);
									A.at(ll,icol) = 0;
									for (unsigned long l = 0; l < row; l++) A.at(ll,l) -= A.at(icol,l)*dum;
									for (unsigned long l = 0; l < column; l++) B.at(ll,l)  -= B.at(icol,l)*dum;
								}
							}
						}
					}
				}
			}
		}

		// This is the end of the main loop over columns of the reduction. It only remains to unscramble
		// the solution in view of the column interchanges. We do this by interchanging pairs of
		// columns in the reverse order that the permutation was built up.
		for (int l = (int)row - 1; l >= 0; l--) 
		{
			if (indxr.at(l,0) != indxc.at(l,0))
			{
				for (unsigned long k = 0; k< row; k++) std::swap(A.at(k,indxr.at(l,0)), A.at(k,indxc.at(l,0)));
			}
		}

		return true;
	};
};
