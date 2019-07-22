#ifndef MARS_MATRIX_HPP
#define MARS_MATRIX_HPP 
#include "mars_tensor_base.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (July 2019)                                                                    ////////
////// We define the Matrix class:                                                                             ////////
////// Matrix<typename T, Integer Rows, Integer Cols>                                                          ////////
////// We want them to be eventually constexpr and static. A proper constructor is then needed                 ////////
////// To build it, we must define a TensorBase class from which Matrix inherit                                ////////
////// In this way we can use Matrix<T,Rows,Cols> as we would normally do, but we can also define:             ////////
////// static constexpr Matrix<T,Rows,Cols> static_mat{T1,T2....};                                             ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace mars {

	template<typename T, Integer Rows_, Integer Cols_>
	class Matrix: public TensorBase<T, std::make_index_sequence<Rows_*Cols_>> 
	{
	public:
		static constexpr Integer Rows=Rows_;
		static constexpr Integer Cols=Cols_;
		using type= Matrix<T,Rows,Cols>;
		using subtype=T;
		using MB = TensorBase<T, std::make_index_sequence<Rows*Cols>>;
		using MB::MB;
		using MB::values;

		inline constexpr static Integer rows() { return Rows; }
		inline constexpr static Integer cols() { return Cols; }

        // access matrix direclty by using I*Col+J index
		inline constexpr T &operator()(const Integer i)
		{
			assert(i < Rows*Cols);
			return values[i];
		}

		inline constexpr const T &operator()(const Integer i)const
		{
			assert(i < Rows*Cols);
			return values[i];
		}

		inline constexpr T &operator()(const Integer i, const Integer j)
		{
			assert(i < Rows);
			assert(j < Cols);
			return values[i*cols() + j];
		}

		inline constexpr const T &operator()(const Integer i, const Integer j) const
		{
			assert(i < Rows);
			assert(j < Cols);
			return values[i*cols() + j];
		}

		operator T &() { 
           // Throw here if size is not 1x1...
			static_assert(Rows_==1&&Cols_==1,"Matrix<T,Rows,Cols> returns T if and only if Rows=Cols=1");
			return (*this)( 0, 0 ); 
		}


		inline constexpr void row(const Integer r, const Vector<T, Cols> &v)
		{
			assert(r < Rows && " row index must be smaller than number of rows");

			for(Integer d = 0; d < Cols; ++d) {
				(*this)(r,d) = v(d);
			}
		}

		inline constexpr void row(const Integer r, const Matrix<T, 1, Cols> &v)
		{
			assert(r < Rows && " row index must be smaller than number of rows");

			for(Integer d = 0; d < Cols; ++d) {
				(*this)(r,d) = v(1,d);
			}
		}


		inline constexpr void col(const Integer c, const Vector<T, Rows> &v)
		{
			assert(c < Cols);

			for(Integer d = 0; d < Rows; ++d) {
				(*this)(d, c) = v(d);
			}
		}


		inline constexpr Vector<T, Cols> get_row(const Integer r) const
		{
			assert(r < Rows && " row index must be smaller than number of rows");
            Vector<T, Cols> v;
			for(Integer d = 0; d < Cols; ++d) {
				v(d)=(*this)(r,d);
			}

			return v;
		}


		inline constexpr void get_row(const Integer r, Vector<T, Cols> &v) const
		{
			assert(r < Rows && " row index must be smaller than number of rows");

			for(Integer d = 0; d < Cols; ++d) {
				v(d)=(*this)(r,d);
			}
		}

		inline constexpr void get_row(const Integer r, Matrix<T, 1,Cols> &v) const
		{
			assert(r < Rows && " row index must be smaller than number of rows");

			for(Integer d = 0; d < Cols; ++d) {
				v(1,d)=(*this)(r,d);
			}
		}


		inline constexpr void get_col(const Integer c, Vector<T, Rows> &v) const
		{
			assert(c < Cols);

			for(Integer d = 0; d < Rows; ++d) {
				v(d) = (*this)(d, c);
			}
		}

		inline constexpr void zero()
		{
			std::fill(begin(values), end(values), 0.);
		}


		void describe(std::ostream &os) const
		{
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					os << (*this)(i, j) << " ";
				}
				os << "\n";
			}

			os << "\n";
		}

		friend std::ostream &operator<<(std::ostream &os, const Matrix &m)
		{
			m.describe(os);
			return os;
		}





		inline constexpr Matrix<T, Rows,Cols>& operator = (const Matrix<T, Rows,Cols> &m) 
		{            
			if(this==&m) return *this;
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) = m(i,j);
					}
				}
			}
			return *this;
		} 


		inline constexpr Matrix<T, Rows,Cols>& operator = (const T &value) 
		{            
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) = value;
					}
				}
			}
			return *this;
		} 

		inline constexpr Matrix<T, Rows,Cols>& operator /= (const Real &alpha) 
		{            
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) /= alpha;
					}
				}
			}
			return *this;
		}       


		inline constexpr Matrix<T, Rows,Cols>& operator *= (const Real &alpha)
		{

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) *= alpha;
					}
				}
			}
			return *this;
		} 


		inline constexpr Matrix<T, Rows,Cols>& operator += (const Matrix<T, Rows,Cols> &mat)
		{        
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) += mat(i,j);
					}
				}
			}
			return *this;
		} 

		inline constexpr Matrix<T, Rows,Cols>& operator -= (const Matrix<T, Rows,Cols> &mat)
		{        
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						(*this)(i, j) -= mat(i,j);
					}
				}
			}
			return *this;
		} 

		inline constexpr Matrix<T, Rows, Cols> operator - ()const
		{
			Matrix<T, Rows, Cols> result;

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					result(i, j) = -(*this)(i, j);
				}
			}

			return result;
		}

		inline constexpr Matrix<T, Rows, Cols> operator + (const Matrix<T, Rows, Cols> &other)const
		{
			Matrix<T, Rows, Cols> ret;

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					ret(i, j) = (*this)(i, j) + other(i, j);
				}
			}

			return ret;
		}

		inline constexpr Matrix<T, Rows, Cols> operator - (const Matrix<T, Rows, Cols> &other)const
		{
			Matrix<T, Rows, Cols> ret;

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					ret(i, j) = (*this)(i, j) - other(i, j);
				}
			}

			return ret;
		}

	    template<Integer OtherCols>
		inline constexpr Matrix<T, Rows, OtherCols> operator * (const Matrix<T, Cols, OtherCols> &other) const
		{
			Matrix<T, Rows, OtherCols> ret;
			ret.zero();

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					for(Integer k = 0; k < OtherCols; ++k) {
						ret(i, k) += (*this)(i, j) * other(j, k);
					}
				}
			}
			return ret;
		}


		inline constexpr Vector<T, Rows> operator * (const Vector<T, Cols> &other) const
		{
			Vector<T, Rows> ret;
			ret.zero();

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						ret[i] += (*this)(i, j) * other(j);
					}
				}
			}

			return ret;
		}

		inline constexpr Matrix<T, Rows,Cols> operator * (const Real &alpha) const
		{
			Matrix<T, Rows, Cols> ret;
			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						ret(i,j) = (*this)(i, j) * alpha;
					}
				}
			}
			return ret;
		}

		inline constexpr Matrix<T, Rows,Cols> operator / (const Real &alpha) const
		{
			Matrix<T, Rows, Cols> ret;

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						ret(i,j) = (*this)(i, j) / alpha;
					}
				}
			}
			return ret;
		}

	};

	template<typename T, Integer N>
	inline constexpr void m_minor(
		const Integer cof_i,
		const Integer cof_j,
		const Matrix<T, N, N> &mat,
		Matrix<T, N-1, N-1> &m)
	{
		Integer i_offset = 0;
		for(Integer i = 0; i < N; ++i) {
			if(i == cof_i) {
				i_offset = -1;
				continue;
			}

			Integer j_offset = 0;

			for(Integer j = 0; j < N; ++j) {
				if(j == cof_j) {
					j_offset = -1;
					continue;
				}

				m(i + i_offset, j + j_offset) = mat(i, j);
			}
		}
	}

	template<typename T, Integer N>
	inline constexpr T det_aux(const Matrix<T, N, N> &m)
	{
		static_assert(N < 7, "max size is 6");

		std::array<Integer, N> nnz;
		std::fill(std::begin(nnz), std::end(nnz), 0);

		for(Integer i = 0; i < N; ++i) {
			for(Integer j = 0; j < N; ++j) {
				nnz[i] += m(i, j) != 0.;
			}
		}

		Integer row = 0;
		for(Integer i = 0; i < N; ++i) {
			if(nnz[row] > nnz[i]) {
				row = i;
			}
		}
		
		assert(row < N);

		if(nnz[row] == 0) return 0.;

		Integer ret = 0.;
		Matrix<T, N-1, N-1> mij;
		for(Integer j = 0; j < N; ++j) {
			const Real coff = m(row, j);
			if(coff == 0.) continue;

			const Real sign = (j % 2) == 0? 1 : -1;
			m_minor(row, j, m, mij);
			ret += sign * coff * det(mij);
		}

		return ret;
	}

	template<typename T, Integer N>
	inline constexpr T det(const Matrix<T, N, N> &m)
	{
		return det_aux(m);
	}

	template<typename T, Integer Rows, Integer Cols>
	inline constexpr T det(const Matrix<T, Rows, Cols> &m)
	{
		if(Rows > Cols) {
			Matrix<T, Rows, Rows> m_square;
			for(Integer i = 0; i < Rows; ++i) {

				for(Integer j = 0; j < Cols; ++j) {
					m_square(i, j) = m(i, j);
				}

				for(Integer j = Cols; j < Rows; ++j) {
					m_square(i, j) = i == j;
				}
			}

			return det(m_square);
		} else {
			Matrix<T, Cols, Cols> m_square;

			for(Integer j = 0; j < Cols; ++j) {
				for(Integer i = 0; i < Rows; ++i) {
					m_square(i, j) = m(i, j);
				}

				for(Integer i = Rows; i < Cols; ++i) {
					m_square(i, j) = i == j;
				}
			}

			return det(m_square);
		}
	}

	template<typename T, Integer Rows, Integer Cols>
	constexpr Matrix<T, Cols, Rows> transpose(const Matrix<T, Rows, Cols> &mat)
	{
		Matrix<T, Cols, Rows> ret;
		for(Integer i = 0; i < Rows; ++i) {
			for(Integer j = 0; j < Cols; ++j) {
				ret(j, i) = mat(i, j);
			}
		}

		return ret;
	}
	
	template<typename T>
	inline constexpr T det(const Matrix<T, 1, 1> &m)
	{
		return m(0, 0);
	}
	
	template<typename T>
	inline constexpr T det(const Matrix<T, 2, 2> &m)
	{
		return m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
	}
	
	template<typename T>
	inline constexpr T det(const Matrix<T, 3, 3> &m)
	{
		return
		m(0,0) * m(1,1) * m(2,2)  +
		m(0,1) * m(1,2) * m(2,0)  +
		m(0,2) * m(1,0) * m(2,1)  -
		m(0,0) * m(1,2) * m(2,1)  -
		m(0,1) * m(1,0) * m(2,2)  -
		m(0,2) * m(1,1) * m(2,0);
	}
	
	template<typename T>
	inline constexpr T det(const Matrix<T, 4, 4> &m)
	{
		const auto m00 = m(0, 0);
		const auto m01 = m(0, 1);
		const auto m02 = m(0, 2);
		const auto m03 = m(0, 3);

		const auto m10 = m(1, 0);
		const auto m11 = m(1, 1);
		const auto m12 = m(1, 2);
		const auto m13 = m(1, 3);

		const auto m20 = m(2, 0);
		const auto m21 = m(2, 1);
		const auto m22 = m(2, 2);
		const auto m23 = m(2, 3);

		const auto m30 = m(3, 0);
		const auto m31 = m(3, 1);
		const auto m32 = m(3, 2);
		const auto m33 = m(3, 3);

		return (m00 == 0. ? 0. : (
			m00 * det(
				Matrix<T, 3, 3>({
					m11, m12, m13,
					m21, m22, m23,
					m31, m32, m33

				}))))
		-
		(m01 == 0. ? 0. :
			m01 * det(
				Matrix<T, 3, 3>({
					m10, m12, m13,
					m20, m22, m23,
					m30, m32, m33
				})))
		+
		(m02 == 0. ? 0. : (
			m02 * det(
				Matrix<T, 3, 3>({
					m10, m11, m13,
					m20, m21, m23,
					m30, m31, m33

				}))))
		-
		(m03 == 0. ? 0. : (
			m03 * det(
				Matrix<T, 3, 3>({
					m10, m11, m12,
					m20, m21, m22,
					m30, m31, m32

				}))));
	}





        template<typename T,Integer Rows,Integer Cols>
	constexpr const Matrix<T,1,1> contraction(const Matrix<T, Rows,Cols> &A,const Matrix<T, Rows,Cols> &B)
	{
            // result must bbe initializable with zero
		Matrix<T,1,1> result=0;
		for(Integer i = 0; i < Rows; ++i) 
			for(Integer j = 0; j < Cols; ++j) 
				result(0,0) += A(i, j) * B(i,j);
			return result;
		};



  //       template<typename T>
		// const T contraction(const T &A,const Matrix<T,1,1> &B)
	 //    {
	 //    	return A*B(0,0);
	 //    }

		// const Real contraction(const Matrix<Real,1,1> &A, const Real &B)
	 //    {
	 //    	return A(0,0)*B;
	 //    }

        template<typename T,Integer Rows,Integer Cols>
		inline constexpr Matrix<T, Rows,Cols> operator * (const Real &alpha,const Matrix<T,Rows,Cols>& mat)
		{
			Matrix<T, Rows, Cols> ret;

			for(Integer i = 0; i < Rows; ++i) {
				for(Integer j = 0; j < Cols; ++j) {
					{
						ret(i,j) = mat(i, j) * alpha;
					}
				}
			}

			return ret;
		}

        template<typename T>
		inline constexpr Matrix<T, 1,1> inverse(const Matrix<T,1,1>& mat)
		{
			Matrix<T, 1,1> inv{1/mat(1,1)};
			return inv;
		}

        template<typename T>
		inline constexpr Matrix<T, 2,2> inverse(const Matrix<T,2,2>& mat)
		{
			Matrix<T, 2,2> inv{mat(1,1), -mat(0,1), -mat(1,0), mat(0,0)};
			inv/=det(mat);
			return inv;
		}

        template<typename T>
		inline constexpr Matrix<T, 3,3> inverse(const Matrix<T,3,3>& mat)
		{
        	Matrix<T,3,3> inv{+(mat(2,2)*mat(1,1)-mat(2,1)*mat(1,2)),  // a00
        		               -(mat(2,2)*mat(0,1)-mat(2,1)*mat(0,2)),  // a01
        		               +(mat(1,2)*mat(0,1)-mat(1,1)*mat(0,2)),  // a02
                               -(mat(2,2)*mat(1,0)-mat(2,0)*mat(1,2)),  // a10 
                               +(mat(2,2)*mat(0,0)-mat(2,0)*mat(0,2)),  // a11
                               -(mat(1,2)*mat(0,0)-mat(1,0)*mat(0,2)),  // a12
							   +(mat(2,1)*mat(1,0)-mat(2,0)*mat(1,1)),  // a20 
		     				   -(mat(2,1)*mat(0,0)-mat(2,0)*mat(0,1)),  // a21
		     				   +(mat(1,1)*mat(0,0)-mat(1,0)*mat(0,1))}; // a22

		     				   inv/=det(mat);
		     				   return inv;
		     				}


      template<typename T, Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
		     	constexpr void assign(Matrix<T, Rows1,Cols1>& mat1, const Matrix<T, Rows2,Cols2>& mat2, Integer II, Integer JJ)

		     				{
		     					static_assert(Rows1>=Rows2," only a smaller matrix can be assigned: Rows1>=Rows2");
		     					static_assert(Cols1>=Cols2," only a smaller matrix can be assigned: Cols1>=Cols2");
		     					assert(Rows2+II<=Rows1 && " Rows2+II<=Rows1");
		     					assert(Cols2+JJ<=Cols1 && " Cols2+JJ<=Cols1");
		     					for(Integer i=0; i<Rows2 ;i++)
		     						for(Integer j=0; j<Cols2 ;j++)
		     							mat1(i+II,j+JJ)=mat2(i,j);

		     					}


		     				}

#endif //MARS_MATRIX_HPP
