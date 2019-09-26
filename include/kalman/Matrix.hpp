// The MIT License (MIT)
//
// Copyright (c) 2015 Markus Herb
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#ifndef KALMAN_MATRIX_HPP_
#define KALMAN_MATRIX_HPP_

#include <cmath>

#include <Eigen/Dense>

#define DEFINE_CONST_GET(T, var, VAR_INDEX) T var() const { return this->internal_[VAR_INDEX]; }
#define DEFINE_REF_GET(T, var, VAR_INDEX) T& var() { return this->internal_[VAR_INDEX]; }

#define KALMAN_VECTOR(NAME, T, N) \
    private: \
        Kalman::Vector<T, N> internal_; \
    public: \
        typedef Eigen::VectorBlock<Kalman::Vector<T, N>> SegmentReturnType; \
        typedef const Eigen::VectorBlock<const Kalman::Vector<T, N>> ConstSegmentReturnType; \
        template<int Size> struct FixedSegmentReturnType { typedef Eigen::VectorBlock<Kalman::Vector<T, N>, Size> Type; }; \
        template<int Size> struct ConstFixedSegmentReturnType { typedef const Eigen::VectorBlock<const Kalman::Vector<T, N>, Size> Type; }; \
        typedef typename Kalman::Vector<T, N>::Scalar Scalar; \
        enum { \
            RowsAtCompileTime = Kalman::Vector<T, N>::RowsAtCompileTime, \
            ColsAtCompileTime = Kalman::Vector<T, N>::ColsAtCompileTime \
        }; \
        NAME(void) {} \
        NAME(const Kalman::Vector<T, N>& internal) : internal_(internal) {} \
        Kalman::Vector<T, N>& get() { return internal_; } \
        Kalman::Vector<T, N> get() const { return internal_; }

namespace Kalman {
    const int Dynamic = Eigen::Dynamic;

    /**
     * @class Kalman::Matrix
     * @brief Template type for matrices
     * @param T The numeric scalar type
     * @param rows The number of rows
     * @param cols The number of columns
     */
    template<typename T, int rows, int cols>
    using Matrix = Eigen::Matrix<T, rows, cols>;
    
    /**
     * @brief Template type for vectors
     * @param T The numeric scalar type
     * @param N The vector dimension
     */
    template<typename T, int N>
    using Vector = Matrix<T, N, 1>;
    
    /**
     * @brief Cholesky square root decomposition of a symmetric positive-definite matrix
     * @param _MatrixType The matrix type
     * @param _UpLo Square root form (Eigen::Lower or Eigen::Upper)
     */
    template<typename _MatrixType, int _UpLo = Eigen::Lower>
    class Cholesky : public Eigen::LLT< _MatrixType, _UpLo >
    {
    public:
        Cholesky() : Eigen::LLT< _MatrixType, _UpLo >() {}
        
        /**
         * @brief Construct cholesky square root decomposition from matrix
         * @param m The matrix to be decomposed
         */
        Cholesky(const _MatrixType& m ) : Eigen::LLT< _MatrixType, _UpLo >(m) {}
        
        /**
         * @brief Set decomposition to identity
         */
        Cholesky& setIdentity()
        {
            this->m_matrix.setIdentity();
            this->m_isInitialized = true;
            return *this;
        }
        
        /**
         * @brief Check whether the decomposed matrix is the identity matrix
         */
        bool isIdentity() const
        {
            eigen_assert(this->m_isInitialized && "LLT is not initialized.");
            return this->m_matrix.isIdentity();
        }
        
        /**
         * @brief Set lower triangular part of the decomposition
         * @param matrix The lower part stored in a full matrix
         */
        template<typename Derived>
        Cholesky& setL(const Eigen::MatrixBase <Derived>& matrix)
        {
            this->m_matrix = matrix.template triangularView<Eigen::Lower>();
            this->m_isInitialized = true;
            return *this;
        }
        
        /**
         * @brief Set upper triangular part of the decomposition
         * @param matrix The upper part stored in a full matrix
         */
        template<typename Derived>
        Cholesky& setU(const Eigen::MatrixBase <Derived>& matrix)
        {
            this->m_matrix = matrix.template triangularView<Eigen::Upper>().adjoint();
            this->m_isInitialized = true;
            return *this;
        }
    };
}

#endif
