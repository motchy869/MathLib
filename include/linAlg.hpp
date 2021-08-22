/**
 * @author motchy (motchy869[at]gmail.com)
 * @brief Linear Algebra library
 */
#ifndef __MATH_LIB__LINALG_HPP__
#define __MATH_LIB__LINALG_HPP__

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <numeric>
#include "common.hpp"
#include "analysis.hpp"

namespace MathLib {
    namespace LinAlg {
        namespace Debug { /* For only debug use. Not intended to be used in field environment. */
            /**
             * @brief Prints real matrix "A" in standard output.
             *
             * @tparam T the number type of the elements of "A"
             * @param[in] m the number of the rows in "A"
             * @param[in] n the number of the columns in "A"
             * @param[in] A the matrix "A"
             * @param[in] format format string which is applied to each element in "A"
             */
            template <typename T>
            void printRealMat(size_t m, size_t n, const T *A, const char *format) {
                assert(m*n >= 1);
                const T *ptr_A = A;
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        printf(format, *ptr_A);
                        if (j < n-1) {
                            printf(", ");
                        }
                        ++ptr_A;
                    }
                    puts("");
                }
            }

            /**
             * @brief Prints real vector in standard output.
             *
             * @tparam T the number type of the elements of the input vector
             * @param[in] m the length of the input vector
             * @param[in] vec the input vector
             * @param[in] format format string which is applied to each element in the input vector
             */
            template <typename T>
            void printRealVec(size_t m, const T *vec, const char *format) {
                assert(m >= 1);
                printRealMat(1, m, vec, format);
            }

            /**
             * @brief Prints complex matrix "A" in standard output.
             *
             * @tparam T the number type of the real and imaginary part of the elements of "A"
             * @param[in] m the number of the rows in "A"
             * @param[in] n the number of the columns in "A"
             * @param[in] A the matrix "A"
             * @param[in] format format string which is applied to each element in "A"
             */
            template <typename T>
            void printComplexMat(size_t m, size_t n, const std::complex<T> *A, const char *format) {
                assert(m*n >= 1);
                const std::complex<T> *ptr_A = A;
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        printf(format, (*ptr_A).real(), (*ptr_A).imag());
                        if (j < n-1) {
                            printf(", ");
                        }
                        ++ptr_A;
                    }
                    puts("");
                }
            }

            /**
             * @brief Prints complex vector in standard output.
             *
             * @tparam T the number type of the real and imaginary part of the elements of the input vector
             * @param[in] m the length of the input vector
             * @param[in] vec the input vector
             * @param[in] format format string which is applied to each element in the input vector
             */
            template <typename T>
            void printComplexVec(size_t m, const std::complex<T> *vec, const char *format) {
                assert(m >= 1);
                printComplexMat(1, m, vec, format);
            }
        }

        /**
         * @brief Checks whether given two matrices "A" and "B" are equal or not.
         * @details "A" and "B" are considered to be equal if and only if the absolute value of all the element of "A-B" is strictly less than a given number, epsilon.
         *
         * @tparam T1 the number type of the elements of A and B
         * @tparam T2 the number type of the elements of epsilon
         * @param[in] m the number of the rows in the input matrices
         * @param[in] n the number of the columns in the input matrices
         * @param[in] A the matrix A
         * @param[in] B the matrix B
         * @param[in] epsilon the threshold epsilon
         * @retval true A and B are equal.
         * @retval false A and B are not equal.
         */
        template <typename T1, typename T2>
        bool isEqualMat(const size_t m, const size_t n, const T1 *const A, const T1 *const B, const T2 epsilon) {
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                if (std::abs(A[i] - B[i]) > epsilon) {return false;}
            }
            return true;
        }

        /**
         * @brief Checks whether given two vectors "a" and "b" are equal or not.
         * @details "a" and "b" are considered to be equal if and only if the absolute value of all the element of "a-b" is strictly less than a given number, epsilon.
         *
         * @tparam T1 the number type of the elements of "a" and "b"
         * @tparam T2 the number type of the elements of epsilon
         * @param[in] m the length of the input vectors
         * @param[in] a the vector "a"
         * @param[in] b the vector "b"
         * @param[in] epsilon the threshold epsilon
         * @retval true "a" and "b" are equal.
         * @retval false "a" and "b" are not equal.
         */
        template <typename T1, typename T2>
        bool isEqualVec(const size_t m, const T1 *const a, const T1 *const b, const T2 epsilon) {
            return isEqualMat(1, m, a, b, epsilon);
        }

        /**
         * @brief Construct complex matrix "A" from real part "A_real" and imaginary part "A_imag".
         *
         * @tparam T the number type of "A_real" and "A_imag"
         * @param[in] m the number of the rows in the input matrices
         * @param[in] n the number of the columns in the input matrices
         * @param[in] A_real the real part of "A"
         * @param[in] A_imag the imaginary part of "A"
         * @param[out] A output buffer for "A"
         */
        template <typename T>
        void complexMat(const size_t m, const size_t n, const T *const A_real, const T *const A_imag, std::complex<T> *const A) {
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                A[i].real(A_real[i]); A[i].imag(A_imag[i]);
            }
        }

        /**
         * @brief Construct complex vector "x" from real part "x_real" and imaginary part "x_imag".
         *
         * @tparam T the number type of "x_real" and "x_imag"
         * @param[in] m the length of the vector "x"
         * @param[in] x_real the real part of "x"
         * @param[in] x_imag the imaginary part of "x"
         * @param[out] x output buffer for "x"
         */
        template <typename T>
        void complexVec(const size_t m, const T *const x_real, const T *const x_imag, std::complex<T> *const x) {
            complexMat(m, 1, x_real, x_imag, x);
        }

        /**
         * @brief Set a given vector "d" to the diagonal entries of a given square matrix "A".
         * The length of "d" must be equal to the number of the rows of "A".
         *
         * @tparam T the number type of the elements of "d, A".
         * @param[in] m the number of the rows in the matrices "A"
         * @param[in] d "d"
         * @param[out] A "A"
         */
        template <typename T>
        void setDiag(const size_t m, const T *const d, T *const A) {
            for (size_t r=0; r<m; ++r) {A[r*m+r] = d[r];}
        }

        /**
         * @brief Add a given vector "d" to the diagonal entries of a given square matrix "A".
         *
         * @tparam T the number type of the elements of "d, A".
         * @param[in] m the number of the rows in the matrices "A"
         * @param[in] d "d"
         * @param[out] A "A"
         */
        template <typename T>
        void addDiag(const size_t m, const T *const d, T *const A) {
            for (size_t r=0; r<m; ++r) {A[r*m+r] += d[r];}
        }

        /**
         * @brief Drop contiguous rows and columns from a given "m"-by-"n" matrix "A" and store the result to "B".
         * The "r1, r1+1, ..., r2"-th rows and "c1, c1+1, ..., c2"-th columns are dropped, where "0<=r1<=r2<=m-1, 0<=c1<=c2<=n-1".
         * Parameter check for "r1, r2, c1, c2" is performed ONLY under bug hunting mode.
         *
         * @tparam T the number type of the elements of "A".
         * @param[in] m the number of the rows in the matrices "A"
         * @param[in] n the number of the columns in the matrices "A"
         * @param[in] r1 "r1"
         * @param[in] r2 "r2"
         * @param[in] c1 "c1"
         * @param[in] c2 "c2"
         * @param[in] A "A"
         * @param[out] B the output buffer, "B"
         */
        template <typename T>
        void dropSubMat(const size_t m, const size_t n, const size_t r1, const size_t r2, const size_t c1, const size_t c2, const T *const A, T *const B) {
            #if MATH_LIB_ENABLE_BUG_HUNTING_MODE
                if (!(0 <= r1 && r1 <= r2 && r2 < m && 0 <= c1 && c1 <= c2 && c2 < n)) {
                    std::cerr << "BUG, FILE: " << __FILE__ << ", LINE: " << __LINE__ << std::endl;
                    exit(EXIT_FAILURE);
                }
            #endif
            const T *A_ptr = A; T *B_ptr = B;
            const size_t rowGap = r2-r1+1, colGap = c2-c1+1;
            for (size_t r=0; r<r1; ++r) {
                for (size_t c=0; c<c1; ++c) {*B_ptr++ = *A_ptr++;} // upper left part
                A_ptr += colGap;
                for (size_t c=c2+1; c<n; ++c) {*B_ptr++ = *A_ptr++;} // upper right part
            }
            A_ptr += rowGap*n;
            for (size_t r=r2+1; r<m; ++r) {
                for (size_t c=0; c<c1; ++c) {*B_ptr++ = *A_ptr++;} // lower left part
                A_ptr += colGap;
                for (size_t c=c2+1; c<n; ++c) {*B_ptr++ = *A_ptr++;} // lower right part
            }
        }

        /**
         * @brief Calculates transpose of a matrix "A" as "B".
         *
         * @tparam T the number type of the elements of input matrix
         * @param[in] m the number of the rows in the input matrix "A"
         * @param[in] n the number of the columns in the input matrix "A"
         * @param[in] A the input matrix
         * @param[out] B the output matrix
         */
        template <typename T>
        void transposeMat(size_t m, size_t n, const T *A, T *B) {
            for (size_t r=0; r<m; ++r) {
                const T *ptr_A = &A[r*n];
                T *ptr_B = B + r;
                for (size_t c=0; c<n; ++c) {
                    ptr_B[c*m] = ptr_A[c];
                }
            }
        }

        /**
         * @brief Calculates conjugate of a matrix "A" as "B".
         *
         * @tparam T the number type of real and imaginary part of the elements of the matrices "A" and "B"'.
         * @param[in] m the number of the rows in the matrices "A" and "B"
         * @param[in] n the number of the columns in the matrices "A" and "B"
         * @param[in] A the matrix "A"
         * @param[out] B the matrix "B"
         */
        template <typename T>
        void conjugateMat(const size_t m, const size_t n, const std::complex<T> *const A, std::complex<T> *const B) {
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                B[i] = std::conj(A[i]);
            }
        }

        /**
         * @brief Calculates conjugate of a vector "x" as "y".
         *
         * @tparam T the number type of real and imaginary part of the elements of the vectors "x" and "y".
         * @param[in] m the length of the vector "x"
         * @param[in] x the vector "x"
         * @param[out] y the vector "y"
         */
        template <typename T>
        void conjugateVec(const size_t m, const std::complex<T> *const x, std::complex<T> *const y) {
            conjugateMat(m, 1, x, y);
        }

        /**
         * @brief Calculates the sum of given two matrices "A" and "B". The result is stored in "A".
         *
         * @tparam T the number type of the elements of the matrices "A" and "B"
         * @param[in] m the number of the rows in the matrices "A" and "B"
         * @param[in] n the number of the columns in the matrices "A" and "B"
         * @param[inout] A the matrix "A"
         * @param[in] B the matrix "B"
         */
        template <typename T>
        void addMat_inplace(const size_t m, const size_t n, T *const A, const T *const B) {
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {A[i] += B[i];}
        }

        /**
         * @brief Calculates the sum of given two vector "x" and "y". The result is stored in "x".
         *
         * @tparam T the number type of the elements of the vectors "x" and "y"
         * @param[in] m the length of the vector "x" and "y"
         * @param[inout] x  the vector "x"
         * @param[in] y the vector "y"
         */
        template <typename T>
        void addVec_inplace(const size_t m, T *const x, const T *const y) {
            addMat_inplace(m, 1, x, y);
        }

        /**
         * @brief Calculates a scaled matrix "cA" as "B" where "A" is a matrix and "c" is a scalar.
         *
         * @tparam T the number type of "c" and the elements of "A"
         * @param[in] m the number of the rows in the input matrix "A"
         * @param[in] n the number of the columns in the input matrix "A"
         * @param[in] c the scalar "c"
         * @param[in] A the matrix "A"
         * @param[out] B the matrix "B"
         */
        template <typename T>
        void scaleMat(const size_t m, const size_t n, const T c, const T *const A, T *const B) {
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {B[i] = c*A[i];}
        }

        /**
         * @brief Multiply each row of a given "m x n" matrix "A", and store result to "B"
         * Multiply c[i] to the i-th row of "A" where "c" is a vector with length "m".
         *
         * @tparam Tc the number type of the elements of "c"
         * @tparam TA the number type of the elements of "A"
         * @tparam TB the number type of the elements of "B"
         * @param[in] m the number of the rows in the matrix "A"
         * @param[in] n the number of the columns in the matrix "A"
         * @param[in] c the vector "m"
         * @param[in] A the matrix "A"
         * @param[out] B the output matrix "B"
         */
        template <typename Tc, typename TA, typename TB>
        void scaleMatEachRow(const size_t m, const size_t n, const Tc *const c, const TA *const A, TB *const B) {
            for (size_t i=0; i<m; ++i) {
                const TA *const ptr_A = &A[i*n];
                TB *const ptr_B = &B[i*n];
                const Tc ci = c[i];
                for (size_t j=0; j<n; ++j) {
                    ptr_B[j] = ci*ptr_A[j];
                }
            }
        }

        /**
         * @brief Calculates a scaled vector "ca" as "b" where "a" is a vector and "c" is a scalar.
         *
         * @tparam T the number type of "c" and the elements of "a"
         * @param[in] m the length of the vector "a"
         * @param[in] c the scalar "c"
         * @param[in] a the vector "a"
         * @param[out] b the vector "b"
         */
        template <typename T>
        void scaleVec(const size_t m, const T c, const T *const a, T *const b) {
            scaleMat(1, m, c, a, b);
        }

        /**
         * @brief Calculate a self outer product of a given real vector "x".
         *
         * @tparam T the number type of the entries of "x"
         * @param[in] M the length of "x"
         * @param[in] x "x"
         * @param[out] X the output buffer
         * @param[in] LUA An option for calculation, defaults to 'A'. Due to the symmetry of "X", there is 3 way to calculate "X":
         * - 'L' only diagonal and lower elements are calculated
         * - 'U' only diagonal and upper elements are calculated
         * - 'A' all elements are calculated
         * - other do nothing
         */
        template <typename T>
        void vecSelfOuterProd(const int M, const T *const x, T *const X, const char LUA='A') {
            static_assert(std::is_floating_point<T>::value, "T must be floating point number type.");
            #define MEM_OFFSET(row,col) ((row)*M+col)
            for (int m=0; m<M; ++m) {X[m*(M+1)] = Analysis::sqAbs(x[m]);} // Calculate diagonal part.

            /* Calculate lower part. */
            if (LUA == 'L' || LUA == 'A') {
                for (int r=1; r<M; ++r) {
                    T *const X_ptr = &X[MEM_OFFSET(r,0)];
                    for (int c=0; c<r; ++c) {
                        X_ptr[c] = x[r]*x[c];
                    }
                }
            }

            /* Calculate upper part. */
            if (LUA == 'U') {
                for (int r=0; r<M-1; ++r) {
                    T *const X_ptr = &X[MEM_OFFSET(r,0)];
                    for (int c=r+1; c<M; ++c) {
                        X_ptr[c] = x[r]*x[c];
                    }
                }
            } else if (LUA == 'A') { // The lower part is already calculated.
                for (int r=0; r<M-1; ++r) {
                    T *const X_ptr = &X[MEM_OFFSET(r,0)];
                    const T *const X_T_ptr = &X[MEM_OFFSET(0,r)];
                    for (int c=r+1; c<M; ++c) {
                        X_ptr[c] = X_T_ptr[c*M];
                    }
                }
            }
            #undef MEM_OFFSET
        }

        /**
         * @brief Calculate a self outer product of a given complex vector "x".
         *
         * @tparam T the number type of the entries of "x"
         * @param[in] M the length of "x"
         * @param[in] x "x"
         * @param[out] X the output buffer
         * @param[inout] workspace a buffer to temporarily store conj(x).
         * @param[in] LUA An option for calculation, defaults to 'A'. Due to the Hermitian symmetry of "X", there is 3 way to calculate "X":
         * - 'L' only diagonal and lower elements are calculated
         * - 'U' only diagonal and upper elements are calculated
         * - 'A' all elements are calculated
         * - other do nothing
         */
        template <typename T>
        void vecSelfOuterProd(const int M, const std::complex<T> *const x, std::complex<T> *const X, std::complex<T> *const workspace, const char LUA='A') {
            #define MEM_OFFSET(row,col) ((row)*M+col)
            for (int m=0; m<M; ++m) {X[m*(M+1)] = Analysis::sqAbs(x[m]);} // Calculate diagonal part.

            /* Pre-calculate conj(x) */
            std::complex<T> *const conj_x = workspace;
            for (int i=0; i<M; ++i) {conj_x[i] = std::conj(x[i]);}

            /* Calculate lower part. */
            if (LUA == 'L' || LUA == 'A') {
                for (int r=1; r<M; ++r) {
                    std::complex<T> *const X_ptr = &X[MEM_OFFSET(r,0)];
                    for (int c=0; c<r; ++c) {
                        X_ptr[c] = x[r]*conj_x[c];
                    }
                }
            }

            /* Calculate upper part. */
            if (LUA == 'U') {
                for (int r=0; r<M-1; ++r) {
                    std::complex<T> *const X_ptr = &X[MEM_OFFSET(r,0)];
                    for (int c=r+1; c<M; ++c) {
                        X_ptr[c] = x[r]*conj_x[c];
                    }
                }
            } else if (LUA == 'A') { // The lower part is already calculated.
                for (int r=0; r<M-1; ++r) {
                    std::complex<T> *const X_ptr = &X[MEM_OFFSET(r,0)];
                    const std::complex<T> *const X_T_ptr = &X[MEM_OFFSET(0,r)];
                    for (int c=r+1; c<M; ++c) {
                        X_ptr[c] = Analysis::conj(X_T_ptr[c*M]);
                    }
                }
            }
            #undef MEM_OFFSET
        }

        /**
         * @brief Calculates matrix multiplication "AB" as "C" where "A", "B", "C" are compatible matrices.
         *
         * @tparam T the number type of the elements of the matrices
         * @param[in] l the number of the rows in the matrix "A"
         * @param[in] m the number of columns in the matrix "A" (= the number of the rows in the matrix "B")
         * @param[in] n the number of columns columns in the matrix "B"
         * @param[in] A the matrix A
         * @param[in] B the matrix B
         * @param[out] C the matrix C
         */
        template <typename T>
        void mulMat(const size_t l, const size_t m, const size_t n, const T *const A, const T *const B, T *const C) {
            T *ptr_C = C;
            for (size_t r=0; r<l; ++r) {
                for (size_t c=0; c<n; ++c) {
                    T sum = static_cast<T>(0);
                    const T *const ptr_A = A + r*m, *const ptr_B = B + c;
                    for (size_t i=0; i<m; ++i) {sum += ptr_A[i]*ptr_B[i*n];}
                    *ptr_C = sum; ++ptr_C;
                }
            }
        }

        /**
         * @brief Calculates Hermitian inner product of given 2 complex vectors "x", "y".
         * Hermitian inner product of "x" and "y" is defined as "<x^*, y>"", where "^*" represents conjugate and "<,>" represents inner product.
         *
         * @tparam T value type of complex number
         * @param[in] N vector length
         * @param[in] vec1 1st input vector, "x"
         * @param[in] vec2 2nd input vector, "y"
         * @param[in] stride1 The sampling interval for input vector "x". "x" is represented as "x = [vec1[0], vec1[stride1], ..., vec1[(N-1)*stride2]]".
         * @param[in] stride2 same as above for "y".
         * @return the Hermitian inner product of 2 input vectors
         */
        template <typename T>
        std::complex<T> hermitianInnerProduct(const size_t N, const std::complex<T> *const vec1, const std::complex<T> *const vec2, const size_t stride1, const size_t stride2) {
            auto sum = std::complex<T>(0, 0);
            for (size_t n=0, n1=0, n2=0; n<N; ++n, n1+=stride1, n2+=stride2) {sum += std::conj(vec1[n1])*vec2[n2];}
            return sum;
        }

        /**
         * @brief Calculates Hermitian inner product of given 2 complex vectors "x", "y".
         * Hermitian inner product of "x" and "y" is defined as "<x^*, y>", where "^*" represents conjugate and "<,>" represents inner product.
         *
         * @tparam T the number type of complex number's real and imaginary part
         * @param[in] N vector length
         * @param[in] vec1 1st input vector, "x"
         * @param[in] vec2 2nd input vector, "y"
         * @param[in] stride The sampling interval for input vector "x". "x" is represented as "x = [vec1[0], vec1[stride], ..., vec1[(N-1)*stride]]". The same applies to "y" too.
         * @return the Hermitian inner product of 2 input vectors
         */
        template <typename T>
        std::complex<T> hermitianInnerProduct(const size_t N, const std::complex<T> *const vec1, const std::complex<T> *const vec2, const size_t stride=1) {
            return hermitianInnerProduct(N, vec1, vec2, stride, stride);
        }

        /**
         * @brief Calculates L-2 norm of a given complex vector "x".
         *
         * @tparam T the number type of complex number's real and imaginary part
         * @param[in] N vector length
         * @param[in] vec input vector
         * @param[in] stride The sampling interval for input vectors. "x" is represented as "x = [vec1[0], vec1[stride], ..., vec1[(N-1)*stride]]".
         * @return the L-2 norm of input vector
         */
        template <typename T>
        T l2Norm(const size_t N, const std::complex<T> *const vec, const size_t stride=1) {
            return std::sqrt(hermitianInnerProduct(N, vec, vec, stride).real());
        }

        /**
         * @brief Calculates the LDL decomposition of a given Hermitian-and-invertible matrix "A".
         * Find a lower-triangle matrix "L" and a diagonal matrix D such that "A = LDL^*".
         * The diagonal entries of "L" are all 1, and the diagonal entries of "D" are all real numbers.
         *
         * @tparam T the number type of complex number's real and imaginary part
         * @param[in] m the number of the rows and columns in the matrix "A"
         * @param[in] A the matrix "A"
         * @param[out] d the diagonal elements of D
         * @param[out] L The matrix "L". The upper part is NOT modified by this function.
         * @param[inout] workspace The pointer to a continuous memory space whose size is "(m-1)*sizeof(std::complex<T>)" bytes. This space is used during calculation.
         * @param[in] epsilon The threshold used for zero-division detection. Zero-division is detected when the absolute value of a divider is smaller than "epsilon".
         * @retval false A zero-division is detected and calculation is aborted. Maybe "A" is non-invertible.
         * @retval true The calculation is successfully done.
         */
        template <typename T>
        bool ldlDecomp(const int m, const std::complex<T> *const A, T *const d, std::complex<T> *const L, std::complex<T> *const workspace, const T epsilon=1e-12) {
            #define MEM_OFFSET(row, col) ((row)*m+col)
            static_assert(std::is_floating_point<T>::value, "T must be floating point number type.");
            constexpr std::complex<T> ONE = 1;
            const std::complex<T> *A_diag_ptr = A;
            for (int i=0; i<m; ++i) {
                T di = A_diag_ptr->real(); A_diag_ptr += m+1; // d[i] <- Re(A[i,i])
                {
                    std::complex<T> *const L_row_ptr = &L[MEM_OFFSET(i,0)];
                    for (int j=0; j<i; ++j) {di -= d[j]*MathLib::Analysis::sqAbs(L_row_ptr[j]);} // d[i] <- d[i] - d[j]*|L[i,j]|^2
                }
                d[i] = di;
                if (std::abs(di) < epsilon) {
                    return false;
                }
                const T inv_di = 1/d[i];
                L[MEM_OFFSET(i,i)] = ONE;

                /* Construct "d[j]*conj(L[i,j]) (j=0,1, ..., i-1)" */
                std::complex<T> *const dL_ptr = workspace;
                {
                    std::complex<T> *const L_row_ptr = &L[MEM_OFFSET(i,0)];
                    for (int j=0; j<i; ++j) {dL_ptr[j] = d[j]*std::conj(L_row_ptr[j]);}
                }

                const std::complex<T> *A_col_ptr = &A[MEM_OFFSET(i+1,i)];
                std::complex<T> *L_col_ptr = &L[MEM_OFFSET(i+1,i)];
                for (int k=i+1; k<m; ++k) {
                    std::complex<T> L_ki = *A_col_ptr; A_col_ptr += m; // L[k,i] <- A[k,i]
                    std::complex<T> *const L_row_ptr = &L[MEM_OFFSET(k,0)];
                    for (int j=0; j<i; ++j) {L_ki -= L_row_ptr[j]*dL_ptr[j];} // L[k,i] <- L[k,i] - L[k,j]d[j]*conj(L[i,j])
                    *L_col_ptr = inv_di*L_ki; L_col_ptr += m; // L[k,i] <- L[k,i]/d[i]
                }
            }
            return true;
            #undef MEM_OFFSET
        }

        /**
         * @brief Calculates the LDL decomposition of a given Hermitian-and-invertible matrix "A".
         * Find a lower-triangle matrix "L" and a diagonal matrix D such that "A = LDL^*".
         * The diagonal entries of "L" are all 1, and the diagonal entries of "D" are all real numbers.
         *
         * @tparam T the number type of matrix "A"
         * @param[in] m the number of the rows and columns in the matrix "A"
         * @param[in] A the matrix "A"
         * @param[out] d the diagonal elements of D
         * @param[out] L The matrix "L". The upper part is NOT modified by this function.
         * @param[inout] workspace The pointer to a continuous memory space whose size is "(m-1)*sizeof(T)" bytes. This space is used during calculation.
         * @param[in] epsilon The threshold used for zero-division detection. Zero-division is detected when the absolute value of a divider is smaller than "epsilon".
         * @retval false A zero-division is detected and calculation is aborted. Maybe "A" is non-invertible.
         * @retval true The calculation is successfully done.
         */
        template <typename T>
        bool ldlDecomp(const int m, const T *const A, T *const d, T *const L, T *const workspace, const T epsilon=1e-12) {
            #define MEM_OFFSET(row, col) ((row)*m+col)
            static_assert(std::is_floating_point<T>::value, "T must be floating point number type.");
            const T *A_diag_ptr = A;
            for (int i=0; i<m; ++i) {
                T di = *A_diag_ptr; A_diag_ptr += m+1; // d[i] <- A[i,i]
                {
                    T *const L_row_ptr = &L[MEM_OFFSET(i,0)];
                    for (int j=0; j<i; ++j) {di -= d[j]*MathLib::Analysis::sqAbs(L_row_ptr[j]);} // d[i] <- d[i] - d[j]*|L[i,j]|^2
                }
                d[i] = di;
                if (std::abs(di) < epsilon) {
                    return false;
                }
                const T inv_di = 1/d[i];
                L[MEM_OFFSET(i,i)] = 1;

                /* Construct "d[j]*conj(L[i,j]) (j=0,1, ..., i-1)" */
                T *const dL_ptr = workspace;
                {
                    T *const L_row_ptr = &L[MEM_OFFSET(i,0)];
                    for (int j=0; j<i; ++j) {dL_ptr[j] = d[j]*L_row_ptr[j];}
                }

                const T *A_col_ptr = &A[MEM_OFFSET(i+1,i)];
                T *L_col_ptr = &L[MEM_OFFSET(i+1,i)];
                for (int k=i+1; k<m; ++k) {
                    T L_ki = *A_col_ptr; A_col_ptr += m; // L[k,i] <- A[k,i]
                    T *const L_row_ptr = &L[MEM_OFFSET(k,0)];
                    for (int j=0; j<i; ++j) {L_ki -= L_row_ptr[j]*dL_ptr[j];} // L[k,i] <- L[k,i] - L[k,j]d[j]*L[i,j]
                    *L_col_ptr = inv_di*L_ki; L_col_ptr += m; // L[k,i] <- L[k,i]/d[i]
                }
            }
            return true;
            #undef MEM_OFFSET
        }

        /**
         * @brief Solves linear equation "Ax = b" by Gaussian elimination, where "A" is a square invertible matrix of size "m", and "b" is a vector of length "m".
         * The matrix "A" must be invertible, otherwise the behavior is undefined.
         * When "A" is known to be Hermitian, use `solveLinEqHermitian` function rather than this function, to achieve higher performance.
         *
         * @tparam T the number type of the elements of the matrix "A", vector "b" and "x"
         * @param[in] m the number of the rows and columns in the matrix "A"
         * @param[in] A the matrix "A"
         * @param[in] b the vector "b"
         * @param[out] x the vector "x"
         * @param[out] workspace The pointer to a continuous memory space whose size is "m*(m+1)*sizeof(T) + m*sizeof(size_t)" bytes. This space is used for the extended-matrix "[A, b]" and the row exchanging table.
         */
        template <typename T>
        void solveLinearEquation(const size_t m, const T *const A, const T *const b, T *const x, char *workspace) {
            #define MEM_OFFSET_E(row, col) ((row)*(m+1)+col)
            constexpr T ZERO = 0.0;
            constexpr T ONE = 1.0;

            /* Creates an extended matrix "E". */
            T *E = reinterpret_cast<T *>(workspace);
            workspace += m*(m+1)*sizeof(T);
            {
                const T *ptr_A = A, *ptr_b = b;
                T *ptr_E = E;
                for (size_t r=0; r<m; ++r) {
                    std::memcpy(ptr_E, ptr_A, m*sizeof(T));
                    ptr_E += m; ptr_A += m;
                    *ptr_E = *ptr_b;
                    ++ptr_E; ++ptr_b;
                }
            }

            /* Initializes pivoting table. */
            size_t *rowExchgTable = reinterpret_cast<size_t *>(workspace);
            std::iota(&rowExchgTable[0], &rowExchgTable[m], 0);

            for (size_t r=0; r<m; ++r) {
                /* Pivoting. Finds a row index r2 that |E[r2,r]| is largest in the set {|E[r,r]|, |E[r,r+1]|, ..., |E[r,m-1]|}. */
                {
                    size_t r2 = r;
                    auto maxAbsVal = std::abs(E[MEM_OFFSET_E(rowExchgTable[r], r)]);
                    for (size_t r3=r+1; r3<m; ++r3) {
                        const auto currentAbsVal = std::abs(E[MEM_OFFSET_E(rowExchgTable[r3], r)]);
                        if (currentAbsVal > maxAbsVal) {
                            r2 = r3;
                            maxAbsVal = currentAbsVal;
                        }
                    }
                    std::swap(rowExchgTable[r], rowExchgTable[r2]);
                }

                /* Normalizes current row. */
                {
                    T *const ptr_E = E + MEM_OFFSET_E(rowExchgTable[r], 0);
                    const T inv_E_rr = ONE/ptr_E[r];
                    ptr_E[r] = ONE;
                    for (size_t c=r+1; c<m+1; ++c) {ptr_E[c] *= inv_E_rr;}
                }

                /* Reduces the rows over/under current row. */
                {
                    for (size_t r2=0; r2<m; ++r2) {
                        if (r2 == r) {
                            continue;
                        }
                        T *const ptr2_E = E + MEM_OFFSET_E(rowExchgTable[r2], 0);
                        const T minus_E_r2r = -ptr2_E[r];
                        ptr2_E[r] = ZERO;
                        const T *const ptr1_E = E + MEM_OFFSET_E(rowExchgTable[r], 0);
                        for (size_t c=r+1; c<m+1; ++c) {
                            ptr2_E[c] += minus_E_r2r*ptr1_E[c];
                        }
                    }
                }
            }

            for (size_t r=0; r<m; ++r) {x[r] = E[MEM_OFFSET_E(rowExchgTable[r], m)];} // Stores the result.
            #undef MEM_OFFSET_E
        }

        /**
         * @brief Solve linear equation "Ax = b" using LDL decomposition, where "A" is a Hermitian and invertible matrix of size "m", and "b" is a vector of length "m".
         * The matrix "A" MUST be invertible.
         *
         * @tparam T the number type of the elements of the matrix "A", vector "b" and "x"
         * @param[in] m the number of the rows and columns in the matrix "A"
         * @param[in] A the matrix "A"
         * @param[in] b the vector "b"
         * @param[out] x the vector "x"
         * @param[out] workspace The pointer to a continuous memory space. Its size is calculated as follows:@n
         * - when A's type is floating-point T: (m*m+2*m-1)*sizeof(T)
         * - when A's type is std::complex<T>: m*sizeof(T) + (m*m+m-1)*sizeof(std::complex<T>)
         * @retval false A zero-division is detected and calculation is aborted. Maybe "A" is non-invertible.
         * @retval true The calculation is successfully done.
         */
        template <typename T>
        bool solveLinEqHermitian(const int m, const T *const A, const T *const b, T *const x, char *workspace) {
            #define MEM_OFFSET(row, col) ((row)*m+col)
            if (m <= 0) {
                return true;
            }

            /* LDL decomposition */
            typedef decltype(std::real(*A)) T_real; // floating-point type "A" is ok
            T_real *const d = reinterpret_cast<T_real *>(workspace); workspace += m*sizeof(T_real);
            T *const L = reinterpret_cast<T *>(workspace); workspace += m*m*sizeof(T);
            const bool noZeroDiv = ldlDecomp(m, A, d, L, reinterpret_cast<T *>(workspace));
            if (!noZeroDiv) {
                return false;
            }

            /* forward substitution (Ly == b) */
            T *const y = x; // Use "x" as a temporary buffer.
            for (int i=0; i<m; ++i) {
                T yi = b[i];
                T *const L_ptr = &L[MEM_OFFSET(i,0)];
                for (int j=0; j<i; ++j) {yi -= L_ptr[j]*y[j];}
                y[i] = yi;
            }

            /* backward substitution (DL^*x == y) */
            for (int i=m-1; i>=0; --i) {
                T xi = y[i]/d[i];
                for (int j=i+1; j<m; ++j) {
                    xi -= MathLib::Analysis::conj(L[MEM_OFFSET(j,i)])*x[j];
                }
                x[i] = xi;
            }

            return true;
            #undef MEM_OFFSET
        }
    }
}

#endif // __MATH_LIB__LINALG_HPP__
