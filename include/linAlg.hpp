#ifndef __MOTCHY_MATH_LIB__LINALG_HPP__
#define __MOTCHY_MATH_LIB__LINALG_HPP__

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <numeric>

namespace MotchyMathLib {
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
         * It is assumed that the matrices' elements are aligned on memory in row-oriented order.
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
        bool isEqualMat(size_t m, size_t n, const T1 *A, const T1 *B, T2 epsilon) {
            const T1 *ptr_A = A, *ptr_B = B;
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                if (std::abs(*ptr_A - *ptr_B) > epsilon) {return false;}
                ++ptr_A; ++ptr_B;
            }
            return true;
        }

        /**
         * @brief Checks whether given two vectors "a" and "b" are equal or not.
         * @details "a" and "b" are considered to be equal if and only if the absolute value of all the element of "a-b" is strictly less than a given number, epsilon.
         * It is assumed that the vector's elements are aligned on memory in row-oriented order.
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
        bool isEqualVec(size_t m, const T1 *a, const T1 *b, T2 epsilon) {
            return isEqualMat(1, m, a, b, epsilon);
        }

        /**
         * @brief Calculates transpose of a matrix "A" as "B".
         * @details It is assumed that the matrix's elements are aligned on memory in row-oriented order.
         *
         * @tparam T the number type of the elements of input matrix
         * @param[in] m the number of the rows in the input matrix "A"
         * @param[in] n the number of the columns in the input matrix "A"
         * @param[in] A the input matrix
         * @param[out] B the output matrix
         */
        template <typename T>
        void transposeMat(size_t m, size_t n, const T *A, T *B) {
            const T *ptr_A = A;
            T *ptr_B = B;
            for (size_t r=0; r<m; ++r) {
                ptr_B = B + r;
                for (size_t c=0; c<n; ++c) {
                    *ptr_B = *ptr_A;
                    ++ptr_A;
                    ptr_B += m;
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
        void conjugateMat(size_t m, size_t n, const std::complex<T> *A, std::complex<T> *B) {
            const std::complex<T> *ptr_A = A;
            std::complex<T> *ptr_B = B;
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                *ptr_B = std::conj(*ptr_A);
                ++ptr_A; ++ptr_B;
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
        void conjugateVec(size_t m, const std::complex<T> *x, std::complex<T> *y) {
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
        void addMat_inplace(size_t m, size_t n, T *A, const T *B) {
            T *ptr_A = A;
            const T *ptr_B = B;
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                *ptr_A += *ptr_B;
                ++ptr_A; ++ptr_B;
            }
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
        void addVec_inplace(size_t m, T *x, const T *y) {
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
        void scaleMat(size_t m, size_t n, T c, const T *A, T *B) {
            const T *ptr_A = A;
            T *ptr_B = B;
            const size_t L = m*n;
            for (size_t i=0; i<L; ++i) {
                *ptr_B = c*(*ptr_A);
                ++ptr_A; ++ptr_B;
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
        void scaleVec(size_t m, T c, const T *a, T *b) {
            scaleMat(1, m, c, a, b);
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
        void mulMat(size_t l, size_t m, size_t n, const T *A, const T *B, T *C) {
            T *ptr_C = C;
            for (size_t r=0; r<l; ++r) {
                for (size_t c=0; c<n; ++c) {
                    T sum = static_cast<T>(0);
                    const T *ptr_A = A + r*m, *ptr_B = B + c;
                    for (size_t i=0; i<m; ++i) {
                        sum += (*ptr_A)*(*ptr_B);
                        ++ptr_A;
                        ptr_B += n;
                    }
                    *ptr_C = sum;
                    ++ptr_C;
                }
            }
        }

        /**
         * @brief Calculates Hermitian inner product of given 2 vectors "x", "y".
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
        std::complex<T> hermitianInnerProduct(size_t N, const std::complex<T> *vec1, const std::complex<T> *vec2, size_t stride1, size_t stride2) {
            auto sum = std::complex<T>(0, 0);
            const std::complex<T> *ptr1 = vec1, *ptr2 = vec2;

            for (size_t n=0; n<N; ++n) {
                sum += std::conj((*ptr1))*(*ptr2);
                ptr1 += stride1;
                ptr2 += stride2;
            }

            return sum;
        }

        /**
         * @brief Calculates Hermitian inner product of given 2 vectors "x", "y".
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
        std::complex<T> hermitianInnerProduct(size_t N, const std::complex<T> *vec1, const std::complex<T> *vec2, size_t stride=1) {
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
        T l2Norm(size_t N, const std::complex<T> *vec, size_t stride=1) {
            return std::sqrt(hermitianInnerProduct(N, vec, vec, stride).real());
        }

        /**
         * @brief Solves linear equation "Ax = b" by Gaussian elimination, where "A" is a square invertible matrix of size "m", and "b" is a vector of length "m".
         * @details The matrix "A" must be invertible, otherwise the behavior is undefined.
         *
         * @tparam T the number type of the elements of the matrix "A", vector "b" and "x"
         * @param[in] m the number of the rows and columns in the matrix "A"
         * @param[in] A the matrix "A"
         * @param[in] b the vector "b"
         * @param[out] x the vector "x"
         * @param[out] workspace The pointer to a continuous memory space whose size is "m*(m+1)*sizeof(T) + m*sizeof(size_t)" bytes. This space is used for the extended-matrix "[A, b]" and the row exchanging table.
         */
        template <typename T>
        void solveLinearEquation(size_t m, const T *A, const T *b, T *x, char *workspace) {
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
                    T *ptr_E = E + MEM_OFFSET_E(rowExchgTable[r], r);
                    const T E_rr = *ptr_E;
                    *ptr_E = ONE;
                    ++ptr_E;
                    for (size_t c=r+1; c<m+1; ++c) {
                        *ptr_E /= E_rr;
                        ++ptr_E;
                    }
                }

                /* Reduces the rows over/under current row. */
                {
                    for (size_t r2=0; r2<m; ++r2) {
                        if (r2 == r) {
                            continue;
                        }
                        T *ptr2_E = E + MEM_OFFSET_E(rowExchgTable[r2], r);
                        const T E_r2r = *ptr2_E;
                        *ptr2_E = ZERO;
                        ++ptr2_E;
                        const T *ptr1_E = E + MEM_OFFSET_E(rowExchgTable[r], r+1);
                        for (size_t c=r+1; c<m+1; ++c, ++ptr2_E, ++ptr1_E) {
                            *ptr2_E -= E_r2r*(*ptr1_E);
                        }
                    }
                }
            }

            /* Stores the result. */
            {
                T *ptr_x = x;
                for (size_t r=0; r<m; ++r, ++ptr_x) {
                    *ptr_x = E[MEM_OFFSET_E(rowExchgTable[r], m)];
                }
            }

            #undef MEM_OFFSET_E
        }
    }
}

#endif // __MOTCHY_MATH_LIB__LINALG_HPP__
