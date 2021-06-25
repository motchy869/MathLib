#ifndef __SIGPROC_HPP__
#define __SIGPROC_HPP__

#include <algorithm>
#include <cstddef>

namespace Utils {
    namespace Math {
        namespace SigProc {
            /**
             * @brief Calculates the first "c(>=0)" elements of the convolution of 2 vectors "x1" and "x2" as "y", with stride "d(>=1)".
             * @details Let "N1(>=1)", "N2(>=1)" be the length of "x1", "x2" respectively, and consider "x2[n] (n<=-1, n>=N2)" to be 0 virtually.
             * The length of the result vector "y" becomes "Ly = 1+floor((N1+N2-2)/d)" when "c<0", otherwise "Ly = c".
             * Each element of "y" is calculated as: "y[n] = x1[0]*x2[n*d] + x1[1]*x2[n*d-1] + ... + x1[N1-1]*x2[n*d-N1+1] (n=0,1, ..., Ly-1)".
             * When "min(N1,N2,d)==0", this function does nothing and returns immediately.
             *
             * @tparam T the number type of the elements of the vector "x1", "x2"
             * @param[in] N1 the length of "x1"
             * @param[in] N2 the length of "x2"
             * @param[in] x1 "x1"
             * @param[in] x2 "x2"
             * @param[out] y "y"
             * @param[in] c "c"
             * @param[in] d "d"
             */
            template <typename T>
            void convolve(size_t N1, size_t N2, const T *x1, const T *x2, T *y, int c, size_t d) {
                if (0 == N1 || 0 == N2 || 0 == d) {
                    return;
                }

                constexpr T ZERO = static_cast<T>(0);

                /* Assigns shorter to "v1", longer to "v2". */
                const T *v1 = x1, *v2 = x2;
                if (N1 > N2) {
                    v1 = x2; v2 = x1;
                    std::swap(N1,N2);
                }
                // Now, 1 <= length(v1) == N1 <= N2 == length(v2).

                const size_t Ly = (c<0) ? 1 + (N1+N2-2)/d : c;
                const T *const v1_end = v1+N1-1;
                const T *const v2_end = v2+N2-1;
                const T *const y_end = y+Ly-1;
                T *ptr_y = y;
                size_t i = 0;

                /* stage 1 */
                for (const size_t i_end = (N2-1)/d; i<=i_end && ptr_y<=y_end; ++i) {
                    const T *ptr_v1 = v1, *ptr_v2 = v2+i*d;
                    T sum = ZERO;
                    while (ptr_v1 <= v1_end && ptr_v2 >= v2) {
                        sum += (*ptr_v1)*(*ptr_v2);
                        ++ptr_v1; --ptr_v2;
                    }
                    *ptr_y = sum;
                    ++ptr_y;
                }

                /* stage 2*/
                for (const size_t i_end = (N1+N2-2)/d; i<=i_end && ptr_y<=y_end; ++i) {
                    const T *ptr_v1 = (v1+i*d+1)-N2, *ptr_v2 = v2_end;
                    T sum = ZERO;
                    while (ptr_v1 <= v1_end) {
                        sum += (*ptr_v1)*(*ptr_v2);
                        ++ptr_v1; --ptr_v2;
                    }
                    *ptr_y = sum;
                    ++ptr_y;
                }
            }

            /**
             * @brief Calculates the convolution of 2 vectors "x1" and "x2" as "y".
             * @details Let "N1", "N2" be the length of "x1", "x2" respectively. The length of "y" becomes "N1+N2-1", and "y[0] = x1[0]x2[0]", "y[N1+N2-2] = x1[N1-1]x2[N2-1]": That is to say, the elements out of range are considered to be 0.
             *
             * @tparam T the number type of the elements of the vector "x1", "x2"
             * @param[in] N1 the length of "x1"
             * @param[in] N2 the length of "x2"
             * @param[in] x1 "x1"
             * @param[in] x2 "x2"
             * @param[out] y "y"
             */
            template <typename T>
            void convolve(size_t N1, size_t N2, const T *x1, const T *x2, T *y) {
                #if true
                    convolve(N1, N2, x1, x2, y, -1, 1);
                #else // Replaced by `void convolve(size_t N1, size_t N2, const T *x1, const T *x2, T *y, int c, size_t d)`.
                    constexpr T ZERO = static_cast<T>(0);

                    /* Assigns shorter to "v1", longer to "v2". */
                    const T *v1 = x1, *v2 = x2;
                    if (N1 > N2) {
                        v1 = x2; v2 = x1;
                        std::swap(N1,N2);
                    }
                    // Now, length(v1) == N1 <= N2 == length(v2).

                    if (0 == N1) {
                        *y = ZERO;
                        return;
                    }

                    const T *const v1_end = v1+N1;
                    const T *const v2_end = v2+N2;
                    T *ptr_y = y;

                    /* stage 1: y[0] ... y[N2-1] */
                    for (size_t i=0; i<N2; ++i) {
                        const T *ptr_v1 = v1, *ptr_v2 = v2+i;
                        T sum = ZERO;
                        while (ptr_v1 < v1_end && ptr_v2 >= v2) {
                            sum += (*ptr_v1)*(*ptr_v2);
                            ++ptr_v1; --ptr_v2;
                        }
                        *ptr_y = sum;
                        ++ptr_y;
                    }

                    /* stage 2  y[N2] ... y[N1+N2-2] if N1 >= 2 */
                    for (size_t i=1; i<N1; ++i) {
                        const T *ptr_v1 = v1+i, *ptr_v2 = v2_end-1;
                        T sum = ZERO;
                        while (ptr_v1 < v1_end) {
                            sum += (*ptr_v1)*(*ptr_v2);
                            ++ptr_v1; --ptr_v2;
                        }
                        *ptr_y = sum;
                        ++ptr_y;
                    }
                #endif
            }

            /**
             * @brief Calculates the convolution of 2 vectors "x1" and "x2" as "y".
             * @details Let "L1, N" be positive integers. The output is calculated as: "y[n] = \sum_{k=0}^{L1-1} x1[k]*x2[n*d-k] (n=0,1,...,N-1)".
             * @tparam T the number type of the elements of the vector "x1", "x2"
             * @param[in] L1 "L1"
             * @param[in] x1 "x1"
             * @param[in] x2 "x2"
             * @param[out] y "y"
             * @param[in] N "N"
             * @param[in] d "d"
             */
            template <typename T>
            void convolve_type2(size_t L1, const T *x1, const T *x2, T *y, size_t N, size_t d=1) {
                constexpr T ZERO = static_cast<T>(0);

                for (size_t n=0; n<N; ++n) {
                    T sum = ZERO;
                    for (size_t k=0; k<L1; ++k) {
                        const size_t nd = n*d;
                        sum += x1[k]*x2[nd-k];
                    }
                    y[n] = sum;
                }
            }
        }
    }
}

#endif // __SIGPROC_HPP__
