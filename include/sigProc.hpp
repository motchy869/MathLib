/**
 * @author motchy (motchy869[at]gmail.com)
 * @brief Signal Processing library
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include "common.hpp"
#include "analysis.hpp"

namespace MathLib {
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
        void convolve(size_t N1, size_t N2, const T *const x1, const T *const x2, T *const y, const int c, const size_t d) {
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

            const size_t Ly = (c<0 ? 1 + (N1+N2-2)/d : c);
            const int i_v1_end = N1-1;
            const T *const y_end = y+Ly-1;
            T *ptr_y = y;
            size_t i = 0;

            /* stage 1 */
            for (const size_t i_end = (N2-1)/d; i<=i_end && ptr_y<=y_end; ++i) {
                T sum = ZERO;
                const int i_v2_end=0;
                for (int i_v1=0, i_v2 = i*d; i_v1<=i_v1_end && i_v2>=i_v2_end; ++i_v1, --i_v2) {
                    Analysis::addProd(v1[i_v1], v2[i_v2], sum);
                }
                *ptr_y = sum; ++ptr_y;
            }

            /* stage 2*/
            for (const size_t i_end = (N1+N2-2)/d; i<=i_end && ptr_y<=y_end; ++i) {
                T sum = ZERO;
                sum = ZERO;
                for (int i_v1=i*d+1-N2, i_v2=N2-1; i_v1<=i_v1_end; ++i_v1, --i_v2) {
                    Analysis::addProd(v1[i_v1], v2[i_v2], sum);
                }
                *ptr_y = sum; ++ptr_y;
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
            convolve(N1, N2, x1, x2, y, -1, 1);
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
                const T *const x2_ptr = &x2[n*d];
                for (size_t k=0; k<L1; ++k) {Analysis::addProd(x1[k], x2_ptr[-k], sum);}
                y[n] = sum;
            }
        }

        /**
         * @brief Exponential smoothing filter.
         * For input "x[n]", Output "y[n]" is calculated with a smoothing coefficient "lambda" as follows:
         *   y[n] = lambda*x[n] + (1-lambda)*y[n-1]
         * where:
         *   y[-1] := static_cast<T_signal>(0)
         *   0 < lambda < 1.
         *
         * @details "y" can be set to any value in arbitrary time by calling a method `setState`.
         * One typically set "y" to "x[0]" before start inputting the sequence "x[0], x[1], ... " so that the 1st output from the filter begins with "x[0]" rather than 0.
         *
         * @tparam T_lambda data type of lambda
         * @tparam T_signal data type of x[n]
         */
        template <typename T_lambda, typename T_signal>
        class ExpSmoothingFilter {
            private:
                const T_lambda m_lambda, m_co_lambda;
                T_signal m_y;
            public:
                /**
                 * @brief Construct a new filter object
                 *
                 * @param[in] lambda a smoothing factor, described in class docstring.
                 */
                ExpSmoothingFilter(T_lambda lambda) : m_lambda(lambda), m_co_lambda(static_cast<T_lambda>(1) - lambda), m_y(static_cast<T_signal>(0)) {
                    static_assert(std::is_floating_point<T_lambda>::value, "`lambda` must be a floating point number.");
                    assert(lambda > 0 && lambda < 1);
                }

                /**
                 * @brief Set the internal state "y" to a given value.
                 *
                 * @param[in] y the value to which "y" is set
                 */
                void setState(T_signal y) {m_y = y;}

                /**
                 * @brief Input a value to filter, and get output.
                 *
                 * @param[in] x input value, "x[n]"
                 * @return filter output, "y[n]"
                 */
                T_signal apply(T_signal x) {
                    m_y = m_lambda*x + m_co_lambda*m_y;
                    return m_y;
                }
        };

        /**
         * @brief Stable recursive moving average filter.
         * This filter calculates the moving "N"-backward average of the input sequence RECURSIVELY, where "N" is a positive integer.
         * This filter is stable against computational errors in floating point numbers.
         * @details The system equation is:
         *   sum[n] = alpha*sum[n-1] + u[n] - (alpha^N)*u[n-N]
         *   y[n] = sum[n]/N
         * The system transfer function is:
         *    H(z) = (1 - (alpha^N)*z^(-N))/(1 - alpha*z^(-1))
         * The block diagram of this filter is:
         *   docs/sigProc/StableRecursiveMovAvgFilterA/block-diagram.png
         * Where "0 < alpha < 1", "u" is an input sequence, "y" is an output sequence, and "sum" is an internal state.
         * Usually, "alpha" is chosen to be close to 1.
         * The calculation error that occurs in every moment is converged to 0 in the process of circulating the loop.
         * @tparam T_alpha the data type of "alpha"
         * @tparam T_signal the data type of the input signal
         */
        template <typename T_alpha, typename T_signal>
        class StableRecursiveMovAvgFilterA {
            private:
                static constexpr T_alpha ONE_ALPHA = static_cast<T_alpha>(1);
                static constexpr T_signal ZERO_SIGNAL = static_cast<T_signal>(0);
                const int m_N;
                const T_alpha m_inv_N;
                const T_alpha m_alpha, m_alpha_pow_N;
                T_signal *m_buffer;
                int m_idx;
                T_signal m_sum;

            public:
                /**
                 * @brief Construct a new filter object
                 *
                 * @param[in] N the averaging window size
                 * @param[in] alpha the stability factor, described in class docstring.
                 */
                StableRecursiveMovAvgFilterA(const int N, const T_alpha alpha) : m_N(N), m_inv_N(ONE_ALPHA/N), m_alpha(alpha), m_alpha_pow_N(std::pow(alpha, N)) {
                    static_assert(std::is_floating_point<T_alpha>::value, "`alpha` must be a floating point number.");
                    assert(N > 0);
                    assert(0 < alpha && alpha < 1);

                    m_buffer = new T_signal[N];
                    reset();
                }

                /**
                 * @brief Set the delay line and internal state to 0.
                 */
                void reset() {
                    std::fill_n(m_buffer, m_N, ZERO_SIGNAL);
                    m_idx = 0;
                    m_sum = ZERO_SIGNAL;
                }

                /**
                 * @brief Update the internal state with new input and get the output.
                 * @param[in] u the input value
                 * @return the output value
                 */
                T_signal step(const T_signal u) {
                    m_sum = m_alpha*m_sum - m_alpha_pow_N*m_buffer[m_idx] + u;
                    m_buffer[m_idx] = u;
                    m_idx += 1;
                    if (m_idx == m_N) {m_idx = 0;}
                    return m_inv_N*m_sum;
                }

                ~StableRecursiveMovAvgFilterA() {delete[] m_buffer;}
        };
    }
}
