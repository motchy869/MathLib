/**
 * @author motchy (motchy869[at]gmail.com)
 * @brief Analysis library
 */
#pragma once

#include <cmath>
#include <complex>
#include <type_traits>
#include <vector>
#include "common.hpp"

namespace MathLib {
    namespace Analysis {
        /**
         * @brief Take the conjugate of a complex number "x"
         *
         * @tparam T the number type of the real part of "x"
         * @param[in] x "x"
         * @return the conjugate of "x"
         */
        template <typename T>
        inline static std::complex<T> __attribute__((always_inline)) conj(std::complex<T> x) {
            return std::move(std::complex<T>(x.real(), -x.imag()));
        }

        /**
         * @brief Take the conjugate of a REAL number "x"
         *
         * @tparam T the number type of "x"
         * @param[in] x "x"
         * @return x
         */
        template <typename T>
        inline T conj(T x) {return x;}

        /**
         * @brief Calculate "|x|^2 = Re(x.conj()*x)" for complex number "x".
         *
         * @tparam T the number type of the real part of "x"
         * @param[in] x "x"
         */
        template <typename T>
        inline static T __attribute__((always_inline)) sqAbs(std::complex<T> x) {
            return x.real()*x.real() + x.imag()*x.imag();
        }

        /**
         * @brief Calculate "|x|^2" for real number "x".
         *
         * @tparam T the number type of "x"
         * @param[in] x "x"
         */
        template <typename T>
        inline static T __attribute__((always_inline)) sqAbs(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            return x*x;
        }

        /**
         * @brief Calculate the product of 2 floating point numbers "x1" and "x2".
         * Calling this function is completely equivalent to just writing "x1*x2".
         *
         * @tparam T the number type of "x1" and "x2"
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @return "x1*x2"
         */
        template <typename T>
        inline static T __attribute__((always_inline)) prod(const T x1, const T x2) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            return x1*x2;
        }

        /**
         * @brief Calculate the product of a floating point number "x1" and a complex number "x2".
         * This function is expected to work about 7 times faster than normal operation such as "y = x1*x2".@n
         * CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB@n
         * N = 1e8@n
         * normal: 350 ms@n
         * this method: 51 ms
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @return "x1*x2"
         */
        template <typename T>
        inline static std::complex<T> __attribute__((always_inline)) prod(const T x1, const std::complex<T> x2) {
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            return std::move(std::complex<T>(x1*x2_vec[0], x1*x2_vec[1]));
        }

        /**
         * @brief Calculate the product of 2 complex numbers "x1" and "x2".
         * This function is expected to work about 7 times faster than normal operation such as "y = x1*x2".@n
         * CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB@n
         * N = 1e8@n
         * normal: 350 ms@n
         * this method: 51 ms
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @return "x1*x2"
         */
        template <typename T>
        inline static std::complex<T> __attribute__((always_inline)) prod(const std::complex<T> x1, const std::complex<T> x2) {
            const T *const x1_vec = reinterpret_cast<const T *>(&x1);
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            const T real = x1_vec[0]*x2_vec[0] - x1_vec[1]*x2_vec[1];
            const T imag = x1_vec[0]*x2_vec[1] + x1_vec[1]*x2_vec[0];
            return std::move(std::complex<T>(real, imag));
        }

        /**
         * @brief Add the product of a floating point number "x1" and a complex number "x2".
         * This function is expected to work much faster than normal operation such as "y += x1*x2".
         * (see also: addProd(const std::complex<T> &x1, const std::complex<T> &x2, std::complex<T> &y))
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) addProd(const T x1, const std::complex<T> x2, std::complex<T> &y) {
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            T *const y_vec = reinterpret_cast<T *>(&y);
            y_vec[0] += x1*x2_vec[0]; // real part
            y_vec[1] += x1*x2_vec[1]; // imaginary part
        }

        /**
         * @brief Add the product of two complex numbers "x1" and "x2" to "y".
         * This function is expected to work about 20% faster than normal operation such as "y += x1*x2".@n
         * CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB@n
         * N = 1e8@n
         * normal: 411 ms@n
         * this method: 330 ms
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) addProd(const std::complex<T> x1, const std::complex<T> x2, std::complex<T> &y) {
            const T *const x1_vec = reinterpret_cast<const T *>(&x1);
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            T *const y_vec = reinterpret_cast<T *>(&y);
            y_vec[0] += x1_vec[0]*x2_vec[0] - x1_vec[1]*x2_vec[1]; // real part
            y_vec[1] += x1_vec[0]*x2_vec[1] + x1_vec[1]*x2_vec[0]; // imaginary part
        }

        /**
         * @brief Add the product of two floating point numbers "x1" and "x2" to "y".
         * Calling this function is completely equivalent to just writing "y += x1*x2".
         *
         * @tparam T the number type
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) addProd(const T x1, const T x2, T &y) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            y += x1*x2;
        }

        /**
         * @brief Add the product of two complex numbers "conj(x1)" and "x2" to "y".
         * This function is expected to work about 20% faster than normal operation such as "y += std::conj(x1)*x2".@n
         * (see also: addProd(const std::complex<T> &x1, const std::complex<T> &x2, std::complex<T> &y))
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) addConjProd(const std::complex<T> x1, const std::complex<T> x2, std::complex<T> &y) {
            const T *const x1_vec = reinterpret_cast<const T *>(&x1);
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            T *const y_vec = reinterpret_cast<T *>(&y);
            y_vec[0] += x1_vec[0]*x2_vec[0] + x1_vec[1]*x2_vec[1]; // real part
            y_vec[1] += x1_vec[0]*x2_vec[1] - x1_vec[1]*x2_vec[0]; // imaginary part
        }

        /**
         * @brief Subtract the product of two complex numbers "x1" and "x2" from "y"
         * This function is expected to work about 20% faster than normal operation such as "y -= x1*x2" (see also: addProd).
         *
         * @tparam T the number type of real and imaginary parts
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) subtractProd(const std::complex<T> x1, const std::complex<T> x2, std::complex<T> &y) {
            const T *const x1_vec = reinterpret_cast<const T *>(&x1);
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            T *const y_vec = reinterpret_cast<T *>(&y);
            y_vec[0] -= x1_vec[0]*x2_vec[0] - x1_vec[1]*x2_vec[1]; // real part
            y_vec[1] -= x1_vec[0]*x2_vec[1] + x1_vec[1]*x2_vec[0]; // imaginary part
        }

        /**
         * @brief Subtract the product of a floating point number "x1" and a complex number "x2" from "y"
         * This function is expected to work much faster than normal operation such as "y -= x1*x2" (see also: addProd).
         *
         * @tparam T the number type of "x1" and real and imaginary parts of "x2"
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) subtractProd(const T x1, const std::complex<T> x2, std::complex<T> &y) {
            const T *const x2_vec = reinterpret_cast<const T *>(&x2);
            T *const y_vec = reinterpret_cast<T *>(&y);
            y_vec[0] -=  x1*x2_vec[0];
            y_vec[1] -=  x1*x2_vec[1];
        }

        /**
         * @brief Subtract the product of two non-complex numbers "x1" and "x2" from "y".
         * Calling this function is completely equivalent to just writing "y -= x1*x2".
         *
         * @tparam T the number type
         * @param[in] x1 "x1"
         * @param[in] x2 "x2"
         * @param[inout] y "y"
         */
        template <typename T>
        inline static void __attribute__((always_inline)) subtractProd(const T x1, const T x2, T &y) {
            y -= x1*x2;
        }

        /**
         * @brief table-based fast square root calculation
         * @details Calculate the square root of "x" (sqrt(x)) using table lookup and linear interpolation.
         * This algorithm takes 2 parameters:
         * - max: the maximum value of "x"
         * - N: the number of the table entries
         * The 0-th table entry corresponds to 0, and the (N-1)-th corresponds to sqrt(x).
         * The approximation always be under estimation.
         * The error gets larger as "x" gets close to 0, so this method can be used only when error caused by small values of "x" is acceptable.
         * @note Modern CPUs in personal computers can perform floating point number division very fast, and `std::sqrt` is faster than this method.
         * This method wins against `std::sqrt` on only a few embedded CPUs.
         * Following is performance report:
         * - condition: max=3, N=100
         * - Intel Core i5-1035G7
         *   - iteration: 1e6
         *   - `std::sqrt`:       9 ms
         *   - `SqrtTable.calc`: 41 ms
         * - TMS320C6748
         *   - iteration: 1e4
         *   - `std::sqrt`:      10.1 ms
         *   - `SqrtTable.calc`: 4.83 ms
         *
         * @tparam T the number type of the input
         */
        template<typename T>
        class SqrtTable {
            static_assert(std::is_floating_point<T>::value, "Type parameter must be floating point number.");

            private:
                const T max, m_inv_max;
                const int m_N, m_Nm1;
                const T m_Nm1_T;
                std::vector<T> m_table;

            public:
                /**
                 * @brief constructor
                 *
                 * @param[in] max the maximum value of "x"
                 * @param[in] N the number of the table entries
                 */
                SqrtTable(const T max, const int N) : max(max), m_inv_max(1/max), m_N(N), m_Nm1(N-1), m_Nm1_T(N-1), m_table(N) {
                    #if MATH_LIB_ENABLE_CANARY_MODE
                        if (max <= 0 || N < 2) {
                            std::cerr << "BUG, FILE: " << __FILE__ << ", LINE: " << __LINE__ << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    #endif

                    const T sqrt_max = std::sqrt(max);
                    const T inv_Nm1 = 1/static_cast<T>(N-1);
                    for (int i=0; i<N; ++i) {
                        m_table[i] = sqrt_max*sqrt(static_cast<T>(i)*inv_Nm1); // approximates sqrt(max) * sqrt(x/max)
                    }
                }

                /**
                * @brief Calculate the square root of "x" (sqrt(x)) using table lookup.
                *
                * @param[in] x the input value
                * @return approximated sqrt(x)
                */
                T calc(const T x) {
                    const T idx_f = x*m_Nm1_T*m_inv_max;
                    const int idx = static_cast<int>(idx_f);
                    const int idx2 = idx + 1;
                    const bool cond1 = (0 <= idx), cond2 = (idx2 <= static_cast<int>(m_Nm1));
                    if (cond1 && cond2) {
                        const T w2 = idx_f - idx;
                        const T w1 = 1 - w2;
                        return w1*m_table[idx] + w2*m_table[idx2];
                    }
                    /* saturation */
                    else if (!cond1) {
                        return 0;
                    } else { // !cond2
                        return m_table[m_N-1];
                    }
                }
        };

        /**
         * @brief Calculates sin(x) using 9 degree-polynomial approximation.
         * "x" must be a floating point real number.
         * @details The 5 coefficients a1,a3,...,a9 were calculated as they minimize the cost function f(a1,a3,...,a9) := \int_0^1 (a1*x + a3*x^3 + ... + a9*x^9 - sin(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than 2.0*10^(-6)
         *
         * Performance meas result:
         *   CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB
         *   iteration = 1e7 (1e7 random samples from [-10pi, 10pi])
         *   std::sin: 352 ms
         *   sin_polyApprox: 130 ms
         *   max abs error: 1.89102e-06
         *
         * @tparam T the number type of the input value, x
         * @param[in] x input value
         * @return sin(x)
         */
        template <typename T>
        T sin_polyApprox(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T ONE = static_cast<T>(1);
            constexpr T PI = static_cast<T>(3.141592653589793);
            constexpr T INV_PI = ONE/PI;
            constexpr T HALF_PI = PI/2;
            constexpr T a1 = 1.0;
            constexpr T a3 = -0.166667;
            constexpr T a5 = 8.33296e-3;
            constexpr T a7 = -1.98048e-4;
            constexpr T a9 = 2.59811e-6;

            const T sign1 = std::copysign(ONE, x); // When x<0, negate the final result.
            const T x2 = std::abs(x);
            const int n = static_cast<int>(INV_PI*x2);
            const T sign2 = 1 - 2*(n&0b1); // When |x| is in odd pi-length interval, negate the final result.
            const T x3 = x2 - n*PI; // x3 in [0,pi)
            const T x4 = x3 <= HALF_PI ? x3 : PI - x3; // x4 in [0,pi/2]

            const T y = x4*x4;
            return sign1*sign2*x4*(a1 + y*(a3 + y*(a5 + y*(a7 + a9*y))));
        }

        /**
         * @brief Calculates cos(x) using sin_polyApprox(x+pi/2).
         * "x" must be a floating point real number.
         *
         * @tparam T the number type of the input value, x
         * @param[in] x input value
         * @return cos(x)
         */
        template <typename T>
        inline static T __attribute__((always_inline)) cos_polyApprox(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T HALF_PI = static_cast<T>(0.5*3.141592653589793);
            return sin_polyApprox(x + HALF_PI);
        }

        /**
         * @brief Calculates arc tangent of input value "x", using 7 degree-polynomial approximation.
         * "x" must be in the range [-1, 1], otherwise the calculation error increases.
         * @details The coefficients a1,a3,...,a7 were calculated as they minimize the cost function f(a1,a3,...,a7) := \int_0^1 (a1*x + a3*x^3 + ... + a7*x^7 - atan(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than 1.8*10^(-4) when 0<=x<=1.
         *
         * @tparam T the number type of the input value
         * @param[in] x input value
         * @return arc tangent of "x"
         */
        template <typename T>
        #if MATH_LIB_INLINE_AGGRESSIVELY
        inline static T __attribute__((always_inline))
        #else
        T
        #endif
        atan_polyApprox_deg7(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T a1 = 0.999298;
            constexpr T a3 = -0.322084;
            constexpr T a5 = 0.148508;
            constexpr T a7 = -0.0404899;
            const T y = x*x;
            return x*(a1 + y*(a3 + y*(a5 + a7*y)));
        }

        /**
         * @brief Calculates arc tangent of input value "x", using 9 degree-polynomial approximation.
         * "x" must be in the range [-1, 1], otherwise the calculation error increases.
         * @details The 5 coefficients a1,a3,...,a9 were calculated as they minimize the cost function f(a1,a3,...,a9) := \int_0^1 (a1*x + a3*x^3 + ... + a9*x^9 - atan(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than 2.5*10^(-5) when 0<=x<=1.
         *
         * Performance meas result:
         *   CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB
         *   iteration = 1e7 (1e7 random samples from unit disc (not include (0,0)))
         *   std::atan2: 861 ms
         *   atan2_polyApprox_deg9: 317 ms
         *   max abs error: 2.32363e-05
         *
         * @tparam T the number type of the input value
         * @param[in] x input value
         * @return arc tangent of "x"
         */
        template <typename T>
        #if MATH_LIB_INLINE_AGGRESSIVELY
        inline static T __attribute__((always_inline))
        #else
        T
        #endif
        atan_polyApprox_deg9(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T a1 = 0.99988;
            constexpr T a3 = -0.330534;
            constexpr T a5 = 0.181156;
            constexpr T a7 = -0.0866971;
            constexpr T a9 = 0.0216165;
            const T y = x*x;
            return x*(a1 + y*(a3 + y*(a5 + y*(a7 + a9*y))));
        }

        /* template to call `atan_polyApprox_deg7`, `atan_polyApprox_deg9` */
        template <typename T, int deg>
        class functor_atan_polyApprox;

        template <typename T>
        class functor_atan_polyApprox<T, 7> {
            public:
                inline static T __attribute__((always_inline)) doCalc(const T x) {
                    return atan_polyApprox_deg7(x);
                }
        };

        template <typename T>
        class functor_atan_polyApprox<T, 9> {
            public:
                inline static T __attribute__((always_inline)) doCalc(const T x) {
                    return atan_polyApprox_deg9(x);
                }
        };

        /**
         * @brief @brief Calculates arc tangent of input value "x", using n(= 7 or 9) degree-polynomial approximation.
         * "x" must be in the range [-1, 1], otherwise the calculation error increases.
         * @details The 1+(n-1)/2 coefficients a1,a3,..., an were calculated as they minimize the cost function f(a1,a3,...,an) := \int_0^1 (a1*x + a3*x^3 + ... + an*x^n - atan(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than E when 0<=x<=1 where E = 1.8*10^(-4) (n=7), 2.5*10^(-5) (n=9)
         *
         * Performance measurement result for n=9:
         *   CPU: Core 2 Quad Q9650, RAM: DDR2 800MHz 8GiB
         *   iteration = 1e7 (1e7 random samples from unit disc (not include (0,0)))
         *   std::atan2: 861 ms
         *   atan2_polyApprox_deg9: 317 ms
         *   max abs error: 2.32363e-05
         *
         * @tparam T the number type of the input value
         * @tparam deg The degree of polynomial used in approximation. this must be 7 or 9.
         * @return arc tangent of "x"
         */
        template <typename T, int deg>
        inline static T __attribute__((always_inline)) atan_polyApprox(const T x) {
            return functor_atan_polyApprox<T, deg>::doCalc(x);
        }

        /**
         * @brief Calculates "atan2(y,x)" using polynomial approximation (internally calls `atan_polyApprox` function).
         *
         * @tparam T the number type of the input value
         * @tparam deg The degree of polynomial used in approximation. this must be 7 or 9.
         * @param y "y"
         * @param x "x"
         * @return "atan2(y,x)"
         */
        template <typename T, int deg>
        T atan2_polyApprox(const T y, const T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            #if true
                constexpr T POSITIVE_ZERO = +0.0;
                constexpr T NEGATIVE_ZERO = -0.0;
                constexpr T pi = 3.14159265358979323846;
                constexpr T half_pi = pi/2;

                if (x == NEGATIVE_ZERO || x == POSITIVE_ZERO) {
                    if (y >= POSITIVE_ZERO) {
                        return half_pi;
                    }
                    if (y <= NEGATIVE_ZERO) {
                        return -half_pi;
                    }
                }
                if (x > POSITIVE_ZERO) {
                    if (-x <= y && y <= x) {
                        return atan_polyApprox<T, deg>(y/x);
                    }
                    if (y > x) {
                        return half_pi - atan_polyApprox<T, deg>(x/y);
                    } else { // y < -x
                        return -half_pi - atan_polyApprox<T, deg>(x/y);
                    }
                } else { // x < NEGATIVE_ZERO
                    if (POSITIVE_ZERO <= y && y <= -x) {
                        return pi + atan_polyApprox<T, deg>(y/x);
                    }
                    if (-x < y) {
                        return half_pi - atan_polyApprox<T, deg>(x/y);
                    }
                    if (x <= y && y <= NEGATIVE_ZERO) {
                        return atan_polyApprox<T, deg>(y/x) - pi;
                    } else { // y < x
                        return -half_pi - atan_polyApprox<T, deg>(x/y);
                    }
                }
            #else // Simpler but slower than the above code.
                static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
                constexpr T pi = 3.14159265358979323846;
                constexpr T half_pi = pi/2;
                const T abs_x = std::abs(x);
                const T abs_y = std::abs(y);

                if (x>=0 && abs_y<=x) { // region 1
                    return atan_polyApprox<T, deg>(y/x);
                }
                if (y>=0 && abs_x<=y) { // region 2
                    return half_pi + atan_polyApprox<T, deg>(-x/y);
                }
                if (x<=0 && abs_y<=abs_x) { // region 3
                    return atan_polyApprox<T, deg>(y/x) + (y>0 ? pi : -pi);
                }
                return atan_polyApprox<T, deg>(-x/y) - half_pi; // region 4
            #endif
        }
    }
}
