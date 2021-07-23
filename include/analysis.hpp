#ifndef __MOTCHY_MATH_LIB__ANALYSIS_HPP__
#define __MOTCHY_MATH_LIB__ANALYSIS_HPP__

#include <cmath>
#include <type_traits>

namespace MotchyMathLib {
    namespace Analysis {
        /**
         * @brief Calculates sin(x) using 9 degree-polynomial approximation.
         * "x" must be a floating point real number.
         * @details The 5 coefficients a1,a3,...,a9 were calculated as they minimize the cost function f(a1,a3,...,a9) := \int_0^1 (a1*x + a3*x^3 + ... + a9*x^9 - sin(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than 2.0*10^(-6)
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
            constexpr T HALF_PI = PI/2;
            constexpr T a1 = 1.0;
            constexpr T a3 = -0.166667;
            constexpr T a5 = 8.33296e-3;
            constexpr T a7 = -1.98048e-4;
            constexpr T a9 = 2.59811e-6;

            const T sign1 = std::copysign(ONE, x);
            const T x2 = std::abs(x);
            const int n = static_cast<int>(x2/PI);
            const T sign2 = 1 - 2*(n&0b1);
            const T x3 = x2 - n*PI;
            const T x4 = x3 <= HALF_PI ? x3 : PI - x3;

            const T y = x4*x4;
            return sign1*sign2*x4*(a1 + y*(a3 + y*(a5 + y*(a7 + a9*y))));
        }

        /**
         * @brief Calculates arc tangent of input value "x", using 9 degree-polynomial approximation.
         * "x" must be in the range [-pi/4, pi/4], otherwise the calculation error increases.
         * @details The 5 coefficients a1,a3,...,a9 were calculated as they minimize the cost function f(a1,a3,...,a9) := \int_0^1 (a1*x + a3*x^3 + ... + a9*x^9 - atan(x))^2 \mathrm{d}x.
         * The maximal absolute error is less than 2.5*10^(-5) when 0<=x<=1.
         *
         * @tparam T the number type of the input value
         * @param[in] x input value
         * @return arc tangent of "x"
         */
        template <typename T>
        T atan_polyApprox(T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T a1 = 0.99988;
            constexpr T a3 = -0.330534;
            constexpr T a5 = 0.181156;
            constexpr T a7 = -0.0866971;
            constexpr T a9 = 0.0216165;
            const T y = x*x;
            return x*(a1 + y*(a3 + y*(a5 + y*(a7 + a9*y))));
        }

        /**
         * @brief Calculates "atan2(y,x)" using polynomial approximation (internally calls `atan_polyApprox` function).
         *
         * @tparam T the number type of the input value
         * @param y "y"
         * @param x "x"
         * @return "atan2(y,x)"
         */
        template <typename T>
        T atan2_polyApprox(T y, T x) {
            static_assert(std::is_floating_point<T>::value, "argument type must be floating point number.");
            constexpr T POSITIVE_ZERO = +0.0;
            constexpr T NEGATIVE_ZERO = -0.0;
            constexpr T pi = 3.14159265358979323846;

            if (x == NEGATIVE_ZERO || x == POSITIVE_ZERO) {
                if (y >= POSITIVE_ZERO) {
                    return pi/2;
                }
                if (y <= NEGATIVE_ZERO) {
                    return -pi/2;
                }
            }
            if (x > POSITIVE_ZERO) {
                if (-x <= y && y <= x) {
                    return atan_polyApprox(y/x);
                }
                if (y > x) {
                    return pi/2 - atan_polyApprox(x/y);
                } else { // y < -x
                    return -pi/2 - atan_polyApprox(x/y);
                }
            } else { // x < NEGATIVE_ZERO
                if (POSITIVE_ZERO <= y && y <= -x) {
                    return pi + atan_polyApprox(y/x);
                }
                if (-x < y) {
                    return pi/2 - atan_polyApprox(x/y);
                }
                if (x <= y && y <= NEGATIVE_ZERO) {
                    return atan_polyApprox(y/x) - pi;
                } else { // y < x
                    return -pi/2 - atan_polyApprox(x/y);
                }
            }
        }
    }
}

#endif // __MOTCHY_MATH_LIB__ANALYSIS_HPP__
