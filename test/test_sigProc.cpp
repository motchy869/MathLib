#include <complex>
#include <gtest/gtest.h>
#include "../include/linAlg.hpp"
#include "../include/sigProc.hpp"

namespace {
    class MySigProcLibTest : public ::testing::Test{};

    TEST_F(MySigProcLibTest, convolve_d_1_case1) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        const std::complex<double> x1[5] = {-0.88, 0.61, 2.56, -2.57, 1.84};
        const std::complex<double> x2[10] = {2.11, -1.46, -4.71, -0.48, -1.31, 2.97, 2.81, 1.89, -0.86, -0.24};
        std::complex<double> y[5+10-1];
        const std::complex<double> y_ans[5+10-1] = {-1.8568, 2.5719, 8.6558, -11.611, -3.563, 4.7768, -11.4475, 10.1376, -0.94, 2.7681, -2.0349, 5.0734, -0.9656, -0.4416};
        Utils::Math::SigProc::convolve(5, 10, &x1[0], &x2[0], &y[0]);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_d_1_case2) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        const std::complex<double> x1[10] = {4.03, -4.44, -4.65, -3.77, -3.87, -3.62, -4.5, 0.5, -3.46, -2.63};
        const std::complex<double> x2[5] = {-1.36, 3.74, -3.97, 4.64, -3.88};
        std::complex<double> y[10+5-1];
        const std::complex<double> y_ans[10+5-1] = {-5.4808, 21.1106, -26.2807, 24.0622, -26.6141, 1.0675, 8.4943, -6.4678, 22.6594, -18.183, 23.68, -7.5533, 1.2216, 10.2044};
        Utils::Math::SigProc::convolve(5, 10, &x1[0], &x2[0], &y[0]);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_c_m1_d_2_case1) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        constexpr size_t d=2;
        constexpr size_t N1 = 3, N2 = 5;
        constexpr size_t Ly = 1+(N1+N2-2)/d;
        const std::complex<double> x1[N1] = {1,2,3};
        const std::complex<double> x2[N2] = {11,12,13,14,15};
        std::complex<double> y[Ly];
        const std::complex<double> y_ans[Ly] = {11,70,82,45};
        Utils::Math::SigProc::convolve(N1, N2, &x1[0], &x2[0], &y[0], -1, d);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(Ly, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_c_m1_d_2_case2) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        constexpr size_t d=2;
        constexpr size_t N1 = 5, N2 = 3;
        constexpr size_t Ly = 1+(N1+N2-2)/d;
        const std::complex<double> x1[N1] = {11,12,13,14,15};
        const std::complex<double> x2[N2] = {1,2,3};
        std::complex<double> y[Ly];
        const std::complex<double> y_ans[Ly] = {11,70,82,45};
        Utils::Math::SigProc::convolve(N1, N2, &x1[0], &x2[0], &y[0], -1, d);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(Ly, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_c_m1_d_2_case3) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        constexpr size_t d=2;
        constexpr size_t N1 = 5, N2 = 1;
        constexpr size_t Ly = 1+(N1+N2-2)/d;
        const std::complex<double> x1[N1] = {11,12,13,14,15};
        const std::complex<double> x2[N2] = {1};
        std::complex<double> y[Ly];
        const std::complex<double> y_ans[Ly] = {11,13,15};
        Utils::Math::SigProc::convolve(N1, N2, &x1[0], &x2[0], &y[0], -1, d);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(Ly, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_c_2_d_2_case4) {
        /* These test case numeric values are created using Julia, and tested by DSP.jl package. */
        constexpr size_t c=2, d=2;
        constexpr size_t N1 = 5, N2 = 1;
        constexpr size_t Ly = c;
        const std::complex<double> x1[N1] = {11,12,13,14,15};
        const std::complex<double> x2[N2] = {1};
        std::complex<double> y[Ly];
        const std::complex<double> y_ans[Ly] = {11,13};
        Utils::Math::SigProc::convolve(N1, N2, &x1[0], &x2[0], &y[0], c, d);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(Ly, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MySigProcLibTest, convolve_type2) {
        constexpr size_t d = 2;
        constexpr size_t L1 = 3;
        constexpr size_t N = 2;
        const std::complex<double> x1[L1] = {1,2,3};
        const std::complex<double> x2[] = {11,12,13,14,15};
        std::complex<double> y[N] = {};
        const std::complex<double> y_ans[N] = {70, 82};
        Utils::Math::SigProc::convolve_type2(L1, x1, &x2[2], y, N, d);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(N, &y[0], &y_ans[0], 1.0e-7l));
    }
}
