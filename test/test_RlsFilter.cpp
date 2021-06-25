#include <array>
#include <iostream>
#include <gtest/gtest.h>
#include "../include/RlsFilter.hpp"

namespace {
    template <typename T>
    static void printComplexArray(const std::complex<T> *array, int N, const char *format) {
        assert(N >= 1);
        for (int n=0; n<N; ++n) {
            printf(format, array[n].real(), array[n].imag());
            if (n == N-1) {
                printf("\n");
            } else {
                printf(", ");
            }
        }
    }

    class MyRlsFilterTest : public ::testing::Test{};

    TEST_F(MyRlsFilterTest, case1_thinningRate_1) {
        const std::array<std::complex<double>, 10> vec_x = {{
            {0.06412935, 0.88332087}, {0.99244617, 0.12382741}, {0.64271018, 0.45058685},
            {0.99968758, 0.67847759}, {0.94431826, 0.48793753}, {0.4698219 , 0.01647779},
            {0.93984455, 0.37076489}, {0.01225688, 0.67889914}, {0.43928084, 0.22816118},
            {0.66261685, 0.13262789}
        }};

        const std::array<std::complex<double>, 10> vec_y = {{
            {0.06412935, 0.88332087}, {1.00847851, 0.34465762}, {0.89883789, 0.59195881},
            {1.2844209 , 0.80660273}, {1.27457893, 0.71388028}, {0.83086241, 0.22327187},
            {1.17533981, 0.43587653}, {0.30594575, 0.77365009}, {0.55982563, 0.44423158},
            {0.77396917, 0.27453058},
        }};

        constexpr int p = 2;
        constexpr int thinningRate = 1;
        RlsFilter<std::complex<double>> rlsFilter(p, thinningRate);

        std::array<std::complex<double>, p+1> w_opt;
        rlsFilter.trainHard(vec_x.data(), vec_y.data(), vec_x.size(), w_opt.data());

        const std::array<std::complex<double>, p+1> w_opt_ans = {{
            {1.0036556, -0.0025984}, {-0.2384301, 0.0000261}, {-0.0427164, 0.0031466}
        }};
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(w_opt.size(), w_opt.data(), w_opt_ans.data(), 1.0e-7l));
    }
}