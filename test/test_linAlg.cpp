#include <complex>
#include <gtest/gtest.h>
#include "linAlg.hpp"

namespace {
    class LinAlgLibTest : public ::testing::Test{};

    TEST_F(LinAlgLibTest, isEqualMat__equal_case) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, isEqualMat__not_equal_case) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[2][3] = {
            {1,2,3},
            {4,5,7},
        };

        EXPECT_EQ(false, MotchyMathLib::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, isEqualVec__equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,3};

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, isEqualVec__not_equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,4};

        EXPECT_EQ(false, MotchyMathLib::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, transposeMat) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        std::complex<double> B[3][2];

        const std::complex<double> B_ans[3][2] = {
            {1,4},
            {2,5},
            {3,6},
        };

        MotchyMathLib::LinAlg::transposeMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(3, 2, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, conjugateMat) {
        const std::complex<double> A[2][3] = {
            {{1, 0.6}, {2, 0.5}, {3, 0.4}},
            {{4, -0.3}, {5, -0.2}, {6, -0.1}},
        };

        std::complex<double> B[2][3];

        const std::complex<double> B_ans[2][3] = {
            {{1, -0.6}, {2, -0.5}, {3, -0.4}},
            {{4, 0.3}, {5, 0.2}, {6, 0.1}},
        };

        MotchyMathLib::LinAlg::conjugateMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(2, 3, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, conjugateVec) {
        const std::complex<double> x[3] = {{1, 0.6}, {2, 0.5}, {3, 0.4}};
        std::complex<double> y[3];
        const std::complex<double> y_ans[3] = {{1, -0.6}, {2, -0.5}, {3, -0.4}};
        MotchyMathLib::LinAlg::conjugateVec(3, &x[0], &y[0]);
        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(3, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, addMat_inplace) {
        std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[2][3] = {
            {0.1, 0.2, 0.3},
            {-0.4, -0.5, -0.6},
        };

        const std::complex<double> A_ans[2][3] = {
            {1.1, 2.2, 3.3},
            {3.6, 4.5, 5.4},
        };

        MotchyMathLib::LinAlg::addMat_inplace(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(3, 2, &A[0][0], &A_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, addVec_inplace) {
        std::complex<double> x[3] = {1,2,3};
        const std::complex<double> y[3] = {0.1, 0.2, 0.3};
        const std::complex<double> y_ans[3] = {1.1, 2.2, 3.3};
        MotchyMathLib::LinAlg::addVec_inplace(3, &x[0], &y[0]);
        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(3, &x[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, scaleMat) {
        std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> a(1, 0.5);

        std::complex<double> B[2][3];

        const std::complex<double> B_ans[2][3] = {
            {{1, 0.5}, {2, 1}, {3, 1.5}},
            {{4, 2}, {5, 2.5}, {6, 3}},
        };

        MotchyMathLib::LinAlg::scaleMat(2, 3, a, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(3, 2, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, scaleVec) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> c(1, -0.5);

        std::complex<double> b[3];
        const std::complex<double> b_ans[3] = {{1, -0.5}, {2, -1}, {3, -1.5}};

        MotchyMathLib::LinAlg::scaleVec(3, c, &a[0], &b[0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(3, &b[0], &b_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, mulMat) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[3][2] = {
            {1,2},
            {3,-4},
            {-5,-6},
        };

        std::complex<double> C[2][2];

        const std::complex<double> C_ans[2][2] = {
            {-8, -24},
            {-11, -48},
        };

        MotchyMathLib::LinAlg::mulMat(2, 3, 2, &A[0][0], &B[0][0], &C[0][0]);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualMat(2, 2, &C[0][0], &C_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, hermitianInnerProduct_2_other_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[3] = {11, -10.5, 9};
        const std::complex<double> c = MotchyMathLib::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2, 1);
        const std::complex<double> c_ans = 17;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(LinAlgLibTest, hermitianInnerProduct_2_same_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[6] = {11, -10.5, 10, -9.5, 9, -8.5};
        const std::complex<double> c = MotchyMathLib::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2);
        const std::complex<double> c_ans = 58;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(LinAlgLibTest, l2Norm) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> c = MotchyMathLib::LinAlg::l2Norm(3, &a[0], 2);
        const std::complex<double> c_ans = 3.7416573867739413;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(LinAlgLibTest, solveLinearEquation) {
        const std::complex<double> A[3][3] = {
            {5, -4, 4},
            {5, 1, 4},
            {3, 4, 5},
        };

        const std::complex<double> b[3] = {25, 15, 10};
        std::complex<double> x[3];
        const std::complex<double> x_ans[3] = {1, -2, 3};
        char workspace[3*(3+1)*sizeof(std::complex<double>) + 3*sizeof(size_t)];

        MotchyMathLib::LinAlg::solveLinearEquation(3, &A[0][0], &b[0], &x[0], workspace);

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(3, &x[0], &x_ans[0], 1.0e-7l));
    }
}
