#include <complex>
#include <gtest/gtest.h>
#include "../include/linAlg.hpp"

namespace {
    class MyLinAlgLibTest : public ::testing::Test{};

    TEST_F(MyLinAlgLibTest, isEqualMat__equal_case) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, isEqualMat__not_equal_case) {
        const std::complex<double> A[2][3] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> B[2][3] = {
            {1,2,3},
            {4,5,7},
        };

        EXPECT_EQ(false, Utils::Math::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, isEqualVec__equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,3};

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, isEqualVec__not_equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,4};

        EXPECT_EQ(false, Utils::Math::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, transposeMat) {
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

        Utils::Math::LinAlg::transposeMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(3, 2, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, conjugateMat) {
        const std::complex<double> A[2][3] = {
            {{1, 0.6}, {2, 0.5}, {3, 0.4}},
            {{4, -0.3}, {5, -0.2}, {6, -0.1}},
        };

        std::complex<double> B[2][3];

        const std::complex<double> B_ans[2][3] = {
            {{1, -0.6}, {2, -0.5}, {3, -0.4}},
            {{4, 0.3}, {5, 0.2}, {6, 0.1}},
        };

        Utils::Math::LinAlg::conjugateMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(2, 3, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, conjugateVec) {
        const std::complex<double> x[3] = {{1, 0.6}, {2, 0.5}, {3, 0.4}};
        std::complex<double> y[3];
        const std::complex<double> y_ans[3] = {{1, -0.6}, {2, -0.5}, {3, -0.4}};
        Utils::Math::LinAlg::conjugateVec(3, &x[0], &y[0]);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, addMat_inplace) {
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

        Utils::Math::LinAlg::addMat_inplace(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(3, 2, &A[0][0], &A_ans[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, addVec_inplace) {
        std::complex<double> x[3] = {1,2,3};
        const std::complex<double> y[3] = {0.1, 0.2, 0.3};
        const std::complex<double> y_ans[3] = {1.1, 2.2, 3.3};
        Utils::Math::LinAlg::addVec_inplace(3, &x[0], &y[0]);
        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &x[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, scaleMat) {
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

        Utils::Math::LinAlg::scaleMat(2, 3, a, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(3, 2, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, scaleVec) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> c(1, -0.5);

        std::complex<double> b[3];
        const std::complex<double> b_ans[3] = {{1, -0.5}, {2, -1}, {3, -1.5}};

        Utils::Math::LinAlg::scaleVec(3, c, &a[0], &b[0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &b[0], &b_ans[0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, mulMat) {
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

        Utils::Math::LinAlg::mulMat(2, 3, 2, &A[0][0], &B[0][0], &C[0][0]);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualMat(2, 2, &C[0][0], &C_ans[0][0], 1.0e-7l));
    }

    TEST_F(MyLinAlgLibTest, hermitianInnerProduct_2_other_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[3] = {11, -10.5, 9};
        const std::complex<double> c = Utils::Math::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2, 1);
        const std::complex<double> c_ans = 17;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(MyLinAlgLibTest, hermitianInnerProduct_2_same_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[6] = {11, -10.5, 10, -9.5, 9, -8.5};
        const std::complex<double> c = Utils::Math::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2);
        const std::complex<double> c_ans = 58;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(MyLinAlgLibTest, l2Norm) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> c = Utils::Math::LinAlg::l2Norm(3, &a[0], 2);
        const std::complex<double> c_ans = 3.7416573867739413;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(MyLinAlgLibTest, solveLinearEquation) {
        const std::complex<double> A[3][3] = {
            {5, -4, 4},
            {5, 1, 4},
            {3, 4, 5},
        };

        const std::complex<double> b[3] = {25, 15, 10};
        std::complex<double> x[3];
        const std::complex<double> x_ans[3] = {1, -2, 3};
        char workspace[3*(3+1)*sizeof(std::complex<double>) + 3*sizeof(size_t)];

        Utils::Math::LinAlg::solveLinearEquation(3, &A[0][0], &b[0], &x[0], workspace);

        EXPECT_EQ(true, Utils::Math::LinAlg::isEqualVec(3, &x[0], &x_ans[0], 1.0e-7l));
    }
}
