#include <complex>
#include <gtest/gtest.h>
#include "../include/linAlg.hpp"

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

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
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

        EXPECT_EQ(false, MathLib::LinAlg::isEqualMat(2, 3, &A[0][0], &B[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, isEqualVec__equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,3};

        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, isEqualVec__not_equal_case) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> b[3] = {1,2,4};

        EXPECT_EQ(false, MathLib::LinAlg::isEqualVec(3, &a[0], &b[0], 1.0e-7l));
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

        MathLib::LinAlg::transposeMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(3, 2, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, conjugateMat_inplace) {
        std::complex<double> A[2][3] = {
            {{1, 0.6}, {2, 0.5}, {3, 0.4}},
            {{4, -0.3}, {5, -0.2}, {6, -0.1}},
        };

        const std::complex<double> conj_A_ans[2][3] = {
            {{1, -0.6}, {2, -0.5}, {3, -0.4}},
            {{4, 0.3}, {5, 0.2}, {6, 0.1}},
        };

        MathLib::LinAlg::conjugateMat(2, 3, &A[0][0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(2, 3, &A[0][0], &conj_A_ans[0][0], 1.0e-7l));
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

        MathLib::LinAlg::conjugateMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(2, 3, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, conjugateVec_inplace) {
        std::complex<double> x[3] = {{1, 0.6}, {2, 0.5}, {3, 0.4}};
        const std::complex<double> conj_x_ans[3] = {{1, -0.6}, {2, -0.5}, {3, -0.4}};
        MathLib::LinAlg::conjugateVec(3, &x[0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &x[0], &conj_x_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, conjugateVec) {
        const std::complex<double> x[3] = {{1, 0.6}, {2, 0.5}, {3, 0.4}};
        std::complex<double> y[3];
        const std::complex<double> y_ans[3] = {{1, -0.6}, {2, -0.5}, {3, -0.4}};
        MathLib::LinAlg::conjugateVec(3, &x[0], &y[0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &y[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, setDiag) {
        constexpr size_t m = 3;
        int A[m][m] = {}; // initialize with 0
        const int d[m] = {1,2,3};
        const int A_true[m][m] = {
            {1,0,0},
            {0,2,0},
            {0,0,3}
        };
        MathLib::LinAlg::setDiag(m, d, &A[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, m, &A[0][0], &A_true[0][0], 0));
    }

    TEST_F(LinAlgLibTest, fillLowTri_lowSubDiag) {
        constexpr size_t m = 5;
        int A[m][m] = {}; // initialize with 0
        const int A_true[m][m] = {
            {0,0,0,0,0},
            {0,0,0,0,0},
            {7,0,0,0,0},
            {7,7,0,0,0},
            {7,7,7,0,0},
        };
        MathLib::LinAlg::fillLowTri(m, &A[0][0], 7, -2);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, m, &A[0][0], &A_true[0][0], 0));
    }

    TEST_F(LinAlgLibTest, fillLowTri_upSubDiag) {
        constexpr size_t m = 5;
        int A[m][m] = {}; // initialize with 0
        const int A_true[m][m] = {
            {7,7,7,0,0},
            {7,7,7,7,0},
            {7,7,7,7,7},
            {7,7,7,7,7},
            {7,7,7,7,7},
        };
        MathLib::LinAlg::fillLowTri(m, &A[0][0], 7, 2);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, m, &A[0][0], &A_true[0][0], 0));
    }

    TEST_F(LinAlgLibTest, addDiag) {
        constexpr size_t m = 3;
        int A[m][m] = {
            {1,2,3},
            {4,5,6},
            {7,8,9}
        };
        const int d[m] = {1,-2,3};
        const int A_true[m][m] = {
            {2,2,3},
            {4,3,6},
            {7,8,12}
        };
        MathLib::LinAlg::addDiag(m, d, &A[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, m, &A[0][0], &A_true[0][0], 0));
    }

    TEST_F(LinAlgLibTest, dropSubMat) {
        constexpr size_t m = 3, n = 4;
        const int A[m][n] = {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9,10,11,12},
        };
        constexpr size_t r1=1, r2=1, c1=1, c2=2;
        constexpr size_t m2=m-1, n2=n-2;
        const int B_true[m2][n2] = {
            {1,4},
            {9,12},
        };
        int B[m2][n2];
        MathLib::LinAlg::dropSubMat(m, n, r1, r2, c1, c2, &A[0][0], &B[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m2, n2, &B[0][0], &B_true[0][0], 0));
    }

    TEST_F(LinAlgLibTest, addMat) {
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

        MathLib::LinAlg::addMat(2, 3, &A[0][0], &B[0][0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(3, 2, &A[0][0], &A_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, addVec) {
        std::complex<double> x[3] = {1,2,3};
        const std::complex<double> y[3] = {0.1, 0.2, 0.3};
        const std::complex<double> y_ans[3] = {1.1, 2.2, 3.3};
        MathLib::LinAlg::addVec(3, &x[0], &y[0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &x[0], &y_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, scaleMat) {
        const size_t m=2, n=3;
        const std::complex<double> A[m][n] = {
            {1,2,3},
            {4,5,6},
        };

        const std::complex<double> a(1, 0.5);

        std::complex<double> B[m][n];

        const std::complex<double> B_ans[m][n] = {
            {{1, 0.5}, {2, 1}, {3, 1.5}},
            {{4, 2}, {5, 2.5}, {6, 3}},
        };

        MathLib::LinAlg::scaleMat(m, n, a, &A[0][0], &B[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, n, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, scaleMatEachRow) {
        const size_t m=2, n=3;
        const float A[m][n] = {
            {1,2,3},
            {4,5,6},
        };
        const float c[m] = {2,-1};
        float B[m][n];
        const float B_ans[m][n] = {
            {2,4,6},
            {-4,-5,-6},
        };
        MathLib::LinAlg::scaleMatEachRow(2, 3, c, &A[0][0], &B[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(m, n, &B[0][0], &B_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, scaleVec) {
        const std::complex<double> a[3] = {1,2,3,};

        const std::complex<double> c(1, -0.5);

        std::complex<double> b[3];
        const std::complex<double> b_ans[3] = {{1, -0.5}, {2, -1}, {3, -1.5}};

        MathLib::LinAlg::scaleVec(3, c, &a[0], &b[0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &b[0], &b_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, vecSelfOuterProd_real) {
        const size_t M = 4;
        const float x[M] = {1, -2, 3, -4};
        float X[M][M];
        const float X_true[M][M] = {
            { 1,  -2,    3,   -4},
            {-2,   4,   -6,    8},
            { 3,  -6,    9,  -12},
            {-4,   8,  -12,   16},
        };
        MathLib::LinAlg::vecSelfOuterProd(M, x, &X[0][0]);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(M, M, &X[0][0], &X_true[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, vecSelfOuterProd_complex) {
        const size_t M = 4;
        const std::complex<float> x[M] = {{1, 11}, {-2, -12}, {3, 13}, {-4, -14}};
        std::complex<float> X[M][M];
        const std::complex<float> X_true[M][M] = {
            {{ 122,  +0}, {-134, -10}, { 146, +20}, {-158, -30}},
            {{-134, +10}, { 148,  +0}, {-162, -10}, { 176, +20}},
            {{ 146, -20}, {-162, +10}, { 178,  +0}, {-194, -10}},
            {{-158, +30}, { 176, -20}, {-194, +10}, { 212,  +0}},
        };
        std::complex<float> workspace[M];
        MathLib::LinAlg::vecSelfOuterProd(M, x, &X[0][0], workspace);
        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(M, M, &X[0][0], &X_true[0][0], 1.0e-7l));
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

        MathLib::LinAlg::mulMat(2, 3, 2, &A[0][0], &B[0][0], &C[0][0]);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualMat(2, 2, &C[0][0], &C_ans[0][0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, innerProd) {
        const size_t m = 5;
        const std::complex<float> x[m]= {{0.293476, -0.211899}, {0.0108039, -0.478889}, {0.350097, -0.40511}, {0.0391371, 0.156917}, {-0.494805, 0.383847}};
        const std::complex<float> y[m]= {{-0.00120354, -0.407946}, {-0.113012, 0.285927}, {0.35409, -0.182355}, {-0.332509, -0.0444678}, {0.035946, 0.0887727}};
        const std::complex<float> ip_ans = std::complex<float>(0.04110481688090317, -0.35358968160291226);
        const std::complex<float> ip = MathLib::LinAlg::innerProd(m, x, y);
        EXPECT_EQ(true, std::abs(ip-ip_ans) < 1.0e-6f);
    }

    TEST_F(LinAlgLibTest, hermitianInnerProduct_2_other_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[3] = {11, -10.5, 9};
        const std::complex<double> c = MathLib::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2, 1);
        const std::complex<double> c_ans = 17;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(LinAlgLibTest, hermitianInnerProduct_2_same_strides) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> b[6] = {11, -10.5, 10, -9.5, 9, -8.5};
        const std::complex<double> c = MathLib::LinAlg::hermitianInnerProduct(3, &a[0], &b[0], 2);
        const std::complex<double> c_ans = 58;
        EXPECT_EQ(true, std::abs(c_ans-c) < 1.0e-7l);
    }

    TEST_F(LinAlgLibTest, l2Norm) {
        const std::complex<double> a[6] = {1, -1, 2, -1.5, 3, -2};
        const std::complex<double> c = MathLib::LinAlg::l2Norm(3, &a[0], 2);
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

        MathLib::LinAlg::solveLinearEquation(3, &A[0][0], &b[0], &x[0], workspace);

        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(3, &x[0], &x_ans[0], 1.0e-7l));
    }

    TEST_F(LinAlgLibTest, ldlDecomp_real) {
        constexpr size_t m = 5;
        const float A[m][m] = {
            { 0.59198756,  0.5907925 , -0.12247961,  0.0801843 , -0.40894139},
            { 0.5907925 , -0.86382306, -0.62571074,  0.11681897, -0.44449101},
            {-0.12247961, -0.62571074,  0.97233429,  0.42309137,  0.90466214},
            { 0.0801843 ,  0.11681897,  0.42309137,  0.54580854, -0.33315596},
            {-0.40894139, -0.44449101,  0.90466214, -0.33315596, -0.20125677},
        };
        float workspace[m];

        /* tested on Numpy */
        const float d_true[m] = {0.59198756, -1.45342291,  1.12140311,  0.37333903, -2.05156301};
        const float L_true[m][m] = {
            { 1.        ,  0.        ,  0.        ,  0.        ,  0.        },
            { 0.99798128,  1.        ,  0.        ,  0.        ,  0.        },
            {-0.20689557,  0.34640873,  1.        ,  0.        ,  0.        },
            { 0.13544929, -0.02531716,  0.38071457,  1.        ,  0.        },
            {-0.69079389,  0.02502723,  0.74251145, -1.59557335,  1.        },
        };

        float d[m] = {};
        float L[m][m] = {};
        const bool noZeroDiv = MathLib::LinAlg::ldlDecomp(m, &A[0][0], d, &L[0][0], workspace);
        const bool isDiagOk = MathLib::LinAlg::isEqualVec(m, d, d_true, 1e-5);
        const bool isLOk = MathLib::LinAlg::isEqualMat(m, m, &L[0][0], &L_true[0][0], 1e-5);

        EXPECT_EQ(true, noZeroDiv && isDiagOk && isLOk);
    }

    TEST_F(LinAlgLibTest, ldlDecomp_complex) {
        constexpr size_t m = 5;
        const float A_real[m][m] = {
            {-0.14583693, -0.07703409,  0.41802944,  0.12720957,  0.7996055 },
            {-0.07703409, -0.95275531,  0.05454494, -0.13649376,  0.66307892},
            { 0.41802944,  0.05454494,  0.67408461, -0.17179905, -0.16748154},
            { 0.12720957, -0.13649376, -0.17179905, -0.16333842,  0.31039082},
            { 0.7996055 ,  0.66307892, -0.16748154,  0.31039082, -0.31746929},
        };
        const float A_imag[m][m] = {
            { 0.        , -0.55608743, -0.91554274, -0.01458005,  0.00506915},
            { 0.55608743,  0.        , -0.2074604 , -0.20955483,  0.63859046},
            { 0.91554274,  0.2074604 ,  0.        ,  0.19562549, -0.51800175},
            { 0.01458005,  0.20955483, -0.19562549,  0.        , -0.26862068},
            {-0.00506915, -0.63859046,  0.51800175,  0.26862068,  0.        },
        };
        std::complex<float> workspace[m];

        /* tested on Numpy */
        const float d_true[m] = {-0.14583693, 1.20833987, -4.42254581, -0.1202024 , 3.13005036};
        const float L_true_real[m][m] = {
            { 1.        ,  0.        ,  0.        ,  0.        ,  0.        },
            { 0.5282207 ,  1.        ,  0.        ,  0.        ,  0.        },
            {-2.86641689,  2.75151731,  1.        ,  0.        ,  0.        },
            {-0.87227267, -0.12255959, -0.05732815,  1.        ,  0.        },
            {-5.48287396,  0.18321145,  0.95407657,  0.03473893,  1.        },
        };
        const float L_true_imag[m][m] = {
            { 0.        ,  0.        ,  0.        ,  0.        ,  0.        },
            {-3.81307695,  0.        ,  0.        ,  0.        ,  0.        },
            {-6.27785256, -1.54768265,  0.        ,  0.        ,  0.        },
            {-0.099975  , -0.23437654, -0.01266549,  0.        ,  0.        },
            { 0.034759  , -3.04953122, -1.19388732, -0.30690583,  0.        },
        };
        std::complex<float> L_true[m][m];
        MathLib::LinAlg::complexMat(m, m, &L_true_real[0][0], &L_true_imag[0][0], &L_true[0][0]);

        std::complex<float> A[m][m];
        MathLib::LinAlg::complexMat(m, m, &A_real[0][0], &A_imag[0][0], &A[0][0]);

        float d[m] = {};
        std::complex<float> L[m][m] = {};
        const bool noZeroDiv = MathLib::LinAlg::ldlDecomp(m, &A[0][0], d, &L[0][0], workspace);
        const bool isDiagOk = MathLib::LinAlg::isEqualVec(m, d, d_true, 1e-5);
        const bool isLOk = MathLib::LinAlg::isEqualMat(m, m, &L[0][0], &L_true[0][0], 1e-5);

        EXPECT_EQ(true, noZeroDiv && isDiagOk && isLOk);
    }

    TEST_F(LinAlgLibTest, solveLinEqHermitian_real) {
        constexpr size_t m = 5;

        const float A[m][m] = {
            { 0.15582773, -0.21069345,  0.15042193, -0.46729612,  0.00596341},
            {-0.21069345, -0.53092562, -0.044495  ,  0.39701956, -0.19034964},
            { 0.15042193, -0.044495  , -0.53684812,  0.07480983,  0.2489673 },
            {-0.46729612,  0.39701956,  0.07480983,  0.83141263,  0.04988819},
            { 0.00596341, -0.19034964,  0.2489673 ,  0.04988819, -0.91536905},
        };
        const float b[m] = {0.72965851, 0.30325472, 0.40853002, 0.51387289, 0.19281199};
        const float x_true[m] = {-74.32784309, 6.42575663, -33.99366676, -40.35903688, -13.47647763}; // tested on Numpy

        float x[m];
        char workSpace[(m*m+2*m-1)*sizeof(float)];

        const bool noZeroDiv = MathLib::LinAlg::solveLinEqHermitian(m, &A[0][0], b, x, workSpace);
        const bool isSolutionOk = MathLib::LinAlg::isEqualVec(m, x, x_true, 1e-3);

        EXPECT_EQ(true, noZeroDiv && isSolutionOk);
    }

    TEST_F(LinAlgLibTest, solveLinEqHermitian_complex) {
        constexpr size_t m = 5;

        const float A_real[m][m] = {
            { 0.38934498,  0.3323428 , -0.17064305, -0.93412886,  0.32600936},
            { 0.3323428 , -0.27048964,  0.29548155,  0.35787475, -0.13053154},
            {-0.17064305,  0.29548155,  0.2378897 ,  0.51834123,  0.29144506},
            {-0.93412886,  0.35787475,  0.51834123,  0.08597033,  0.04272254},
            { 0.32600936, -0.13053154,  0.29144506,  0.04272254, -0.28545381},
        };
        const float A_imag[m][m] = {
            { 0.        ,  0.60230465, -0.47143342,  0.08291764,  0.72430698},
            {-0.60230465,  0.        ,  0.20668524, -0.72295811, -0.05486585},
            { 0.47143342, -0.20668524,  0.        , -0.14381023,  0.7543405 },
            {-0.08291764,  0.72295811,  0.14381023,  0.        , -0.10504604},
            {-0.72430698,  0.05486585, -0.7543405 ,  0.10504604,  0.        },
        };
        std::complex<float> A[m][m];
        MathLib::LinAlg::complexMat(m, m, &A_real[0][0], &A_imag[0][0], &A[0][0]);

        const float b_real[m] = {0.46009212, -0.17479708, 0.062278, 0.35397522, 0.46410087};
        const float b_imag[m] = {0.29789619, 0.29789619, 0.29789619, 0.29789619, 0.29789619};
        std::complex<float> b[m];
        MathLib::LinAlg::complexVec(m, b_real, b_imag, b);

        /* tested on Numpy */
        const float x_true_real[m] = {-0.12165307, 0.21805922, 0.06650101, 0.08355123, 0.41089338};
        const float x_true_imag[m] = {0.12005935, -0.24947381, 0.69593837, -0.00897475, 0.12241331};
        std::complex<float> x_true[m];
        MathLib::LinAlg::complexVec(m, x_true_real, x_true_imag, x_true);

        std::complex<float> x[m];
        char workSpace[m*sizeof(float) + (m*m+m-1)*sizeof(std::complex<float>)];

        const bool noZeroDiv = MathLib::LinAlg::solveLinEqHermitian(m, &A[0][0], b, x, workSpace);
        const bool isSolutionOk = MathLib::LinAlg::isEqualVec(m, x, x_true, 1e-3);

        EXPECT_EQ(true, noZeroDiv && isSolutionOk);
    }
}
