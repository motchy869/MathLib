#include <algorithm>
#include <gtest/gtest.h>
#include "../include/analysis.hpp"
#include "../include/linAlg.hpp"

namespace {
    class AnalysisLibTest : public ::testing::Test{};

    TEST_F(AnalysisLibTest, addProd) {
        /* test case produced using Julia */
        constexpr int m=5;
        const std::complex<float> vec_x[m] = {{0.734619, 0.686623}, {0.683544, 0.388953}, {0.101859, 0.298419}, {0.288473, 0.0104381}, {0.952811, 0.546041}};
        const std::complex<float> sum_ans = {0.48406428, 1.2532825};

        std::complex<float> sum = 0;
        for (int i=0; i<m-1; ++i) {
            MathLib::Analysis::addProd(vec_x[i], vec_x[i+1], sum);
        }

        EXPECT_EQ(true, std::abs(sum - sum_ans) < 1e-6);
    }

    TEST_F(AnalysisLibTest, subtractProd_complex_complex) {
        /* test case produced using Julia */
        constexpr int m=5;
        const std::complex<float> vec_x[m] = {{0.259803, 0.960717}, {0.203663, 0.538169}, {0.358519, 0.981518}, {0.951527, 0.551323}, {0.0274041, 0.460809}};
        const std::complex<float> minus_sum_ans = {1.3472931, -2.3135047};
        std::complex<float> minus_sum = 0;
        for (int i=0; i<m-1; ++i) {
            MathLib::Analysis::subtractProd(vec_x[i], vec_x[i+1], minus_sum);
        }
        EXPECT_EQ(true, std::abs(minus_sum - minus_sum_ans) < 1e-6);
    }

    TEST_F(AnalysisLibTest, subtractProd_float_complex) {
        constexpr int m=5;
        const float vec_x1[m] = {0.309286, 0.660707, 0.620851, 0.448405, 0.886714};
        const std::complex<float> vec_x2[m] = {{0.259803, 0.960717}, {0.203663, 0.538169}, {0.358519, 0.981518}, {0.951527, 0.551323}, {0.0274041, 0.460809}};
        const std::complex<float> minus_sum_ans = {-0.8884709436304, -1.9179065568040001};
        std::complex<float> minus_sum = 0;
        for (int i=0; i<m; ++i) {
            MathLib::Analysis::subtractProd(vec_x1[i], vec_x2[i], minus_sum);
        }
        EXPECT_EQ(true, std::abs(minus_sum - minus_sum_ans) < 1e-6);
    }

    TEST_F(AnalysisLibTest, SqrtTable) {
        const float max = 3;
        const int N = 100;
        MathLib::Analysis::SqrtTable<float> sqrtTable(max, N);
        EXPECT_EQ(0, sqrtTable.calc(-1));
        EXPECT_EQ(true, (sqrtTable.calc(max+1) - sqrt(max)) < 1e-6);

        const int N2 = 1000;
        std::vector<float> errors(N2);
        for (int i=1; i<N2; ++i) {
            const float x = max*static_cast<float>(i)/(N2-1);
            errors[i] = sqrtTable.calc(x) - std::sqrt(x);
        }
        EXPECT_EQ(true, abs(*std::max_element(errors.begin(), errors.end(), [](const float &a, const float &b) { return std::abs(a) < std::abs(b); })) < 0.05);
    }

    TEST_F(AnalysisLibTest, sin_polyApprox) {
        constexpr size_t N = 20;
        const double vec_x[N] = {-5.422284533254419, -8.368253325506895, 7.230370357282717, 8.540266743128377, 8.796032571059138, -9.0058423693436, -4.1458323100147165, 6.6625543111462315, 4.714898106104982, -6.991907637267804, -4.039460886872965, -8.222201344331197, -6.5703506066193365, 1.5968050614449092, -2.7632828820450843, 6.233630363227799, -1.805797535337036, 4.069744312298307, -8.991169132560007, -9.276178027269587};
        const double ans[N] = {0.7584299540206195, -0.8706512087514109, 0.8117748754722436, 0.7736053539586267, 0.5881305349111817, -0.4067883225567331, 0.843754111510663, 0.3703344119084007, -0.9999968521457114, -0.6508643030941429, 0.7820000033037139, -0.9329696480010395, -0.2832347454042438, 0.999661791926789, -0.3693502849318244, -0.04953466449150098, -0.972514059931373, -0.800513570965066, -0.4201483819628514, -0.14805364120001024};
        double result[N];
        for (size_t i=0; i<N; ++i) {
            result[i] = MathLib::Analysis::sin_polyApprox(vec_x[i]);
        }

        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(N, &result[0], &ans[0], 2.0e-6l));
    }

    TEST_F(AnalysisLibTest, cos_polyApprox) {
        constexpr size_t N = 20;
        const double vec_x[N] = {-5.422284533254419, -8.368253325506895, 7.230370357282717, 8.540266743128377, 8.796032571059138, -9.0058423693436, -4.1458323100147165, 6.6625543111462315, 4.714898106104982, -6.991907637267804, -4.039460886872965, -8.222201344331197, -6.5703506066193365, 1.5968050614449092, -2.7632828820450843, 6.233630363227799, -1.805797535337036, 4.069744312298307, -8.991169132560007, -9.276178027269587};
        const double ans[N] = {0.651754558744533, -0.49190087690479595, 0.5839705057209854, -0.633667701817402, -0.8087660192571688, -0.9135224467036808, -0.5367299128135599, 0.9288985000291792, 0.002509123087504593, 0.7591940851704363, -0.6232784248094839, -0.3599550470667362, 0.9590506133650054, -0.026005802462551894, -0.9292902490722584, 0.99877240501213, -0.2328441608368093, -0.5993146274710456, -0.9074554188135073, -0.9889793321032643};
        double result[N];
        for (size_t i=0; i<N; ++i) {
            result[i] = MathLib::Analysis::cos_polyApprox(vec_x[i]);
        }

        EXPECT_EQ(true, MathLib::LinAlg::isEqualVec(N, &result[0], &ans[0], 2.0e-6l));
    }

    TEST_F(AnalysisLibTest, atan_polyApprox) {
        constexpr size_t N = 11;
        const float vec_x[N] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        float result_deg7[N];
        float result_deg9[N];
        float ans[N];
        for (size_t i=0; i<N; ++i) {
            const float x = vec_x[i];
            ans[i] = std::atan(x);
            result_deg7[i] = MathLib::Analysis::atan_polyApprox<float, 7>(x);
            result_deg9[i] = MathLib::Analysis::atan_polyApprox<float, 9>(x);
        }

        const bool isOk_deg7 = MathLib::LinAlg::isEqualVec(N, &result_deg7[0], &ans[0], 1.8e-4f);
        const bool isOk_deg9 = MathLib::LinAlg::isEqualVec(N, &result_deg9[0], &ans[0], 2.5e-5f);
        EXPECT_EQ(true, isOk_deg7 && isOk_deg9);
    }

    TEST_F(AnalysisLibTest, atan2_polyApprox) {
        constexpr size_t N = 10;
        const float vec_x[N] = {0.00124446, 0.531533, -0.953014, 0.563438, 0.0766272, -0.0719563, 0.255811, 0.06043, 0.660943, -0.824122};
        const float vec_y[N] = {-0.298056, 0.727663, 0.523291, -0.880508, 0.657289, -0.0940576, 0.130962, -0.0782994, 0.960371, 0.276083};
        const float ans[N] = {-1.56662, 0.939915, 2.63945, -1.00155, 1.45474, -2.22384, 0.473159, -0.913501, 0.968023, 2.81834};
        float result_deg7[N];
        float result_deg9[N];

        for (size_t i=0; i<N; ++i) {
            result_deg7[i] = MathLib::Analysis::atan2_polyApprox<float, 7>(vec_y[i], vec_x[i]);
            result_deg9[i] = MathLib::Analysis::atan2_polyApprox<float, 9>(vec_y[i], vec_x[i]);
        }

        const bool isOk_deg7 = MathLib::LinAlg::isEqualVec(N, &result_deg7[0], &ans[0], 1.8e-4f);
        const bool isOk_deg9 = MathLib::LinAlg::isEqualVec(N, &result_deg9[0], &ans[0], 2.5e-5f);
        EXPECT_EQ(true, isOk_deg7 && isOk_deg9);
    }
}
