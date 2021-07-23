#include <gtest/gtest.h>
#include "../include/analysis.hpp"
#include "../include/linAlg.hpp"

namespace {
    class AnalysisLibTest : public ::testing::Test{};

    TEST_F(AnalysisLibTest, sin_polyApprox) {
        constexpr size_t N = 20;
        const double vec_x[N] = {-5.422284533254419, -8.368253325506895, 7.230370357282717, 8.540266743128377, 8.796032571059138, -9.0058423693436, -4.1458323100147165, 6.6625543111462315, 4.714898106104982, -6.991907637267804, -4.039460886872965, -8.222201344331197, -6.5703506066193365, 1.5968050614449092, -2.7632828820450843, 6.233630363227799, -1.805797535337036, 4.069744312298307, -8.991169132560007, -9.276178027269587};
        const double ans[N] = {0.7584299540206195, -0.8706512087514109, 0.8117748754722436, 0.7736053539586267, 0.5881305349111817, -0.4067883225567331, 0.843754111510663, 0.3703344119084007, -0.9999968521457114, -0.6508643030941429, 0.7820000033037139, -0.9329696480010395, -0.2832347454042438, 0.999661791926789, -0.3693502849318244, -0.04953466449150098, -0.972514059931373, -0.800513570965066, -0.4201483819628514, -0.14805364120001024};
        double result[N];
        for (size_t i=0; i<N; ++i) {
            result[i] = MotchyMathLib::Analysis::sin_polyApprox(vec_x[i]);
        }
        //MotchyMathLib::LinAlg::Debug::printRealVec(N, result, "%g");

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(N, &result[0], &ans[0], 2.0e-6l));
    }

    TEST_F(AnalysisLibTest, atan2_polyApprox) {
        constexpr size_t N = 10;
        const float vec_x[N] = {0.00124446, 0.531533, -0.953014, 0.563438, 0.0766272, -0.0719563, 0.255811, 0.06043, 0.660943, -0.824122};
        const float vec_y[N] = {-0.298056, 0.727663, 0.523291, -0.880508, 0.657289, -0.0940576, 0.130962, -0.0782994, 0.960371, 0.276083};
        const float ans[N] = {-1.56662, 0.939915, 2.63945, -1.00155, 1.45474, -2.22384, 0.473159, -0.913501, 0.968023, 2.81834};

        float result[N];
        for (size_t i=0; i<N; ++i) {
            result[i] = MotchyMathLib::Analysis::atan2_polyApprox(vec_y[i], vec_x[i]);
        }

        EXPECT_EQ(true, MotchyMathLib::LinAlg::isEqualVec(N, &result[0], &ans[0], 2.5e-5f));
    }
}
