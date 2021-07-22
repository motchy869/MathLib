#include <gtest/gtest.h>
#include "../include/analysis.hpp"
#include "../include/linAlg.hpp"

namespace {
    class AnalysisLibTest : public ::testing::Test{};

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
