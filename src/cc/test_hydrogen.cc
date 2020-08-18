#include <gtest/gtest.h>
#include "hydrogen.hpp"

// Given the function, return the graidents 
// w.r.t. the coordinates
template<typename T, typename Fn> 
Eigen::Matrix<T, 1, 3> gradient(Eigen::Matrix<T, 1, 3> p, Fn func) {
    const T eps = 1e-6;
    // +2h +h -h -2h, and divide by 12
    T coeff[] = {-1, 8, -8, 1};
    T ds[] = {+2, +1, -1, -2};

    auto dfunc = [&](Eigen::Matrix<T, 1, 3> dp){
        T ret = 0.0;
        for(int i=0; i<4; i++) {
            ret += coeff[i]*func(p+ds[i]*dp);
        }
        return ret/12.0;
    };

    T delta = 1e-6
    Eigen::Matrix<T, 1, 3> dx;
    dx << delta << 0.0 << 0.0;
    Eigen::Matrix<T, 1, 3> dy;
    dy << 0.0 << delta << 0.0;
    Eigen::Matrix<T, 1, 3> dz;
    dz << 0.0 << 0.0 << delta;

    Matrix<T, 1, 3> ret;
    ret << dfunc(dx) << dfunc(dy) << dfunc(dz);
    return ret;
}

TEST(AtomicWaveFnTEST, Gradient) {
    // ASSERT_DOUBLE_EQ()
}

TEST(AtomicWaveFnTEST, Laplace) {

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 