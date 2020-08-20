#include <gtest/gtest.h>
#include "hydrogen.hpp"
//#include <unsupported/Eigen/MatrixFunctions>

// Given the function, return the graidents 
// w.r.t. the coordinates
template<typename T, typename Fn> 
Eigen::Matrix<T, 1, 3> gradient(Eigen::Matrix<T, 1, 3> p, Fn func) {
    const T eps = 1e-6;
    // +2h +h -h -2h, and divide by 12
    T coeff[] = {-1, 8, -8, 1};
    T ds[] = {+2, +1, -1, -2};

    auto dfunc = [&](Eigen::Matrix<T, 1, 3>& dp){
        T ret = 0.0;
        for(int i=0; i<4; i++) {
            ret += coeff[i]*func(p+ds[i]*dp);
        }
        return ret/(12.0*eps);
    };

    Eigen::Matrix<T, 1, 3> dx;
    dx << eps, 0.0, 0.0;
    Eigen::Matrix<T, 1, 3> dy;
    dy << 0.0, eps, 0.0;
    Eigen::Matrix<T, 1, 3> dz;
    dz << 0.0, 0.0, eps;

    Eigen::Matrix<T, 1, 3> ret;
    ret << dfunc(dx), dfunc(dy), dfunc(dz);
    return ret;
}

TEST(SimpleGradientTEST, Gradient) {
    auto func = [](Eigen::Matrix<double, 1, 3>&& p) {
        return p.array().pow(2).sum();
    };

    auto dfunc = [](Eigen::Matrix<double, 1, 3>& p) {
        auto ret = 2.0*p;
        return ret;
    };

    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random();
    // std::cout<<p<<std::endl;
    auto nderv = gradient(p, func);
    // std::cout<<nderv<<std::endl;
    auto derv = dfunc(p);
    // std::cout<<derv<<std::endl;
    ASSERT_NEAR((nderv - derv).norm(), 0.0, 1e-6);
}

//TEST(AtomicWaveFnTEST, Gradient) {
    // ASSERT_DOUBLE_EQ()
//}

//TEST(AtomicWaveFnTEST, Laplace) {
//
//}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 