#include <gtest/gtest.h>
#include "hydrogen.hpp"
//#include <unsupported/Eigen/MatrixFunctions>

// Given the function, return the graidents 
// w.r.t. the coordinates
template<typename T, typename Fn> 
Eigen::Matrix<T, 1, 3> gradient(const Eigen::Matrix<T, 1, 3>& p, Fn func) {
    const T eps = 1e-6;
    // +2h +h -h -2h, and divide by 12
    T coeff[] = {-1, 8, -8, 1};
    T ds[] = {+2, +1, -1, -2};

    auto dfunc = [&](const Eigen::Matrix<T, 1, 3>& dp){
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

// Given the function, return the laplacian
// w.r.t. the coordinates
template<typename T, typename Fn>
T laplace(const Eigen::Matrix<T, 1, 3>& p, Fn func) {
    const T eps = 1e-6;

    Eigen::Matrix<T, 1, 3> dx;
    dx << eps, 0.0, 0.0;
    Eigen::Matrix<T, 1, 3> dy;
    dy << 0.0, eps, 0.0;
    Eigen::Matrix<T, 1, 3> dz;
    dz << 0.0, 0.0, eps;

    return (func(p+dx) + func(p-dx) + 
            func(p+dy) + func(p-dy) + 
            func(p+dz) + func(p-dz) - 6*func(p))/(eps*eps);
}

TEST(SimpleGradientTEST, Gradient) {
    auto func = [](const Eigen::Matrix<double, 1, 3>& p) {
        return p.array().pow(2).sum();
    };

    auto dfunc = [](const Eigen::Matrix<double, 1, 3>& p) {
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

TEST(SimpleLaplacianTEST, Laplacian) {
    auto func = [](const Eigen::Matrix<double, 1, 3>& p) {
        return p.array().pow(2).sum();
    };

    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random();
    double nderv2 = laplace(p, func);
    double derv2 = 6.0;

    ASSERT_NEAR(nderv2, derv2, 1e-2);
}

TEST(AtomicWaveFn, Gradient) {
    auto wfn = new AtomicWaveFn<double>(0.5, 1.0);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
        return wfn->density(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv = gradient(p, func);
    auto derv = wfn->grad(p);
    delete wfn;
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