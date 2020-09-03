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
        return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv = gradient(p, func);
    auto derv = wfn->grad(p);
    ASSERT_NEAR((nderv - derv).norm(), 0.0, 1e-6);
    delete wfn;
}

TEST(AtomicWaveFn, Laplacian) {
    auto wfn = new AtomicWaveFn<double>(0.5, 1.0);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
        return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv2 = laplace(p, func);
    auto derv2 = wfn->laplace(p);
    ASSERT_NEAR(nderv2, derv2, 1e-2);
    delete wfn;
}

TEST(VBWaveFn, Gradient) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new VBWaveFn<double>(0.5, 1.0, r1, r2);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv = gradient(p, func);
    auto derv = wfn->grad(p);
    ASSERT_NEAR((nderv - derv).norm(), 0.0, 1e-6);
    delete wfn;
}

TEST(VBWaveFn, Laplacian) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new VBWaveFn<double>(0.5, 1.0, r1, r2);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv2 = laplace(p, func);
    auto derv2 = wfn->laplace(p);
    ASSERT_NEAR(nderv2, derv2, 1e-2);
    delete wfn;
}

TEST(MOWaveFn, Gradient) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new MOWaveFn<double>(0.5, 1.0, r1, r2);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv = gradient(p, func);
    auto derv = wfn->grad(p);
    ASSERT_NEAR((nderv - derv).norm(), 0.0, 1e-6);
    delete wfn;
}

TEST(MOWaveFn, Laplacian) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new MOWaveFn<double>(0.5, 1.0, r1, r2);
    auto func = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p);
    };
    Eigen::Matrix<double, 1, 3> p = Eigen::Matrix<double, 1, 3>::Random(); 
    auto nderv2 = laplace(p, func);
    auto derv2 = wfn->laplace(p);
    ASSERT_NEAR(nderv2, derv2, 1e-2);
    delete wfn;
}

TEST(JastrowWfn, Gradient) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new JastrowWfn<double>(2.0);
    auto func1 = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p, r2);
    };
    auto func2 = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(r1, p);
    };
    Eigen::Matrix<double, 1, 3> derv_f1, derv_f2;
    auto nderv_f1 = gradient(r1, func1);
    auto nderv_f2 = gradient(r2, func2);
    std::tie(derv_f1, derv_f2) = wfn->grad(r1, r2);
    ASSERT_NEAR((nderv_f1 - derv_f1).norm(), 0.0, 1e-6);
    ASSERT_NEAR((nderv_f2 - derv_f2).norm(), 0.0, 1e-6);
    delete wfn;
}

TEST(JastrowWfn, Laplacian) {
    Eigen::Matrix<double, 1, 3> r1 = Eigen::Matrix<double, 1, 3>::Random(); 
    Eigen::Matrix<double, 1, 3> r2 = Eigen::Matrix<double, 1, 3>::Random(); 
    auto wfn = new JastrowWfn<double>(2.0);
    auto func1 = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(p, r2);
    };
    auto func2 = [&](const Eigen::Matrix<double, 1, 3>& p) {
         return wfn->value(r1, p);
    };
    auto nderv2_f1 = laplace(r1, func1);
    auto nderv2_f2 = laplace(r2, func2);
    double derv2_f1, derv2_f2;
    std::tie(derv2_f1, derv2_f2) = wfn->laplace(r1, r2);
    ASSERT_NEAR(nderv2_f1, derv2_f1, 1e-2);
    ASSERT_NEAR(nderv2_f2, derv2_f2, 1e-2);
    delete wfn;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 