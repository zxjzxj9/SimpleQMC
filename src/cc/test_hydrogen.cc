#include <gtest/gtest.h>
#include "hydrogen.hpp"

// Given the function, return the graidents 
// w.r.t. the coordinates
template<typename T, typename Fn> 
Eigen::Matrix<T, 1, 3> gradient(Eigen::Matrix<T, 1, 3> p, Fn func) {
    const T eps = 1e-6;
    // +2h +h -h -2h
    T coeff[] = {-1, 8, -8, 1};
}

TEST(AtomicWaveFnTEST, Gradient) {
    ASSERT_DOUBLE_EQ()
}

TEST(AtomicWaveFnTEST, Laplace) {

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 