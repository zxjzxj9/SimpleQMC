#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>

template<typename T>
class H2MolQMC {
public:
    H2MolQMC() {
    };

    std::pair<T, T> sample(int maxstep=10000) {
        return {0, 0};
    }

private:
    Eigen::Matrix<T, 2, 3> coord;    
};