#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>

template<typename T>
class WaveFn {
    // return the value of wave function
    virtual T operator()(Eigen::Matrix<T, 1, 3>) = 0;

    // return three components of grad wave function
    virtual Eigen::Matrix<T, 1, 3> grad(Eigen::Matrix<T, 1, 3>) = 0;
    
    // return laplacian of the wave function
    virtual T laplace(Eigen::Matrix<T, 1, 3>) = 0;

};

template<typename T, typename wfn>
class H2MolQMC {
public:
    H2MolQMC() {
    };

    std::pair<T, T> sample(int maxstep=10000) {
        return {0, 0};
    }

private:
    Eigen::Matrix<T, 2, 3> coord;

    // T c, alpha;
};