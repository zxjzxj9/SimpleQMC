#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>

template<typename T>
class WaveFn {
public:
    // return the value of wave function
    virtual T operator()(Eigen::Matrix<T, 1, 3>) = 0;

    // return three components of grad wave function
    virtual Eigen::Matrix<T, 1, 3> grad(Eigen::Matrix<T, 1, 3>) = 0;
    
    // return laplacian of the wave function
    virtual T laplace(Eigen::Matrix<T, 1, 3>) = 0;

};

// Simple Wave function $$\phi(r) = (1+cr)e^{-\alpha r}$$
template <typename T>
class AtomicWaveFn: public WaveFn<T> {
public:
    AtomicWaveFn(T c, T alpha): c(c), alpha(alpha) {};

    virtual T operator()(Eigen::Matrix<T, 3, 1> coord) {
        auto r = coord.norm();
        return (1+c*r)*std::exp(-alpha*r);
    }

    virtual Eigen::Matrix<T, 1, 3> grad(Eigen::Matrix<T, 1, 3> coord)  {
        auto r = coord.norm();
        auto coeff = (1-alpha*(c*r+1))*std::exp(-alpha*r);
        return coeff*coord.normalized();
    }

private:
    T c, alpha;
};

template<typename T> //, typename wfn>
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