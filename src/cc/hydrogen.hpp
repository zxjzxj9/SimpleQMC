#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>

template<typename T>
class WaveFn {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    // return the value of wave function
    virtual T operator()(const PCoord&) = 0;

    // return three components of grad wave function
    virtual PCoord grad(const PCoord&) = 0;
    
    // return laplacian of the wave function
    virtual T laplace(const PCoord&) = 0;

};

// Simple Wave function $$\phi(r) = (1+cr)e^{-\alpha r}$$
template <typename T>
class AtomicWaveFn: public WaveFn<T> {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    AtomicWaveFn(T c, T alpha): c(c), alpha(alpha) {};

    T operator()(const PCoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        return (1+c*r)*std::exp(-ar);
    }

    PCoord grad(const PCoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        auto coeff = (1-alpha*(c*r+1))*std::exp(-ar);
        return coeff*coord.normalized();
    }

    T laplace(const Pcoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        return (c*(ar*(ar-4)+2) + alpha*(ar-2))*std::exp(-ar)/r;
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