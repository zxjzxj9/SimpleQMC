#pragma once

#include <Eigen/Dense>
#include <random>

template<typename T>
using PCoord = Eigen::Matrix<T, 1, 3>;

int factorial(int n) {
    if(n==0) return 1;
    int ret = 1;
    for(int i=1; i<=n; i++) ret *= i;
    return ret;
}

template<typename T>
class WaveFn {
public:
    virtual ~WaveFn(){};
    // return the value of wave function
    virtual T value(const PCoord<T>&)=0;

    // return three components of grad wave function
    virtual PCoord<T> grad(const PCoord<T>&)=0;
    
    // return laplacian of the wave function
    virtual T laplace(const PCoord<T>&)=0;
};

// n define the quantum number n of wave function
template<typename T, int N>
class SlaterWaveFn {
public:
    const int n = N;
    const int nterm = 7;
    SlaterWaveFn(const PCoord<T>& r0): r0(r0){
        for(int i=0; i<nterm; i++) {
            // recaluate zeta_val accorindg to orbital
            zeta_val[i] = zeta_val[i]/slap[n-1];
            norm_const[i] = std::pow(2*zeta_val[i], pnu_val[i]-0.5)/
                            std::sqrt(factorial(2*pnu_val[i]));
        }
    }
    ~SlaterWaveFn(){};

    T value(const PCoord<T>& r) {
        auto dist = (r - r0).norm();
        T ret = 0.0;
        for(int i=0; i<nterm; i++) {
            ret += norm_const[i]*phi_val[n-1][i]*
                   std::pow(dist, pnu_val[i]-1)*std::exp(-dist*zeta_val[i]);
        }
        return ret;
    }

    PCoord<T> grad(const PCoord<T>& r) {
        auto dist = (r - r0).norm();
        auto nvec = (r - r0).normalized();
        auto scalar_derv = 0;
        for(int i=0; i<nterm; i++) {
            scalar_derv += norm_const[i]*phi_val[n-1][i]*(
                                (pnu_val[i] - 1)*pow(dist, pnu_val[i] - 2) -zeta_val[i]*pow(dist, pnu_val[i] - 1))*
                                std::exp(-dist*zeta_val[i]);
        }
        return scalar_derv*nvec;
    }

    T laplace(const PCoord<T>& r) {
        auto dist = (r - r0).norm();
        auto scalar_derv = 0;
        auto scalar_dderv = 0;
        for(int i=0; i<nterm; i++) {
            scalar_derv += norm_const[i]*phi_val[n-1][i]*(
                                (pnu_val[i] - 1)*pow(dist, pnu_val[i] - 2) -zeta_val[i]*pow(dist, pnu_val[i] - 1))*
                                std::exp(-dist*zeta_val[i]);
            // TODOs: formula for the second derivative
            scalar_dderv = 0.0;
        }

    }
    
private:
    PCoord<T> r0;

    T phi_val[2][7] = {
        {-0.12220686, 1.11273225, 0.04125378, 0.09306499, -0.10260021, -0.00034191, 0.00021963},
        {0.47750469, 0.11140449,-1.25954273, -0.18475003, -0.02736293, -0.00025064, 0.00057962}
    };
    T zeta_val[7] = {
        0.72089388, 2.61691643, 0.69257443, 1.37137558, 3.97864549, 13.52900016, 19.30801440
    };
    T norm_const[7];
    int pnu_val[7] = {
        1, 1, 2, 2, 2, 2, 3 // 1s, 1s, 2s, 2s, 2s, 2s, 3s
    };
    T slap[2] = {1.00, 0.95};
};