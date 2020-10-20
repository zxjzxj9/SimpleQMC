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
// for lithium atom only (n=1, 2)
template<typename T, int N>
class SlaterWaveFn: public WaveFn<T> {
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
                                (pnu_val[i] - 1)*pow(dist, pnu_val[i] - 2) - zeta_val[i]*pow(dist, pnu_val[i] - 1))*
                                std::exp(-dist*zeta_val[i]);
        }
        return scalar_derv*nvec;
    }

    T laplace(const PCoord<T>& r) {
        auto dist = (r - r0).norm();
        auto scalar_derv = 0;
        auto scalar_dderv = 0;
        for(int i=0; i<nterm; i++) {
            // (ab)'' = a''b + 2a'b'+ ab''
            scalar_derv += norm_const[i]*phi_val[n-1][i]*(
                                (pnu_val[i] - 1)*pow(dist, pnu_val[i] - 2) - zeta_val[i]*pow(dist, pnu_val[i] - 1))*
                                std::exp(-dist*zeta_val[i]);
            scalar_dderv += norm_const[i]*phi_val[n-1][i]*(
                                (pnu_val[i] - 1)*(pnu_val[i] - 2)*pow(dist, pnu_val[i] - 3) -
                                2*zeta_val[i]*(pnu_val[i] - 1)*pow(dist, pnu_val[i] - 2) +
                                zeta_val[i]*zeta_val[i]*pow(dist, pnu_val[i] - 1))* 
                                std::exp(-dist*zeta_val[i]);
        }
        scalar_dderv = scalar_dderv + 2*scalar_derv/dist; // why?
        return  scalar_derv;
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

// For Jastrow term of wave function
template<typename T>
class Jastrow: {

};

// For lithium atom + Jastrow wave function
// Not a child class of wave function
template<typename T>
class SlaterDet: {
public:
    // Radom initialization, r0 is the position of nucleus
    SlaterDet(const PCoord<T> r0): r0(r0) {
        r1 = PCoord<T>::Random();
        r2 = PCoord<T>::Random();
        r3 = PCoord<T>::Random();
    }

    // Given initial coordinate
    SlaterDet(const PCoord<T> r0, const PCoord<T> r1, const PCoord<T> r2, const PCoord<T> r3)
        :r0(r0), r1(r1), r2(r2), r3(r3) {}

    ~SlaterDate() {
        delete s1, s2, s3;
    }

    T eval() {
        // For initialize evaluation and validation
        sdet << s1->value(r1), s2 -> value(r1), 0.0,
                s1->value(r2), s2 -> value(r2), 0.0,
                0.0, 0.0, s3 -> value(r3); 
        inv_sdet = sdet.inverse();
        return sdet.determinant();
    }

    T update(PCoord<T> r, int i) {
        // update the slater det matrix and the inverse matrix
        // i = 1, .., 3
        r1_t = r1_t + PCoord<T>::Random();
        PCoord<T> svec;
        switch (i) {
            case 1:
                svec << s1 -> value(r1_t), s2 -> value(r1_t), 0.0;
                break;
            case 2:
                svec << s1 -> value(r1_t), s2 -> value(r1_t), 0.0;
                break;
            case 3;
                svec << 0.0, 0.0, s3->value(r1_t);
                break;
            default:
                break;
        }

        // apply Sherman–Morrison formula

    }

    T value() {

    }

    PCoord<T> value(int i) {

    }

    T laplace(int i) {

    }

private:
    PCoord<T> r0;
    PCoord<T> r1, r2, r3;
    SlaterWaveFn<T, 1>* s1;
    SlaterWaveFn<T, 1>* s2; // 1s alpha, 1s beta
    SlaterWaveFn<T, 2>* s3; // 2s alpha
    SlaterDet<T, 3, 3> sdet, inv_sdet; // Slater determinant
};