#pragma once

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
    const int n = n;
    SlaterWaveFn(){

    }
    ~SlaterWaveFn(){};
private:
    T phi_val[2][6] = {
        {-0.12220686, 1.11273225, 0.04125378, 0.09306499, -0.10260021, -0.00034191, 0.00021963},
        {0.47750469, 0.11140449,-1.25954273, -0.18475003, -0.02736293, -0.00025064, 0.00057962}
    };
    T zeta_val[6] = {
        0.72089388, 2.61691643, 0.69257443, 1.37137558, 3.97864549, 13.52900016, 19.30801440
    };

    T slap[2] = {1.00, 0.95};
};