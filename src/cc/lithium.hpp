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

template<typename T>
class SlaterWaveFn {
public:
    SlaterWaveFn(){

    }
    ~SlaterWaveFn(){};
private:

};