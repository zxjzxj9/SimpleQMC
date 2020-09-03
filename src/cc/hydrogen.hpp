#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>

enum JastrowType {
    SIMPLE_JASTROW,
};

enum AtomicWfnType {
    VB,
    MO,
};

template<typename T>
class WaveFn {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    virtual ~WaveFn(){};
    // return the value of wave function
    virtual T value(const PCoord&)=0;

    // return three components of grad wave function
    virtual PCoord grad(const PCoord&)=0;
    
    // return laplacian of the wave function
    virtual T laplace(const PCoord&)=0;
};

// Simple Wave function $$\phi(r) = (1+cr)e^{-\alpha r}$$
template <typename T>
class AtomicWaveFn: public WaveFn<T> {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    AtomicWaveFn(T c, T alpha): c(c), alpha(alpha) {};
    ~AtomicWaveFn(){};

    T value(const PCoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        return (1+c*r)*std::exp(-ar);
    }

    PCoord grad(const PCoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        auto coeff = (c-alpha*(c*r+1))*std::exp(-ar);
        return coeff*coord.normalized();
    }

    T laplace(const PCoord& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        return (c*(ar*(ar-4)+2) + alpha*(ar-2))*std::exp(-ar)/r;
    }

private:
    T c, alpha;
};

// Composed Wave functions: VB
template<typename T>
class VBWaveFn: public WaveFn<T> {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    VBWaveFn(T c, T alpha, const PCoord& R1, const PCoord& R2): 
        c(c), alpha(alpha), R1(R1), R2(R2) {
        phi1 = new AtomicWaveFn<T>(c, alpha);
        phi2 = new AtomicWaveFn<T>(c, alpha);
    }

    ~VBWaveFn() {
        delete phi1, phi2;
    }

    T value(const PCoord& r) {
        return phi1->value(r-R1)*phi2->value(r-R2);
    }
    
    PCoord grad(const PCoord& r) {
        return phi2->value(r-R2)*phi1->grad(r-R1) + 
               phi1->value(r-R1)*phi2->grad(r-R2);
    }

    T laplace(const PCoord& r) {
        return phi2->value(r-R2)*phi1->laplace(r-R1) + 
               phi1->value(r-R1)*phi2->laplace(r-R2) +
               2*(phi1->grad(r-R1)).dot(phi2->grad(r-R2));
    }

private:
    T c, alpha;
    PCoord R1, R2;
    AtomicWaveFn<T>* phi1;
    AtomicWaveFn<T>* phi2;
};

// Composed Wave functions: MO
template <typename T>
class MOWaveFn: public WaveFn<T> {
public:
    using PCoord = Eigen::Matrix<T, 1, 3>;
    MOWaveFn(T c, T alpha, const PCoord& R1, const PCoord& R2): 
        c(c), alpha(alpha), R1(R1), R2(R2) {
        phi1 = new AtomicWaveFn<T>(c, alpha);
        phi2 = new AtomicWaveFn<T>(c, alpha);
    }

    T value(const PCoord& r) {
        return phi1->value(r-R1)+phi2->value(r-R2);
    }
    
    PCoord grad(const PCoord& r) {
        return phi1->grad(r-R1) + phi2->grad(r-R2);
    }

    T laplace(const PCoord& r) {
        return phi1->laplace(r-R1) + phi2->laplace(r-R2);
    }

private:
    T c, alpha;
    PCoord R1, R2;
    AtomicWaveFn<T>* phi1;
    AtomicWaveFn<T>* phi2;
};

// Define Jastrow wavefunction
template<typename T>
// class JastrowWfn: public WaveFn<T> {
class JastrowWfn {
    using PCoord = Eigen::Matrix<T, 1, 3>;
public:
    JastrowWfn(T factor): factor(factor) {}

    ~JastrowWfn() {}

    T value(const PCoord& r1, const PCoord& r2) {
        auto r12 = (r1 - r2).norm();
        return factor/(2 + 2*r12/factor);
    }

    PCoord grad(const PCoord& r1, const PCoord& r2) {
        auto r12 = r1 - r2;
        auto r = r12.norm();
        return -r12/(2*r*std::pow(1+r/factor, 2));
    }

    T laplace(const PCoord& r1, const PCoord& r2) {
        auto r12 = (r1 - r2).norm();
        return -1/(2*r12*std::pow(1+r12/factor, 3)); 
    }
private:
    T factor;
};

template<typename T, JastrowType Jastrow, AtomicWfnType AtomicWfn>
class H2Mol {

public:
    const JastrowType jastrow_type = Jastrow;
    const AtomicWfnType atomicwfn_type = AtomicWfn;
    using PCoord = Eigen::Matrix<T, 1, 3>;

    H2Mol(T factor, T c, T alpha, PCoord& R1, PCoord& R2) {
        switch (Jastrow) {
            case JastrowType::SIMPLE_JASTROW:
                jastrow = new JastrowWfn<T>(factor);
                break;
            default:
                throw std::runtime_error("Invalid jastrow function.");
                break;
        }

        switch(AtomicWfn) {
            case AtomicWfnType::MO:
                atomicwfn = new MOWaveFn<T>(c, alpha, R1, R2);
                break;
            case AtomicWfnType::VB:
                atomicwfn = new VBWaveFn<T>(c, alpha, R1, R2);
                break;
            default:
                throw std::runtime_error("Invalid atomic wave function");
                break;
        }
    }

    ~H2Mol() {
        delete jastrow;
        delete atomicwfn;
    }

    // Notice this calculate \psi^* \psi
    T density(PCoord& r1, PCoord& r2) {
        return std::pow(jastrow->value(r1, r2)* \
               atomicwfn->value(r1)*atomicwfn->value(r2), 2);
    }

    // Calculate the energy with current density
    T energy(PCoord& r1, PCoord& r2) {
        T ret = 0.0;
        return ret;
    }

private:
    WaveFn<T>* jastrow = nullptr;
    WaveFn<T>* atomicwfn = nullptr;
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