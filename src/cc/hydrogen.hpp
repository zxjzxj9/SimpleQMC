#pragma once
#include <fmt/core.h>
#include <Eigen/Dense>
#include <random>

template<typename T>
using PCoord = Eigen::Matrix<T, 1, 3>;

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
    virtual ~WaveFn(){};
    // return the value of wave function
    virtual T value(const PCoord<T>&)=0;

    // return three components of grad wave function
    virtual PCoord<T> grad(const PCoord<T>&)=0;
    
    // return laplacian of the wave function
    virtual T laplace(const PCoord<T>&)=0;
};

// Simple Wave function $$\phi(r) = (1+cr)e^{-\alpha r}$$
template <typename T>
class AtomicWaveFn: public WaveFn<T> {
public:
    AtomicWaveFn(T c, T alpha): c(c), alpha(alpha) {};
    ~AtomicWaveFn(){};

    T value(const PCoord<T>& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        return (1+c*r)*std::exp(-ar);
    }

    PCoord<T> grad(const PCoord<T>& coord) {
        auto r = coord.norm();
        auto ar = alpha*r;
        auto coeff = (c-alpha*(c*r+1))*std::exp(-ar);
        return coeff*coord.normalized();
    }

    T laplace(const PCoord<T>& coord) {
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
    VBWaveFn(T c, T alpha, const PCoord<T>& R1, const PCoord<T>& R2): 
        c(c), alpha(alpha), R1(R1), R2(R2) {
        phi1 = new AtomicWaveFn<T>(c, alpha);
        phi2 = new AtomicWaveFn<T>(c, alpha);
    }

    ~VBWaveFn() {
        delete phi1, phi2;
    }

    T value(const PCoord<T>& r) {
        return phi1->value(r-R1)*phi2->value(r-R2);
    }
    
    PCoord<T> grad(const PCoord<T>& r) {
        return phi2->value(r-R2)*phi1->grad(r-R1) + 
               phi1->value(r-R1)*phi2->grad(r-R2);
    }

    T laplace(const PCoord<T>& r) {
        return phi2->value(r-R2)*phi1->laplace(r-R1) + 
               phi1->value(r-R1)*phi2->laplace(r-R2) +
               2*(phi1->grad(r-R1)).dot(phi2->grad(r-R2));
    }

private:
    T c, alpha;
    PCoord<T> R1, R2;
    AtomicWaveFn<T>* phi1;
    AtomicWaveFn<T>* phi2;
};

// Composed Wave functions: MO
template <typename T>
class MOWaveFn: public WaveFn<T> {
public:
    MOWaveFn(T c, T alpha, const PCoord<T>& R1, const PCoord<T>& R2): 
        c(c), alpha(alpha), R1(R1), R2(R2) {
        phi1 = new AtomicWaveFn<T>(c, alpha);
        phi2 = new AtomicWaveFn<T>(c, alpha);
    }

    T value(const PCoord<T>& r) {
        return phi1->value(r-R1)+phi2->value(r-R2);
    }
    
    PCoord<T> grad(const PCoord<T>& r) {
        return phi1->grad(r-R1) + phi2->grad(r-R2);
    }

    T laplace(const PCoord<T>& r) {
        return phi1->laplace(r-R1) + phi2->laplace(r-R2);
    }

private:
    T c, alpha;
    PCoord<T> R1, R2;
    AtomicWaveFn<T>* phi1;
    AtomicWaveFn<T>* phi2;
};

// Define Jastrow wavefunction
template<typename T>
// class JastrowWfn: public WaveFn<T> {
class JastrowWfn {
public:
    JastrowWfn(T factor): factor(factor) {}

    ~JastrowWfn() {}

    T value(const PCoord<T>& r1, const PCoord<T>& r2) {
        auto r12 = (r1 - r2).norm();
        return factor/(2 + 2*r12/factor);
    }

    // The gradient is symmetric w.r.t r1 and r2
    std::pair<PCoord<T>, PCoord<T>> 
    grad(const PCoord<T>& r1, const PCoord<T>& r2) {
        auto r12 = r1 - r2;
        auto r = r12.norm();
        auto v1 = -r12/(2*r*std::pow(1+r/factor, 2));
        return {v1, -v1};
    }

    std::pair<T, T> 
    laplace(const PCoord<T>& r1, const PCoord<T>& r2) {
        auto r = (r1 - r2).norm();
        auto val = -1/(r*std::pow(1+r/factor, 3)); 
        return {val, val};
    }
private:
    T factor;
};

template<typename T, JastrowType Jastrow, AtomicWfnType AtomicWfn>
class H2Mol {

public:
    const JastrowType jastrow_type = Jastrow;
    const AtomicWfnType atomicwfn_type = AtomicWfn;

    H2Mol(T factor, T c, T alpha, const PCoord<T>& R1, const PCoord<T>& R2):
        R1(R1), R2(R2) {
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
    T density(const PCoord<T>& r1, const PCoord<T>& r2) {
        return std::pow(jastrow->value(r1, r2)* \
               atomicwfn->value(r1)*atomicwfn->value(r2), 2);
    }

    // Calculate the energy with current density
    T energy(const PCoord<T>& r1, const PCoord<T>& r2) {
        T ret = 0.0;
        PCoord<T> dJdr1, dJdr2;
        T d2Jdr1, d2Jdr2;
        auto jval = jastrow->value(r1, r2);
        std::tie(d2Jdr1, d2Jdr2) = jastrow->laplace(r1, r2);
        std::tie(dJdr1, dJdr2) = jastrow->grad(r1, r2);
        ret += (d2Jdr1+d2Jdr2)/jval;
        auto val1 = atomicwfn->value(r1);
        auto val2 = atomicwfn->value(r2);
        ret += atomicwfn->laplace(r1)/val1;
        ret += atomicwfn->laplace(r2)/val2;
        
        ret += 2*dJdr1.dot(atomicwfn->grad(r1))/(jval*val1);
        ret += 2*dJdr2.dot(atomicwfn->grad(r2))/(jval*val2);

        ret += -1.0/(r1-R1).norm();
        ret += -1.0/(r1-R2).norm();
        ret += -1.0/(r2-R1).norm();
        ret += -1.0/(r2-R2).norm();
        ret += 1.0/(R1-R2).norm();

        return ret;
    }

private:
    JastrowWfn<T>* jastrow = nullptr;
    WaveFn<T>* atomicwfn = nullptr;
    PCoord<T> R1, R2;
};


template<typename T> //, typename wfn>
class H2MolQMC {
public:
    H2MolQMC(T factor, T c, T alpha, const PCoord<T>& R1, const PCoord<T>& R2, T dr): dr(dr){
        mol = new H2Mol<T, JastrowType::SIMPLE_JASTROW, AtomicWfnType::VB>(factor, c, alpha, R1, R2);
        
    }

    ~H2MolQMC() {
        delete mol;
    }

    std::pair<T, T> sample(int maxstep=10000) {
        
        PCoord<T> r1 = PCoord<T>::Random(); 
        PCoord<T> r2 = PCoord<T>::Random();
        T energy_old = mol->energy(r1, r2);
        T density_old = mol->density(r1, r2);
        T energy_new = 0.0;
        T density_new = 0.0;

        //energy_old = mol->energy();
        T energy_tot = 0.0;
        T energy_sq_tot = 0.0;
        for(int i=0; i<maxstep; i++) {
            r1 += 2*dr*(PCoord<T>::Random()-1.0);
            r2 += 2*dr*(PCoord<T>::Random()-1.0);
            energy_new = mol->energy(r1, r2);
            density_new = mol->density(r1, r2);

            if(density_new/density_old > 
                std::uniform_real_distribution<T>(rgen)) {
                energy_old = energy_new;
                density_old = density_new;
            }
            energy_tot += energy_old;
            energy_sq_tot += energy_old*energy_old;
        }
        auto energy_avg = energy_tot/maxstep;
        auto energy_std = std::sqrt(energy_sq_tot/maxstep - energy_avg*energy_avg);
        return {energy_avg, energy_std};
    }

private:
    H2Mol<T, JastrowType::SIMPLE_JASTROW, AtomicWfnType::VB>* mol;
    PCoord<T> R1, R2;
    T dr;
    std::mt19937 rgen(std::random_device());
    // T c, alpha;
};