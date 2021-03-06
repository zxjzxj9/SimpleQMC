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

    inline T value(const PCoord<T>& r1, const PCoord<T>& r2) {
        return std::exp(-ufunc(r1, r2));
    }

    // The gradient is symmetric w.r.t r1 and r2
    std::pair<PCoord<T>, PCoord<T>> 
    grad(const PCoord<T>& r1, const PCoord<T>& r2) {
        auto val = value(r1, r2);
        auto r12 = r1 - r2;
        auto r = r12.norm();
        auto dv = -r12/(2*r*std::pow(1+r/factor, 2));
        return {-dv*val, dv*val};
    }

    std::pair<T, T> 
    laplace(const PCoord<T>& r1, const PCoord<T>& r2) {
        auto val = value(r1, r2);
        auto r = (r1 - r2).norm();
        auto dv = -1/(2*std::pow(1+r/factor, 2));
        auto d2v = -1/(r*std::pow(1+r/factor, 3));
        auto ret = (dv*dv - d2v)*val;
        return {ret, ret};
    }
private:
    T factor;
    inline T ufunc(const PCoord<T>& r1, const PCoord<T> r2) {
        auto r = (r1 - r2).norm();
        return factor/(2 + 2*r/factor);
    }
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

        ret = -0.5*ret; // kintetic energy

        ret += -1.0/(r1-R1).norm();
        ret += -1.0/(r1-R2).norm();
        ret += -1.0/(r2-R1).norm();
        ret += -1.0/(r2-R2).norm();
        ret += 1.0/(r1-r2).norm();
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
        mol = new H2Mol<T, JastrowType::SIMPLE_JASTROW, AtomicWfnType::MO>(factor, c, alpha, R1, R2);
    }

    ~H2MolQMC() {
        delete mol;
    }

    std::pair<T, T> sample(int maxstep=10000) {
        
        // Thsi method return the number between [-1, 1]
        PCoord<T> r1 = PCoord<T>::Random(); 
        PCoord<T> r2 = PCoord<T>::Random();

        PCoord<T> r1_new, r2_new;
        T energy = mol->energy(r1, r2);
        T density = mol->density(r1, r2);
        std::uniform_real_distribution<T> rnum(0, 1);

        //energy_old = mol->energy();
        T energy_tot = 0.0;
        T energy_sq_tot = 0.0;
        int accept = 0;
        for(int i=0; i<maxstep; i++) {

            dr = std::max(dr, 0.1);
            dr = std::min(dr, 10.0);

            r1_new = r1 + dr*PCoord<T>::Random();
            r2_new = r2 + dr*PCoord<T>::Random();

            auto dtmp = mol->density(r1_new, r2_new);
            auto ratio = dtmp/density;
            
            if(ratio > 1 || ratio > rnum(rgen)) {
                r1 = r1_new;
                r2 = r2_new;
                density = dtmp;
                energy = mol->energy(r1, r2);
                accept += 1;
            }

            if(static_cast<T>(accept)/(i+1) > 0.5) dr*=scale;
            else dr/=scale;

            energy_tot += energy;
            energy_sq_tot += energy*energy;
        }
        auto energy_avg = energy_tot/maxstep;
        auto energy_std = std::sqrt(energy_sq_tot/maxstep - energy_avg*energy_avg);
        fmt::print("Accept ratio: {}\n", (double) accept/maxstep);
        return {energy_avg, energy_std};
    }

private:
    H2Mol<T, JastrowType::SIMPLE_JASTROW, AtomicWfnType::MO>* mol;
    PCoord<T> R1, R2;
    T dr;
    std::random_device rd;
    std::mt19937 rgen{rd()};
    const T scale = 1.01;
    // T c, alpha;
};