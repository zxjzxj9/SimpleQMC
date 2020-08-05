#pragma once
#include <cmath>
#include <random>
#include <fmt/core.h>

template <typename T>
class NaiveQMC {

public:
    NaiveQMC(T c, T alpha, T dr): c(c), alpha(alpha), dr(dr) {
    }

    // return rho1 / rho2
    T inline rho_ratio(T r1, T r2) {
        // r1 = std::abs(r1);
        // r2 = std::abs(r2);
        // if(r1 < 1e-6) return 0;
        return std::pow(c*r1+1, 2)/std::pow(c*r2+1, 2)*
            std::exp(2*alpha*(r2-r1));
    }

    T inline energy_func(T r) {
        // r = std::abs(r);
        return (-c*(r*(alpha*(alpha*r-4)+2)+2)+
            alpha*(2-alpha*r)-2)/(2*r*(c*r+1));
    }

    std::pair<T, T> sample(int maxstep=10000) {
        std::uniform_real_distribution<T> dist(-1.0, 1.0);
        std::uniform_real_distribution<T> rnum(0, 1);
        T rnew, rold = 1.0;
        T energy = energy_func(rold);
        // fmt::print("{:f}\n",rold);
        // exit(0);
        T etot = 0;
        T etot_sq = 0;
        int accept = 0;
        for(int i=0; i<maxstep; i++) {
            dr = std::min(dr, rold);
            rnew = rold + dr*dist(rgen);
            fmt::print("{}\n", rold);
            auto ratio = rho_ratio(rnew, rold);
            if(ratio > 1 || ratio > rnum(rgen)) {
                rold = rnew;
                energy = energy_func(rold);
                accept++;
            }

            // if(static_cast<T>(accept)/(i+1) > 0.5) dr*=scale;
            // else dr/=scale;

            etot += energy;
            etot_sq += std::pow(energy, 2);
        }
        auto mean = etot/static_cast<T>(maxstep);
        // auto std = etot_sq/static_cast<T>(maxstep);
        fmt::print("Accept ration: {:.2f}, step size: {:6.4f}, last place: {:6.4f}\n", static_cast<T>(accept)/maxstep, dr, rold);
        auto std = std::sqrt(etot_sq/static_cast<T>(maxstep) - std::pow(mean, 2));
        return std::make_pair(mean, std);
    }

private:
    T c, alpha;
    T dr;
    std::random_device rd;
    std::mt19937 rgen{rd()};
    const T scale = 1.01;
};