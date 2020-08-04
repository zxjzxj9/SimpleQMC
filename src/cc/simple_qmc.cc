#include <functional>
#include <fmt/core.h>
#include <cxxopts.hpp>
#include "simple_qmc.hpp"

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

int main(int argc, char** argv) {
    cxxopts::Options options("SimpleQMC", "Simple Quantum Monte Carlo Program");
    options.add_options()
        ("r,range", "Use range exploration", cxxopts::value<bool>()->default_value("false"))
        ("c", "Trial wave function parameter c", cxxopts::value<double>()->default_value("0.0"))
        ("alpha", "Trial wave function parameter alpha", cxxopts::value<double>()->default_value("1.0"))
        ("s,step", "Monte Carlo step size", cxxopts::value<double>()->default_value("1.0"))
        ("n,nstep", "Monte Carlo step size", cxxopts::value<int>()->default_value("1000000"))
        ("h,help", "Print usage")
        ("cmin", "Minimum parameter c", cxxopts::value<double>())
        ("cmax", "Maximum parameter c", cxxopts::value<double>())
        ("alphamin", "Minimum parameter alpha", cxxopts::value<double>())
        ("alphamax", "Maximum parameter alpha", cxxopts::value<double>())
        ("ngridc", "Number of grid in c", cxxopts::value<int>())
        ("ngridalpha", "Number of grid alpha", cxxopts::value<int>())
    ;
    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        fmt::print("%s\n", options.help());
        exit(0);
    }
    std::string banner =                                                
        R"(                               ____                  )" "\n" 
        R"(        ,----..              ,'  , `.   ,----..      )" "\n" 
        R"(       /   /   \          ,-+-,.' _ |  /   /   \     )" "\n" 
        R"(      /   .     :      ,-+-. ;   , || |   :     :    )" "\n" 
        R"(     .   /   ;.  \    ,--.'|'   |  ;| .   |  ;. /    )" "\n" 
        R"(    .   ;   /  ` ;   |   |  ,', |  ': .   ; /--`     )" "\n" 
        R"(    ;   |  ; \ ; |   |   | /  | |  || ;   | ;        )" "\n" 
        R"(    |   :  | ; | '   '   | :  | :  |, |   : |        )" "\n" 
        R"(    .   |  ' ' ' :   ;   . |  ; |--'  .   | '___     )" "\n" 
        R"(    '   ;  \; /  |   |   : |  | ,     '   ; : .'|    )" "\n" 
        R"(     \   \  ',  . \  |   : '  |/      '   | '/  :    )" "\n" 
        R"(      ;   :      ; | ;   | |`-'       |   :    /     )" "\n" 
        R"(       \   \ .'`--"  |   ;/            \   \ .'      )" "\n" 
        R"(        `---`        '---'              `---`        )" "\n" ;

    fmt::print(banner);
    fmt::print("Simple Quantum Monte Carlo Program\n");
    if(result["range"].as<bool>()) {
        fmt::print("Start range exploration...\n");
        fmt::print("{:>20s}\t{:>20s}\t{:>20s}\t{:>20s}\n", "c", "alpha", "mean", "std");
        auto cmin = result["cmin"].as<double>();
        auto cmax = result["cmax"].as<double>();
        auto ngridc = result["ngridc"].as<int>();
        auto alphamin = result["alphamin"].as<double>();
        auto alphamax = result["alphamax"].as<double>();
        auto ngridalpha = result["ngridalpha"].as<int>();
        double dr = result["step"].as<double>();
        int nstep = result["nstep"].as<int>();

        for(auto c: linspace(cmin, cmax, ngridc)) {
            for(auto alpha: linspace(alphamin, alphamax, ngridalpha)) {
                NaiveQMC<double> sampler(c, alpha, dr);
                double mean, std;
                //do {
                std::tie(mean, std) = sampler.sample(nstep);
                //} while(std > 0.01);
                fmt::print("{:>20.6f}\t{:>20.6f}\t{:>20.6f}\t{:>20.6f}\n", c, alpha, mean, std);
            }
        }

    } else {
        fmt::print("Single point calculation...\n");
        fmt::print("{:>20s}\t{:>20s}\t{:>20s}\t{:>20s}\n", "c", "alpha", "mean", "std");
        double c = result["c"].as<double>();
        double alpha = result["alpha"].as<double>();
        double dr = result["step"].as<double>();
        int nstep = result["nstep"].as<int>();
        NaiveQMC<double> sampler(c, alpha, dr);
        double mean, std;
        std::tie(mean, std) = sampler.sample(nstep);
        fmt::print("{:>20.6f}\t{:>20.6f}\t{:>20.6f}\t{:>20.6f}\n", c, alpha, mean, std);
    }

    return 0;
}