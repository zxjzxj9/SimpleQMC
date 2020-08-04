#include <functional>
#include <fmt/core.h>
#include <cxxopts.hpp>
#include "simple_qmc.hpp"

int main(int argc, char** argv) {
    cxxopts::Options options("SimpleQMC", "Simple Quantum Monte Carlo Program");
    options.add_options()
        ("r,range", "Use range exploration", cxxopts::value<bool>()->default_value("false"))
        ("c", "Trial wave function parameter c", cxxopts::value<double>()->default_value("0.0"))
        ("alpha", "Trial wave function parameter alpha", cxxopts::value<double>()->default_value("1.0"))
        ("s,step", "Monte Carlo step size", cxxopts::value<double>()->default_value("0.001"))
        ("n,nstep", "Monte Carlo step size", cxxopts::value<int>()->default_value("10000"))
        ("h,help", "Print usage")
        ("minc", "Minimum parameter c")
        ("maxc", "Maximum parameter c")
        ("minalpha", "Minimum parameter alpha")
        ("maxalpha", "Maximum parameter alpha")
        ("ngridc", "Number of grid in c")
        ("ngridalpha", "Number of grid alpha")
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