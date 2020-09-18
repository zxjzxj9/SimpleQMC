#include <fmt/core.h>
#include <cxxopts.hpp>
#include <vector>
#include "hydrogen.hpp"

const double Hartree = 27.21138602;

int main(int argc, char** argv) {
    // H2MolQMC<double> h2qmc;
    // r1 << 0.40, 0.0, 0.0;
    // r2 << -0.40, 0.0, 0.0;
    // H2MolQMC<double> h2qmc(1.0, 0.5, 1.0, r1, r2, 0.1);
    cxxopts::Options options("HydrogenQMC", "Quantum Monte Carlo Program for Hydrogen Molecule");
    options.add_options()
        ("F", "Jastrow factor parameter", cxxopts::value<double>()->default_value("1.0"))
        ("c", "Trial wave function parameter c", cxxopts::value<double>()->default_value("0.0"))
        ("alpha", "Trial wave function parameter alpha", cxxopts::value<double>()->default_value("1.0"))
        ("s,step", "Monte Carlo step size", cxxopts::value<double>()->default_value("1.0"))
        ("n,nstep", "Monte Carlo step size", cxxopts::value<int>()->default_value("1000000"))
        ("r1", "First atom coordinateds", cxxopts::value<std::vector<double>>()->default_value("0.5 0.0 0.0"))
        ("r2", "Second atom coordinateds", cxxopts::value<std::vector<double>>()->default_value("-0.5 0.0 0.0"))
        ("h,help", "Print usage")
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
    fmt::print("Simple Quantum Monte Carlo Program for H2\n");

    PCoord<double> r1; 
    PCoord<double> r2;
    auto r1_opt = result["r1"].as<std::vector<double>>();
    auto r2_opt = result["r2"].as<std::vector<double>>();
    r1<<r1_opt[0],r1_opt[1],r1_opt[2];
    r2<<r2_opt[0],r2_opt[1],r2_opt[2];
    auto F = result["F"].as<double>();
    auto c = result["c"].as<double>();
    auto alpha = result["alpha"].as<double>();
    auto s = result["step"].as<double>();
    auto nstep = result["nstep"].as<int>();
    H2MolQMC<double> h2qmc(F, c, alpha, r1, r2, s);
    double energy, energy_std;
    std::tie(energy, energy_std) = h2qmc.sample(nstep);
    fmt::print("{:>20s}\t{:>20s}\n", "Energy (eV rel. 2H)", "Energy Std");
    fmt::print("{:>20.8f}\t{:>20.8f}\n", (energy+1)*Hartree, energy_std*Hartree);
    return 0;
}