#include <fmt/core.h>
#include <cxxopts.hpp>
#include "hydrogen.hpp"

int main(int argc, char** argv) {
    // H2MolQMC<double> h2qmc;
    PCoord<double> r1; 
    PCoord<double> r2;
    r1 << 0.40, 0.0, 0.0;
    r2 << -0.40, 0.0, 0.0;
    H2MolQMC<double> h2qmc(1.0, 0.5, 1.0, r1, r2, 0.1);
    cxxopts::Options options("HydrogenQMC", "Quantum Monte Carlo Program for Hydrogen Molecule");
    options.add_options()
        ("r,range", "Use range exploration", cxxopts::value<bool>()->default_value("false"))
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
    fmt::print("Simple Quantum Monte Carlo Program\n");
    return 0;
}