#include <fmt/core.h>
#include "hydrogen.hpp"

int main() {
    // H2MolQMC<double> h2qmc;
    PCoord<double> r1; 
    PCoord<double> r2;
    r1 << 0.40, 0.0, 0.0;
    r2 << -0.40, 0.0, 0.0;
    H2MolQMC<double> h2qmc(1.0, 0.5, 1.0, r1, r2, 0.1);
    return 0;
}