# Quantum Monte Carlo

## 1. Simple QMC for the ground state of hydrogen atom

### 1.1 Theoretical backgrounds:
If we use the following trial function:
$$\psi(r) = (c r+1) e^{-\alpha r}$$

The the ration of the probabilty density at two different sites r1 and r2 becomes:
$$\frac{\rho(r_1)}{\rho(r_2)} = \frac{(c r_1+1)^2 e^{2 \alpha  (r_2-r_1)}}{(c r_2+1)^2}$$

We have the expectation of the energy:
$$\langle \psi|\hat{H}|\psi\rangle = \underset{r \sim \rho(r)}{E} \frac{\hat{H}\psi(r)}{\psi(r)}$$

Given the trial wavefunction, we have the energy expression as follows,

$$\frac{\hat{H}\psi(r)}{\psi(r)} = \frac{-c (r (\alpha  (\alpha  r-4)+2)+2)+\alpha  (2-\alpha  r)-2}{2 r (c r+1)}$$