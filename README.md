# Orbit
Solving orbital problems numerically.
1. [Basic Stormer-Verlet for n-body problem](https://github.com/mintDan/MonteCarlo#mc-integration)
2. [Path through minimization of Action](https://github.com/mintDan/MonteCarlo#mc-calculate-pi)

## Basic Stormer-Verlet for n-body problem
Simulation of Earth, Sun, Mars, Mercury, Jupiter and Venus.  
Currently the time step is loaded from a fortran module. The fortran module is converted to be callable within Python from f2py.
In the future the numerical steps will be made in fortran to compare speed against NumPy.  
Based on real parameters(mass, orbit speed) from Wikipedia.  
###Simulation
![Orbit.png](https://github.com/Bootlegg/Orbit/blob/master/Orbit.png)
###Virial Theorem
Testing the virial theorem
![Virial.png](https://github.com/Bootlegg/Orbit/blob/master/Virial.png)

##Path through minimization of Action
The true path travelled by an object will minimize the Action S, which truely is a phenomenal and jaw-dropping insight. The Lagrangian L is given by the potential energy subtracted from the kinetic energy.

![Lapprox.png](https://github.com/Bootlegg/Orbit/blob/master/Lapprox.png)

The action is given by the integral

![Sint.png](https://github.com/Bootlegg/Orbit/blob/master/Sint.png)

Which is approximated with discretization to

![Sapprox.png](https://github.com/Bootlegg/Orbit/blob/master/Sapprox.png)

Varying the Action with respect to the path in x and y direction gives a system of non-linear equations

![Svary.png](https://github.com/Bootlegg/Orbit/blob/master/Svary.png)