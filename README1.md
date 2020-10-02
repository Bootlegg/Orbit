# Orbit
Solving orbital problems numerically.  

1. [Basic Stormer-Verlet for n-body problem](https://github.com/mintDan/Orbit#basic-stormer-verlet-for-n-body-problem)
2. [Barnes-Hut algorithm](https://github.com/mintDan/Orbit#barnes-hut-algorithm)
3. [Path through minimization of Action](https://github.com/mintDan/Orbit#path-through-minimization-of-action)

## Basic Stormer-Verlet for n-body problem
Simulation of Earth, Sun, Mars, Mercury, Jupiter and Venus.  
Currently the time step is loaded from a fortran module. The fortran module is converted to be callable within Python from f2py.
In the future the numerical steps will be made in fortran to compare speed against NumPy.  
Based on real parameters(mass, orbit speed) from Wikipedia.  
The force is calculated directly and thus takes 1/2*N*(N-1) evaluations. Symmetry, by Newton's 3rd law, reduces the number of direct evaluations required.  

### Simulation
![Orbit.png](https://github.com/mintDan/Orbit/blob/master/figs/Orbit.png)
### Virial Theorem and Total Energy
Testing the virial theorem  and conservation of total energy.  
![Virial.png](https://github.com/mintDan/Orbit/blob/master/figs/Virial.png)
### Total Angular Momentum
Testing conservation of total angular momentum
![AngularMomentum.png](https://github.com/mintDan/Orbit/blob/master/figs/AngularMomentum.png)
### Orbital periods with DFT
Seeing if we can observe the orbital periods of the planets through Discrete Fourier Transform. Looking here at speed of the Sun.  
The first big peak corresponds to Jupiter period of around 12 years, and gives as known the biggest impact on the sun, ie. the center of mass of the Solar System is mostly between Sun-Jupiter
The peak at around k = 410 gives a period of 1 year, which is the pull from Earth on the sun.  
A small peak at around k = 490 gives 0.6 years of orbital period, which is from Venus.  
The peaks for Mars and Mercury seem to be displaced to greater extent than the rest, maybe due to poor initial conditions,
giving an orbit period of 1.44 for Mars, where expected orbit period is around 1.88.
![DFT.png](https://github.com/mintDan/Orbit/blob/master/figs/DFT.png)

## Barnes-Hut algorithm
Uses the Barnes-Hut algo made with dictionaries. Uses a Quadtree for O(NlogN) efficiency, instead of direct summation O(N^2).
If a planet/particle is sufficiently far away, the force from distant masses can be approximated by their center of mass, using a multipole expansion.

![BH.png](https://github.com/mintDan/Orbit/blob/master/figs/BH.png)

## Path through minimization of Action
The true path travelled by an object will minimize the Action S, which truely is a phenomenal and jaw-dropping insight. The Lagrangian L is given by the potential energy subtracted from the kinetic energy.

![Lapprox.png](https://github.com/mintDan/Orbit/blob/master/figs/Lapprox.png)

The action is given by the integral

![Sint.png](https://github.com/mintDan/Orbit/blob/master/figs/Sint.png)

Which is approximated with discretization to

![Sapprox.png](https://github.com/mintDan/Orbit/blob/master/figs/Sapprox.png)

Varying the Action with respect to the path in x and y direction gives a system of non-linear equations

![Svary.png](https://github.com/mintDan/Orbit/blob/master/figs/Svary.png)

These systems of equations are solved by the Jacobi Method, or multidimensional Newton-Raphson method.

The result is seen below for Earth's orbit around the Sun.

![ActioNOrbit.png](https://github.com/Bootlegg/Orbit/blob/master/ActionOrbit.png)

The initial guessed path is a straight line from start to end position. As seen, the curved true path has a lower total action S.