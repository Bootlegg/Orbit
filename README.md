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
###Simulation
![Orbit.png](https://github.com/mintDan/Orbit/blob/master/figs/Orbit.png)
###Virial Theorem
Testing the virial theorem
![Virial.png](https://github.com/mintDan/Orbit/blob/master/figs/Virial.png)

## Barnes-Hut algorithm
Uses the Barnes-Hut algo made with dictionaries. Uses a Quadtree for O(NlogN) efficiency, instead of direct summation O(N^2).
If a planet/particle is sufficiently far away, the force from distant masses can be approximated by their center of mass, using a multipole expansion.

![BH.png](https://github.com/mintDan/Orbit/blob/master/figs/BH.png)

##Path through minimization of Action
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