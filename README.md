# Orbit
Solving orbital problems numerically.
## nbody problem.
Simulation of Earth, Sun, Mars, Mercury, Jupiter and Venus.  
Currently the time step is loaded from a fortran module. The fortran module is converted to be callable within Python from f2py.
In the future the numerical steps will be made in fortran to compare speed against NumPy.  
Based on real parameters(mass, orbit speed) from Wikipedia.
###Simulation
![Orbit.png](https://github.com/Bootlegg/Orbit/blob/master/Orbit.png)
###Virial Theorem
Testing the virial theorem
![Virial.png](https://github.com/Bootlegg/Orbit/blob/master/Virial.png)