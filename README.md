
# Galaxy Collision Simulation
## Overview
This project simulates simplified models of galaxy collisions, inspired by the Toomre model. While it’s a basic model compared to full-scale galaxy interaction simulations, this program allows for the recreation of morphological features observed in real galactic collisions. The approach focuses on using gravitational forces to study the dynamics between galaxies and stars during collisions.

## Theory & Approach
### Gravitational N-Body Problem
The gravitational N-body problem calculates gravitational interactions between multiple stars around galaxy cores. Each particle is treated as a point mass, experiencing forces from all other particles. A more detailed overview of the underlying math behind this code can be found in the galaxy_sim_overview.pdf file. 

### Numerical Solution
Using finite difference approximation, the problem is discretized over time to compute particle positions at each time step. Initial conditions are set to provide realistic starting conditions for galaxy and star positions and velocities.

### Galaxy Model (Toomre Model)
The Toomre model represents galaxies with a central core surrounded by massless stars in circular orbits, simplifying calculations by omitting star-to-star gravitational influences.

## Functions
nbodyaccn.m: Calculates gravitational acceleration for each particle in the galaxy.

tnbodyaccn.m: Tests gravitational acceleration and simulates circular orbits for two particles.

twobodysim.m: Simulates two-body interactions using the level and initial conditions.

t2bodyconverge.m: Checks convergence of the solution across multiple levels.

fastnbodyaccn.m: An optimized acceleration calculator that uses matrix operations for efficiency.

randcirclepts.m: Generates random star positions within a defined range around the core.

## Simulations
### Single Galaxy Simulation
The simulation for a single galaxy uses onegalaxysim.m, initializing random star positions around a central core to form a stable orbit. Galaxy motion as a whole is achieved by shifting all particle positions by a constant initial velocity.

### Two Galaxy Collision Simulation
twogalaxysim.m simulates two galaxies in motion toward each other, calculating gravitational interactions and visualizing collisions.

### Convergence Testing
The script t2bodyconverge.m uses multiple levels of time-step resolution to assess solution stability and accuracy. Differences between simulation levels are plotted to ensure convergence.

## Results
Detailed results, including figures from the simulations and convergence testing, showcase the morphological changes during galaxy collisions, offering insight into the dynamics of interacting galaxies.

## Dependencies
The simulation was created using MATLAB, and requires:

MATLAB (R2021a or later)
Basic linear algebra and plotting toolboxes
