# Coursera, Spacecraft Dynamics

### Summary

This repository contains all of the functions and scripts I have written to-date while completing coursework for the [Kinematics: Describing the Motions of Spacecraft](https://www.coursera.org/learn/spacecraft-dynamics-kinematics) course on Coursera, taught by Prof. Hanspeter Schaub of the University of Colorado at Boulder. The course is divided into four modules, with the relevant files sorted into folders by the same name.

1. **Module 1 (COMPLETED):** Provides an introduction to the mathematics describing the motion of particles and rigid bodies, touching on angular velocities, rotating reference frames, and relating relative motion between reference frames. Most of the work for this module was done on paper

2. **Module 2 (COMPLETED):** Module 2 goes into greater depth on rigid body kinematics, specifically focusing on rotating reference frames via Direction Cosine matrices. For this module, I wrote numerous functions leveraging the SymPy Mechanics library, such as the following which:
    - Determine the equivalent Euler angles for any arbitrary rotation sequence, given another arbitrary rotation sequence and a set of angles for that initial sequence
    - Compute values using the *Spherical* Laws of Sine and Cosine
    - Add or subtract consecutive (symmetric) rotations to return a resultant equivalent

3. **Module 3 (IN PROGRESS):** Module 3 covers principle rotation vectors, quaternions, Rodrigues parameters, and other modern attitude descriptions. This module is currently in-work, and I have written functions to:
    - Return the Principle Rotation Vector given a DCM or Euler angles
    - Add or subtract principle rotation vectors to return a resultant equivalent

4. **Module 4:** [From Coursera] "This module covers how to take an instantaneous set of observations (sun heading, magnetic field direction, star direction, etc.) and compute a corresponding 3D attitude measure. The attitude determination methods covered include the TRIAD method, Devenport's q-method, QUEST as well as OLAE. The benefits and computation challenges are reviewed for each algorithm."
