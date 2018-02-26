# Simple_MonteCarlo

This simplified Monte Carlo code was developed as student project for the "Monte Carlo method in particle transport simulations" course held at Aalto University. 

Its only intention was to be a learning experience, and should NOT under any circumstances be used for anything even close to real world applications! The results suffer from extreme simplifications done to make the project fit into the schedule of a course.

The code is presented in the presentation slides. The geometry part of the code is found in the Input, Universe, subspace, and Surface cpp-files. The nuclear physics calculations is in the NData.cpp file and the auxillary functions are in their respective files. Currently the way to input which simulations and calculations to run is hard codded in the main function, some examples provided in comment bracket. The nuclear cross sections are not included in the source code! It should be provided to the code in the following format:

```
SYM Z A AW T 
NNU 
E1 NU1 
E2 NU2 
... 
MT1 Q1 NE1 
E1 XS1 
E2 XS2 
...
MT2 Q2 NE2 
... 
```
