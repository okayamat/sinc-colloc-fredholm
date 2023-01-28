# sinc-colloc-fredholm
Numerical solvers for Fredholm integral equations of the second kind by Sinc-collocation methods

## Overview
These programs solve four examples of Fredholm integral equations of the
second kind, conducted in [2, Example 1--4].

Those problems are solved by means of the following 4 methods:
* Original SE-Sinc-collocation method [1]
* New SE-Sinc-collocation method [2]
* Original DE-Sinc-collocation method [1]
* New DE-Sinc-collocation method [2]

The name of the program denotes the method and example number. For
example, SE_orig_ex2.c denotes the original SE-Sinc-collocation method
for example 2, and DE_new_ex3.c denotes the new DE-Sinc-collocation method
for example 3. LAPACK in Apple's Accelerate framework is used for
computation of the system of linear equations. If you want to use another
LAPACK library, modify make files according to your installation.

Each program solves those problems increasing N as N = 5, 10, 15, 20, ...,
and outputs N and maximum error over the target interval. Especially in
Example 3, in addition to N and maximum error, those programs also output
computation time.

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs for creating graphs are also stored in the directory.

computation environment:

OS: Mac OS X 10.12.6  
CPU: 1.7 GHz Intel Core i7  
Memory: 8 GB 1600 MHz DDR3  
Compiler: Apple LLVM version 9.0.0  
Library: Apple's Accelerate framework

## References
1. T. Okayama, T. Matsuo, M. Sugihara:
 Improvement of a Sinc-collocation method for Fredholm integral equations
 of the second kind, BIT Numer. Math., Vol. 51 (2011), pp. 339--366.
2. T. Okayama: Sinc-collocation methods with consistent collocation points
 for Fredholm integral equations of the second kind, arXiv, submitted.
