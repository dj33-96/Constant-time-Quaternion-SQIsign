# SQIsign implementation using a naice constant-time qlapoti



## Origin of code

This implementation is based on an existing [Qlapoti C implementation](https://github.com/KULeuven-COSIC/Qlapoti) by Giacomo Borin, Maria Corte-Real Santos, Jonathan Komada Eriksen, Riccardo Invernizzi, Marzio Mula, Sina Schaeffler, and Frederik Vercauteren. 


We mostly modified the file src/quaternion/ref/generic/qlapoti.c which cntains the norm equation solving, and renamed their README to Qlapoti_README.md. Some precomputation files might also have been modified by our experiments.

## Requirements and how to run

For instructions on how to run this, look at the SQIsign_README.md.

## Constant-time subroutines used or not used

This code is only meant as a proof-of-concept for making the qlapoti algrithm naively constant-time, essentially by iterationg until the bound given in our baseline implementation (which is not the theorectical bound from the Qlapoti paper, but far worse due to additional constraints). 

In particular, we did not use any lower-level constant-tme subroutines, such as HNF, lattice reduction, cornacchia's algorithm for prime-norm inputs, or even constant-time arithmetic. We did make sure however to only use algorithms for which known constant-time alternatives exist.

As a consequence, the resulting implementation is far from constant-time, despite its name.