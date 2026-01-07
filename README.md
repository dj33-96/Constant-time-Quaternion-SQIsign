# Constant-time-Quaternion

This repository provides constant-time implementations of quaternion-based subroutines used in the SQIsign digital signature scheme, a compact and post-quantum secure signature protocol based on isogenies. These benchmarks were obtained using SQIsign's GMP-based integer arithmetic (not constant-time) and non-constant-time LLL.

## Overview


This repository is organized into three primary components:
* Python implementation: A reference implementation of the constant-time Represent Integer algorithm (algorithm-4 & 9) in Python.
* C implementation: A constant-time version of our quaternion arithmetic routines integrated into the SQIsign Round-2 submission.
* Qlapoti_C: A proof-of-concept implementation of a naive constant-time adaptation of the new Qlapoti algorithm (https://eprint.iacr.org/2025/1604)

```text
.
├── python/
│   └── ct_represent-integer.py         # Python reference implementation of ct-Represent Integer (q=1)
|   └── ct_represent-integer2.py        # Python reference implementation of ct-Represent Integer (q!=1)
├── Qlapoti-C/                          # C implementation, also includes Qlapoti
│      
└── cc/
    └── ct-quaternion-NIST2-main/       # Modified SQIsign Round 2 source with our constant-time quaternion functions
```

## For Python Implementation

1. The Python files (with .py extension) were run in Spyder version-6 with Python version-3.9.21
2. Please note that the desired output of the represent integer algorithm (on input value specified under "working example" in the code) is demonstrated in the comments for Python files


## Getting Started with the C Implementation

To compile and run the constant-time version of SQIsign: 

1. Navigate to the SQIsign source directory: 
    
    cd cc/ct-quaternion-NIST2-main

2. Follow the build instructions as described in the original SQIsign README.md

This version integrates constant-time logic only in the quaternion arithmetic layer, and for the functions described in the paper. 
All other parts of the SQIsign (including the integers used) are unchanged.

## Example

Build and run the benchmark for the quaternion module for SQIsign round 2 version using:

    cd build
    make
    ./src/quaternion/ref/generic/bench/sqisign_bench_quaternion 100 260 1

Where 100 is the number of runs for the benchmark, 260 being the average size of the input to the functions for NIST lvl1 we used for the benchmark.

## Benchmark Results

1. This section compares the performance of the old (SQIsign NIST-v2 code) and our constant-time (ct) implementations of the ideal_generator, hnf_mod_core and represent-integer functions in the C implementation.
2. Benchmarks are run over 100 iterations for each of the three NIST security levels, with both median and average clock cycles reported.
3. All these results were obtained using GMP integer functions.
4. Versions 'old' refers to the non-constant time SQIsign NIST-v2 code whereas 'ct' refers to our proposed constant-time quaternion functions.
5. Be aware that including constant-time integer arithmetic likely increases the runtime additionally.
6. Results for ct-represent-integer (q = 1 and q = 5) have been reported for the experimental and theoretical bound of the inner loops (value 'l' in Step-3, Algorithm-4 and details in sec.5.4 of the paper) for comparison.
7. Results for ct-represent-integer with q = 5 (one of the six values reported in SQIsign Round-2 submission document) is for the case of NIST lvl-1 and 3.
8. TimeCop for constant-time verification: We used TimeCop(https://crocs-muni.github.io/ct-tools/tutorials/timecop/) for verifying the constant-time behavior of the proposed quaternion algorithms. Apart from the (expected) leakages due to the use of underlying GMP-based integer arithmetic (Level-1 as per to Section-1 of the paper), TimeCop did not find any data-dependent branches or memory-accesses in our implementation of the quaternion-based algorithms(Level-2 as per Section-1 of the paper).


| Function               | Version | Median lvl 1  | Median lvl 3  | Median lvl 5  | Average lvl 1 | Average lvl 3 | Average lvl 5 |
|------------------------|---------|---------------|---------------|---------------|---------------|---------------|---------------|
| `ideal_generator`      | old     | 28,980        | 32,112        | 30,092        | 60,326        | 63,646        | 72,896        |
|                        | ct      | 1,567         | 1,588         | 1,573         | 1,579         | 1,654         | 1,622         |
| `hnf_mod_core`         | old     | 286,757       | 317,063       | 352,199       | 290,043       | 331,192       | 387,030       |
|                        | ct      | 543,034       | 546,829       | 554,583       | 542,698       | 546,304       | 555,377       |
| `represent-integer q=1`| old     | 66,047,343    | 77,561,202    | 159,838,719   | 97,571,124    | 145,479,164   | 211,963,034   |
| `q=1, 300 loops`       | ct      | 80,118,223    | 137,064,247   | 257,256,202   | 80,940,708    | 138,613,351   | 260,160,987   |
| `q=1, 6000 loops`      | ct      | 1,803,208,021 | 3,122,012,985 | 5,836,999,077 | 1,840,203,914 | 3,189,369,482 | 5,962,446,500 |
| `represent-integer q=5`| old     | 131,426,966   | 205,108,076   | 255,437,093   | 247,064,890   | 311,030,797   | 368,036,369   |
| `q=5, 300*q loops`     | ct      | 972,900,769   | 1,235,153,780 | -             | 993,622,687   | 1,247,785,628 | -             |


<img width="3000" height="1800" alt="loops_vs_cycles_RI" src="https://github.com/user-attachments/assets/bdf68091-85fa-4bf1-b636-3a12aa841174" />
This graph shows an almost linear increase in cycle counts with increasing number of inner loops (value 'l' in Step-3, Algorithm-4) of the ct-Represent Integer algorithm with q=1 (in case the graph is not diplayed, refer to the file loops_vs_cycles_RI.png).

We also include the theoretical failure rate curve as a function of the number of loops, as estimated in our paper, shown in the graph below.
![Failure_vs_loops](https://github.com/user-attachments/assets/db4ff2ad-db32-4f9d-bc7e-6443c34531ec)


These measurements were done on a laptop with the following specifications:
- Device: LENOVO_MT_20X1_BU_Think_FM_ThinkPad L14 Gen 2
- CPU: 11th Gen Intel(R) Core(TM) i5-1135G7 @ 2.40GHz
- Architecture: Intel x86-64 (64-bit)
- Memory: 16 GiB
- Operating System: Ubuntu 24.10 (Linux-based)
- Turbo Boost and Hyperthreading disabled
