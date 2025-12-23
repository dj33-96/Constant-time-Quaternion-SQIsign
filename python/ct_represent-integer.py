# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 16:12:16 2025

@author: David1
"""
"""
This file contains a proof-of-concept of constant-time Represent Integer function (Algorithm-5). 

	We define the following helper functions:
	1. ct_mod_exp is a function used for fast exponentiation. 
    These exponentiation functions are not constant-time and hence are not used in the constant-time C version of RepresentInteger.
	2. ct_cmov is a constant-time conditional move.
	3. c_hgcd is the constant-time halfgcd algorithm (Algorithm-4). 
	4. z_k and t_k are the functions to generate the sequence of z and t as per Table-1.
    5. There can be some issue with the sqrt function for large integer (line 202)

This was done on Python 3.9.
"""
import math
import random

#Prime selection
p1 = 5 * 2**(248) - 1  #sqisign prime lvl1
p3 = 65 * 2**(376) - 1 #sqisign prime lvl3

# constant-time half-gcd and helper functions

def Rednorm(g,r,p):
    norm = g[0]**2 + g[1]**2 + p*(g[2]**2 + g[3]**2)
    rednorm = norm //(r**2)
    return rednorm

def nbits(a):
    if a == 0:
        return 0
    else:
        s = int(math.log2(a))
        return s+1

def c_hgcd(a,b):
    """
    (near) Constant time half gcd algorithm.

    Parameters
    ----------
    a (int): input
    b (int): modular square root

    Returns
    -------
    a0 (int): a0^2 = 0 mod a

    """
    d = max(nbits(abs(a)), nbits(abs(b)))
    n = d + 2
    a0 = a
    b0 = b
    l = int(math.sqrt(a))
    while n>0:
        c = max(a0,b0)
        if (l-c)>0:
            sign = 0
        else:
            sign = 1        
        d = min(a0,b0)
        d2 = d * sign
        b0 = d2
        ka = nbits(c)
        kb = nbits(d)
        if ka-kb-1>0:
            alpha = ka-kb-1
        else:
            alpha = 0
        a0 = c-(d2<<alpha)
        n -= 1
    return a0


def ct_mod_exp(a, b, m):
    """
    Constant-time modular exponentiation: computes a^b mod m.

    Parameters:
    a (int): base
    b (int): exponent
    m (int): modulus

    Returns:
    int: a^b mod m
    """
    bit_length = b.bit_length()

    result = 1
    a = a % m
    for i in reversed(range(bit_length)):
        result = (result * result) % m
        mask = -((b >> i) & 1)  # -1 if bit is 1, 0 if bit is 0
        result = (result * a) % m if mask == -1 else result
    return result

def ct_cmov(x, y, condition):
    """
    Constant-time conditional move:
    Return y if condition != 1, else return x.

    Parameters:
    x (int): value to return if condition == 1
    y (int): value to return if condition != 1
    condition (int): condition to check

    Returns:
    int: y if condition != 1, else x
    """
    mask = -int(condition == 1)
    return (x & mask) | (y & ~mask)    

"""
	Functions to generate the sequences of z and t
"""  
def z_k(d3,k):
    term1 = 12*k
    term2 = 8 * (d3) *(d3-2) // (1*(-1)) + 72*12
    # print("z =",term1+term2)
    return term1 + term2

def t_k(d3,k):
    term1 = (12*k)+1
    term2 = 8 * (d3-1) *(d3) // (1*2)
    return term1 + term2

# Select SQIsign level
p = p1
#For random inputs
size = 280 # Set the max range
high = 2**size
low = p # Set the min range
m = random.randint(low, high)
M = 2*m + 1

# A working example, should output (for 30*20 loops) : 
# (222417499386461302576201484589313118050171, 239303315499048742612826190720423116111146, 648, 421)
# M = 27021525655500012275513012332961039414091864278697710519289531008301815110339137753

print("starting represent integer :")

def ct_represent_integer(M):
    d12 = M%12
    xf,yf,zf,tf = 0,0,0,0
    gamma = [0,0,0,0]
    # Two loop version, can also work with one loop.
    # Here we use 600 = 30 * 20 loops, can be changed
    for k1 in range(1, 30):
        z = z_k(d12,k1)
        for k2 in range(1, 20):
            t = t_k(d12,k2)
            N = 4*M - p*(z**2 + t**2)
            if N < 0:
                print("N is negative at iteration of k =", k1,k2) # For debug, in case of failure due to N becoming negative => increase input size
            if N%3!=2 or N%4!=1:
                print("error z,t, N=",N,"mod3, mod4",N%3,N%4) # For debug, NOT supposed to happen
            n_f = (N - 1)//4
            u0 = ct_mod_exp(3,n_f,N) # u0 = 3^(N-1)/4
            x = c_hgcd(N,u0)
            y = math.isqrt(N - x**2)
            test = 4*M - x**2-y**2-p*(z**2+t**2)
            # swap x and y
            mask = (x ^ t) & 1   
            tmp = mask * (x ^ y)
            x ^= tmp
            y ^= tmp
            cond = (test == 0) and (y == int(y)) and ((x-t)%4 == 2) and ((y-z)%4 == 2)
            gamma[0] = ct_cmov(x, gamma[0], cond)
            gamma[1] = ct_cmov(y, gamma[1], cond)
            gamma[2] = ct_cmov(z, gamma[2], cond)
            gamma[3] = ct_cmov(t, gamma[3], cond)
    return gamma

gamma = ct_represent_integer(M)
N = Rednorm(gamma,2,p) #reduce norm check here as r = 2
print("M =",M)
print("gamma =", gamma,"norm :",N,", is M ?",N==M)
