# -*- coding: utf-8 -*-
"""
This file contains a proof-of-concept of constant-time Represent Integer function for q !=1 (Algorithm-6). 

	We define the following helper functions:
	1. ct_mod_exp is a function used for fast exponentiation. 
    These exponentiation functions are not constant-time and hence are not used in the constant-time C version of RepresentInteger.
	2. ct_cmov is a constant-time conditional move.
	3. ct_ts is the constant-time Tonelli-Shanks algorithm (Algorithm-2).
	4. c_hgcd is the constant-time halfgcd algorithm (Algorithm-4). 
	5. z_k and t_k are the functions to generate the sequence of z and t as per Table-1.
    6. There can be some issue with the isqrt function for large integer

This version only works for q in lvl1 and lvl3.

This was done on Python 3.9.
"""
import math
import random

p1 = 5 * 2**(248) - 1 #sqisign prime lvl1
p3 = 65 * 2**(376) - 1 #sqisign prime lvl3

def Rednorm(g,r,p,q):
    norm = g[0]**2 + q*g[1]**2 + p*(g[2]**2 + q*g[3]**2)
    rednorm = norm //(r**2)
    return rednorm

def nbits(a):
    if a == 0:
        return 0
    elif a>0:
        s = int(math.log2(int(a)))
        return s+1
    else:
        print("a is negative")

# Exponentiation functions
def safe_pow(base, exponent, modulus, chunk_bits=40):
    result = 1
    base_mod = base % modulus
    power_chunk = base_mod
    exp_bits = nbits(exponent)
    num_chunks = (exp_bits + chunk_bits - 1) // chunk_bits
    for i in range(num_chunks):
        lower = i * chunk_bits
        upper = min(exp_bits, (i + 1) * chunk_bits)
        chunk = (exponent >> lower) & ((1 << (upper - lower)) - 1)
        if chunk != 0:
            result = (result * pow(power_chunk, chunk, modulus)) % modulus
        power_chunk = pow(power_chunk, 1 << chunk_bits, modulus)
    return result

def ct_mod_exp(a, b, m, bit_length=None):
    """
    Constant-time modular exponentiation: computes a^b mod m.

    Parameters:
    a (int): base
    b (int): exponent
    m (int): modulus
    bit_length (int): fixed bit-length of the exponent for constant-time behavior

    Returns:
    int: a^b mod m
    """
    if bit_length is None:
        bit_length = nbits(b)

    result = 1
    a = a % m

    for i in reversed(range(bit_length)):
        result = (result * result) % m
        # Use bit masking instead of conditional branching
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
    is_different = condition ^ 1 #xor sign changes for python and sage
    mask = -int(bool(is_different))
    return (x & ~mask) | (y & mask)

def ct_ts(n,p):
    """
    constant-time Tonelli-Shanks used to find the square root of 'n' modulo 'p'. 
    This function is customised for the case when 2^e|(p-1) only for e = {2, 3}
    
    """
    x1 = 0
    d = p%16
    # f = ((d-1)//4+1)&1
    # print(f)
    e = 2
    q = (p-1)//2**e
    # print("ts value :",q,p)
    z = safe_pow(3, q, p)
    exp2 = (p - 2**e-1)//2**(e+1)
    y = ct_mod_exp(n, exp2, p)
    s = (n * y)%p
    t = (s * y)%p
    for k in range(e,1,-1):
        b = t
        for i in range(1,k-1):
            b = (b**2)%p
        s = ct_cmov(s,(s*z)%p,b)
        z = (z**2)%p
        t = ct_cmov(t,(t*z)%p,b)
    x2 = s    
    return (x1,x2)

def c_hgcd(a,b,printf=0):
    """
    (Near) constant-time halfgcd algorithm with debug variable 'printf'. 
     set printf == 1 in the debug mode to print all intermediate values.
    
    """
    d = max((nbits(abs(a))), (nbits(abs(b))))
    n = d + 2
    if printf == 1:
        print("nb of iterations : ",n)
    a0 = a
    b0 = b
    l = int(math.sqrt(a))
    while n>0:
        if printf == 1:
            print("input =",a0,b0)
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
        if printf == 1:
            print("sign =",sign,l-c)
            print("man,min =",c,d)
        if ka-kb-1>0:
            alpha = ka-kb-1
        else:
            alpha = 0
        if printf == 1:
            print("alpha =",alpha)
        a0 = c-(d2<<alpha)
        if printf == 1:
            print("output =",a0,b0)
        n -= 1
    return a0

def z_k(d,q,q8,k):
    result = 12*k + 6 * (q8 - 1)//(4)
    result = result + 4 * (d - 1)*(d - 2)//(1*2)
    result = result + 4 * (d)*(d - 2)//(1*(-1)) * (q-2)//(-1)
    return result

def t_k(d,k):
    result = 12*k + 1
    result = result + 8
    result = result + 4 * (d)*(d - 2)//(1*(-1))
    return result
# 4 * M - p(z²+qt²)

def ct_represent_integer(M,q,printf = 0):
    """
    Constant-time represent integer
    
    """
    xf = 0
    yf = 0
    zf = 0
    tf = 0
    found = 0
    d = M%3
    q3 = q%3
    q8 = q%8
    count = 0
    n1,n2 = 300,10*q
    for k1 in range(1, n1):
        if printf == 1:
            print("entering iteration", k1)
        z = z_k(d,q3,q8,k1)
        for k2 in range(1, n2):
            t = t_k(d, k2)
            N = 4*M - p1*(z**2 + q*t**2)
            if(N<0):
                print("N is negative at iteration of:", n1,n2) # For debug, in case of failure due to N becoming negative
            u0 = ct_ts(N-q, N)[1]
            n_f = (N - 1)//4
            # u0 = ct_mod_exp(3,n_f,N)
            x = c_hgcd(N,u0,0) 
            y_sq = (N - x**2)//q
            y = math.isqrt(y_sq)
            s = (4*M - (x**2 + q*y**2 + p1*(z**2 + q*t**2)))
            if printf == 1:
                print("s =", s)
            mask = 1 ^ ((s != 0) * 1)
            condition = (s == 0)
            xf = ct_cmov(x, xf, condition)
            if printf == 1:
                print("xf=", xf)
            yf = mask * y + (1 - mask) * yf 
            zf = ct_cmov(z, zf, condition)
            tf = ct_cmov(t, tf, condition)
            found = ct_cmov(1, found, condition)
    return (xf, yf, zf, tf)



# Select which q
Q = [5, 17, 37, 41, 53, 97] # q values for NIST-1
q = Q[0]
# Select lvl
p = p1
high = 2**(32) * p
m = random.randint(0, high)
M = 8*m + 7

# A working example, it should output: 
# (928366428836047150845353931280403482181253540177017627688, 234946754292183649317535424853135364927041239980638635551, 3586, 513)
# M = 284466028237988494763970178932217686760481835540054313610030844452355852304850148319078969945563537596721366979647
# p = p1 and q = 5

gamma = ct_represent_integer(M,q)
x,y,z,t = gamma[0],gamma[1],gamma[2],gamma[3]
print("M =", M)
print("gamma=", gamma)
N = Rednorm(gamma,2,p,q) #reduce norm check here as r = 2
print("norm =", N)
print("Valid rednorm of gamma :",N == M)
