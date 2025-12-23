/** @file
 *
 * @authors Sina Schaeffler
 *
 * @brief Declarations for big integer functions only used in quaternion functions
 */

#ifndef INTBIG_INTERNAL_H
#define INTBIG_INTERNAL_H

#include "intbig.h"

/** @internal
 * @ingroup quat_helpers
 * @defgroup ibz_helper Internal integer functions (gmp-based)
 * @{
 */

/********************************************************************/

/**
 *  @brief r = a ^ b
 */
void ibz_xor(ibz_t *r, const ibz_t *a, const ibz_t *b);

/**
 *  @brief r = a | b
 */
void ibz_or(ibz_t *r, const ibz_t *a, const ibz_t *b);
/**
 *  @brief sign of an integer
 */
int8_t ibz_sgn(const ibz_t *a);

/**
 *  @brief 1 if op negative, 0 if positive
 */
int is_negative(const ibz_t *op);

/**
 *  @brief r = a & b
 */
void ibz_and(ibz_t *r, const ibz_t *a, const ibz_t *b);

/** @brief prod=a*2^exp
 */
void ibz_mul_2exp(ibz_t *prod, const ibz_t *a, uint64_t exp);

/**
 *  @brief Return 1 if a > 0, 0 if a = 0, and -1 if a < 0
 */
int8_t ibz_sgn(const ibz_t *a);

/** @brief Euclidean division of a by b
 *
 * Computes quotient, remainder so that remainder+quotient*b = a where 0<=|remainder|<|b|
 * The quotient is rounded towards minus infinity.
 */
void ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d);

/** @brief generate random value in [a, b]
 *  assumed that a >= 0, b >= 0 and a < b
 * @returns 1 on success, 0 on failure
 */
int ibz_rand_interval_i(ibz_t *rand, int32_t a, int32_t b);

/** @brief generate random value in [-2^m, 2^m]
 *  assumed that m > 0 and bitlength of m < 32 bit
 * @returns 1 on success, 0 on failure
 */
int ibz_rand_interval_bits(ibz_t *rand, uint32_t m);

/** @brief generate random value in in [-2^m, 2^m]
 */
void rand_exact_m_bits(ibz_t *result, gmp_randstate_t state, uint32_t m);

/** @brief set str to a string containing the representation of i in base
 *
 * Base should be 10 or 16
 *
 * str should be an array of length enough to store the representation of in
 * in base, which can be obtained by ibz_sizeinbase(i, base) + 2, where the 2
 * is for the sign and the null terminator
 *
 * Case for base 16 does not matter
 *
 * @returns  1 if the integer could be converted to a string, 0 otherwise
 */
int ibz_convert_to_str(const ibz_t *i, char *str, int base);

/** @brief print num in base to stdout
 *
 * Base should be 10 or 16
 */
void ibz_print(const ibz_t *num, int base);

/** @brief set i to integer contained in string when read as number in base
 *
 * Base should be 10 or 16, and the number should be written without ponctuation or whitespaces
 *
 * Case for base 16 does not matter
 *
 * @returns  1 if the string could be converted to an integer, 0 otherwise
 */
int ibz_set_from_str(ibz_t *i, const char *str, int base);

/**
 * @brief Probabilistic primality test
 *
 * @param n The number to test
 * @param reps Number of Miller-Rabin repetitions. The more, the slower and the less likely are
 * false positives
 * @return 1 if probably prime, 0 if certainly not prime, 2 if certainly prime
 *
 * Using GMP's implementation:
 *
 * From GMP's documentation: "This function performs some trial divisions, a Baillie-PSW probable
 * prime test, then reps-24 Miller-Rabin probabilistic primality tests."
 */
int ibz_probab_prime(const ibz_t *n, int reps);

/**
 * @brief GCD and Bézout coefficients u, v such that ua + bv = gcd
 *
 * @param gcd Output: Set to the gcd of a and b
 * @param u Output: integer such that ua+bv=gcd
 * @param v Output: Integer such that ua+bv=gcd
 * @param a
 * @param b
 */
void ibz_xgcd(ibz_t *gcd,
              ibz_t *u,
              ibz_t *v,
              const ibz_t *a,
              const ibz_t *b); // integers, dim4, test/integers, test/dim4

/**
 * @brief Square root modulo a prime
 *
 * @returns 1 if square root of a mod p exists and was computed, 0 otherwise
 * @param sqrt Output: Set to a square root of a mod p if any exist
 * @param a number of which a square root mod p is searched
 * @param p assumed prime
 */
int ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p);

/**
 * @brief Integer square root of a perfect square
 *
 * @returns 1 if an integer square root of a exists and was computed, 0 otherwise
 * @param sqrt Output: Set to a integer square root of a if any exist
 * @param a number of which an integer square root is searched
 */
int ibz_sqrt(ibz_t *sqrt, const ibz_t *a);

/**
 * @brief Legendre symbol of a mod p
 *
 * @returns Legendre symbol of a mod p
 * @param a
 * @param p assumed prime
 *
 * Uses GMP's implementation
 *
 * If output is 1, a is a square mod p, if -1, not. If 0, it is divisible by p
 */
int ibz_legendre(const ibz_t *a, const ibz_t *p);

/** @}
 */

/** @internal
 * @ingroup ibz_helper
 * @defgroup ibz_cond_moves Conditional moves
 * @{
 */

/** @brief if c then x else y
 *
 * @param res Output: if c, res = x, else res = y
 * @param x
 * @param y
 * @param c condition: must be 0 or 1
 */
void ibz_conditional_assign(ibz_t *res, const ibz_t *x, const ibz_t *y, int c);

/** @brief if c then exchange a and b, else do nothing
 *
 * @param a
 * @param b
 * @param c condition: must be 0 or 1
 */
void ibz_conditional_swap(ibz_t *a, ibz_t *b, int c);

/** @}
 */

// end of ibz_all
/** @}
 */
#endif
