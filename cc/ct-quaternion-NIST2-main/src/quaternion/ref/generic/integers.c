#include <quaternion.h>
#include "internal.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Random prime generation for tests
int
ibz_generate_random_prime(ibz_t *p, int is3mod4, int bitsize, int probability_test_iterations)
{
    assert(bitsize != 0);
    int found = 0;
    ibz_t two_pow, two_powp;

    ibz_init(&two_pow);
    ibz_init(&two_powp);
    ibz_pow(&two_pow, &ibz_const_two, (bitsize - 1) - (0 != is3mod4));
    ibz_pow(&two_powp, &ibz_const_two, bitsize - (0 != is3mod4));

    int cnt = 0;
    while (!found) {
        cnt++;
        if (cnt % 100000 == 0) {
            printf("Random prime generation is still running after %d attempts, this is not "
                   "normal! The expected number of attempts is %d \n",
                   cnt,
                   bitsize);
        }
        ibz_rand_interval(p, &two_pow, &two_powp);
        ibz_add(p, p, p);
        if (is3mod4) {
            ibz_add(p, p, p);
            ibz_add(p, &ibz_const_two, p);
        }
        ibz_add(p, &ibz_const_one, p);

        found = ibz_probab_prime(p, probability_test_iterations);
    }
    ibz_finalize(&two_pow);
    ibz_finalize(&two_powp);
    return found;
}

// solves x^2 + n y^2 == p for positive integers x, y
// assumes that p is prime and -n mod p is a square
int
ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p)
{
    ibz_t r0, r1, r2, a, prod;
    ibz_init(&r0);
    ibz_init(&r1);
    ibz_init(&r2);
    ibz_init(&a);
    ibz_init(&prod);

    int res = 0;

    // manage case p = 2 separately
    if (!ibz_cmp(p, &ibz_const_two)) {
        if (ibz_is_one(n)) {
            ibz_set(x, 1);
            ibz_set(y, 1);
            res = 1;
        }
        goto done;
    }
    // manage case p = n separately
    if (!ibz_cmp(p, n)) {
        ibz_set(x, 0);
        ibz_set(y, 1);
        res = 1;
        goto done;
    }

    // test coprimality (should always be ok in our cases)
    ibz_gcd(&r2, p, n);
    if (!ibz_is_one(&r2))
        goto done;

    // get sqrt of -n mod p
    ibz_neg(&r2, n);
    if (!ibz_sqrt_mod_p(&r2, &r2, p))
        goto done;

    // run loop
    ibz_copy(&prod, p);
    ibz_copy(&r1, p);
    ibz_copy(&r0, p);
    while (ibz_cmp(&prod, p) >= 0) {
        ibz_div(&a, &r0, &r2, &r1);
        ibz_mul(&prod, &r0, &r0);
        ibz_copy(&r2, &r1);
        ibz_copy(&r1, &r0);
    }
    // test if result is solution
    ibz_sub(&a, p, &prod);
    ibz_div(&a, &r2, &a, n);
    if (!ibz_is_zero(&r2))
        goto done;
    if (!ibz_sqrt(y, &a))
        goto done;

    ibz_copy(x, &r0);
    ibz_mul(&a, y, y);
    ibz_mul(&a, &a, n);
    ibz_add(&prod, &prod, &a);
    res = !ibz_cmp(&prod, p);

done:
    ibz_finalize(&r0);
    ibz_finalize(&r1);
    ibz_finalize(&r2);
    ibz_finalize(&a);
    ibz_finalize(&prod);
    return res;
}

// Constant-time cornacchia with helpers

void
ibz_ct_cmov(ibz_t *out, const ibz_t *in, int8_t sel)
{
    ibz_t mask, not_mask, t1, t2;

    ibz_init(&mask);
    ibz_init(&not_mask);
    ibz_init(&t1);
    ibz_init(&t2);

    // Compute mask = (sel ^ 1) ? -1 : 0
    // i.e., mask = -1 if sel == 0, mask = 0 if sel == 1
    ibz_set(&mask, (sel ^ 1) ? -1 : 0);

    // not_mask = ~mask (bitwise not)
    ibz_neg(&not_mask, &mask);
    ibz_sub(&not_mask, &not_mask, &ibz_const_one); // not_mask = -mask - 1 = ~mask

    // t1 = out & ~mask
    ibz_and(&t1, out, &not_mask);

    // t2 = in & mask
    ibz_and(&t2, in, &mask);

    // out = (out & ~mask) | (in & mask)
    ibz_or(out, &t1, &t2);

    ibz_finalize(&mask);
    ibz_finalize(&not_mask);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
}

void
ibz_max_2(ibz_t *max, const ibz_t *a, const ibz_t *b)
{
    uint8_t s, mova, movb;
    int8_t sp;
    ibz_t c;
    ibz_init(&c);
    ibz_sub(&c, a, b);
    // We need to get the sign of c.
    // s = c->sign;//Version with sign large integers
    sp = ibz_sgn(&c); // Version without sign large integers
    // Change sp to s here -> s needs to be the format s = 1 if c negative, else s = 0
    s = (sp & 1) & ((sp >> 1) & 1);
    // Only if sel of ibz_ct_cmov is for -1
    mova = s - 1;
    movb = -s;
    ibz_ct_cmov(max, a, mova);
    ibz_ct_cmov(max, b, movb);
    ibz_finalize(&c);
}

void
ibz_min_2(ibz_t *min, const ibz_t *a, const ibz_t *b)
{
    uint8_t s, mova, movb;
    int8_t sp;
    ibz_t c;
    ibz_init(&c);
    ibz_sub(&c, a, b);
    // We need to get the sign of c.
    // s = c->sign;//Version with sign large integers
    sp = ibz_sgn(&c); // Version without sign large integers
    // Change sp to s here -> s needs to be the format s = 1 if c negative, else s = 0
    s = (sp & 1) & ((sp >> 1) & 1);
    // Only if sel of ibz_ct_cmov is for -1
    movb = s - 1;
    mova = -s;
    ibz_ct_cmov(min, a, mova);
    ibz_ct_cmov(min, b, movb);
    ibz_finalize(&c);
}
void
ibz_constant_halfgcd(ibz_t *x, const ibz_t *M, const ibz_t *u, const ibz_t *l)
{
    ibz_t c, Mn, un, min, lc, sign_ibz, ll;
    int64_t si;
    uint64_t n, sm;
    int64_t sign;
    int64_t kp, kf, k2, k;
    ibz_init(&c);
    ibz_init(&lc);
    ibz_init(&ll);
    ibz_init(&min);
    ibz_init(&Mn);
    ibz_init(&sign_ibz);
    ibz_init(&un);

    ibz_copy(&Mn, M);
    ibz_copy(&ll, l);
    ibz_copy(&un, u);
    // Bitsize and loop computation
    sm = ibz_bitsize(M);

    n = sm + 2;

    while (n > 0) {

        si = ibz_cmp(&Mn, &un);
        if (si == 1) {
            ibz_copy(&c, &Mn);
            ibz_copy(&min, &un);
        } else {
            ibz_copy(&c, &un);
            ibz_copy(&min, &Mn);
        }
        k = ibz_bitsize(&c);
        ibz_sub(&lc, &ll, &c);
        // Works for the new structure
        sign = ibz_sgn(&lc);

        k2 = ibz_bitsize(&min);
        ibz_set(&sign_ibz, sign);

        ibz_mul(&min, &min, &sign_ibz);
        ibz_copy(&un, &min);

        kp = k - k2 - 1;

        kf = ~(kp >> 63);
        kp = kp & (kf);

        ibz_mul_2exp(&min, &min, kp);
        ibz_sub(&Mn, &c, &min);
        n = n - 1;
    }
    ibz_copy(x, &Mn);
    ibz_finalize(&c);
    ibz_finalize(&lc);
    ibz_finalize(&ll);
    ibz_finalize(&min);
    ibz_finalize(&Mn);
    ibz_finalize(&sign_ibz);
    ibz_finalize(&un);
}