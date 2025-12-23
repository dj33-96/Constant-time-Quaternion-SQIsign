#include "hnf_internal.h"
#include "internal.h"

// Small helper for integers
void
ibz_mod_not_zero(ibz_t *res, const ibz_t *x, const ibz_t *mod)
{
    ibz_t m, t;
    ibz_init(&m);
    ibz_init(&t);
    ibz_mod(&m, x, mod);
    ibz_set(&t, ibz_is_zero(&m));
    ibz_mul(&t, &t, mod);
    ibz_add(res, &m, &t);
    ibz_finalize(&m);
    ibz_finalize(&t);
}

// centered and rather positive then negative
void
ibz_centered_mod(ibz_t *remainder, const ibz_t *a, const ibz_t *mod)
{
    assert(ibz_cmp(mod, &ibz_const_zero) > 0);
    ibz_t tmp, d, t;
    ibz_init(&tmp);
    ibz_init(&d);
    ibz_init(&t);
    ibz_div_floor(&d, &tmp, mod, &ibz_const_two);
    ibz_mod_not_zero(&tmp, a, mod);
    ibz_set(&t, ibz_cmp(&tmp, &d) > 0);
    ibz_mul(&t, &t, mod);
    ibz_sub(remainder, &tmp, &t);
    ibz_finalize(&tmp);
    ibz_finalize(&d);
    ibz_finalize(&t);
}

// mpz_gcdext specification specifies unique outputs used here
// understand how to put the non-zero on the 1st one
void
ibz_xgcd_with_u_not_0(ibz_t *d, ibz_t *u, ibz_t *v, const ibz_t *x, const ibz_t *y)
{
    ibz_t q, x1, y1, t, old_x, old_y;
    int c;
    int c_zeros;
    ibz_init(&q);
    ibz_init(&t);
    ibz_init(&x1);
    ibz_init(&y1);
    ibz_init(&old_x);
    ibz_init(&old_y);
    ibz_copy(&x1, x);
    ibz_copy(&y1, y);
    ibz_copy(&old_x, x);
    ibz_copy(&old_y, y);
    // prepare inputs:
    // set x to 1 if both 0
    c_zeros = ibz_is_zero(&x1) & ibz_is_zero(&y1);
    ibz_set(&t, c_zeros);
    ibz_add(&x1, &x1, &t);
    // xgcd
    ibz_xgcd(d, u, v, &x1, &y1);
    // make sure u!=0 (v can be 0 if needed)
    ibz_copy(&x1, &old_x);
    c = ibz_is_zero(&y1);
    ibz_set(&t, c);
    ibz_add(&y1, &y1, &t);
    ibz_div(&q, &t, &x1, &y1);

    c = ibz_is_zero(u);
    ibz_set(&t, c);
    ibz_add(u, u, &t); // if u is 0, set u to 1 else leave it u -> u+c1
    ibz_sub(&t, &ibz_const_one, &q);
    ibz_conditional_assign(&t, v, &t, ibz_is_zero(&old_x)); // if x is zero, set t to v (v
                                                            // unchanged)
    ibz_copy(&q, v);
    ibz_conditional_assign(v, &t, &q, c); // prepare u-1 put it into u , else not: 1-q

    // try to match sign rules
    c = (ibz_cmp(d, &ibz_const_zero)) < 0;
    ibz_set(&t, c);
    ibz_mul(&t, &t, &ibz_const_two);
    ibz_sub(&t, &ibz_const_one, &t);
    ibz_mul(u, u, &t);
    ibz_mul(d, d, &t);
    ibz_mul(v, v, &t);

    ibz_mul(&t, u, x);
    ibz_mul(&q, v, y);
    ibz_add(&q, &t, &q);

    int test = !c_zeros && !(0 == ibz_cmp(&q, d));
    ibz_neg(&q, v);
    ibz_neg(&t, u);
    ibz_conditional_assign(u, &t, u, test);
    ibz_conditional_assign(v, &q, v, test);

    ibz_finalize(&x1);
    ibz_finalize(&y1);
    ibz_finalize(&old_x);
    ibz_finalize(&old_y);
    ibz_finalize(&q);
    ibz_finalize(&t);
}
