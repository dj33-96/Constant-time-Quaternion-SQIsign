#include <quaternion.h>
#include "internal.h"
#include <stdio.h>

#include <stdlib.h>

void
quat_lattice_O0_set(quat_lattice_t *O0)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(O0->basis[i][j]), 0);
        }
    }
    ibz_set(&(O0->denom), 2);
    ibz_set(&(O0->basis[0][0]), 2);
    ibz_set(&(O0->basis[1][1]), 2);
    ibz_set(&(O0->basis[2][2]), 1);
    ibz_set(&(O0->basis[1][2]), 1);
    ibz_set(&(O0->basis[3][3]), 1);
    ibz_set(&(O0->basis[0][3]), 1);
}

void
quat_lattice_O0_set_extremal(quat_p_extremal_maximal_order_t *O0)
{
    ibz_set(&O0->z.coord[1], 1);
    ibz_set(&O0->t.coord[2], 1);
    ibz_set(&O0->z.denom, 1);
    ibz_set(&O0->t.denom, 1);
    O0->q = 1;
    quat_lattice_O0_set(&(O0->order));
}

void
quat_order_elem_create(quat_alg_elem_t *elem,
                       const quat_p_extremal_maximal_order_t *order,
                       const ibz_vec_4_t *coeffs,
                       const quat_alg_t *Bpoo)
{

    // var dec
    quat_alg_elem_t quat_temp;

    // var init
    quat_alg_elem_init(&quat_temp);

    // elem = x
    quat_alg_scalar(elem, &(*coeffs)[0], &ibz_const_one);

    // quat_temp = i*y
    quat_alg_scalar(&quat_temp, &((*coeffs)[1]), &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->z, &quat_temp, Bpoo);

    // elem = x + i*y
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = z * j
    quat_alg_scalar(&quat_temp, &(*coeffs)[2], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);

    // elem = x + i* + z*j
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = t * j * i
    quat_alg_scalar(&quat_temp, &(*coeffs)[3], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);
    quat_alg_mul(&quat_temp, &quat_temp, &order->z, Bpoo);

    // elem =  x + i*y + j*z + j*i*t
    quat_alg_add(elem, elem, &quat_temp);

    quat_alg_elem_finalize(&quat_temp);
}

int
quat_old_represent_integer(quat_alg_elem_t *gamma,
                           const ibz_t *n_gamma,
                           int non_diag,
                           const quat_represent_integer_params_t *params)
{

    if (ibz_is_even(n_gamma)) {
        return 0;
    }
    // var dec
    int found;
    ibz_t counter;
    ibz_t cornacchia_target;
    ibz_t adjusted_n_gamma, q;
    ibz_t bound, sq_bound, temp;
    ibz_t test;
    ibz_vec_4_t coeffs; // coeffs = [x,y,z,t]
    quat_alg_elem_t quat_temp;

    if (non_diag)
        assert(params->order->q % 4 == 1);

    // var init
    found = 0;
    ibz_init(&counter);
    ibz_init(&bound);
    ibz_init(&test);
    ibz_init(&temp);
    ibz_init(&q);
    ibz_init(&sq_bound);

    ibz_vec_4_init(&coeffs);
    quat_alg_elem_init(&quat_temp);
    ibz_init(&adjusted_n_gamma);
    ibz_init(&cornacchia_target);

    ibz_set(&q, params->order->q);

    // this could be removed in the current state
    int standard_order = (params->order->q == 1);

    // adjusting the norm of gamma (multiplying by 4 to find a solution in an order of odd level)
    if (non_diag || standard_order) {
        ibz_mul(&adjusted_n_gamma, n_gamma, &ibz_const_two);
        ibz_mul(&adjusted_n_gamma, &adjusted_n_gamma, &ibz_const_two);
    } else {
        ibz_copy(&adjusted_n_gamma, n_gamma);
    }
    // computation of the first bound = sqrt (adjust_n_gamma / p - q)
    ibz_div(&sq_bound, &bound, &adjusted_n_gamma, &((params->algebra)->p));
    ibz_set(&temp, params->order->q);
    ibz_sub(&sq_bound, &sq_bound, &temp);
    ibz_sqrt_floor(&bound, &sq_bound);

    // the size of the search space is roughly n_gamma / (p√q)
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_sqrt_floor(&temp, &temp);
    ibz_div(&counter, &temp, &adjusted_n_gamma, &temp);

    // entering the main loop
    while (!found && ibz_cmp(&counter, &ibz_const_zero) != 0) {
        // decreasing the counter
        ibz_sub(&counter, &counter, &ibz_const_one);

        // we start by sampling the first coordinate
        ibz_rand_interval(&coeffs[2], &ibz_const_one, &bound);

        // then, we sample the second coordinate
        // computing the second bound in temp as sqrt( (adjust_n_gamma - p*coeffs[2]²)/qp )
        ibz_mul(&cornacchia_target, &coeffs[2], &coeffs[2]);
        ibz_mul(&temp, &cornacchia_target, &(params->algebra->p));
        ibz_sub(&temp, &adjusted_n_gamma, &temp);
        ibz_mul(&sq_bound, &q, &(params->algebra->p));
        ibz_div(&temp, &sq_bound, &temp, &sq_bound);
        ibz_sqrt_floor(&temp, &temp);

        if (ibz_cmp(&temp, &ibz_const_zero) == 0) {
            continue;
        }
        // sampling the second value
        ibz_rand_interval(&coeffs[3], &ibz_const_one, &temp);

        // compute cornacchia_target = n_gamma - p * (z² + q*t²)
        ibz_mul(&temp, &coeffs[3], &coeffs[3]);
        ibz_mul(&temp, &q, &temp);
        ibz_add(&cornacchia_target, &cornacchia_target, &temp);
        ibz_mul(&cornacchia_target, &cornacchia_target, &((params->algebra)->p));
        ibz_sub(&cornacchia_target, &adjusted_n_gamma, &cornacchia_target);
        assert(ibz_cmp(&cornacchia_target, &ibz_const_zero) > 0);

        // applying cornacchia
        if (ibz_probab_prime(&cornacchia_target, params->primality_test_iterations))
            found = ibz_cornacchia_prime(&(coeffs[0]), &(coeffs[1]), &q, &cornacchia_target);
        else
            found = 0;

        if (found && non_diag && standard_order) {
            // check that we can divide by two at least once
            // the treatmeat depends if the basis contains (1+j)/2 or (1+k)/2
            // we must have x = t mod 2 and y = z mod 2
            // if q=1 we can simply swap x and y
            if (ibz_is_odd(&coeffs[0]) != ibz_is_odd(&coeffs[3])) {
                ibz_swap(&coeffs[1], &coeffs[0]);
            }
            // we further check that (x-t)/2 = 1 mod 2 and (y-z)/2 = 1 mod 2 to ensure that the
            // resulting endomorphism will behave well for dim 2 computations
            found = found && ((ibz_get(&coeffs[0]) - ibz_get(&coeffs[3])) % 4 == 2) &&
                    ((ibz_get(&coeffs[1]) - ibz_get(&coeffs[2])) % 4 == 2);
        }
        if (found) {

#ifndef NDEBUG
            ibz_set(&temp, (params->order->q));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&test, &(coeffs[0]), &(coeffs[0]));
            ibz_add(&temp, &temp, &test);
            assert(0 == ibz_cmp(&temp, &cornacchia_target));

            ibz_mul(&cornacchia_target, &(coeffs[3]), &(coeffs[3]));
            ibz_mul(&cornacchia_target, &cornacchia_target, &(params->algebra->p));
            ibz_mul(&temp, &(coeffs[1]), &(coeffs[1]));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_set(&temp, (params->order->q));
            ibz_mul(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[0]), &coeffs[0]);
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[2]), &coeffs[2]);
            ibz_mul(&temp, &temp, &(params->algebra->p));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            assert(0 == ibz_cmp(&cornacchia_target, &adjusted_n_gamma));
#endif
            // translate x,y,z,t into the quaternion element gamma
            quat_order_elem_create(gamma, (params->order), &coeffs, (params->algebra));
#ifndef NDEBUG
            quat_alg_norm(&temp, &(coeffs[0]), gamma, (params->algebra));
            assert(ibz_is_one(&(coeffs[0])));
            assert(0 == ibz_cmp(&temp, &adjusted_n_gamma));
            assert(quat_lattice_contains(NULL, &((params->order)->order), gamma));
#endif
            // making gamma primitive
            // coeffs contains the coefficients of primitivized gamma in the basis of order
            quat_alg_make_primitive(&coeffs, &temp, gamma, &((params->order)->order));

            if (non_diag || standard_order)
                found = (ibz_cmp(&temp, &ibz_const_two) == 0);
            else
                found = (ibz_cmp(&temp, &ibz_const_one) == 0);
        }
    }

    if (found) {
        // new gamma
        ibz_mat_4x4_eval(&coeffs, &(((params->order)->order).basis), &coeffs);
        ibz_copy(&gamma->coord[0], &coeffs[0]);
        ibz_copy(&gamma->coord[1], &coeffs[1]);
        ibz_copy(&gamma->coord[2], &coeffs[2]);
        ibz_copy(&gamma->coord[3], &coeffs[3]);
        ibz_copy(&gamma->denom, &(((params->order)->order).denom));
    }
    // var finalize
    ibz_finalize(&counter);
    ibz_finalize(&bound);
    ibz_finalize(&temp);
    ibz_finalize(&sq_bound);
    ibz_vec_4_finalize(&coeffs);
    quat_alg_elem_finalize(&quat_temp);
    ibz_finalize(&adjusted_n_gamma);
    ibz_finalize(&cornacchia_target);
    ibz_finalize(&q);
    ibz_finalize(&test);

    return found;
}

int
prod_odd(int d, int a, int n)
{
    int64_t s = 1;
    for (int i = 1; i < n; i = i + 2)
        if (i != a)
            s = s * (d - i);
    return s;
}

int64_t
prod_dif(int a, int n)
{
    int64_t s = 1;
    for (int i = 1; i < n; i = i + 2)
        if (i != a)
            s = s * abs(a - i);
    return s;
}

// tonelli-shanks
void
tsss(ibz_t *x2, ibz_t *n, ibz_t *p)
{
    ibz_t d, q, z, exp2, y, s, t, b, nn, pp;
    ibz_init(&b);
    ibz_init(&d);
    ibz_init(&exp2);
    ibz_init(&nn);
    ibz_init(&pp);
    ibz_init(&q);
    ibz_init(&s);
    ibz_init(&t);
    ibz_init(&y);
    ibz_init(&z);
    ibz_copy(&nn, n);
    ibz_copy(&pp, p);
    int e, ed, f;
    ibz_set(&s, 16);
    ibz_mod(&d, &pp, &s);
    ibz_sub(&t, &d, &ibz_const_one);
    ibz_div_2exp(&d, &t, 2);
    f = (ibz_get(&d) & 1) ^ 1;
    e = 2 + f;
    ibz_sub(&s, &pp, &ibz_const_one);
    ibz_div_2exp(&q, &s, e);
    ibz_set(&s, 3);
    ibz_pow_mod(&z, &s, &q, &pp);
    ed = 1;
    for (int i = 0; i < e; i++) {
        ed = ed * 2;
    }
    ed = ed - 1;

    ibz_set(&s, ed);
    ibz_sub(&d, &pp, &s);

    ibz_div_2exp(&exp2, &d, e + 1);
    ibz_set(&s, e + 1);

    ibz_pow_mod(&y, n, &exp2, &pp);
    ibz_mul(&d, &nn, &y);
    ibz_mod(&s, &d, &pp);
    ibz_mul(&d, &s, &y);
    ibz_mod(&t, &d, &pp);
    for (int k = e; k > 1; k--) {
        ibz_copy(&b, &t);

        for (int i = 1; i < k - 1; i++) {
            ibz_mul(&d, &b, &b);
            ibz_mod(&b, &d, &pp);
        }
        ibz_mul(&d, &s, &z);
        ibz_mod(&y, &d, &pp);
        int8_t con = 2 * ibz_is_one(&b) - 1;
        ibz_ct_cmov(&s, &y, con);
        ibz_mul(&d, &z, &z);
        ibz_mod(&z, &d, &pp);
        ibz_mul(&d, &t, &z);
        ibz_mod(&y, &d, &pp);
        ibz_ct_cmov(&t, &y, con);
    }

    ibz_copy(x2, &s);
    ibz_finalize(&b);
    ibz_finalize(&d);
    ibz_finalize(&exp2);
    ibz_finalize(&nn);
    ibz_finalize(&pp);
    ibz_finalize(&q);
    ibz_finalize(&s);
    ibz_finalize(&t);
    ibz_finalize(&y);
    ibz_finalize(&z);
}

int
quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,
                                         const ibz_t *norm,
                                         int is_prime,
                                         const quat_represent_integer_params_t *params,
                                         const ibz_t *prime_cofactor)
{

    ibz_t n_temp, norm_d;
    ibz_t disc;
    quat_alg_elem_t gen, gen_rerand;
    int found = 0;
    ibz_init(&n_temp);
    ibz_init(&norm_d);
    ibz_init(&disc);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_rerand);

    // when the norm is prime we can be quite efficient
    // by avoiding to run represent integer
    // the first step is to generate one ideal of the correct norm
    if (is_prime) {

        // we find a quaternion element of norm divisible by norm
        while (!found) {
            // generating a trace-zero element at random
            ibz_set(&gen.coord[0], 0);
            ibz_sub(&n_temp, norm, &ibz_const_one);
            for (int i = 1; i < 4; i++)
                ibz_rand_interval(&gen.coord[i], &ibz_const_zero, &n_temp);

            // first, we compute the norm of the gen
            quat_alg_norm(&n_temp, &norm_d, &gen, (params->algebra));
            assert(ibz_is_one(&norm_d));

            // and finally the negation mod norm
            ibz_neg(&disc, &n_temp);
            ibz_mod(&disc, &disc, norm);
            // now we check that -n is a square mod norm
            // and if the square root exists we compute it
            found = ibz_sqrt_mod_p(&gen.coord[0], &disc, norm);
            found = found && !quat_alg_elem_is_zero(&gen);
        }
    } else {
        assert(prime_cofactor != NULL);
        // if it is not prime or we don't know if it is prime, we may just use represent integer
        // and use a precomputed prime as cofactor
        assert(!ibz_is_zero(norm));
        ibz_mul(&n_temp, prime_cofactor, norm);
        // printf("\n input ngamma to repint =");
        // ibz_print(&n_temp,10);
        found = quat_represent_integer(&gen, &n_temp, 0, params);
        // printf("\n the gen after rep int  is  for found =%d",found);
        // quat_alg_elem_print(&gen);
        found = found && !quat_alg_elem_is_zero(&gen);
    }
#ifndef NDEBUG
    if (found) {
        // first, we compute the norm of the gen
        quat_alg_norm(&n_temp, &norm_d, &gen, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_mod(&n_temp, &n_temp, norm);
        assert(ibz_cmp(&n_temp, &ibz_const_zero) == 0);
    }
#endif

    // now we just have to rerandomize the class of the ideal generated by gen
    found = 0;
    while (!found) {
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval(&gen_rerand.coord[i], &ibz_const_one, norm);
        }
        quat_alg_norm(&n_temp, &norm_d, &gen_rerand, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_gcd(&disc, &n_temp, norm);
        found = ibz_is_one(&disc);
        found = found && !quat_alg_elem_is_zero(&gen_rerand);
        //   printf("\n the genrand after while for found=%d ",found);
        // quat_alg_elem_print(&gen_rerand);
    }
    // printf("\n the gen before  mul is ");
    //   quat_alg_elem_print(&gen);
    //   printf("\n the gen_rerand before  mul is ");
    //  quat_alg_elem_print(&gen_rerand);
    quat_alg_mul(&gen, &gen, &gen_rerand, (params->algebra));
    // printf("\n the gen is ");
    // quat_alg_elem_print(&gen);
    // in both cases, whether norm is prime or not prime,
    // gen is not divisible by any integer factor of the target norm
    // therefore the call below will yield an ideal of the correct norm

    quat_lideal_create(lideal, &gen, norm, &((params->order)->order), (params->algebra));

    assert(ibz_cmp(norm, &(lideal->norm)) == 0);

    ibz_finalize(&n_temp);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_rerand);
    ibz_finalize(&norm_d);
    ibz_finalize(&disc);
    return (found);
}

void
quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_copy(&(*vec)[2], &el->coord[2]);
    ibz_add(&(*vec)[2], &(*vec)[2], &(*vec)[2]); // double (not optimal if el->denom is even...)
    ibz_copy(&(*vec)[3], &el->coord[3]);         // double (not optimal if el->denom is even...)
    ibz_add(&(*vec)[3], &(*vec)[3], &(*vec)[3]);
    ibz_sub(&(*vec)[0], &el->coord[0], &el->coord[3]);
    ibz_sub(&(*vec)[1], &el->coord[1], &el->coord[2]);

    assert(ibz_divides(&(*vec)[0], &el->denom));
    assert(ibz_divides(&(*vec)[1], &el->denom));
    assert(ibz_divides(&(*vec)[2], &el->denom));
    assert(ibz_divides(&(*vec)[3], &el->denom));

    ibz_div(&(*vec)[0], &tmp, &(*vec)[0], &el->denom);
    ibz_div(&(*vec)[1], &tmp, &(*vec)[1], &el->denom);
    ibz_div(&(*vec)[2], &tmp, &(*vec)[2], &el->denom);
    ibz_div(&(*vec)[3], &tmp, &(*vec)[3], &el->denom);

    ibz_finalize(&tmp);
}

int64_t
prod_diff(const int64_t *arr, int size, int64_t val, int64_t exclude)
{
    int64_t result = 1;
    for (int i = 0; i < size; i++) {
        if (arr[i] != exclude) {
            result *= (val - arr[i]);
        }
    }
    return result;
}

int64_t
floor_div(int64_t a, int64_t b)
{
    int64_t q = a / b;
    int64_t r = a % b;
    // If remainder is nonzero and signs differ, adjust
    if ((r != 0) && ((b < 0) != (a < 0))) {
        q--;
    }
    return q;
}

int64_t
correction_element(int64_t d,
                   int64_t q,
                   int64_t di,
                   int64_t qi,
                   int64_t coeff,
                   const int64_t *D,
                   int D_size,
                   const int64_t *Q,
                   int Q_size)
{
    int64_t d_num = prod_diff(D, D_size, d, di);
    int64_t d_den = prod_diff(D, D_size, di, di);
    int64_t d_basis = floor_div(d_num, d_den);

    int64_t q_num = prod_diff(Q, Q_size, q, qi);
    int64_t q_den = prod_diff(Q, Q_size, qi, qi);
    int64_t q_basis = floor_div(q_num, q_den);

    int64_t out = d_basis * q_basis * 16 * coeff;
    return out;
}


void z_kq(ibz_t *z, int64_t d, int64_t q, int64_t q8, int64_t k)
{
     int64_t result = 12*k+6*(q8-1)/4;

    result = result + 4 * (d - 1)*(d - 2)/(1*2);

    result = result + 4 * (d)*(d - 2)/(1*(-1)) * (q-2)/(-1);

    ibz_set(z, result);
}

void t_kq(ibz_t *t, int64_t d, int64_t k)
{
  
   int64_t  result = 12*k + 1;

    result = result + 8;

    result = result + 4 * (d)*(d - 2)/(1*(-1));
    ibz_set(t, result);
}
int
quat_ct_represent_integerwithq(quat_alg_elem_t *gamma,
                               const ibz_t *n_gamma,
                               int non_diag,
                               const quat_represent_integer_params_t *params)
{
    // variable dec
    int found, s;
    int64_t q = params->order->q;

    // variable dec
    uint16_t n11, n22;
    ibz_t pp;
    ibz_t x2;
    ibz_t u, u0, l, y2, a0, N;
    ibz_t x, y, z, w, t, t0, t1, t2, z2, dd, mod, mm,m4;
    ibz_t qq, nq;
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    ibz_init(&a0);
    ibz_init(&dd);
    ibz_init(&l);
    ibz_init(&mod);
    ibz_init(&mm);
    ibz_init(&N);
    ibz_init(&nq);
    ibz_init(&pp);
    ibz_init(&qq);
    ibz_init(&t);
    ibz_init(&t0);
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&u);
    ibz_init(&u0);
    ibz_init(&w);
    ibz_init(&x);
    ibz_init(&x2);
    ibz_init(&y);
    ibz_init(&y2);
    ibz_init(&z);
    ibz_init(&z2);
    ibz_set(&qq, q);
    ibz_set(&mod, 3);
    ibz_copy(&mm, n_gamma);
    ibz_t ibz_c4;
    ibz_set(&ibz_c4,4);
    ibz_mul(&m4,&mm,&ibz_c4);
    ibz_copy(&pp, &params->algebra->p);
    found = 0;
    s = 0;
    n11 = 30;
    //for q=5,n2=q*10
    n22 = 10*q;
    ibz_mod(&dd, &mm, &mod);
    int d3 = ibz_get(&dd);
    ibz_init(&dd);
    ibz_mod(&dd, &qq, &mod);
    int q3=ibz_get(&dd);
    ibz_set(&mod, 8);
    ibz_mod(&dd, &qq, &mod);
    int q8=ibz_get(&dd);
    //printf("\n q=%d,d=%d",q3,d3);
     int sq;
    // entering the main loop
   // printf("\n start \n");
    for (int i = 1; i < n11; i++) {
        z_kq(&z, d3, q3,q8, i);
        //printf("\n zz");
        //ibz_print(&z,10);
        for (int j = 1; j < n22; j++) {

          
            t_kq(&t, d3, j);
           // if (i==1){
//printf("\n tt");
     // ibz_print(&t,10);}
            // N=4M-p * (z² + qt²)
            ibz_mul(&z2, &z, &z);

            ibz_mul(&t2, &t, &t);
            ibz_mul(&t2, &t2, &qq);
            ibz_add(&t1, &z2, &t2);

            ibz_mul(&t2, &pp, &t1);
            
            ibz_sub(&N, &m4, &t2);

            // u0=tonelli-shanks(n-q,N)
            ibz_sub(&nq, &N, &qq);
            ibz_set(&u0, 0);
            tsss(&u0, &nq, &N);
// printf("\n u");
 //ibz_print(&u0,10);
            ibz_set(&l, 0);
            ibz_sqrt(&l, &N);

            ibz_constant_halfgcd(&x, &N, &u0, &l);
             //printf("\n x");
// ibz_print(&x,10);
            // y = sqrt(M - xÂ²) new version y = sqrt(N - x²)
            ibz_mul(&w, &x, &x);

            ibz_sub(&t1, &N, &w);
            ibz_set(&u, 0);
            ibz_div(&u, &w, &t1, &qq);
            ibz_sqrt(&y, &u);

//printf("\n y");
// ibz_print(&y,10);
           
           /* ibz_mul(&z2, &z, &z);
            ibz_mul(&t2, &t, &t);
            ibz_mul(&t2, &t2, &qq);
            ibz_add(&w, &z2, &t2);
            ibz_mul(&t1, &pp, &w);*/
           
            ibz_mul(&x2, &x, &x);
            ibz_mul(&y2, &y, &y);
            sq=(ibz_cmp(&y2,&u)==0);
            ibz_mul(&y2, &y2, &qq);
            ibz_add(&t0, &x2, &y2);
            ibz_add(&x2, &t0, &t2);
            
            ibz_sub(&t0, &m4, &x2);
          //  printf("\nk=%d j%d  sq%d sss=",i,j,sq);
           // ibz_print(&t0,10);
            s = (ibz_cmp(&t0, &ibz_const_zero) == 0);
           // printf("\n s=%d",s);
            found = found | s;
            if ((s == 1)&&(sq)) {
         
                found = 1;
                // translate x,y,z,t into the quaternion element gamma
                ibz_copy(&coeffs[0], &x);
                ibz_copy(&coeffs[1], &y);
                ibz_copy(&coeffs[2], &z);
                ibz_copy(&coeffs[3], &t);

    //  printf("\n found ");
       //   ibz_vec_4_print(&coeffs);
               // quat_order_elem_create(gamma, (params->order), &coeffs, (params->algebra));
            }
        }
    }
       if (found){ 
         quat_order_elem_create(gamma, (params->order), &coeffs, (params->algebra));
          quat_alg_make_primitive(&coeffs, &mm, gamma, &((params->order)->order));

       
       
       // Evaluate the solution into O0
        ibz_mat_4x4_eval(&coeffs, &(((params->order)->order).basis), &coeffs);

        // Copy the final solu
        //  printf("\n coeff after  eval");
        //     ibz_vec_4_print(&coeffs);
        // Copy the final solution
        ibz_copy(&gamma->coord[0], &coeffs[0]);
        ibz_copy(&gamma->coord[1], &coeffs[1]);
        ibz_copy(&gamma->coord[2], &coeffs[2]);
        ibz_copy(&gamma->coord[3], &coeffs[3]);
        ibz_copy(&gamma->denom, &(((params->order)->order).denom));
    }
    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&a0);
    ibz_finalize(&dd);
    ibz_finalize(&l);
    ibz_finalize(&mod);
    ibz_finalize(&mm);
    ibz_finalize(&N);
    ibz_finalize(&nq);
    ibz_finalize(&pp);
    ibz_finalize(&qq);
    ibz_finalize(&t);
    ibz_finalize(&t0);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&u);
    ibz_finalize(&u0);
    ibz_finalize(&w);
    ibz_finalize(&x);
    ibz_finalize(&x2);
    ibz_finalize(&y);
    ibz_finalize(&y2);
    ibz_finalize(&z);
    ibz_finalize(&z2);
    // printf("\n found=%d",found);
    return found;
}

// Final version of the RI

// First sequence for the z coordinate
void
z_kff(ibz_t *z, int d, int64_t d12, int k)
{
    int64_t term1 = (12 * k); //+864;
    int64_t term2 = 8 * (d12) * (d12 - 2) * (-1);
    int64_t res = (term1 + term2);
    ibz_set(z, res);
}
// Second sequence for the t coordinate
void
t_kff(ibz_t *t, int d, int64_t d12, int k)
{
    int64_t term1 = (12 * k) + 1;
    int64_t term2 = (8 * (d12 - 1) * d12) / 2;
    int64_t res = (term1 + term2);
    ibz_set(t, res);
}

// latest fast rep_int
int
quat_represent_integer(quat_alg_elem_t *gamma,
                       const ibz_t *n_gamma,
                       int non_diag,
                       const quat_represent_integer_params_t *params)
{
int found=0;
    int s;
     int standard_order=0;
     standard_order= (params->order->q == 1);
    if(standard_order){
    // variable dec and init
    
    ibz_t pp;
    uint16_t n11, n22;
    ibz_t x2;
    ibz_t u, u0, l, y2, a0, N;
    ibz_t x, y, z, w, t, t0, t1, t2, z2, dd, mod, mm, m4, mm4; //,yz;
    ibz_vec_4_t coeffs;
    ibz_t uu;
    ibz_t yyy;
    ibz_t xt;
    ibz_vec_4_init(&coeffs);
    ibz_init(&a0);
    ibz_init(&dd);
    ibz_init(&l);
    ibz_init(&mm);
    ibz_init(&m4);
    ibz_init(&mm4);
    ibz_init(&mod);
    ibz_init(&N);
    ibz_init(&pp);
    ibz_init(&t);
    ibz_init(&t0);
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&u);
    ibz_init(&u0);
    ibz_init(&uu);
    ibz_init(&w);
    ibz_init(&x);
    ibz_init(&xt);
    ibz_init(&x2);
    ibz_init(&y);
    ibz_init(&y2);
    ibz_init(&yyy);
    ibz_init(&z);
    ibz_init(&z2);
    
    ibz_copy(&pp, &params->algebra->p);
    ibz_copy(&mm, n_gamma);
    // printf("\n n-gamma =");
    // ibz_print(n_gamma,10);
    found = 0;
    s = 0;
    ibz_set(&m4, 4);
    n11 = 50;
    n22 = 100;
    ibz_set(&mod, 4);
    ibz_mod(&dd, &mm, &mod);
   
    int d = ibz_get(&dd);
    ibz_set(&mod, 12);
    ibz_set(&dd, 0);
    ibz_mod(&dd, &mm, &mod);
    int d12 = ibz_get(&dd);
    // entering the main loop
    ibz_set(&m4, 4);
    ibz_mul(&mm4, &m4, &mm);
    
   for (int i = 1; i < n11; i++) {
        z_kff(&z, d, d12, i);
        for (int j = 1; j < n22; j++) {
            t_kff(&t, d, d12, j);
            // N=M-p * (z² + t²)
            ibz_mul(&z2, &z, &z);
            ibz_mul(&t2, &t, &t);
            ibz_add(&t1, &z2, &t2);
            // M = 4 * N - p * (z² + w²)
            ibz_mul(&t2, &pp, &t1);
            ibz_sub(&N, &mm4, &t2);
            /*   printf("\n M=");
            //ibz_print(n_gamma,10);
           // printf("\n N=");
            ibz_print(&N,10);
            printf("\n t=");
            ibz_print(&t,10);
            printf("\n z=");
            ibz_print(&z,10);*/
            // printf(" \n N=");
            // ibz_print(&N,10);
            // if(ibz_cmp(&N,&ibz_const_zero)<0)
            // printf(" \n N is negative \n");
            ibz_sub(&u, &N, &ibz_const_one);
            ibz_set(&u0, 0);
            ibz_div_2exp(&uu, &u, 2);
            ibz_pow_mod(&u0, &ibz_const_three, &uu, &N);
            /************* without tonneli ***********/
            ibz_set(&l, 0);
            ibz_sqrt(&l, &N);
            ibz_constant_halfgcd(&x, &N, &u0, &l);
            ibz_mul(&x2, &x, &x);
            ibz_sub(&t1, &N, &x2);
            ibz_sqrt(&y, &t1);
            ibz_mul(&y2, &y, &y);
            int sq = (ibz_cmp(&y2, &t1) == 0);
           // ibz_mul(&z2, &z, &z);
            //ibz_mul(&t2, &t, &t);
            //ibz_add(&w, &z2, &t2);
            //ibz_mul(&t1, &pp, &t1);
           // ibz_mul(&x2, &x, &x);
            //ibz_mul(&y2, &y, &y);
            ibz_add(&t0, &x2, &y2);
            ibz_add(&x2, &t0, &t2);
            ibz_sub(&t0, &mm4, &x2);

            s = (ibz_cmp(&t0, &ibz_const_zero) == 0) & sq;
            /* if(s){
             printf(" \n check 4m ");
             ibz_print(&x2,10);
                printf(" \n ");
                 ibz_copy(&coeffs[0], &x);
                 ibz_copy(&coeffs[1], &y);
                 ibz_copy(&coeffs[2], &z);
                 ibz_copy(&coeffs[3], &t);
                 // ibz_vec_4_print(&coeffs);
                 printf("\n the norm is =");
                 ibz_print(&x2,10);
                 printf("\n");
                 ibz_div(&t0,&t1,&x2,n_gamma);
                 printf("\n the ratio=");
                 ibz_print(&t0,10);
                 printf("\n the remain=");
                 ibz_print(&t1,10);
               // norm = gamma[0]**2 + gamma[1]**2 + p*(gamma[2]**2 + gamma[3]**2);
             }


             ibz_set(&m4,0);
             ibz_set(&mod,4);
             ibz_sub(&m4,&x,&t);
             ibz_mod(&dd, &m4, &mod);
             ibz_set(&yz,0);
             ibz_sub(&m4,&y,&z);
             ibz_mod(&yz, &m4, &mod);*/
            // printf("\n non diag %d , standard order =%d",non_diag,standard_order);
            if (s && non_diag && standard_order) {
                /*  printf("\n m");
                ibz_print(n_gamma,10);
                printf("\n test= ");
                ibz_print(&t0,10);
                printf("\n 4m= ");
                ibz_print(&mm4,10);
                printf("\n x2= ");
                ibz_print(&x2,10);*/
                ibz_sub(&xt, &x, &t);
                // printf("\ncoeff before swap  ");
                // ibz_vec_4_print(&coeffs);
                if (ibz_is_odd(&xt)) {
                    // ibz_swap(&coeffs[0],&coeffs[1]);
                    ibz_swap(&y, &x);
                }
                // printf("\ncoeff  after swap");
                //  ibz_vec_4_print(&coeffs);
                int xtt, yzz;
                ibz_set(&mod, 4);
                ibz_sub(&m4, &x, &t);
                ibz_mod(&dd, &m4, &mod);
                // printf("\n x-t mod 4= %d",ibz_get(&dd));
                xtt = (ibz_cmp(&dd, &ibz_const_two) == 0);
                ibz_sub(&m4, &y, &z);
                ibz_mod(&dd, &m4, &mod);
                // printf("\n y-z mod 4 =%d",ibz_get(&dd));
                yzz = (ibz_cmp(&dd, &ibz_const_two) == 0);
                // ibz_mod(&yz, &m4, &mod);
                if ((yzz) && (xtt)) {
                    // printf(" \n condition xtt, ytt satisfied \n");
                    s = 1;

                } else {
                    s = 0;
                    //  printf(" \n condition not satisfied \n");
                }
            }
            // printf("\n s= %d",s);
            //
            // printf(" \n
            // *************************************************************************************************************found***********************************************************************************************************************************************************************************************
            // \n");
            if (s) {
                //  translate x,y,z,t into the quaternion element gamma
                ibz_copy(&coeffs[0], &x);
                ibz_copy(&coeffs[1], &y);
                ibz_copy(&coeffs[2], &z);
                ibz_copy(&coeffs[3], &t);
                // printf("\ncoeff before make ");
                //  ibz_vec_4_print(&coeffs);

                quat_order_elem_create(gamma, (params->order), &coeffs, (params->algebra));
                quat_alg_make_primitive(&coeffs, &mm, gamma, &((params->order)->order));

                /* printf("\n coeff after ");
                      ibz_vec_4_print(&coeffs);
                      printf("\n"); */
                found = 1;
            }
        }
    }
    if (found) {
        // printf("\n coeff before  eval");
        // ibz_vec_4_print(&coeffs);

        // Evaluate the solution into O0
        ibz_mat_4x4_eval(&coeffs, &(((params->order)->order).basis), &coeffs);

        // Copy the final solu
        //  printf("\n coeff after  eval");
        //     ibz_vec_4_print(&coeffs);
        // Copy the final solution
        ibz_copy(&gamma->coord[0], &coeffs[0]);
        ibz_copy(&gamma->coord[1], &coeffs[1]);
        ibz_copy(&gamma->coord[2], &coeffs[2]);
        ibz_copy(&gamma->coord[3], &coeffs[3]);
        ibz_copy(&gamma->denom, &(((params->order)->order).denom));
    }
     //printf("\n found%d",found);
    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&a0);
    ibz_finalize(&dd);
    ibz_finalize(&l);
    ibz_finalize(&mm);
    ibz_finalize(&mm4);
    ibz_finalize(&mod);
    ibz_finalize(&N);
    ibz_finalize(&pp);
    ibz_finalize(&t);
    ibz_finalize(&t0);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&u);
    ibz_finalize(&u0);
    ibz_finalize(&uu);
    ibz_finalize(&w);
    ibz_finalize(&x);
    ibz_finalize(&x2);
    ibz_finalize(&xt);
    ibz_finalize(&y);
    ibz_finalize(&y2);
    ibz_finalize(&yyy);
    ibz_finalize(&z);
    ibz_finalize(&z2);
    }else{
    quat_ct_represent_integerwithq( gamma,n_gamma, non_diag,params);

     //quat_old_represent_integer( gamma,n_gamma, non_diag,params);

        
    }
    //printf("\n found = %d",found);
    return found;
}
 
