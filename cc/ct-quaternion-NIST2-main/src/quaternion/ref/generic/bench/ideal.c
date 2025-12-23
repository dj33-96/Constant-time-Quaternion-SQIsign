#include "quaternion_bench.h"
#include "internal.h"
#include "lll_internals.h"
#include "poison.h"
int
quat_helper_order_discriminant(ibz_t *disc, const quat_lattice_t *order, const quat_alg_t *alg)
{
    int ok = 0;
    ibz_t det, sqr, div;
    ibz_mat_4x4_t transposed, norm, prod;
    ibz_init(&det);
    ibz_init(&sqr);
    ibz_init(&div);
    ibz_mat_4x4_init(&transposed);
    ibz_mat_4x4_init(&norm);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_transpose(&transposed, &(order->basis));
    // multiply gram matrix by 2 because of reduced trace
    ibz_mul(&(norm[3][3]), &ibz_const_two, &(alg->p));
    ibz_mul(&(norm[2][2]), &ibz_const_two, &(alg->p));
    ibz_copy(&(norm[1][1]), &ibz_const_two);
    ibz_copy(&(norm[0][0]), &ibz_const_two);
    ibz_mat_4x4_mul(&prod, &transposed, &norm);
    ibz_mat_4x4_mul(&prod, &prod, &(order->basis));
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &prod);
    ibz_mul(&div, &(order->denom), &(order->denom));
    ibz_mul(&div, &div, &div);
    ibz_mul(&div, &div, &div);
    ibz_div(&sqr, &div, &det, &div);
    ok = ibz_is_zero(&div);
    ok = ok && ibz_sqrt(disc, &sqr);
    ibz_finalize(&det);
    ibz_finalize(&div);
    ibz_finalize(&sqr);
    ibz_mat_4x4_finalize(&transposed);
    ibz_mat_4x4_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    return (ok);
}

int
quat_bench_helper_max_bitsize_in_lattice(const quat_lattice_t *lat)
{
    int maxi = ibz_bitsize(&(lat->denom));
    int bitsize;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            bitsize = ibz_bitsize(&(lat->basis[i][j]));
            if (maxi < bitsize) {
                maxi = bitsize;
            }
        }
    }
    return (maxi);
}

int
quat_bench_helper_random_ideal_of_order(quat_left_ideal_t *lideal,
                                        const quat_lattice_t *order,
                                        const quat_alg_t *alg,
                                        int intsize)
{
    ibz_t n, g, qn, qd;
    quat_alg_elem_t gen;
    int bitsize = intsize;
    int ind2 = 0;
    int randret = 1;
    ibz_init(&g);
    ibz_init(&n);
    ibz_init(&qn);
    ibz_init(&qd);
    quat_alg_elem_init(&gen);
    // find "random" ideal
    randret = randret && ibz_rand_interval_bits(&n, bitsize);
    if (!randret)
        goto end;
    ibz_abs(&n, &n);
    ibz_add(&n, &n, &ibz_const_one);
    // make random primitive element
    ibz_copy(&(gen.denom), &(order->denom));
    ind2 = 0;
    while (ind2 == 0) {
        randret = randret && ibz_rand_interval_bits(&(gen.coord[0]), bitsize) &&
                  ibz_rand_interval_bits(&(gen.coord[1]), bitsize) &&
                  ibz_rand_interval_bits(&(gen.coord[2]), bitsize) && ibz_rand_interval_bits(&(gen.coord[3]), bitsize);
        if (!randret)
            goto end;
        ibz_vec_4_content(&g, &(gen.coord));
        if (ibz_is_one(&g)) {
            ind2 = 1;
            ibz_mat_4x4_eval(&(gen.coord), &(order->basis), &(gen.coord));
            quat_alg_norm(&qn, &qd, &gen, alg);
            if (!ibz_is_one(&qd)) {
                ind2 = 0;
            }
        }
    }
    // compute ideal
    quat_lideal_create(lideal, &gen, &n, order, alg);
end:;
    ibz_finalize(&g);
    ibz_finalize(&n);
    ibz_finalize(&qn);
    ibz_finalize(&qd);
    quat_alg_elem_finalize(&gen);
    return (randret);
}

int
quat_bench_helper_random_order(quat_lattice_t *order, const quat_alg_t *alg, int intsize)
{
    quat_lattice_t lat;
    quat_left_ideal_t lideal;
    int bitsize = intsize;
    int randret = 1;
    quat_lattice_init(&lat);
    quat_left_ideal_init(&lideal);
    // find random lattice
    quat_lattice_O0_set(&lat);
    randret = quat_bench_helper_random_ideal_of_order(&lideal, &lat, alg, (int)(bitsize / 8));
    if (!randret)
        goto end;
    // compute right order
    quat_lideal_right_order(order, &lideal, alg);
#ifndef NDEBUG
    if (quat_bench_helper_max_bitsize_in_lattice(order) < intsize - 5)
        printf("Warning: Too small order %d %d\n", quat_bench_helper_max_bitsize_in_lattice(order), intsize);
    if (quat_bench_helper_max_bitsize_in_lattice(order) > intsize + 5)
        printf("Warning: Too large order %d %d\n", quat_bench_helper_max_bitsize_in_lattice(order), intsize);
#endif
end:;
    quat_lattice_finalize(&lat);
    quat_left_ideal_finalize(&lideal);
    return (randret);
}

// very slow since guessing
int
quat_bench_helper_random_maximal_order(quat_lattice_t *order, const quat_alg_t *alg, int intsize)
{
    ibz_t d;
    int randret = 1;
    int ind1 = 0;
    quat_lattice_t order0;
    quat_left_ideal_t lideal;
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&order0);
    ibz_init(&d);
    ibz_set(&(order0.basis[0][0]), 2);
    ibz_set(&(order0.basis[1][1]), 2);
    ibz_set(&(order0.basis[2][2]), 1);
    ibz_set(&(order0.basis[3][3]), 1);
    ibz_set(&(order0.basis[0][3]), 1);
    ibz_set(&(order0.basis[1][2]), 1);
    ibz_set(&(order0.denom), 2);
    while (quat_bench_helper_max_bitsize_in_lattice(order) < intsize) {
        ind1 = 0;
        while (ind1 == 0) {
            randret = randret && quat_bench_helper_random_ideal_of_order(&lideal, &order0, alg, 8);
            if (!randret)
                goto end;
            quat_lideal_right_order(order, &lideal, alg);
            quat_helper_order_discriminant(&d, order, alg);
            if (ibz_cmp(&(alg->p), &d) == 0) {
                ind1 = 1;
                ibz_mat_4x4_copy(&(order0.basis), &(order->basis));
                ibz_copy(&(order0.denom), &(order->denom));
            }
        }
    }
#ifndef NDEBUG
    if (quat_bench_helper_max_bitsize_in_lattice(order) < intsize - 10)
        printf("Warning: Too small max order %d %d\n", quat_bench_helper_max_bitsize_in_lattice(order), intsize);
    if (quat_bench_helper_max_bitsize_in_lattice(order) > intsize + 10)
        printf("Warning: Too large max order %d %d\n", quat_bench_helper_max_bitsize_in_lattice(order), intsize);
#endif
end:;
    ibz_finalize(&d);
    quat_lattice_finalize(&order0);
    quat_left_ideal_finalize(&lideal);
    return (randret);
}

// quat_lideal_right_order(quat_order_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg);
int
quat_bench_ideal_right_order(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2;
    int i;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_alg_t alg;
    quat_lattice_t order1, order2, test_order;
    quat_left_ideal_t lideal;
    int bitsize = intsize;
    int randret = 1;
    quat_lattice_init(&order1);
    quat_lattice_init(&order2);
    quat_lattice_init(&test_order);
    quat_left_ideal_init(&lideal);
    // use input prime
    quat_alg_init_set(&alg, prime);

    QUAT_BENCH_CODE_1(iterations)
    // use "random" maximal order
    randret = randret && quat_bench_helper_random_maximal_order(&order1, &alg, bitsize);
    if (!randret)
        goto end;
    randret = randret && quat_bench_helper_random_ideal_of_order(&lideal, &order1, &alg, bitsize);
    if (!randret)
        goto end;

    QUAT_BENCH_CODE_START()
    // compute right order
    quat_lideal_right_order(&order2, &lideal, &alg);
    QUAT_BENCH_CODE_STOP()
    // test if stabilizer
    quat_lattice_mul(&test_order, &(lideal.lattice), &order2, &alg);
    runs = i;
    res = res | !quat_lattice_sublattice(&test_order, &(lideal.lattice));
    QUAT_BENCH_CODE_2("ideal_right_order", res)

end:;
    if (!randret)
        printf("randomness failure in ideal_right_order\n");
    if (res)
        printf("test failure in ideal_right_order\n");
    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order1);
    quat_lattice_finalize(&order2);
    quat_lattice_finalize(&test_order);
    quat_left_ideal_finalize(&lideal);
    free(cycles_list);
    return (res);
}

long
quat_bench_lideal_generator(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking
    int64_t cycles, cycles1, cycles2;
    int i;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_elem_t gen;
    quat_left_ideal_t lideal;
    int bitsize = intsize;
    int randret = 1;
    int out = 0;
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&order);
    // use input prime
    quat_alg_init_set(&alg, prime);
    // Init of the output
    quat_alg_elem_init(&gen);
    QUAT_BENCH_CODE_1(iterations)
    // Set up random lideal by using a "random" maximal order
    randret = randret && quat_bench_helper_random_maximal_order(&order, &alg, bitsize);
    if (!randret)
        goto end;
    // Choose ideal
    randret = randret && quat_bench_helper_random_ideal_of_order(&lideal, &order, &alg, bitsize);
    if (!randret)
        goto end;
    // Start measurements
    QUAT_BENCH_CODE_START()
    //poison(&lideal, sizeof(lideal));
     // poison(&alg, sizeof(alg));
    out = quat_lideal_generator(&gen, &lideal, &alg);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = res | !out;
    QUAT_BENCH_CODE_2("ct ideal_generator", res)
end:;
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&order);
    free(cycles_list);
    return (res);
}

long
quat_bench_lideal_generator_old(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking
    int64_t cycles, cycles1, cycles2;
    int i;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_elem_t gen;
    quat_left_ideal_t lideal;
    int bitsize = intsize;
    int randret = 1;
    int out = 0;
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&order);
    // use input prime
    quat_alg_init_set(&alg, prime);
    // Init of the output
    quat_alg_elem_init(&gen);
    QUAT_BENCH_CODE_1(iterations)
    // Set up random lideal by using a "random" maximal order
    randret = randret && quat_bench_helper_random_maximal_order(&order, &alg, bitsize);
    if (!randret)
        goto end;
    // Choose ideal
    randret = randret && quat_bench_helper_random_ideal_of_order(&lideal, &order, &alg, bitsize);
    if (!randret)
        goto end;
    // Start measurements
    QUAT_BENCH_CODE_START()
    out = quat_lideal_generator_old(&gen, &lideal, &alg);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = res | !out;
    QUAT_BENCH_CODE_2("old ideal_generator", res)
end:;
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&order);
    free(cycles_list);
    return (res);
}
long
quat_bench_ideal(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;

    res = res | quat_bench_lideal_generator_old(iterations, intsize, prime);
    res = res | quat_bench_lideal_generator(iterations, intsize, prime);

    // res = res | quat_bench_ideal_right_order(iterations, intsize, prime);
    return (res);
}
