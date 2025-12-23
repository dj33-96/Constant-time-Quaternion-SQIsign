#include "quaternion_bench.h"

int
quat_lattice_sublattice(const quat_lattice_t *sublat, const quat_lattice_t *lat)
{
    int res = 0;
    quat_lattice_t sum;
    quat_lattice_init(&sum);
    quat_lattice_add(&sum, lat, sublat);
    res = quat_lattice_equal(&sum, lat);
    quat_lattice_finalize(&sum);
    return (res);
}

int
quat_bench_helper_random_lattice(quat_lattice_t *lat, int intsize)
{
    int randret = 1;
    ibz_t det;
    ibz_init(&det);
    ibz_set(&det, 0);
    while (ibz_is_zero(&det)) {
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                randret = randret && ibz_rand_interval_bits(&(lat->basis[a][b]), intsize);
            }
        }
        ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &(lat->basis));
    }
    ibz_set(&(lat->denom), 0);
    while (ibz_is_zero(&(lat->denom)))
        randret = randret && ibz_rand_interval_bits(&(lat->denom), intsize);
    if (!randret)
        goto end;
    quat_lattice_hnf(lat);
end:;
    ibz_finalize(&det);
    return (randret);
}

// void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int
quat_bench_lattice_add(int iterations, int intsize, const ibz_t *prime)
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
    quat_lattice_t sum, lat1, lat2;
    int randret = 1;
    quat_lattice_init(&sum);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    // use input prime
    quat_alg_init_set(&alg, prime);

    QUAT_BENCH_CODE_1(iterations)
    // create random lattices
    randret = randret && quat_bench_helper_random_lattice(&lat1, intsize);
    randret = randret && quat_bench_helper_random_lattice(&lat2, intsize);
    if (!randret)
        break;
    QUAT_BENCH_CODE_START()
    // compute sum
    quat_lattice_add(&sum, &lat1, &lat2);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = res | !quat_lattice_sublattice(&lat1, &sum);
    res = res | !quat_lattice_sublattice(&lat2, &sum);
    QUAT_BENCH_CODE_2("lattice_add", res)

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    quat_alg_finalize(&alg);
    free(cycles_list);
    return (res);
}

// void quat_lattice_mul(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2, const quat_alg_t
// *alg);
int
quat_bench_lattice_mul(int iterations, int intsize, const ibz_t *prime)
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
    quat_lattice_t prod, lat1, lat2;
    int randret = 1;
    quat_lattice_init(&prod);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    // use input prime
    quat_alg_init_set(&alg, prime);

    QUAT_BENCH_CODE_1(iterations)
    // create random lattices
    randret = randret && quat_bench_helper_random_lattice(&lat1, intsize);
    randret = randret && quat_bench_helper_random_lattice(&lat2, intsize);
    if (!randret)
        break;
    QUAT_BENCH_CODE_START()
    // compute sum
    quat_lattice_mul(&prod, &lat1, &lat2, &alg);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    // try to find adequate partial test
    QUAT_BENCH_CODE_2("lattice_mul", res)

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod);
    quat_alg_finalize(&alg);
    free(cycles_list);
    return (res);
}

// void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int
quat_bench_lattice_intersect(int iterations, int intsize, const ibz_t *prime)
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
    quat_lattice_t inter, lat1, lat2;
    int randret = 1;
    quat_lattice_init(&inter);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    // use input prime
    quat_alg_init_set(&alg, prime);

    QUAT_BENCH_CODE_1(iterations)
    // create random lattices
    randret = randret && quat_bench_helper_random_lattice(&lat1, intsize);
    randret = randret && quat_bench_helper_random_lattice(&lat2, intsize);
    if (!randret)
        break;
    QUAT_BENCH_CODE_START()
    // compute sum
    quat_lattice_intersect(&inter, &lat1, &lat2);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = res | !quat_lattice_sublattice(&inter, &lat1);
    res = res | !quat_lattice_sublattice(&inter, &lat2);
    QUAT_BENCH_CODE_2("lattice_intersect", res)

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&inter);
    quat_alg_finalize(&alg);
    free(cycles_list);
    return (res);
}

long
quat_bench_lattice(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    res = res | quat_bench_lattice_add(iterations, intsize, prime);
    res = res | quat_bench_lattice_mul(iterations, intsize, prime);
    res = res | quat_bench_lattice_intersect(iterations, intsize, prime);
    return (res);
}
