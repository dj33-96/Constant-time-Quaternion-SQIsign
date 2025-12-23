#include "quaternion_bench.h"
#include "poison.h"
void
bench_mod_core_lattice_add(ibz_vec_4_t (*generators)[8],
                           ibz_t *mod,
                           const quat_lattice_t *lat1,
                           const quat_lattice_t *lat2)
{
    ibz_mat_4x4_t tmp;
    ibz_t det1, det2, detprod;
    ibz_init(&det1);
    ibz_init(&det2);
    ibz_init(&detprod);
    ibz_mat_4x4_init(&tmp);
    ibz_mat_4x4_scalar_mul(&tmp, &(lat1->denom), &(lat2->basis));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&((*generators)[j][i]), &(tmp[i][j]));
        }
    }
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &det1, &tmp);
    ibz_mat_4x4_scalar_mul(&tmp, &(lat2->denom), &(lat1->basis));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&((*generators)[4 + j][i]), &(tmp[i][j]));
        }
    }
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &det2, &tmp);
    assert(!ibz_is_zero(&det1));
    assert(!ibz_is_zero(&det2));
    ibz_gcd(&detprod, &det1, &det2);
    ibz_copy(mod, &detprod);
    // ibz_mat_4xn_hnf_mod_core(&(res->basis),8,&generators,&detprod)
    ibz_mat_4x4_finalize(&tmp);
    ibz_finalize(&det1);
    ibz_finalize(&det2);
    ibz_finalize(&detprod);
}

// ibz_mat_4xn_hnf_mod_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t (*generators)[generator_number],
// const ibz_t *mod)
int
quat_bench_ibz_mat_4xn_hnf_mod_core(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2, cycles_ave;
    int i;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));
    uint8_t loop_size = 10;
    int64_t table_cycles[4] = { 0, 0, 0, 0 }; // Used to store the four values we want.

    // Define values for hnf_mod_core
    ibz_mat_4x4_t hnf;
    int generator_number = 8;
    ibz_t mod;
    ibz_vec_4_t generators[generator_number];
    int randret = 1;
    ibz_mat_4x4_init(&hnf);
    ibz_init(&mod);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));

    // Define values for lattice add
    quat_alg_t alg;
    quat_lattice_t sum, lat1, lat2;
    quat_lattice_init(&sum);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_alg_init_set(&alg, prime);
    QUAT_BENCH_CODE_1(iterations)
    // create random data for lattice addition
    randret = randret && quat_bench_helper_random_lattice(&lat1, intsize);
    randret = randret && quat_bench_helper_random_lattice(&lat2, intsize);
    if (!randret)
        break;
    // We run the first part of lattice addition until mod_core

    bench_mod_core_lattice_add(&generators, &mod, &lat1, &lat2);
    // We loop 10 times the operation and take the average of the 5-6-7-8 run.
    for (uint8_t i = 0; i < loop_size; i++) {
        QUAT_BENCH_CODE_START()
        // compute hnf_mod_core
       // poison(&(generators[0]), sizeof((generators[0])));
        // poison(&mod, sizeof(mod));
        ibz_mat_4xn_hnf_mod_core(&hnf, generator_number, &(generators[0]), &mod);
        QUAT_BENCH_CODE_STOP()
        QUAT_BENCH_SORTING_SAVE(i)
    }
    QUAT_BENCH_AVE_SET()
    runs = i;
    QUAT_BENCH_CODE_2("ct dim 4 hnf_mod_core", res)

    ibz_mat_4x4_finalize(&hnf);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    quat_alg_finalize(&alg);
    free(cycles_list);
    return (res);
}

int
quat_bench_ibz_mat_4xn_hnf_mod_core_old(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2, cycles_ave;
    int i;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));
    uint8_t loop_size = 10;
    int64_t table_cycles[4] = { 0, 0, 0, 0 }; // Used to store the four values we want.

    // Define values for hnf_mod_core
    ibz_mat_4x4_t hnf;
    int generator_number = 8;
    ibz_t mod;
    ibz_vec_4_t generators[8];
    int randret = 1;
    ibz_mat_4x4_init(&hnf);
    ibz_init(&mod);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));

    // Define values for lattice add
    quat_alg_t alg;
    quat_lattice_t sum, lat1, lat2;
    quat_lattice_init(&sum);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_alg_init_set(&alg, prime);
    QUAT_BENCH_CODE_1(iterations)
    // create random data for lattice addition
    randret = randret && quat_bench_helper_random_lattice(&lat1, intsize);
    randret = randret && quat_bench_helper_random_lattice(&lat2, intsize);
    if (!randret)
        break;
    // We run the first part of lattice addition until mod_core

    bench_mod_core_lattice_add(&generators, &mod, &lat1, &lat2);
    // We loop 10 times the operation and take the average of the 5-6-7-8 run.
    for (uint8_t i = 0; i < loop_size; i++) {
        QUAT_BENCH_CODE_START()
        // compute hnf_mod_core
        ibz_mat_4xn_hnf_mod_core_old(&hnf, generator_number, &(generators[0]), &mod);
        QUAT_BENCH_CODE_STOP()
        QUAT_BENCH_SORTING_SAVE(i)
    }
    QUAT_BENCH_AVE_SET()
    runs = i;
    QUAT_BENCH_CODE_2("old dim 4 hnf_mod_core", res)

    ibz_mat_4x4_finalize(&hnf);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    quat_alg_finalize(&alg);
    free(cycles_list);
    return (res);
}
long
quat_bench_dim4(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    res = res | quat_bench_ibz_mat_4xn_hnf_mod_core_old(iterations, intsize, prime);
    res = res | quat_bench_ibz_mat_4xn_hnf_mod_core(iterations, intsize, prime);

    return (res);
}
