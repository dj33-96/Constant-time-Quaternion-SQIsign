/** @file
 *
 * @authors Sina Schaeffler
 *
 * @brief Declarations for benchmarks of quaternion algebra operations
 */

#ifndef QUATERNION_BENCH_H
#define QUATERNION_BENCH_H

#include <quaternion.h>
#include <stdio.h>
#include <rng.h>
#include <bench.h>
#include "internal.h"

static int
cmpfunc(const void *a, const void *b)
{
    return (1 - 2 * (*(uint64_t **)a < *(uint64_t **)b));
}

#define QUAT_BENCH_CODE_1(r)                                                                                           \
    cycles = 0;                                                                                                        \
    for (i = 0; i < (r); ++i) {

#define QUAT_BENCH_CODE_START(r) cycles1 = cpucycles();

#define QUAT_BENCH_CODE_STOP(r) cycles2 = cpucycles();

#define QUAT_BENCH_SORTING_SAVE(i)                                                                                     \
    if (i == 5)                                                                                                        \
        table_cycles[0] = (cycles2 - cycles1);                                                                         \
    else if (i == 6)                                                                                                   \
        table_cycles[1] = (cycles2 - cycles1);                                                                         \
    else if (i == 7)                                                                                                   \
        table_cycles[2] = (cycles2 - cycles1);                                                                         \
    else if (i == 8)                                                                                                   \
        table_cycles[3] = (cycles2 - cycles1);

#define QUAT_BENCH_AVE_SET(r)                                                                                          \
    cycles_ave = (table_cycles[0] + table_cycles[1] + table_cycles[2] + table_cycles[3]) >> 2;                         \
    cycles2 = cycles_ave;                                                                                              \
    cycles1 = 0;

#define QUAT_BENCH_CODE_2(name, csv)                                                                                   \
    if (i < LIST_SIZE)                                                                                                 \
        cycles_list[i] = (cycles2 - cycles1);                                                                          \
    cycles = cycles + (cycles2 - cycles1);                                                                             \
    }                                                                                                                  \
    qsort(cycles_list, (runs < LIST_SIZE) ? runs : LIST_SIZE, sizeof(uint64_t), cmpfunc);                              \
    if (csv)                                                                                                           \
        printf("(with incorrect outputs)\n");                                                                          \
    else {                                                                                                             \
        printf("  %-20s-> median: %2" PRId64 ", average: %2" PRId64 " ",                                               \
               name,                                                                                                   \
               cycles_list[(runs < LIST_SIZE) ? runs / 2 : LIST_SIZE / 2],                                             \
               (cycles / runs));                                                                                       \
        printf("cycles\n");                                                                                            \
    }

/** @internal
 * @ingroup quat_helpers
 * @defgroup quat_bench Quaternion module benchmark functions
 * @{
 */
/** @internal
 * @defgroup quat_bench_helper Quaternion module benchmark helper functions
 * @{
 */
int quat_lattice_sublattice(const quat_lattice_t *sublat, const quat_lattice_t *lat);
int quat_bench_helper_random_lattice(quat_lattice_t *lat, int intsize);
/** @}
 */
long quat_bench_dim4(int iterations, int intsize, const ibz_t *prime);
long quat_bench_lattice(int iterations, int intsize, const ibz_t *prime);
long quat_bench_ideal(int iterations, int intsize, const ibz_t *prime);

long quat_bench_rep(int iterations, int intsize, const ibz_t *prime);

/** @}
 */

#endif
