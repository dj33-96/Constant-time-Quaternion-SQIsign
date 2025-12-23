#include "quaternion_bench.h"
#include "poison.h"
// Define default values
#define NIST_LVL 1
#define INTSIZE 3000
#define ITERATIONS 100

// run all benchmarks
// should add possibility to select some
int
main(int argc, char *argv[])
{
    int res = 0;
    int intsize = INTSIZE;
    int iterations = ITERATIONS;
    int lvl = NIST_LVL;
    ibz_t prime;
    ibz_init(&prime);
#ifndef NDEBUG
    printf("Compiled in Debug mode, so performences might be affected by test code and assertions\n");
#endif
    if (argc != 4 || !atoi(argv[1]) || !atoi(argv[2]) || !atoi(argv[3]) || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 ||
        ((atoi(argv[3]) != 1) && (atoi(argv[3]) != 3) && (atoi(argv[3]) != 5))) {
        printf("Invalid input: Expect iterations bitsize NIST level for prime choice from set {1,3,5}\n");
        printf("Proceed with default %d %d %d\n", ITERATIONS, INTSIZE, NIST_LVL);
    } else {
        iterations = atoi(argv[1]);
        intsize = atoi(argv[2]);
        lvl = atoi(argv[3]);
        printf("%d iterations, integers of %d bit and a prime for NIST %d\n", iterations, intsize, lvl);
    }
    if (lvl == 1)
        ibz_set_from_str(&prime, "2261564242916331941866620800950935700259179388000792266395655937654553313279", 10);
    else {
        if (lvl == 3)
            ibz_set_from_str(&prime,
                             "73925988009479995968240693512677722569677103469699654884481176219992919948243830045716827"
                             "3992115382691515787640831",
                             10);
        else {
            if (lvl == 5)
                ibz_set_from_str(&prime,
                                 "6005426571433468718907217573637180958155395333430694650046866082961925491649338992503"
                                 "5818409650398945598839482752258919095341032899110775086601744678911",
                                 10);
            else {
                printf("Invalid NIST level %d, should be 1, 3 or 5. Proceed with default: 1\n", lvl);
                ibz_set_from_str(
                    &prime, "2261564242916331941866620800950935700259179388000792266395655937654553313279", 10);
            }
        }
    }

    res = res | quat_bench_dim4(iterations, intsize, &prime); 
    res = res | quat_bench_ideal(iterations, intsize, &prime);
    res = res | quat_bench_rep(iterations, intsize, &prime);
    ibz_finalize(&prime);
    return (res);
}
