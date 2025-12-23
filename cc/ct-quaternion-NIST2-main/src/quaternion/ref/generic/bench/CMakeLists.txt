set(SOURCE_FILES_QUATERNION_GENERIC_REF_BENCHS
    lattice.c
    ideal.c
    dim4.c
    repint.c
    quaternion_bench.c
)
add_executable(sqisign_bench_quaternion ${SOURCE_FILES_QUATERNION_GENERIC_REF_BENCHS})
target_link_libraries(sqisign_bench_quaternion ${LIB_QUATERNION} ${GMP} sqisign_common_test)
target_include_directories(sqisign_bench_quaternion PRIVATE ../internal_quaternion_headers ${INC_COMMON} ${INC_QUATERNION} ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/mini-gmp)

set(BM_BINS ${BM_BINS} sqisign_bench_quaternion CACHE INTERNAL "List of benchmark executables")
