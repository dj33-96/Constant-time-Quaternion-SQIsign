#include "quaternion_bench.h"
#include "internal.h"
#include "poison.h"
long
quat_bench_repint(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2;
    int i;
    ibz_t n_gamma, ga;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    quat_alg_elem_init(&gamma);
    quat_alg_init_set(&Bpoo, prime);

    
    quat_p_extremal_maximal_order_t order;

    quat_represent_integer_params_t params;
    quat_lattice_init(&order.order);
    quat_alg_elem_init(&order.z);
    quat_alg_elem_init(&order.t);
    quat_lattice_O0_set_extremal(&order);
    params.primality_test_iterations = 1;
    params.algebra = &Bpoo;
    params.order = &order;
    // assert(quat_test_special_extremal_setup(params, alg));

    // quat_test_set_params_standard(&params, &order, &Bpoo);
    QUAT_BENCH_CODE_1(iterations);
    ibz_init(&n_gamma);
    ibz_init(&ga);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    rand_exact_m_bits(&n_gamma, state, intsize + 50);
    ibz_mul(&ga, &n_gamma, &ibz_const_two);
    ibz_add(&n_gamma, &ga, &ibz_const_one);
    // ibz_set_from_str(&n_gamma,
    // "203562952599819965249873094651830260029570152201138642387625227520548342261511078168880053295862782462007489033701196404662274340956042194528145046088681493167",
    // 10);
    //  }
    QUAT_BENCH_CODE_START()
     //poison(&n_gamma, sizeof(n_gamma));
     // poison(&params, sizeof(params));
    quat_represent_integer(&gamma, &n_gamma,0, &params);
    QUAT_BENCH_CODE_STOP()

    runs = i;
    res = 0;
    QUAT_BENCH_CODE_2("ct represent_integer", res);

    return 1;
}

long
quat_bench_repintold(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2;
    int i;
    int out;

    ibz_t n_gamma, ga;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_represent_integer_params_t params;

    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    quat_alg_elem_init(&gamma);
    quat_alg_init_set(&Bpoo, prime);
    quat_p_extremal_maximal_order_t order;
    quat_lattice_init(&order.order);
    quat_alg_elem_init(&order.z);
    quat_alg_elem_init(&order.t);
    quat_lattice_O0_set_extremal(&order);
    // quat_alg_elem_print(&order.t);
    out = 1;

    params.primality_test_iterations = 1;
    params.algebra = &Bpoo;
    params.order = &order;
    // assert(quat_test_special_extremal_setup(params, alg));
    // n_gamma= 3573960250147475915502838919373697042179781142194761891096159111893170748566397047
    ibz_init(&n_gamma);
    ibz_init(&ga);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    QUAT_BENCH_CODE_1(iterations);
    rand_exact_m_bits(&n_gamma, state, intsize + 50);
    ibz_mul(&ga, &n_gamma, &ibz_const_two);
    ibz_add(&n_gamma, &ga, &ibz_const_one);
    QUAT_BENCH_CODE_START()
    quat_old_represent_integer(&gamma, &n_gamma, 0, &params);

    QUAT_BENCH_CODE_STOP()

    runs = i;
    res = res | !out;
    QUAT_BENCH_CODE_2("old represent_integer", res);
    return 1;
}

long
quat_bench_repintvarq(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2;
    int i;
    ibz_t n_gamma,ga;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    quat_alg_elem_init(&gamma);
    quat_alg_init_set(&Bpoo, prime);

    // n_gamma= 3573960250147475915502838919373697042179781142194761891096159111893170748566397047

    quat_p_extremal_maximal_order_t order;

    quat_represent_integer_params_t params;
    quat_lattice_init(&order.order);
    quat_alg_elem_init(&order.z);
    quat_alg_elem_init(&order.t);
    params.primality_test_iterations = 1;

    ibz_set_from_str(&order.z.coord[0], "0", 10);
    ibz_set_from_str(&order.z.coord[1], "-1", 10);
    ibz_set(&order.z.coord[2], 0);
    ibz_set_from_str(&order.z.coord[3], "-1", 10);
    ibz_set_from_str(&order.z.denom, "21267647932558653966460912964485513216", 10);
    ibz_set(&order.t.coord[2], 1);
    ibz_set(&order.t.denom, 1);
    order.q = 5;
    ibz_set_from_str(&(order.order.basis[0][0]), "21267647932558653966460912964485513216", 10);
    ibz_set_from_str(&(order.order.basis[1][1]), "-1", 10);
    ibz_set_from_str(&(order.order.basis[2][2]), "10633823966279326983230456482242756608", 10);
    ibz_set(&(order.order.basis[3][1]), -1);
    ibz_set_from_str(&(order.order.basis[0][2]), "10633823966279326983230456482242756608", 10);
    ibz_set_from_str(&(order.order.basis[0][3]), "0", 10);
    ibz_set_from_str(&(order.order.basis[1][2]), "0", 10);
    ibz_set_from_str(
        &(order.order.basis[1][3]), "-226156424291633194186662080095093570025917938800079226639565593765455331328", 10);

    ibz_set_from_str(&(order.order.denom), "21267647932558653966460912964485513216", 10);

    params.algebra = &Bpoo;
    params.order = &order;
    params.algebra = &Bpoo;
    params.order = &order;
    ibz_init(&n_gamma);
ibz_init(&ga);
    
    
    ibz_set_from_str(  &n_gamma, "3573960250147475915502838919373697042179781142194761891096159111893170748566397047", 10);
    //ibz_set_from_str(   &n_gamma,"284466028237988494763970178932217686760481835540054313610030844452355852304850148319078969945563537596721366979647", 10);
        //"3573960250147475915502838919373697042179781142194761891096159111893170748566397047",
        
        /*ibz_rand_interval_bits(&n_gamma,intsize+50);
   //rand_exact_m_bits(&n_gamma,state,intsize+50);
       ibz_mul(&ga,&n_gamma,&ibz_const_two);
       ibz_add(&n_gamma,&ga,&ibz_const_one);*/
      /* if(intsize==350){ 
         ibz_init(&n_gamma);
           ibz_init(&ga);
    ibz_rand_interval_bits(&n_gamma,intsize+50);
   //rand_exact_m_bits(&n_gamma,state,intsize+50);
       ibz_mul(&ga,&n_gamma,&ibz_const_two);
       ibz_add(&n_gamma,&ga,&ibz_const_one); }
*/
       
       
    QUAT_BENCH_CODE_1(iterations);

    QUAT_BENCH_CODE_START()
    quat_ct_represent_integerwithq(&gamma, &n_gamma, 0, &params);
    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = 0;
    QUAT_BENCH_CODE_2("ct represent_integer with varying q", res);

    return 1;
}

long
quat_bench_repintoldwithq(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;
    // Initilialize for benchmarking macros
    int64_t cycles, cycles1, cycles2;
    int i;
    
     
    ibz_t n_gamma,ga;//, ga;
    int LIST_SIZE = iterations;
    int64_t *cycles_list;
    int runs = 0;
    cycles_list = (int64_t *)malloc(iterations * sizeof(int64_t));

    quat_represent_integer_params_t params;

    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    quat_alg_elem_init(&gamma);
    quat_alg_init_set(&Bpoo, prime);
    quat_p_extremal_maximal_order_t order;
    quat_lattice_init(&order.order);
    quat_alg_elem_init(&order.z);
    quat_alg_elem_init(&order.t);
    params.primality_test_iterations = 1;
if(intsize==260){
    ibz_set_from_str(&order.z.coord[0], "0", 10);
    ibz_set_from_str(&order.z.coord[1], "-1", 10);
    ibz_set(&order.z.coord[2], 0);
    ibz_set_from_str(&order.z.coord[3], "-1", 10);
    ibz_set_from_str(&order.z.denom, "21267647932558653966460912964485513216", 10);
    ibz_set(&order.t.coord[2], 1);
    ibz_set(&order.t.denom, 1);
    order.q = 5;
    ibz_set_from_str(&(order.order.basis[0][0]), "21267647932558653966460912964485513216", 10);
    ibz_set_from_str(&(order.order.basis[1][1]), "-1", 10);
    ibz_set_from_str(&(order.order.basis[2][2]), "10633823966279326983230456482242756608", 10);
    ibz_set(&(order.order.basis[3][1]), -1);
    ibz_set_from_str(&(order.order.basis[0][2]), "10633823966279326983230456482242756608", 10);
    ibz_set_from_str(&(order.order.basis[0][3]), "0", 10);
    ibz_set_from_str(&(order.order.basis[1][2]), "0", 10);
    ibz_set_from_str(
        &(order.order.basis[1][3]), "-226156424291633194186662080095093570025917938800079226639565593765455331328", 10);

    ibz_set_from_str(&(order.order.denom), "21267647932558653966460912964485513216", 10);

    params.algebra = &Bpoo;
    params.order = &order;
    
   
    ibz_set_from_str(&n_gamma, "284466028237988494763970178932217686760481835540054313610030844452355852304850148319078969945563537596721366979647", 10);
}else{
    ibz_set_from_str(&order.z.coord[1], "214764116738303679745780048598183569015", 10);
    ibz_set(&order.z.coord[2], 0);
    ibz_set_from_str(&order.z.coord[3], "1", 10);
    ibz_set_from_str(&order.z.denom, "73403150722045993989123427738005336972", 10);
    ibz_set(&order.t.coord[2], 1);
    ibz_set(&order.t.denom, 1);
    order.q = 13;
    ibz_set_from_str(&(order.order.basis[0][0]), "73403150722045993989123427738005336972", 10);
    ibz_set_from_str(&(order.order.basis[1][1]),
                     "2694011267961700664357934052637599020390646823337886018360381743577635064392",
                     10);
    ibz_set_from_str(&(order.order.basis[2][2]), "36701575361022996994561713869002668486", 10);
    ibz_set(&(order.order.basis[3][3]), 1);
    ibz_set_from_str(&(order.order.basis[0][2]), "36701575361022996994561713869002668486", 10);
    ibz_set_from_str(&(order.order.basis[0][3]), "0", 10);
    ibz_set_from_str(&(order.order.basis[1][2]), "0", 10);
    ibz_set_from_str(&(order.order.basis[1][3]), "214764116738303679745780048598183569015", 10);
    ibz_set_from_str(&(order.order.denom), "73403150722045993989123427738005336972", 10);

    params.algebra = &Bpoo;
    params.order = &order;
      ibz_init(&n_gamma);
  ibz_init(&ga);
    // quat_test_set_params_standard(&params, &order, &Bpoo);
  
    ibz_rand_interval_bits(&n_gamma,intsize+50);
   //rand_exact_m_bits(&n_gamma,state,intsize+50);
       ibz_mul(&ga,&n_gamma,&ibz_const_two);
       ibz_add(&n_gamma,&ga,&ibz_const_one); 
    //ibz_set_from_str(&n_gamma, "284466028237988494763970178932217686760481835540054313610030844452355852304850148319078969945563537596721366979647", 10);

}

   QUAT_BENCH_CODE_1(iterations);    
    QUAT_BENCH_CODE_START()
     
    quat_old_represent_integer(&gamma, &n_gamma, 0, &params);

    QUAT_BENCH_CODE_STOP()
    runs = i;
    res = 0;
    QUAT_BENCH_CODE_2("old represent_integer with q=5 :", res);
    return 1;
}

long
quat_bench_rep(int iterations, int intsize, const ibz_t *prime)
{
    int res = 0;

  res = res | quat_bench_repintold(iterations,intsize,prime);
 // 
     res = res | quat_bench_repint(iterations,intsize,prime); 
    // res = res | quat_bench_repintoldwithq(iterations, intsize, prime);
 //   res = res | quat_bench_repintvarq(iterations, intsize, prime);
 
  
    return (res);
}


