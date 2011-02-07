#ifndef __OMP_H__
#define __OMP_H__

typedef int omp_int_t;
static inline omp_int_t omp_get_thread_num() { return 0;}
static inline omp_int_t omp_get_max_threads() { return 1;}

#endif /* omp.h */
