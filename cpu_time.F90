subroutine my_cpu_time(t1)
  real :: t1
#ifdef _OPENMP
    real :: omp_get_wtime
#endif

#ifdef _OPENMP
  t1=omp_get_wtime()
#else
  call cpu_time(t1)
#endif

end  subroutine my_cpu_time
