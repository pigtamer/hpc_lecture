#include <cstdio>
#include <omp.h>

int main() {
  int n = 16;
  double dx = 1. / n;
  double pi = 0;

# pragma omp parallel for 
  for (int i=0; i<n; i++) {
    double x = (i + 0.5) * dx;

#pragma omp atomic update
    pi += 4.0 / (1.0 + x * x) * dx;

    int tid = omp_get_thread_num();
    printf("Thread %d\n", tid);
  }
  printf("%17.15f\n",pi);
}
