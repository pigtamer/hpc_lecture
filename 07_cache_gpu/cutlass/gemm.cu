#include <iostream>
#include <typeinfo>
#include <random>
#include <stdint.h>
#include <cublas_v2.h>
#define DEBUG
#include <gemm/dispatch.h>
#include <gemm/epilogue_function.h>
#include "util/timer.h"

using namespace cutlass;


void sync_device(float *M, float * dM, int m, int n)
{
    size_t bytes = m * n * sizeof(float);
    CUDA_PERROR_EXIT(cudaMemcpy(dM, &M[0], bytes, cudaMemcpyHostToDevice));
}


void sync_host(float *M, float * dM, int m, int n)
{
    size_t bytes = m * n * sizeof(float);
    CUDA_PERROR_EXIT(cudaMemcpy(&M[0], dM, bytes, cudaMemcpyDeviceToHost));
}

void fill_ramp(
  float *M,
  float rs,
  float cs,
  int m, int n)
{
  for (int x = 0; x < n; x++)
  {
      for (int y = 0; y < m; y++)
      {
          (&M[0])[y + (x * m)] = float((y * rs) + (x * cs));
      }
  }
}

void random(float *M,
  int m, int n) {
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      (&M[0])[i + j * m] = float(drand48());
    }
  }
}

int main(int argc, const char **argv) {
  int m = 10240;
  int k = 4096;
  int n = 4096;
  float alpha = 1.0;
  float beta = 0.0;
  static const matrix_transform_t::kind_t TransformA = matrix_transform_t::NonTranspose;
  static const matrix_transform_t::kind_t TransformB = matrix_transform_t::NonTranspose;
  typedef float value_t;
  typedef float accum_t;
  int g_timing_iterations = 10;
  cudaStream_t stream = 0;
  // matrix<value_t> A(m, k);
  // matrix<value_t> B(k, n);
  // matrix<accum_t> C(m, n);
  // matrix<accum_t> C2(m, n);
  float A[m*k], B[k*n], C[m*n], C2[m*n];
  float *d_A, *d_B, *d_C, *d_C2;

  CUDA_PERROR_EXIT(cudaMalloc((void **) &d_A, m*k * sizeof(float)));

  CUDA_PERROR_EXIT(cudaMalloc((void **) &d_B, k * n * sizeof(float)));
  CUDA_PERROR_EXIT(cudaMalloc((void **) &d_C, m * n * sizeof(float)));
  CUDA_PERROR_EXIT(cudaMalloc((void **) &d_C2, m * n * sizeof(float)));

  random(A, m, k);
  random(B, k ,n);

  fill_ramp(C, 0,0, m,n);
  fill_ramp(C2, 0,0,m,n);

  sync_device(A, d_A, m,k);
  sync_device(B, d_B, k,n);
  sync_device(C, d_C, m,n);
  sync_device(C2, d_C2, m,n);

  cublasHandle_t g_cublas_handle;
  cublasCreate(&g_cublas_handle);
  gpu_timer timer;
  for (int i = 0; i < g_timing_iterations+2; i++) {
    if (i == 2) timer.start();
    CUDA_PERROR(cublasSgemm(
                            g_cublas_handle,
                            (cublasOperation_t) TransformA,
                            (cublasOperation_t) TransformB,
                            m,
                            n,
                            k,
                            &alpha,
                            &d_A[0],
                            m,
                            &d_B[0],
                            k,
                            &beta,
                            &d_C[0],
                            m));
  }
  timer.stop();
  int64_t num_flops = (2 * int64_t(m) * int64_t(n) * int64_t(k)) + (2 * int64_t(m) * int64_t(n));
  double tcublas = timer.elapsed_millis() / g_timing_iterations;
  double cublas_flops = double(num_flops) / tcublas / 1.0e6;
  typedef gemm::blas_scaled_epilogue<float, float, float> epilogue_op_t;
  epilogue_op_t epilogue(alpha, beta);
  for (int i = 0; i < g_timing_iterations+2; i++) {
    if (i == 2) timer.start();
    gemm::dispatch<epilogue_op_t>(
        m,
        n,
        k,
        alpha,
        beta,
        d_A,
        d_B,
        d_C,
        stream,
        false);
  }
  timer.stop();
  double tcutlass = timer.elapsed_millis() / g_timing_iterations;
  double cutlass_flops = double(num_flops) / tcutlass / 1.0e6;
  printf("CUBLAS: %.2f Gflops, CUTLASS: %.2f Gflops\n", cublas_flops, cutlass_flops);
  sync_host(C, d_C, m,n);
  sync_host(C2, d_C2, m,n);
  // C2.sync_host();
  double err = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      err += fabs(C[i+j*n] - C2[i+j*n]);
    }
  }
  printf("error: %lf\n", err/n/m);
  cublasDestroy(g_cublas_handle);
}
