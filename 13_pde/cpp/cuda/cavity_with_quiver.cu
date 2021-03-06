/*
COMPILE WITH:
 nvcc cavity.cu -I/usr/include/python2.7 -lpython2.7&& ./a.out

 PLOTTING FUNCTION FROM:

 https://github.com/lava/matplotlib-cpp
 
*/

#include "../matplotlibcpp.h"
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <math.h>
#include <string>

namespace plt = matplotlibcpp;
using namespace std;

float lb = 0;
float ub = 2;
const int nx = 41;
const int ny = 41;
int nt = 500;
int nit = 50;
int c = 1;
float dx = ub / float(nx - 1);
float dy = ub / float(ny - 1);

float rho = 1;
float nu = 0.1;
float dt = 0.001;

void npprint(float *u, int dimx = ny, int dimy = nx, string msg = "OUT: ") {
  printf("%s\n", msg.c_str());
  printf("x-------------------------------x\n");
  printf("[\n");
  for (int i = 0; i < dimx; i++) {
    printf("[");

    for (int k = 0; k < dimy; k++)
      printf("%3.4f, ", u[i * dimy + k]);
    printf("],\n");
  }
  printf("]\n");

  printf("x-------------------------------x\n");
}

void npprint(vector<vector<float>> u, string msg = "OUT: ") {
  printf("%s\n", msg.c_str());
  printf("x-------------------------------x\n");
  printf("[\n");
  for (int i = 0; i < u.size(); i++) {
    printf("[");
    for (int k = 0; k < u[0].size(); k++)
      cout << u[i][k] << ", ";
    printf("],\n");
  }
  printf("]\n");
  printf("x-------------------------------x\n");
}

void linspace(float *x, int lb, int ub, int num) {
  for (int k = 0; k < num; k++) {
    x[k] = k * (ub - lb) / float(num - 1) + lb;
  }
}

__device__ void copy(float *lhs, float *rhs, int dimx, int dimy) {
  for (int k = 0; k < dimy * dimx; k++)
    rhs[k] = lhs[k];
}

void fill(float *x, float fillnum, int dimx, int dimy) {
  for (int k = 0; k < dimy * dimx; k++)
    x[k] = fillnum;
}

void meshgrid(float *x, float *y, float *X, float *Y) {
  int dimx = ny;
  int dimy = ny;

  for (int k = 0; k < dimy; k++)
    for (int i = 0; i < dimx; i++) {
      X[k * dimx + i] = x[i];
      Y[k * dimx + i] = y[k];
    }
}

__global__ void build_up_b(float *b, float rho, float dt, float *u, float *v,
                           float dx, float dy, int nx, int ny) {
  /*

def build_up_b(b, rho, dt, u, v, dx, dy):
  return b
*/
  int xlim = ny;
  int ylim = nx;

  int idx = blockIdx.x;
  int idy = threadIdx.x;
  assert(nx = blockDim.x);
  if (idx >= xlim - 2)
    return;
  if (idy >= ylim - 2)
    return;
  b[(idx + 1) * nx + idy + 1] =
      (rho *
       (1 / dt *
            ((u[(idx + 1) * nx + idy + 2] - u[(idx + 1) * nx + idy]) /
                 (2 * dx) +
             (v[(idx + 2) * nx + idy + 1] - v[(idx)*nx + idy + 1]) / (2 * dy)) -
        pow(((u[(idx + 1) * nx + idy + 2] - u[(idx + 1) * nx + idy]) /
             (2 * dx)),
            2) -
        2 * ((u[(idx + 2) * nx + idy + 1] - u[(idx)*nx + idy + 1]) / (2 * dy) *
             (v[(idx + 1) * nx + idy + 2] - v[(idx + 1) * nx + idy]) /
             (2 * dx)) -
        pow((v[(idx + 2) * nx + idy + 1] - v[(idx)*nx + idy + 1]) / (2 * dy),
            2)));
}

__global__ void pmargin(float *p, int nx, int ny) {
  int xlim = ny;
  int ylim = nx;
  int idx = blockIdx.x;
  int idy = threadIdx.x;

  p[0 * ylim + idy] = p[1 * ylim + idy]; // dp/dy = 0 at y = 0
  p[(xlim - 1) * ylim + idy] = 0;        // p = 0 at y = 2
  __syncthreads();

  p[idx * ylim + ylim - 1] = p[idx * ylim + ylim - 2]; // dp/dx = 0 at x =2
  p[idx * ylim + 0] = p[idx * ylim + 1];               // dp/dx = 0 at x =0
  __syncthreads();
}

__global__ void pupdate(float *p, float *pn, float dx, float dy, float *b,
                        int nx, int ny, int nit) {
  int xlim = ny;
  int ylim = nx;
  int idx = blockIdx.x;
  int idy = threadIdx.x;
  if (idx < xlim - 2 && idy < ylim - 2) {
    copy(p, pn, ny, nx);

    __syncthreads();
    p[(idx + 1) * nx + idy + 1] =

        (((pn[(idx + 1) * nx + idy + 2] + pn[(idx + 1) * nx + idy]) *
              pow(dy, 2) +
          (pn[(idx + 2) * nx + idy + 1] + pn[(idx)*nx + idy + 1]) *
              pow(dx, 2)) /
             (2 * (pow(dx, 2) + pow(dy, 2))) -
         pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
             b[(idx + 1) * nx + idy + 1]);
  }
  __syncthreads();
}

void pressure_poisson(float *p, float *pn, float dx, float dy, float *b, int nx,
                      int ny, int nit) {
  /*
def pressure_poisson(p, dx, dy, b):

*/
  for (int q = 0; q < nit; q++) {
    pupdate<<<ny, nx>>>(p, pn, dx, dy, b, nx, ny, nit);
    cudaDeviceSynchronize();
    pmargin<<<ny, nx>>>(p, nx, ny);
    cudaDeviceSynchronize();
  }
}
__global__ void cupdate(int nt, float *u, float *v, float *un, float *vn,
                        float dt, float dx, float dy, float *p, float rho,
                        float nu) {
  int idx = blockIdx.x;
  int idy = threadIdx.x;
  copy(u, un, ny, nx);
  copy(v, vn, ny, nx);
  __syncthreads();
  u[(idx + 1) * nx + idy + 1] =
      (un[(idx + 1) * nx + idy + 1] -
       un[(idx + 1) * nx + idy + 1] * dt / dx *
           (un[(idx + 1) * nx + idy + 1] - un[(idx + 1) * nx + idy]) -
       vn[(idx + 1) * nx + idy + 1] * dt / dy *
           (un[(idx + 1) * nx + idy + 1] - un[(idx)*nx + idy + 1]) -
       dt / (2 * rho * dx) *
           (p[(idx + 1) * nx + idy + 2] - p[(idx + 1) * nx + idy]) +
       nu * (dt / pow(dx, 2) *
                 (un[(idx + 1) * nx + idy + 2] -
                  2 * un[(idx + 1) * nx + idy + 1] + un[(idx + 1) * nx + idy]) +
             dt / pow(dy, 2) *
                 (un[(idx + 2) * nx + idy + 1] -
                  2 * un[(idx + 1) * nx + idy + 1] + un[(idx)*nx + idy + 1])));
  __syncthreads();

  v[(idx + 1) * nx + idy + 1] =
      (vn[(idx + 1) * nx + idy + 1] -
       un[(idx + 1) * nx + idy + 1] * dt / dx *
           (vn[(idx + 1) * nx + idy + 1] - vn[(idx + 1) * nx + idy]) -
       vn[(idx + 1) * nx + idy + 1] * dt / dy *
           (vn[(idx + 1) * nx + idy + 1] - vn[(idx)*nx + idy + 1]) -
       dt / (2 * rho * dy) *
           (p[(idx + 2) * nx + idy + 1] - p[(idx)*nx + idy + 1]) +
       nu * (dt / pow(dx, 2) *
                 (vn[(idx + 1) * nx + idy + 2] -
                  2 * vn[(idx + 1) * nx + idy + 1] + vn[(idx + 1) * nx + idy]) +
             dt / pow(dy, 2) *
                 (vn[(idx + 2) * nx + idy + 1] -
                  2 * vn[(idx + 1) * nx + idy + 1] + vn[(idx)*nx + idy + 1])));
  __syncthreads();
}

__global__ void cmargin(float *u, float *v, int nx, int ny) {
  int xlim = ny;
  int ylim = nx;
  int idx = blockIdx.x;
  int idy = threadIdx.x;
  u[0 * ylim + idy] = 0;
  u[(xlim - 1) * ylim + idy] = 1;
  v[0 * ylim + idy] = 0;
  v[(xlim - 1) * ylim + idy] = 0;
  __syncthreads();
  u[idx * ylim + 0] = 0;
  v[idx * ylim + 0] = 0;
  u[idx * ylim + ylim - 1] = 0;
  v[idx * ylim + ylim - 1] = 0;
  __syncthreads();
}
void cavity_flow(int nt, float *u, float *v, float *un, float *vn, float dt,
                 float dx, float dy, float *p, float *pn, float rho, float nu) {
  /*
  def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
*/
  float *b;
  cudaMallocManaged(&b, nx * ny * sizeof(float));
  fill(b, 0, nx, ny);
  for (int n = 0; n < nt; n++) {

    build_up_b<<<nx, ny>>>(b, rho, dt, u, v, dx, dy, nx, ny);
    cudaDeviceSynchronize();

    pressure_poisson(p, pn, dx, dy, b, nx, ny, nit);
    cudaDeviceSynchronize();

    cupdate<<<nx, ny>>>(nt, u, v, un, vn, dt, dx, dy, p, rho, nu);
    cudaDeviceSynchronize();
    cmargin<<<nx, ny>>>(u, v, nx, ny);
    cudaDeviceSynchronize();
  }
  cudaFree(b);
}

void reshapeToVector(float *x, vector<float> &vx) {
  int ny = vx.size();
  for (int idx = 0; idx < ny; idx++)
    vx[idx] = x[idx];
}

int main() {
  float x[nx];
  float y[ny];
  float X[ny * nx];
  float Y[ny * nx];
  linspace(x, lb, ub, nx);
  linspace(y, lb, ub, ny);
  meshgrid(x, y, X, Y);

  float *u, *un;
  float *v, *vn;
  float *p, *pn;
  cudaMallocManaged(&p, nx * ny * sizeof(float));
  cudaMallocManaged(&pn, nx * ny * sizeof(float));
  cudaMallocManaged(&u, nx * ny * sizeof(float));
  cudaMallocManaged(&v, nx * ny * sizeof(float));
  cudaMallocManaged(&un, nx * ny * sizeof(float));
  cudaMallocManaged(&vn, nx * ny * sizeof(float));
  fill(u, 0, nx, ny);
  fill(v, 0, nx, ny);
  fill(p, 0, nx, ny);
  fill(pn, 0, nx, ny);
  fill(vn, 0, nx, ny);
  fill(un, 0, nx, ny);

  auto start = std::chrono::high_resolution_clock::now();
  cavity_flow(nt, u, v, un, vn, dt, dx, dy, p, pn, rho, nu);
  cudaDeviceSynchronize();

  auto finish = std::chrono::high_resolution_clock::now();
  npprint(u, ny, nx, "U");
  npprint(v, ny, nx, "V");
  npprint(p, ny, nx, "P");
  std::chrono::duration<double> elapsed = finish - start;
  printf("GPU Elapsed time: %3.3f s\n", elapsed.count());

  vector<float> vec_X(nx * ny);
  vector<float> vec_Y(ny * nx);
  vector<float> vec_u(ny * nx);
  vector<float> vec_v(ny * nx);
  vector<float> vec_p(ny * nx);

  reshapeToVector(u, vec_u);
  reshapeToVector(v, vec_v);
  reshapeToVector(p, vec_p);
  reshapeToVector(X, vec_X);
  reshapeToVector(Y, vec_Y);

  plt::figure();
  plt::contour(vec_X, vec_Y, vec_u);
  plt::quiver(vec_X, vec_Y, vec_u, vec_v);
  plt::show();

  cudaFree(p);
  cudaFree(pn);
  cudaFree(vn);
  cudaFree(un);
  cudaFree(u);
  cudaFree(v);
}