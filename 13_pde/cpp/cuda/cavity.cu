#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

float lb = 0;
float ub = 2;
const int nx = 11;
const int ny = 11;
int nt = 10;
int nit = 5;
int c = 1;
float dx = ub / float(nx - 1);
float dy = ub / float(ny - 1);

float rho = 1;
float nu = 0.1;
float dt = 0.001;

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

void indexfill(float *x, int dimx, int dimy) {
  for (int k = 0; k < dimy * dimx; k++)
    x[k] = pow(k, 2);
}
void meshgrid(float *x, float *y, float *X, float *Y) {
  int dimx = ny;
  int dimy = ny;

  for (int k = 0; k < dimy; k++)
    for (int i = 0; i < dimx; i++) {
      X[k * dimx + i] = x[k];
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

__global__ void pressure_poisson(float *p, float *pn, float dx, float dy,
                                 float *b, int nx, int ny, int nit) {
  /*
def pressure_poisson(p, dx, dy, b):

*/
  int xlim = ny;
  int ylim = nx;
  int idx = blockIdx.x;
  int idy = threadIdx.x;
  // if(idx >= xlim - 2) return;
  // if (idy >= ylim - 2) return;

//   for (int q = 0; q < nit; q++) 
  {
    printf("%d %d; ", idx, idy);

    if (idx < xlim - 2 && idy < ylim - 2) {
        copy(p, pn, ny, nx);
        __syncthreads();
      p[(idx + 1) * nx + idy + 1] =
          //   pn[(idx + 1) * nx + idy+2] + b[(idx + 1) * nx + idy + 1];
          (
            //   (pn[(idx + 1) * nx + idy + 2] + 
            // (
        //   pn[(idx + 1) * nx + idy]) +
            (pn[(idx + 2) * nx + idy + 1] + 
            pn[(idx)*nx + idy + 1]) -
               b[(idx + 1) * nx + idy + 1]);
        
    //   (((pn[(idx + 1) * nx + idy + 2] + pn[(idx + 1) * nx + idy]) * pow(dy, 2) +
    //     (pn[(idx + 2) * nx + idy + 1] + pn[(idx)*nx + idy + 1]) * pow(dx, 2)) /
    //        (2 * (pow(dx, 2) + pow(dy, 2))) -
    //    pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
    //        b[(idx + 1) * nx + idy + 1]);
    }

    // p[0 * ylim + idy] = p[1 * ylim + idy]; // dp/dy = 0 at y = 0
    // __syncthreads();
    // p[(xlim - 1) * ylim + idy] = 0;        // p = 0 at y = 2
    // __syncthreads();

    // p[idx * ylim + ylim - 1] = p[idx * ylim + ylim - 2]; // dp/dx = 0 at x =
    // 2
    // __syncthreads();

    // p[idx * ylim + 0] = p[idx * ylim + 1];               // dp/dx = 0 at x =
    // 0
    // __syncthreads();
  }
}

// __global__ void cavity_flow(int nt, float *u, float *v,
//                  float dt, float dx, float dy, float *p, float rho, float nu)
// {
//     /*
//     def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
// */
//     int xlim = ny;
//     int ylim = nx;
//     float un[ny * nx]; // copy !!!
//     float vn[ny * nx];
//     float b[ny * nx];

//     for (int n = 0; n < nt; n++)
//     {
//         copy(u, un, ny, nx);
//         copy(v, vn, ny, nx);
//         build_up_b(b, rho, dt, u, v, dx, dy);
//         pressure_poisson(p, dx, dy, b);

//         for (int idx = 0; idx < xlim - 2; idx++)
//             for (int idy = 0; idy < ylim - 2; idy++)
//             {
//                 u[(idx + 1) * nx + idy + 1] = (un[(idx + 1) * nx + idy + 1] -
//                                                un[(idx + 1) * nx + idy + 1] *
//                                                dt / dx *
//                                                    (un[(idx + 1) * nx + idy +
//                                                    1] - un[(idx + 1) * nx +
//                                                    idy]) -
//                                                vn[(idx + 1) * nx + idy + 1] *
//                                                dt / dy *
//                                                    (un[(idx + 1) * nx + idy +
//                                                    1] - un[(idx)*nx + idy +
//                                                    1]) -
//                                                dt / (2 * rho * dx) * (p[(idx
//                                                + 1) * nx + idy + 2] - p[(idx
//                                                + 1) * nx + idy]) + nu * (dt /
//                                                pow(dx, 2) *
//                                                          (un[(idx + 1) * nx +
//                                                          idy + 2] - 2 *
//                                                          un[(idx + 1) * nx +
//                                                          idy + 1] + un[(idx +
//                                                          1) * nx + idy]) +
//                                                      dt / pow(dy, 2) *
//                                                          (un[(idx + 2) * nx +
//                                                          idy + 1] - 2 *
//                                                          un[(idx + 1) * nx +
//                                                          idy + 1] +
//                                                          un[(idx)*nx + idy +
//                                                          1])));

//                 v[(idx + 1) * nx + idy + 1] = (vn[(idx + 1) * nx + idy + 1] -
//                                                un[(idx + 1) * nx + idy + 1] *
//                                                dt / dx *
//                                                    (vn[(idx + 1) * nx + idy +
//                                                    1] - vn[(idx + 1) * nx +
//                                                    idy]) -
//                                                vn[(idx + 1) * nx + idy + 1] *
//                                                dt / dy *
//                                                    (vn[(idx + 1) * nx + idy +
//                                                    1] - vn[(idx)*nx + idy +
//                                                    1]) -
//                                                dt / (2 * rho * dy) * (p[(idx
//                                                + 2) * nx + idy + 1] -
//                                                p[(idx)*nx + idy + 1]) + nu *
//                                                (dt / pow(dx, 2) *
//                                                          (vn[(idx + 1) * nx +
//                                                          idy + 2] - 2 *
//                                                          vn[(idx + 1) * nx +
//                                                          idy + 1] + vn[(idx +
//                                                          1) * nx + idy]) +
//                                                      dt / pow(dy, 2) *
//                                                          (vn[(idx + 2) * nx +
//                                                          idy + 1] - 2 *
//                                                          vn[(idx + 1) * nx +
//                                                          idy + 1] +
//                                                          vn[(idx)*nx + idy +
//                                                          1])));
//             }
//         for (int k = 0; k < ylim; k++)
//         {
//             u[0 * ylim + k] = 0;
//             u[(xlim - 1) * ylim + k] = 10;
//             v[0 * ylim + k] = 0;
//             v[(xlim - 1) * ylim + k] = 0;
//         }
//         for (int idx = 0; idx < xlim; idx++)
//         {
//             u[idx * ylim + 0] = 0;
//             v[idx * ylim + 0] = 0;
//             u[idx * ylim + ylim - 1] = 0;
//             v[idx * ylim + ylim - 1] = 0;
//         }
//     }
// }

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
__global__ void npp(float *a, int N) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i > N)
    return;
  a[i] = i;
}
int main() {
  float x[nx];
  float y[ny];
  float X[ny * nx];
  float Y[ny * nx];
  linspace(x, lb, ub, nx);
  // npprint(x, 1, nx);
  linspace(y, lb, ub, ny);
  meshgrid(x, y, X, Y);
  // npprint(X, ny, nx);
  // npprint(Y);

  // float u[ny * nx]; fill(u, 1, nx, ny);
  // float v[ny * nx]; fill(v, 2, nx, ny);
  // float b[ny * nx]; fill(b, 0, nx, ny);
  // float p[ny * nx]; fill(p, 0, nx, ny);

  float *u;
  float *v;
  float *p, *pn;
  float *b;
  cudaMallocManaged(&b, nx * ny * sizeof(float));
  cudaMallocManaged(&p, nx * ny * sizeof(float));
  cudaMallocManaged(&pn, nx * ny * sizeof(float));
  cudaMallocManaged(&u, nx * ny * sizeof(float));
  cudaMallocManaged(&v, nx * ny * sizeof(float));
  fill(u, 1, nx, ny);
  fill(v, 2, nx, ny);
  indexfill(v, nx, ny);
  fill(b, 0, nx, ny);
  fill(p, 0, nx, ny);
  fill(pn, 0, nx, ny);

  auto start = std::chrono::high_resolution_clock::now();
  // cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);
  auto finish = std::chrono::high_resolution_clock::now();
  // npprint(u);
  // npprint(v);
  // npprint(p);
  // npprint(b);
  std::chrono::duration<double> elapsed = finish - start;
  printf("Elapsed time: %3.3f s\n", elapsed.count());

  // float *pg;
  // float *bg;
  // cudaMallocManaged(&pg, nx*ny*sizeof(float));
  // cudaMallocManaged(&bg, nx*ny*sizeof(float));
  // npp<<<nx, ny>>>(pg, nx*ny);
  build_up_b<<<nx, ny>>>(b, rho, dt, u, v, dx, dy, nx, ny);
  // fill(b, 1, nx, ny);
  pressure_poisson<<<nx, ny>>>(p, pn, dx, dy, b, nx, ny, nit);
  cudaDeviceSynchronize();
  // cudaMemcpy(bg, b, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);

  //   npprint(b); // confirmed
  npprint(p);
  // cout << &bg;
  cudaFree(p);
  cudaFree(b);
  cudaFree(u);
  cudaFree(v);
}