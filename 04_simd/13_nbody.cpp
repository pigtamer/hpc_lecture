#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main()
{
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], idx[N], mm[N];
  for (int i = 0; i < N; i++)
  {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    mm[i] = m[i];
    idx[i] = i;
    fx[i] = fy[i] = 0;
  }

  // Original implementation
  // for (int i = 0; i < N; i++)
  // {
  //   for (int j = 0; j < N; j++)
  //   {
  //     if (i != j)
  //     {
  //       float rx = x[i] - x[j];
  //       float ry = y[i] - y[j];
  //       float r = std::sqrt(rx * rx + ry * ry);
  //       fx[i] -= rx * m[j] / (r * r * r);
  //       fy[i] -= ry * m[j] / (r * r * r);
  //     }
  //   }
  //   printf("%d %g %g\n", i, fx[i], fy[i]);
  // }

  // Assignment:
  __m256 fxvec = _mm256_load_ps(fx);
  __m256 fyvec = _mm256_load_ps(fy);
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);

  for (int i = 0; i < N; i++)
  {
    __m256 mvec = _mm256_load_ps(m);
    __m256 idxvec = _mm256_load_ps(idx);
    __m256 j = _mm256_set1_ps(i);
    __m256 mask = _mm256_cmp_ps(idxvec, j, _CMP_EQ_OQ);

    mvec = _mm256_blendv_ps(mvec, _mm256_setzero_ps(), mask);
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);

    __m256 rx = _mm256_sub_ps(xi, xvec);
    __m256 ry = _mm256_sub_ps(yi, yvec);

    __m256 r = _mm256_rsqrt_ps(
        _mm256_add_ps(
            _mm256_mul_ps(rx, rx),
            _mm256_mul_ps(ry, ry)));
    r = _mm256_blendv_ps(r, _mm256_setzero_ps(), mask);

    __m256 fxi, fyi;
    // sum up j-indexed terms for each i
    fxi = _mm256_permute2f128_ps(rx * mvec * (r * r * r), rx * mvec * (r * r * r), 1);
    fxi = _mm256_hadd_ps(fxi, rx * mvec * (r * r * r));
    fxi = _mm256_hadd_ps(fxi, fxi);
    fxi = _mm256_hadd_ps(fxi, fxi);
    fxi = _mm256_blendv_ps(_mm256_setzero_ps(), fxi, mask);

    fyi = _mm256_permute2f128_ps(ry * mvec * (r * r * r), ry * mvec * (r * r * r), 1);
    fyi = _mm256_hadd_ps(fyi, ry * mvec * (r * r * r));
    fyi = _mm256_hadd_ps(fyi, fyi);
    fyi = _mm256_hadd_ps(fyi, fyi);
    fyi = _mm256_blendv_ps(_mm256_setzero_ps(), fyi, mask);

    fxvec = _mm256_sub_ps(fxvec, fxi);
    fyvec = _mm256_sub_ps(fyvec, fyi);

    _mm256_store_ps(fx, fxvec);
    _mm256_store_ps(fy, fyvec);

    //for debug
    // _mm256_store_ps(mm, mvec);
    // for (int im = 0; im < N; im++)
    //   printf("%3.2f\t", fx[im]);
    
    printf("%d %g %g\n", i, fx[i], fy[i]);
  }
}