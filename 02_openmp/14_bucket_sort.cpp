#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

template <class T>
bool check(std::vector<T> &vec, size_t n)
{
  bool res = true;
  for (size_t k = 0; k < n - 1; k++)
  {
    res = res & (vec[k] <= vec[k + 1]);
  }
  printf("%d\n", res);
  return res;
}

int main()
{
  int n = 100;
  int range = 128;
  omp_set_num_threads(range); // significant
  std::vector<int> key(n), a(range), b(range);
  std::vector<int> bucket(range);
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < n; i++)
    {
      key[i] = rand() % range;
      printf("%d ", key[i]);
    }
    // printf("\n");

#pragma omp for
    for (int i = 0; i < n; i++)
    {
#pragma omp atomic update
      bucket[key[i]]++;
    }

// build a global index for the buckets' starting locations
#pragma omp single
    {
      a = bucket;
      b = bucket;
    }

    for (int j = 1; j < range; j <<= 1)
    {
#pragma omp for
      for (int i = 0; i < range; i++)
        b[i] = a[i];
#pragma omp for
      for (int i = j; i < range; i++)
        a[i] += b[i - j];
    }

#pragma omp for
    for (int i = 0; i < range; i++)
    {
      int j = 0; // use a local index inside the bucket instead
      for (; bucket[i] > 0; bucket[i]--)
      {
        key[j + a[omp_get_thread_num() - 1]] = i;
        j++;
      }
    }
  }

  printf("\n");
  for (int i = 0; i < n; i++)
  {
    printf("%d ", key[i]);
  }
  printf("\n");

  check(key, n);
}