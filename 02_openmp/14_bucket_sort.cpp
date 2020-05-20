#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

template <class T>
bool check(std::vector<T> &vec, size_t n){
  bool res = true;
  for(size_t k=0; k< n-1; k++){
    res = res&(vec[k] <= vec[k+1]);
  }
  printf("%d\n", res);
  return res;
}

int main()
{
  int n = 100;
  int range = 128;
  std::vector<int> key(n);
#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    key[i] = rand() % range;
    printf("%d ", key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range, 0);
#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
#pragma omp atomic update
    bucket[key[i]]++;
  }

  // build a global index for the buckets' starting locations
  std::vector<int> cumsum = bucket;
  for (int i = 1; i < range; i++)
  {
    cumsum[i] += cumsum[i - 1];
  }
  for (int i = range - 1; i > 0; i--)
  {
    cumsum[i] = cumsum[i - 1];
  }cumsum[0] = 0;

  int j = 0; // use a local index inside the bucket instead

  #pragma omp parallel
  for (int i = 0; i < range; i++)
  {
    
    for (; bucket[i] > 0; bucket[i]--)
    {
      key[j++] = i;
      printf("%d\n", omp_get_thread_num());
    }
  }

  for (int i = 0; i < n; i++)
  {
    printf("%d ", key[i]);
  }
  printf("\n");

  // check(key, n);
}
