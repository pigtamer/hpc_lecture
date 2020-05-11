#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

int main()
{
  int n = 50;
  int range = 100;
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
  {
    for (int i = 0; i < n; i++)
    {
      bucket[key[i]]++;
      printf("# %d ", omp_get_thread_num());
    }
  } // put value to the corresponding position
  printf("\n");

  int j = 0;
  for (int i = 0; i < range; i++)
  {
    for (; bucket[i] > 0; bucket[i]--)
    {
      j++;
      key[j] = i; // write back

      printf("# %d ", omp_get_thread_num());
    }
  }

  printf("\n");
  for (int i = 0; i < n; i++)
  {
    printf("%d ", key[i]);
  }
  printf("\n");
}
