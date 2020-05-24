#include <cstdio>
#include <cstdlib>
#include <vector>

bool check(int*vec, size_t n)
{
  bool res = true;
  for (size_t k = 0; k < n - 1; k++)
  {
    res = res & (vec[k] <= vec[k + 1]);
  }
  printf("%d\n", res);
  return res;
}

__global__ void bucketSort(int *key, int *bucket, int *a, int *b, int n, int range)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  bucket[i] = 0;
  for (int j=0; j<n; j++){
    if(i==key[j])
    {bucket[i]+=1;
    
    }
__syncthreads();
    
    a[i] = bucket[i];

b[i] = bucket[i];
__syncthreads();

}

  for (int j = 1; j < range; j <<= 1)
  {
    b[i] = a[i];
    __syncthreads();
    if(i>=j)
    a[i] += b[i-j];
    __syncthreads();
  }
  
  printf(" - %d %d \n", a[i], i);
  
  for (int j=0; bucket[i] > 0; bucket[i]--)
  {

    key[j + a[i - 1]] = i;
    j++;
    __syncthreads();
  }


}

int main()
{
  int n = 1000;
  int range = 1024; // seem to be an upper limit for a single block
  
  int *key, *a, *b, *bucket;
  cudaMallocManaged(&key, n * sizeof(int));
  cudaMallocManaged(&a, range * sizeof(int));
  cudaMallocManaged(&b, range * sizeof(int));
  cudaMallocManaged(&bucket, range * sizeof(int));
  
  for (int i = 0; i < n; i++)
  {
    key[i] = rand() % range;
    printf("%d ", key[i]);
  }
  printf("\n");

  bucketSort<<<1, range>>>(key, bucket, a, b, n, range);

  cudaDeviceSynchronize();
  printf("\n");
  for (int i = 0; i < n; i++)
  {
    printf("%d ", key[i]);
  }
  printf("\n");

  check(key, n);

  printf("\n");
  cudaFree(a);
  cudaFree(b);
  cudaFree(bucket);
  cudaFree(key);
}
