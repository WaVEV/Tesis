#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "matrix_helper.h"
#include "helper_cuda.h" // check errors

#if __CUDA_ARCH__ < 600 
__device__ double atomicAdd(double* address, double val) { 
  unsigned long long int* address_as_ull = (unsigned long long int*)address; 
  unsigned long long int old = *address_as_ull, assumed;
  do { 
    assumed = old; 
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); 
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
  } while (assumed != old); return __longlong_as_double(old); } 
#endif


__global__ void integracion(double rr1, double w1, double *f, double *x, double *w, double* Sps)
{
  unsigned int nb = L_INTERVALS + KORD - 3; // tamaño de la base //

  int m = threadIdx.x;
  int n = threadIdx.y;
  int l = blockIdx.y;
  int k = blockIdx.x + KORD;
  unsigned int im = k - KORD + m - 1;
  unsigned int in = k - KORD + n - 1;
  int basek2 = k - KORD;

  double rr2 = eval_xi(basek2, l, x);
  double w2 = eval_wi(basek2, l, w);
  double *Sp = &Sps[idx_dense(basek2, l, INT_G, KORD)];

  if(im < nb)
  {
    if(in < nb)
    {
      if(rr2 <= rr1){
        atomicAdd(&f[idx2(im, in, nb)], Sp[m] * Sp[n] * w2/rr1);
      }else{
        atomicAdd(&f[idx2(im, in, nb)], Sp[m] * Sp[n] * w2/rr2);
      }
    }
  }
}

__global__ void sumatoria_vef(double w1, double *f, double *Sps, int basek, int i, int j, double* Vef, double* strace){
  double *Sp = &Sps[idx_dense(basek, j, INT_G, KORD)];
  int m = threadIdx.x;
  int mp = threadIdx.y;
  int im = i - KORD + m - 1;
  int imp = i - KORD + mp - 1;
  int n = blockIdx.x;
  int np = n - KORD + threadIdx.z;
  if(np >= 0 && im >= 0 && imp >= 0){
    double term = Sp[m]*Sp[mp]*w1* (f[idx(n, np, nb)]) / sqrt(strace[n] * strace[np]);
    atomicAdd(&Vef[idxVef(im, n, imp, np)], term);
  }

}


__host__ void call_integracion(double rr1, double w1, double *f, double *x, double *w, double* Sps, size_t f_size)
{
  checkCudaErrors(cudaMemset(f, 0, f_size));
  dim3 threadsPerBlock(KORD, KORD);
  dim3 BlockNums(L_INTERVALS, INT_G);
  //printf("%d %d %d\n", INT_G, KORD, KORD);
  integracion <<< BlockNums, threadsPerBlock >>> (rr1, w1, f, x, w, Sps);
  checkCudaErrors(cudaPeekAtLastError());
  cudaDeviceSynchronize();
}


__host__ void call_sumatoria_vef(double w1, double *f, double* Sps, int basek, int i, int j, double* Vef, double * strace){
                      // m,    mp,   np
  dim3 threadsPerBlock(KORD, KORD, 2*KORD);
  sumatoria_vef <<< L_INTERVALS + KORD - 3, threadsPerBlock >>> (w1, f, Sps, basek, i, j, Vef, strace);
}