#ifndef MATRIX_HELPER_POZO_H_
#define MATRIX_HELPER_POZO_H_

#ifndef R_MIN
#define R_MIN 0.0 // R minimo donde empieza el intervalo para la integracion //
#endif

#ifndef R_MAX
#define R_MAX 10.0 // R maximo donde termina el intervalo para la integracion //
#endif

#ifndef L_INTERVALS
#define L_INTERVALS 200 // numero de intervalos en el que divido al intervalo [R_MIN, R_MAX] //
#endif

#ifndef KORD
#define KORD 5 // orden de los B-splines, el grado es kord-1 //
#endif

#ifndef RADIO_1
#define RADIO_1 5.0 // radio interio del pozo //
#endif

#ifndef RADIO_2
#define RADIO_2 10.0 // radio exterior del pozo //
#endif

#ifndef ME
#define ME 1.0 // masa de la particula //
#endif

#ifndef INT_G
#define INT_G 100 // grado de integracion por cuadratura //
#endif

#ifndef NEV
#define NEV 10 // numero de autovalores que vamos a calcular //
#endif

#ifndef L_MAX
#define L_MAX 0 // momento angular que vamos a usar //
#endif

#ifndef LAMBDA_IN
#define LAMBDA_IN 0.0 // lambda inicial para el pozo //
#endif

#ifndef LAMBDA_FIN
#define LAMBDA_FIN 20.0 // lambda final para el pozo //
#endif

#ifndef NUMEROS_PUNTO_LAMBDA
#define NUMEROS_PUNTO_LAMBDA 200 // numero de puntos para calcular //
#endif

#ifndef BASE_KORD
#define BASE_KORD 0
#endif

#ifndef JMAX
#define JMAX 100
#endif

#ifndef ETAI
#define ETAI 0.1
#endif

#ifndef ETAF
#define ETAF 1.0
#endif

#ifndef NUM_PUNTOS_ETA
#define NUM_PUNTOS_ETA 10
#endif

#include <cuda_runtime.h>
#include <assert.h>

static const unsigned int nb = L_INTERVALS + KORD - 3; // tamaño de la base //
static const unsigned int nv = (L_INTERVALS+2*KORD-3) * 2 * (KORD);


__host__ __device__ static int idx2(const unsigned int y, const unsigned int x, const unsigned int numcolumns){
  int i=y, j=x;
  if(i - j >= KORD || j - i >= KORD) {
    return 0;
  }
  i = KORD + ( j - i ) - 1;
  j += y > x ? y - x : x - y;
  assert((L_INTERVALS+2*KORD-3) * (2*KORD) > i*numcolumns + j);
  return i*numcolumns + j;
}

__host__ __device__ static int idx(const unsigned int y, const unsigned int x, const unsigned int numcolumns){
  int i, j;
  j = x > y ? x:y;
  i = x > y ? y:x;
  if(j - i >= KORD) return 0;
  i = KORD - (j - i) - 1;

  assert((L_INTERVALS+2*KORD-3) * (KORD) > i*numcolumns + j);
  return i*numcolumns + j;
}

__host__ __device__ static int idxVef(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int w)
{
  return idx2(x, z, nb) * nv + idx2(y, w, nb);
}

__host__ __device__ static int idx_dense(const unsigned int y, const unsigned int x, const unsigned int n, const unsigned int m){
  return y * n * m + x * m;
}


__host__ __device__ static void setxparameters(const int i, double *xm, double *xl){
    const double dr = (R_MAX-R_MIN)/L_INTERVALS;
    const double x1 = R_MIN + i * dr, 
           x2 = x1 + dr;
    *xm = 0.5*(x2 + x1);
    *xl = 0.5*(x2 - x1);
}

__host__ __device__ static double eval_xi(const int i, const int j, const double * x){
    double xm, xl;
    setxparameters(i, &xm, &xl);
    if(j >= (INT_G + 1) / 2){
        return xm + x[j] * xl;
    }else{
        return xm - x[j] * xl;
    }
}

__host__ __device__ static double eval_wi(const int i, const int j, const double * w){
    double xm, xl;
    setxparameters(i, &xm, &xl);
    return 2.0*xl / w[j];
}


__global__ void integracion(double rr1, double w1, double *f, double *x, double *w, double* Sps);

__host__ void call_integracion(double rr1, double w1, double *f, double *x, double *w, double* Sps, size_t f_size);
__host__ void call_sumatoria_vef(double w1, double *f, double* Sps, int basek, int i, int j, double* Vef, double * strace);

#endif