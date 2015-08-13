/* Programa en C que resuelve el problema de una particula 
 * en un pozo de potencial usando B-splines.
 * El potencial es de la forma V(r) = -lambda si r1<r<r2 y cero fuera.
 * El programa guarda los autovalores en funcion de lambda.
 * Usa el metodo variacional de Rayleight-Ritz usando como base del 
 * espacio los B-splines, como esta no es una base ortonormal el 
 * problema de autovalores queda de la siguiente forma:
 *
 *              H |psi> = e S |psi>
 *
 * donde H es la matriz del hamiltoniano y S es la matriz de solapamiento
 * de los B-splines.
 *
 * Usa la funcion gauleg() del Numerical Recipies C para calcular los 
 * puntos de evaluacion y los pesos para la cuadratura.
 * La funcion KNOTS_PESOS() esta hecha por mi guiandome de la version de 
 * fortran, calcula los knots de los B-splines con una distribucion 
 * uniforme solamente, en el caso de querer otra distribucion es solo 
 * cuestion de modificar el codigo para usar la distribucion que uno quiera.
 * Para el calculo de los autovalores usa la funcion dsygvx_() de lapack.
 * Las funciones para evaluar los B-splines y las derivadas son bsplvb() y 
 * bder() respectivamente, son versiones hechas por mi a partir de las 
 * versiones de fortran.
 *
 * Una observacion importante: este codigo anda solo para l>=25, 
 * para valores mas chicos hace cosas raras no se bien por que.
 *
 * Al final del archivo se encuentran las diferencias con la version 
 * anterior, leerlas para tener en cuenta.
*/

#include <stdio.h> // prinft
#include <stdlib.h> // malloc y free
#include <malloc.h>
#include <math.h> // posibles operaciones matematicas
#include <assert.h>
#include <string.h>
#include <omp.h> //omp_get_wtime()

// defino algunas constantes para el programa //
#define EPS 3.0e-14  // EPS precision relativa para gauleg //

// defino algunos parametros del programa //

#ifndef R_MIN
#define R_MIN 0.0 // R minimo donde empieza el intervalo para la integracion //
#endif

#ifndef R_MAX
#define R_MAX 50.0 // R maximo donde termina el intervalo para la integracion //
#endif

#ifndef L_INTERVALS
#define L_INTERVALS 510 // numero de intervalos en el que divido al intervalo [R_MIN, R_MAX] //
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
#define INT_G 500 // grado de integracion por cuadratura //
#endif

#ifndef NEV
#define NEV 15 // numero de autovalores que vamos a calcular //
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

const unsigned int nk = L_INTERVALS+2*KORD-1;
double x[INT_G], w[INT_G];

const unsigned int nb = L_INTERVALS + KORD-3; // tamaño de la base //
    
double  s[ (L_INTERVALS+2*KORD-3) * (KORD)],
        v0[(L_INTERVALS+2*KORD-3) * (KORD)],
        ke[(L_INTERVALS+2*KORD-3) * (KORD)];

// escribo las funciones del programa //
int dsygvx_(int *itype, char *jobz, char *range, char * uplo, 
    int *n, double *a, int *lda, double *b, int *ldb, 
    double *vl, double *vu, int *il, int *iu, double *abstol, 
    int *m, double *w, double *z__, int *ldz, double *work, 
    int *lwork, int *iwork, int *ifail, int *info);


#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

int idx(const unsigned int y, const unsigned int x, const unsigned int numcolumns){
    int i, j;
    j = max(x,y);
    i = min(x,y);
    if(j - i >= KORD) return 0;
    i = KORD - (j - i) - 1;

    assert((L_INTERVALS+2*KORD-3) * (KORD) > i*numcolumns + j);
    return i*numcolumns + j;
}

int cleari(const unsigned int N, int * __restrict__ vec) {
    
    for(unsigned int i = 0; i<N; ++i) 
        vec[i] = 0;

    return 0;
}

int cleard(const unsigned int N, double * __restrict__ vec) {
    memset(vec, 0, sizeof(double) * N);
    //for(unsigned int i = 0; i<N; ++i)
    //    vec[i] = 0.f;

    return 0;
}

void setxparameters(const int i, double *xm, double *xl){
    double dr = (R_MAX-R_MIN)/L_INTERVALS;
    double x1 = R_MIN + i * dr, 
           x2 = x1 + dr;
    *xm = 0.5*(x2 + x1);
    *xl = 0.5*(x2 - x1);
}

double eval_xi(const int i, const int j, const double * x){
    double xm, xl;
    setxparameters(i, &xm, &xl);
    if(j >= (INT_G + 1) / 2){
        return xm + x[j] * xl;
    }else{
        return xm - x[j] * xl;
    }
}

double eval_wi(const int i, const int j, const double * w){
    double xm, xl;
    setxparameters(i, &xm, &xl);
    return 2.0*xl / w[j];
}

void gaulegm(double * __restrict__ x, double * __restrict__ w, int n) {
/* Given the lower and upper limits of integration x1 and x2, 
 * and given n, this routine returns arrays x[1..n] and w[1..n]
 * of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
*/
    int m, j, i;
    double z1, z, pp, p3, p2, p1;

    m = (n+1)/2;

    for (i = 0; i<m; i++) {
        z = cos(M_PI*((i+1)-0.25)/(n+0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1;j<=n;j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp = n*(z*p1-p2)/(z*z-1.0);
            z1 = z;
            z = z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i] = z;
        x[n-i-1] = z;
        w[i] = (1.0-z*z)*pp*pp;
        w[n-i-1] = (1.0-z*z)*pp*pp;
    }
}


double ti(int i){
    double dr = (R_MAX-R_MIN)/L_INTERVALS;
    int pos = (i - KORD + 1);
    return R_MIN + (KORD - 1 < i) * dr * (pos - (pos - L_INTERVALS) * (KORD + L_INTERVALS - 1 < i));
}

int bsplvb(unsigned int jhigh, double rr, int left, double * __restrict__ biatx, double * __restrict__ biatx_1) {
  
    unsigned int j;
    double saved, term;
    

    double deltar[JMAX];
    double deltal[JMAX];

    biatx[0] = 1.0;


    for(j=0; j<jhigh-2; ++j) {
            
        deltar[j] = ti(left+j+1) - rr;
        deltal[j] = rr-ti(left-j);

        saved = 0.0;
        for(unsigned int i = 0; i<j+1; ++i) {
            term = biatx[i]/(deltar[i]+deltal[j-i]);
            biatx[i] = saved + deltar[i]*term;
            saved = deltal[j-i]*term;
        }
        biatx[j+1] = saved;

    }
    memcpy(biatx_1, biatx, sizeof(double) * KORD);

    deltar[j] = ti(left+j+1) - rr;
    deltal[j] = rr-ti(left-j);

    saved = 0.0;
    for(unsigned int i = 0; i<j+1; ++i) {
        term = biatx[i]/(deltar[i]+deltal[j-i]);
        biatx[i] = saved + deltar[i]*term;
        saved = deltal[j-i]*term;
    }
    biatx[j+1] = saved;

    return 0; 
}


double bder(unsigned int indexm, unsigned int left, double * __restrict__ Sp) {

    unsigned int i;
    double dm = 0;
    //assert(ti(0)<rr && rr<ti(nk-1));
    
    //if(ti(0)<rr && rr<ti(nk-1)) {
        //int cnt = 0;
        //if(abs(rr-ti(nk-1))<1.e-10) {
          //  cnt |= 1;

            //dm = ((indexm == nk - KORD)*2-1) * (KORD-1)/(ti(nk-1)-ti(nk-1-KORD));
            /*if(indexm == nk - KORD) {
                dm = (KORD-1)/(ti(nk-1)-ti(nk-1-KORD));
            }
            else {
                dm = -(KORD-1)/(ti(nk-1)-ti(nk-1-KORD));
            }*/
        //}

        /*else */if(indexm-left+KORD>=1) {
            i = indexm-left+KORD;
            if(1==i) {
                dm = (KORD-1)*(-Sp[i-1]/(ti(indexm+KORD)-ti(indexm+1)));
            }
            else if(KORD==i) {
                dm = (KORD-1)*(Sp[i-1-1]/(ti(indexm+KORD-1)-ti(indexm)));
            }
            else {
                dm = (KORD-1)*(Sp[i-1-1]/(ti(indexm+KORD-1)-ti(indexm))
                    - Sp[i-1]/(ti(indexm+KORD)-ti(indexm+1)));
            }
        }

        //if(cnt == 2) printf("entre solo al segundo \n");
    //}


    return dm;
}

void calculo_matrices(const double * __restrict__ const x, const double * __restrict__ const w,
              double * __restrict__ s, double * __restrict__ v0,
              double * __restrict__ ke) {

    double ma;
    double Sp[KORD];
    double Sp_1[KORD];

    ma = 0.5*L_MAX*(L_MAX+1);
    
    double rr, _rr2;
    

    double bders[KORD];
    for(unsigned int basek=1; basek<L_INTERVALS; ++basek) {
        unsigned int i = basek - 1 + KORD;
        for(unsigned int j = 0; j<INT_G; ++j) {
            rr = eval_xi(basek, j, x);
            _rr2= 1.0/(rr*rr);

            bsplvb(KORD, rr, i, Sp, Sp_1);

            for(unsigned int m = i-KORD+1; m<=i && m<nb ; ++m) {
                bders[m - (i-KORD+1)] = bder(m, i, Sp_1);
            }

            for(int k=0 ; k<KORD ; k++){
                for(unsigned int m = i-KORD+1, n = m + k; n<=i && n<nb ; ++m, ++n) {
                    double  bm = bders[m - (i-KORD+1)];
                    double bn = bders[n - (i-KORD+1)];
                    ke[idx(m-1, n-1, nb)] += 0.5*eval_wi(basek, j, w)*bm*bn/ME;
                }
            }

            double wikj = eval_wi(basek, j, w);

            for(int k=0 ; k<KORD ; k++){
                for(unsigned int m = 0, im = i-KORD, n = k, in = i-KORD + k; 
                    in<nb && n < KORD;
                    ++m, ++im, ++n, ++in) {
                    
                    s[idx(im, in, nb)] += Sp[m] * Sp[n] * wikj;

                    ke[idx(im, in, nb)] += ma*Sp[m] * Sp[n] * wikj * _rr2;

                    if(RADIO_1<rr && rr<RADIO_2) v0[idx(im, in, nb)] += Sp[m] * Sp[n] * wikj;
                }
            }
        }
    }

    for(unsigned int j = 0; j<INT_G; ++j) {
        rr = eval_xi(BASE_KORD, j, x);
        _rr2= 1.0/(rr*rr);

        bsplvb(KORD, rr, KORD-1, Sp, Sp_1);

        for(unsigned int m = 1; m<=KORD-1; ++m) {
            bders[m] = bder(m, KORD-1, Sp_1);
        }

        for(unsigned int m = 1; m<=KORD-1; ++m) {
            for(unsigned int n = m; n<=KORD-1; ++n) {

                double  bm = bders[m],
                        bn = bders[n];

                ke[idx(m-1, n-1, nb)] += 0.5* eval_wi(BASE_KORD, j, w)*bm*bn/ME;
            }
        }

        double wikj = eval_wi(0, j, w);

        for(int k=0 ; k<KORD ; k++){
            for(unsigned int m = 1, im = 0,  n = 1 + k, in = k; 
                n < KORD && in<nb ;
                ++m, ++im, ++n, ++in) {

                s[idx(im, in, nb)] += Sp[m] * Sp[n] * wikj;

                ke[idx(im, in, nb)] += ma*Sp[m] * Sp[n] * wikj * _rr2;

                if(RADIO_1<rr && rr<RADIO_2) v0[idx(im, in, nb)] += Sp[m] * Sp[n] * wikj;
            }
        }
    }
}


/*
void eigenvalues(int n, int m, double * __restrict__ a,
         double * __restrict__ b, double * __restrict__ w, 
         double * __restrict__ z) {
    int itype, lda, ldb, ldz;
    int il, iu, lwork, info;
    char jobz, range, uplo;
    double abstol, vl, vu;
    double * work;
    int * iwork, * ifail;
    // le doy los valores correspondientes a las distintas variables //
    itype = 1; lda = n; ldb = n; ldz = n;
    vl = 0.0; vu = 0.0;
    jobz = 'V'; range = 'I'; uplo = 'U';
    il = 1; iu = m; lwork = 9*n;
    // le doy memoria a las matrices que neceista dsygvx //
    work = (double *) malloc(lwork*sizeof(double));
    iwork = (int *) malloc(5*n*sizeof(int));
    ifail = (int *) malloc(n*sizeof(int));
    dsygvx_( &itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, 
         &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, 
         &lwork, iwork, ifail, &info);
    free(work); 
    free(iwork); 
    free(ifail);
} 
void hamiltoniano_autovalores(unsigned int nb, double * __restrict__ s,
                  double * __restrict__ v0, double * __restrict__ ke,
                  FILE * archivo) {
    double * h, * auval, * auvec, * s_copy;
    double lambda, delta;
    
    // doy memoria a las matrices y vectores //
    h = (double *) calloc( nb*nb, sizeof(double));
    s_copy = (double *) calloc( nb*nb, sizeof(double));
    auval = (double *) calloc( NEV, sizeof(double));
    auvec = (double *) calloc( NEV*nb, sizeof(double));
    
    lambda = 1.0;
    delta = (LAMBDA_FIN-LAMBDA_IN)/NUMEROS_PUNTO_LAMBDA;
    // abro el archivo para guardar los datos //
    fprintf(archivo, "# Autovalores calculados\n");
    fprintf(archivo, "# Lambda  auval[0]   auval[1] ....\n");
    for(unsigned int j = 0; j<=NUMEROS_PUNTO_LAMBDA; ++j) {
        lambda = LAMBDA_IN + delta*j;
    
        for(unsigned int m = 0; m<nb; ++m) {
            for(unsigned int n = 0; n<nb; ++n){
                h[idx(n, m, nb)] = ke[idx(n, m, nb)] - lambda*v0[idx(n, m, nb)];
                s_copy[idx(n, m, nb)] = s[idx(n, m, nb)];
            }
        }
        eigenvalues( nb, NEV, h, s_copy, auval, auvec );
        fprintf(archivo, "%.5f   ", lambda);
        for(unsigned int i = 0; i<NEV; ++i) {
            fprintf(archivo, "%.15f   ", auval[i]);
        }
        fprintf(archivo, "\n");
    }
    free(h); free(s_copy); free(auval); free(auvec);
}*/

int main(void) {

    // defino algunas variables que voy a usar //
    double t_in, t_fin, t_n;
    FILE * archivo;

    // controlo algunos parametros //
    assert(INT_G>KORD);
    assert(NEV>0);

    archivo = fopen("autovalores_pozo-1.dat", "w");
    // imprimo los parametros //
    fprintf(archivo, "# Rmin=%.12f y R_MAX=%.12f\n", R_MIN, R_MAX);
    fprintf(archivo, "# Numero de intervalos l=%i\n", L_INTERVALS);
    fprintf(archivo, "# Orden los B-splines kord=%i\n", KORD);
    fprintf(archivo, "# Radios del pozo RADIO_1=%.12f y RADIO_2=%.12f\n", RADIO_1, RADIO_2);
    fprintf(archivo, "# Masa de la particula me=%.12f\n", ME);
    fprintf(archivo, "# Grado de integracion de la cuadratura INT_G=%i\n", INT_G);
    fprintf(archivo, "# Numero de autovalores NEV=%i\n", NEV);
    fprintf(archivo, "# Momento angular que usamos L_MAX=%i\n", L_MAX);
    fprintf(archivo, "# Numero de knots nk=%i\n", nk);
    fprintf(archivo, "# Tamaño de la base nb=%i\n", nb);
    fprintf(archivo, "# Valores inicial y final del pozo, LAMBDA_IN=%.12f, LAMBDA_FIN=%.12f\n", LAMBDA_IN, LAMBDA_FIN);
    fprintf(archivo, "# Numero de puntos lambda = %i\n", NUMEROS_PUNTO_LAMBDA);

   
    cleard(INT_G, x);
    cleard(INT_G, w);
  /*  cleard(nb * (2 * KORD - 1), s);
    cleard(nb * (2 * KORD - 1), v0);
    cleard(nb * (2 * KORD - 1), ke);*/

    // primero calculos los knost y los pesos para hacer la cuadratura //

    t_in = omp_get_wtime();
    gaulegm(x, w, INT_G);
    calculo_matrices(x, w, s, v0, ke);
    t_fin = omp_get_wtime();

    // armo el hamiltoniano y calculo los autovalores //
    // hamiltoniano_autovalores(nb, s, v0, ke, archivo);
    
    t_n = (t_fin-t_in)/(nb*nb*L_INTERVALS*INT_G),
    printf("%i     %i     %i     %.12f  %.12f\n", L_INTERVALS, nk, nb, t_fin-t_in, t_n); 
    FILE *file;
    file = fopen("matrices.dat", "w");
    for(unsigned int i = 0; i<nb; ++i)
        for(unsigned int j=0; j<nb; ++j)
            fprintf(file, "%i\t%i\t%.12f\t%.12f\t%.12f\n", i, j, s[idx(i,j,nb)], v0[idx(i,j,nb)], ke[idx(i,j,nb)]);
    fclose(file);

    // libero la memoria //
    fclose(archivo);

    return 0;
}



