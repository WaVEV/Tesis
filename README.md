Pozo-2
===================

Programa en C++ que resuelve el problema de una particula 
en un pozo de potencial usando B-splines.
El potencial es de la forma V(r) = -lambda si r entre r1 y r2 y cero fuera.
El programa guarda los autovalores en funcion de lambda.
Usa el metodo variacional de Rayleight-Ritz usando como base del 
espacio los B-splines, como esta no es una base ortonormal el 
problema de autovalores queda de la siguiente forma:
            `H |psi> = e S |psi>`
donde H es la matriz del hamiltoniano y S es la matriz de solapamiento
de los B-splines.
Usa la funcion gauleg() del Numerical Recipies C para calcular los 
puntos de evaluacion y los pesos para la cuadratura.
La funcion KNOTS_PESOS() esta hecha por mi guiandome de la version de 
fortran, calcula los knots de los B-splines con una distribucion 
uniforme solamente, en el caso de querer otra distribucion es solo 
cuestion de modificar el codigo para usar la distribucion que uno quiera.
Para el calculo de los autovalores usa la funcion dsygvx_() de lapack.
Las funciones para evaluar los B-splines y las derivadas son bsplvb() y 
bder() respectivamente, son versiones hechas por mi a partir de las 
versiones de fortran.

Una observacion importante: este codigo anda solo para l>=25.

Al final del archivo se encuentran las diferencias con la version
anterior, leerlas para tener en cuenta.

## Install ##

Para instalar primero clonamos el repositorio de manera recursiva (para descargar arpack)

    git clone --recursive https://github.com/WaVEV/Tesis

Ahora compilaremos arpack e instalaremos los complementos para matrices ralas.

    cd Tesis
    sh instalar_dependencias.sh

(tarda un poco en mostrar progreso, sea paciente)


## Correr el programa ##

Para compilar solo hay que hacer make

    make

Luego para correr:

    ./pozo-2

