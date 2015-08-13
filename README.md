# Tesis

    23/06/2015

        Tareas hechas:
            1) Se cambiaron muchos nombres de las variables globales
            2) Se definió como estatica el arreglo Sp
            3) Se hicieron variables globales con los arreglos que se definian en el main

        Problemas:
            1) Los arreglos del punto 3 de tarea se pasaron a global porque los 
                declaraba en el main y se le daba memoria en tiempo de ejecución, 
                y es mas, daba sigsfaul. tenes alguna idea de que puede ser

        Tareas para la proxima:
            1) Hacer que las variables globales sean globales y no se esten pasando por argunemnto
            2) Hacer el dibujo de como opera la matriz
            3) yapa: perf

    28/06/2015
        Tareas hechas:
            1) Se calculo de manera O(1) (con operaciones simpples) los arreglos t y k
            2) Se realizo el pasaje (o no pasaje) de parametros globales
            3) Se empezo a hacer el dibujo de la operacion matricial

        Problemas:
            --

        Tareas para la próxima:
            1) Terminar el dibujo de la addicion en la matriz
            2) Implemetar la primera idea de paralelización
            3) Usar los pragmas OPM

    09/06/2015
        Tareas hechas:
            1) Se modifico el orden de calculo para minimizar cuentas
            2) Se profundizo en el flujo de la funcion bder

        Problemas:
            1) La funcion bder hace cosas "de mas", preguntar el argumento teorico de esas cosas

        Tareas para la próxima
            1) Hacer el blocke omp con una slide window

    07/08/2015
        Tareas hechas:
            1) Se redujo la cantidad de memoria (de O(n²) aprox a O(n * kord) aprox) (donde kord es menos de 10 y mas 1)
            2) Se acomodo los ciclos para tener mas cache hit, por la nueva reestructuracion de la memoria
            3) Se cambio los ciclos para que solamente se calcule la parse superior de la matriz, luego se refleja

        Problemas:
            --

        Tareas para la próxima
            --

    Toma del pulso de mejoras:
        Compilando con el flag de optimizacion O3 y en "mini", con los siguientes parametros:

            EPS= 3.0e-14
            R_MIN= 0.0 
            R_MAX= 50.0
            L_INTERVALS= 510
            KORD= 5 
            RADIO_1= 5.0 
            RADIO_2= 10.0 
            ME= 1.0 
            INT_G= 500 
            NEV= 15 
            L_MAX= 0 
            LAMBDA_IN= 0.0 
            LAMBDA_FIN= 20.0 
            NUMEROS_PUNTO_LAMBDA= 200 
            BASE_KORD= 0

        Performance:
            Sólo se toma el tiempo que tarda el proceso calculo_matrices
            Código base:                        2.727035049000 segundos
            Con cacheo de la funcion bsplvb:    2.177257117000 segundos
            Reemplazando los arrelgos t y k:    1.811051266995 segundos
            Factorizacion de código:            1.676901567989 segundos
            Calculo de Knot_pesos by need       0.185952243977 segundos 
			Mejora de memoria, y calculo de     0.097898888867 segundos :')
                la parte superior de la matriz
                solamente
			Eliminado checkeo redundantes		0.078862830997 segundos
            Refactoriazacion del codigo         0.071713954909 segundos

            

    Apéndice:
            Lo primero que se hará para intentar paralelizar es usar algo asi como una "ventana corredisa" (en verdad esto fue lo primero que se me ocurrió) para operar en el circulo KNOTS para que no haya condicion de carrera cuando varios cores intenten escribir en el mismo lugar.


    PD: le pasaré un corrector, estoy haciendo esto mas para mi que para presentar, y tambien se que
    que le sirve al lector.
