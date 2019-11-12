# Carpeta de códigos

Esta carpeta contiene los códigos de simulación. Hay 3 ficheros además del ejecutable.

1. Fichero 1: brownian_dynamics.c
 2. Código en C# para la dinámica del sistema. No se cambia de una simulación a otra, lee los parámetros de input.dat.
 2. Con cada nueva simulación se crea una nueva carpeta (salvo que ya exista una carpeta para los mismos parámetros) en ~/output_1DSPPdynamics/ que a su vez se crea con la primera simulación.

1. Fichero 2: input.dat
 2. Contiene todos los parámetros necesarios para correr una simulación. Aquí es donde se modifica antes de lanzar una simulación, no en brownian_dynamics.c.

1. Fichero 3: input_HOWTO.dat
 2. Contiene instrucciones para pasar parámetros en input.dat. Se explica qué parámetro corresponde a cada línea del fichero y cómo se debe interpretar e indicar dicho parámetro.


**IMPORTANTE**: Al menos de momento, las dos carpetas del repositorio deben estar a distancia ../.. de home para que todo funcione correctamente. El repositorio debe colocarse en home por lo tanto.
