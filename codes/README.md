# Carpeta de códigos

Esta carpeta contiene los códigos de simulación. Hay 3 ficheros además del ejecutable.

* *Fichero 1*: __brownian_dynamics.c__
  * Código en __C#__ para la dinámica del sistema. No se cambia de una simulación a otra, lee los parámetros de input.dat.
  * Con cada nueva simulación se crea una nueva carpeta (salvo que ya exista una carpeta para los mismos parámetros) en __~/output_1DsppBD/__ que a su vez se crea con la primera simulación.

* *Fichero 2*: __input.dat__
  * Contiene todos los parámetros necesarios para correr una simulación. Aquí es donde se modifica antes de lanzar una simulación, no en __brownian_dynamics.c__.

* *Fichero 3*: __input_HOWTO.dat__
  * Contiene instrucciones para pasar parámetros en __input.dat__. Se explica qué parámetro corresponde a cada línea del fichero y cómo se debe interpretar e indicar dicho parámetro.

**IMPORTANTE**: Al menos de momento, las dos carpetas del repositorio deben estar a distancia __../..__ de home para que todo funcione correctamente. El repositorio debe colocarse en home por lo tanto.
