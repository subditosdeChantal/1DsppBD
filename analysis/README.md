# Carpeta de análisis

Esta carpeta contiene códigos de análisis de datos: plots, etc. Se corren desde aquí y funcionan con todas las simulaciones.

1. Se pueden correr individualmente con gnuplot, ejemplo: __gnuplot -c plot_nclusters sim*[opciones de selección]__

1. O se pueden correr todos juntos: __bash run_analysis.sh sim*[opciones de selección]__

Donde __[opciones de selección]__ son especificaciones en el nombre de las carpetas a analizar. Por ejemplo, con __bash run_analysis.sh sim*f0.500*__ correrá un análisis de todas las simulaciones para fi=0.5. Para analizar todos los outputs basta con no indicar opciones de selección: __bash run_analysis.sh sim*__
