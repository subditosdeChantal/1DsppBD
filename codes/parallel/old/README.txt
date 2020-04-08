As it is the jobgen_run.sh script just creates a folder "jobs" and puts inside the compiled simulation script and different copies of the input file each one with different values for the parameters specefied at the beguining of jobgen_run.sh.
This is intended to be the job generation for later launching this jobs in SLURM as array jobs or by other means.

If we uncomment the lines marked with "# parallel" then jobgen_run.sh also runs the jobs that were created by using GNU-parallel ("parallel" command at the last line of the script), it launches each job in a different thread.

A brief explanation of the script is provided in spanish.





	#!/bin/bash

Estos son el nombre del script en c de la simulacion que vas a utlizar. Tienes que llamarlo ${scriptName}".c". En este caso se deberia llamar "para_brownian_dynamics.c".
EL sweepName da nombre al ejecutable que se va a compilar y a los jobs que se generaran.

	scriptName="para_brownian_dynamics"
	sweepName="prueba_paramod"
	#outputDir="output"
	execName=$scriptName"_"$sweepName

Estos son los avlores de las variables que van a tener los distintos archivos de input que genere este script.

	alpha=(0.001 0.010)
	phi=(0.1 0.2)
	eps=(0 1)

Numero the jobs simultaneos, si el ordenador donde lo lancemos tiene x threads es recomendable poner threads=x-1 para que no se sobrecargue.

	# threads=4 # parallel

"jobs" es la carpeta en la que saldran las copias del archivo de input, cada una con valores distintos de las variables anteriores.

	rm -rf jobs
	mkdir jobs
	# mkdir jobs/$outputDir

Compilamos el codigo de la simulacion.

	gcc -o jobs/$execName $scriptName".c" -lm

	idx=($(seq -w 000 719))
	ii=0

Generamos las copias de los archivos de input y sustituimos en cada una de ellas los valores correspondientes.

	# Generate jobs.
	for i in $(seq -w 0 $(( ${#alpha[@]} - 1)) )
	do
	        for j in $(seq -w 0 $(( ${#phi[@]} - 1)) )
	        do
	                for k in $(seq -w 0 $(( ${#eps[@]} - 1)) )
	                do

	                        jobName=$sweepName"_"${idx[${ii}]}

	                        cp "input.dat" jobs/$jobName".input"

Metemos los comandos en el archivo paralaunch.txt que al final del script daremos de comer al comando parallel.

	                        # echo "./"$execName" < "$jobName".input" >> jobs/paralaunch.txt # parallel

Realizamos las sustituciones.

	                        sed -i "1s/.*/Simulation name: "$sweepName"/" jobs/$jobName".input"

	                        # sed -i "2s/.*/Output dir: "$outputDir"/" jobs/$jobName".input"

	                        sed -i "5s/.*/Tumbling rate (alpha): "${alpha[${i#0}]}"/" jobs/$jobName".input"

	                        sed -i "6s/.*/Density (phi): "${phi[${j#0}]}"/" jobs/$jobName".input"

	                        sed -i "10s/.*/LJ pot. intensity (epsilon): "${eps[${k#0}]}"/" jobs/$jobName".input"

	                        ((ii++))
	                done
	        done
	done

Lanzamos con el comando parallel.

	# Launch jobs.
	# cd jobs # parallel
	# parallel -j $threads < paralaunch.txt # parallel

