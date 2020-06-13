The simSet.sh script just creates a folder "jobs" and puts inside the compiled simulation script and different copies of the input file each one with different values for the parameters specefied at the beguining of simSet.sh.
This is intended to be the job generation for later launching this jobs in SLURM as array jobs or by other means (gnu-parallel).

Now simSet.sh also generates a script called simCont.sh for continuating the simulation from the last step of the current simulation. When run, simCont.sh does the following:

	- Creates a new simulation directory with the same name as the current simulation folder plus "_cont". And copies the current simulation scripts to the new folder.

	- It modifies the simulation names in the scripts adding them the "_cont". And sets the initial state of the new simulation to 4 (enabling the reading of the last frames).

	- It runs simSet.sh for the new simulation.

	- It extracts the last frames from the current simulation and puts them in the jobs folder of the new one. 

So we just need to run simRun.sh in the new (*_cont) simulation and it will resume from the last frame of the current simulation.

Now simSet.sh DOESN'T run the simulation, this is done by simRun.sh which also launches the analysis scripts.



A brief explanation of the simSet.sh script is provided in spanish.


	#!/bin/bash

Estos son el nombre del script en c de la simulacion que vas a utlizar. Tienes que llamarlo ${scriptName}".c". En este caso se deberia llamar "BDp_051.c".
EL sweepName da nombre al ejecutable que se va a compilar y a los jobs que se generaran y debe ser el mismo nombre de la carpeta de simulacion.

	scriptName="BDp_051"
	sweepName="init"
	#outputDir="output"
	execName=$scriptName"_"$sweepName

Estos son los avlores de las variables que van a tener los distintos archivos de input que genere este script.

	alpha=(0.1 0.5)
	phi=(0.1 0.8)
	eps=(0 24)
	dt=(0.01)
	# tsim=(100)
	# tmeas=(1000 10000 100000)
	v=(1)
	beta=(10)

Numero the jobs simultaneos, si el ordenador donde lo lancemos tiene x threads es recomendable poner threads=x-1 para que no se sobrecargue.

	threads=8

"jobs" es la carpeta en la que saldran las copias del archivo de input, cada una con valores distintos de las variables anteriores.

	rm -rf jobs

	mkdir jobs
	# mkdir jobs/$outputDir

	echo $threads > jobs/procfile

Compilamos el codigo de la simulacion.

	# Hay que poner lo de -std por las declaraciones de variable dentro de bucles.
	gcc -o jobs/$execName -std=gnu99 $scriptName".c" -lm

Generamos las copias de los archivos de input y sustituimos en cada una de ellas los valores correspondientes.

	idx=($(seq -w 0 $(( ${#alpha[@]}*${#phi[@]}*${#eps[@]}*${#dt[@]}*${#v[@]}*${#beta[@]} - 1 )) ))
	ii=0

	# Generate jobs.
	for i in $(seq -w 0 $(( ${#alpha[@]} - 1)) )
	do
		for j in $(seq -w 0 $(( ${#phi[@]} - 1)) )
		do
			for k in $(seq -w 0 $(( ${#eps[@]} - 1)) )
			do
	            for l in $(seq -w 0 $(( ${#dt[@]} - 1)) )
	            do
	                for m in $(seq -w 0 $(( ${#v[@]} - 1)) )
	                do
	                	for o in $(seq -w 0 $(( ${#beta[@]} - 1)) )
	                	do
							jobName=$sweepName"_"${idx[${ii}]}
							
							# jobName=$sweepName"_a"${alpha[${i}]}"_f"${phi[${j}]}"_e"${eps[${k}]}"_dt"${dt[${l}]}"_v"${v[${m}]}"_b"${beta[${o}]}

							cp "input.dat" jobs/$jobName".input"

Metemos los comandos en el archivo paralaunch.txt que al final del script daremos de comer al comando parallel.

							echo "./"$execName" < "$jobName".input > "$jobName".out" >> jobs/paralaunch.txt

Realizamos las sustituciones.

							sed -i "1s/.*/Simulation name: "$sweepName"/" jobs/$jobName".input"

							# sed -i "2s/.*/Output dir: "$outputDir"/" jobs/$jobName".input"

							sed -i "5s/.*/Tumbling rate (alpha): "${alpha[${i#0}]}"/" jobs/$jobName".input"

							sed -i "6s/.*/Density (phi): "${phi[${j#0}]}"/" jobs/$jobName".input"

		                  	# sed -i "7s/.*/Simulation time (time units): "${tsim[${l#0}]}"/" jobs/$jobName".input"

							sed -i "10s/.*/LJ pot. intensity (epsilon): "${eps[${k#0}]}"/" jobs/$jobName".input"

					        sed -i "12s/.*/Propulsion force (v, Fp): "${v[${m#0}]}"/" jobs/$jobName".input"

					   		sed -i "14s/.*/Temperature (beta): "${beta[${o#0}]}"/" jobs/$jobName".input"

							sed -i "18s/.*/Time Step: "${dt[${l#0}]}"/" jobs/$jobName".input"

							((ii++))
						done
					done
				done
			done
		done
	done

Creamos el script simCont.sh

	# Create script for continuating simulation from the last step.
	rm -f simCont.sh
	cat > simCont.sh << EOF
	#!/bin/bash

	mkdir ../../${sweepName}_cont
	mkdir ../../${sweepName}_cont/scripts
	cp ${scriptName}.c ../../${sweepName}_cont/scripts
	cp simRun.sh simSet.sh input.dat ../../${sweepName}_cont/scripts
	cp -r analysis ../../${sweepName}_cont/scripts

	sed -i "1s/.*/Simulation name: "${sweepName}"_cont/" "../../"${sweepName}"_cont/scripts/input.dat"
	sed -i "17s/.*/Initial state: 4/" "../../"${sweepName}"_cont/scripts/input.dat"
	sed -i "4s/.*/sweepName="${sweepName}"_cont/" "../../"${sweepName}"_cont/scripts/simSet.sh"

	cd ../../${sweepName}_cont/scripts
	./simSet.sh
	cd -

	cd ../output_1DsppBD
	dirnames=\$(ls | grep sim)
	for dir in \$dirnames
	do
		if [[ -d \$dir ]];
		then
			cd \$dir
			tail -n 1 system_ss.dat > ../../../${sweepName}_cont/scripts/jobs/frame0_\${dir:4:14}\${dir:30:40}.dat
			cd ..
		fi
	done

	EOF

Y le damos permisos de ejecucion.

	chmod u+x simCont.sh