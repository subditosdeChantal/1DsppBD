#!/bin/bash

scriptName="BDp_02"
sweepName="L2000"
#outputDir="output"
execName=$scriptName"_"$sweepName

alpha=(0.1)
phi=(0.5)
eps=(0 24)
dt=(0.01)
# tsim=(100)
# tmeas=(1000 10000 100000)
v=(1)
beta=(1000 10)

# threads=14

rm -rf jobs

mkdir jobs
# mkdir jobs/$outputDir

# Hay que poner lo de -std por las declaraciones de variable dentro de bucles.
gcc -o jobs/$execName -std=gnu99 $scriptName".c" -lm

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

						# echo "./"$execName" < "$jobName".input > "$jobName".out" >> jobs/paralaunch.txt

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

# Launch jobs.
# cd jobs
# parallel -j $threads < paralaunch.txt
