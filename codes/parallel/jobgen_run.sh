#!/bin/bash

scriptName="para_brownian_dynamics"
sweepName="prueba_paramod"
#outputDir="output"
execName=$scriptName"_"$sweepName

alpha=(0.001 0.010)
phi=(0.1 0.2)
eps=(0 1)

# threads=4

rm -rf jobs

mkdir jobs
# mkdir jobs/$outputDir

gcc -o jobs/$execName $scriptName".c" -lm

idx=($(seq -w 000 719))
ii=0

# Generate jobs.
for i in $(seq -w 0 $(( ${#alpha[@]} - 1)) )
do
	for j in $(seq -w 0 $(( ${#phi[@]} - 1)) )
	do
		for k in $(seq -w 0 $(( ${#eps[@]} - 1)) )
		do
			jobName=$sweepName"_"${idx[${ii}]}

			cp "input.dat" jobs/$jobName".input"

			# echo "./"$execName" < "$jobName".input" >> jobs/paralaunch.txt

			sed -i "1s/.*/Simulation name: "$sweepName"/" jobs/$jobName".input"

			# sed -i "2s/.*/Output dir: "$outputDir"/" jobs/$jobName".input"

			sed -i "5s/.*/Tumbling rate (alpha): "${alpha[${i#0}]}"/" jobs/$jobName".input"

			sed -i "6s/.*/Density (phi): "${phi[${j#0}]}"/" jobs/$jobName".input"

			sed -i "10s/.*/LJ pot. intensity (epsilon): "${eps[${k#0}]}"/" jobs/$jobName".input"

			((ii++))
		done
	done
done

# Launch jobs.
# cd jobs
# parallel -j $threads < paralaunch.txt
