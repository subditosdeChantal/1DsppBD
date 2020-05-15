#!/bin/bash

scriptName="BDp_00"
sweepName="sweepdt"
#outputDir="output"
execName=$scriptName"_"$sweepName

alpha=(0.05)
phi=(0.2)
eps=(0 100)
dt=(1 0.1 0.01)
# tsim=(100)
# tmeas=(1000 10000 100000)
v=(0 1)

threads=12

rm -rf jobs

mkdir jobs
# mkdir jobs/$outputDir

gcc -o jobs/$execName $scriptName".c" -lm

idx=($(seq -w 0 $(( ${#alpha[@]}*${#phi[@]}*${#eps[@]}*${#dt[@]}*${#v[@]} - 1 )) ))
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

					jobName=$sweepName"_"${idx[${ii}]}

					cp "input.dat" jobs/$jobName".input"

					echo "./"$execName" < "$jobName".input" >> jobs/paralaunch.txt

					sed -i "1s/.*/Simulation name: "$sweepName"/" jobs/$jobName".input"

# 					sed -i "2s/.*/Output dir: "$outputDir"/" jobs/$jobName".input"

					sed -i "5s/.*/Tumbling rate (alpha): "${alpha[${i#0}]}"/" jobs/$jobName".input"

					sed -i "6s/.*/Density (phi): "${phi[${j#0}]}"/" jobs/$jobName".input"

#                    			sed -i "7s/.*/Simulation time (time units): "${tsim[${l#0}]}"/" jobs/$jobName".input"

					sed -i "10s/.*/LJ pot. intensity (epsilon): "${eps[${k#0}]}"/" jobs/$jobName".input"

			                sed -i "12s/.*/Propulsion force (v, Fp): "${v[${m#0}]}"/" jobs/$jobName".input"

					sed -i "18s/.*/Time Step: "${dt[${l#0}]}"/" jobs/$jobName".input"

					((ii++))
				done

			done
		done
	done
done

# Launch jobs.
cd jobs
parallel -j $threads < paralaunch.txt
