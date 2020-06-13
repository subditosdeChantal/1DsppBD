#!/bin/bash

scriptName="BDp_051"
sweepName="init"
#outputDir="output"
execName=$scriptName"_"$sweepName

alpha=(0.1 0.5)
phi=(0.1 0.8)
eps=(0 24)
dt=(0.01)
# tsim=(100)
# tmeas=(1000 10000 100000)
v=(1)
beta=(10)

threads=8

rm -rf jobs

mkdir jobs
# mkdir jobs/$outputDir

echo $threads > jobs/procfile

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

						echo "./"$execName" < "$jobName".input > "$jobName".out" >> jobs/paralaunch.txt

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

chmod u+x simCont.sh