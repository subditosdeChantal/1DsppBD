#!/bin/bash

echo "Launching simulation loop..."
cd codes/
./BD
echo "Simulation's over"
cd ../analysis/
echo ""
echo "Launching analysis..."
echo ""
bash run_analysis.sh sim*mu*
cd

