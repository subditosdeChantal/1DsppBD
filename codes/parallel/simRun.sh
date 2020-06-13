#!/bin/bash

# Launch jobs.
cd jobs
parallel -j procfile < paralaunch.txt

# Run analysis
cd ../analysis
bash run_analysis.sh sim*
cd ..