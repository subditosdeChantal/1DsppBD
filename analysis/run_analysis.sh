#!/bin/bash

echo ""
echo "Plotting cluster size distributions..."
gnuplot -c plot_size_distr $1
gnuplot -c plot_size_distr_norm $1

echo ""
echo "Plotting direction distributions..."
gnuplot -c plot_dir_distr $1

echo ""
echo "Plotting total number of clusters evolution..."
gnuplot -c plot_nclusters $1

echo ""
echo "Computing correlations..."
python3 correlations.py $1

echo ""
echo "Plotting correlations..."
gnuplot -c plot_correlations $1

echo ""
echo "Plotting system snapshots..."
gnuplot -c plot_snapshot_beg $1
gnuplot -c plot_snapshot_ss $1

echo ""
echo "Done!"
