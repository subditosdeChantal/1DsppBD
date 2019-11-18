#!/bin/bash

echo "Plotting cluster size distributions..."
gnuplot plot_size_distr
gnuplot plot_size_distr_norm

echo "Plotting direction distributions..."
gnuplot plot_dir_distr

echo "Plotting total number of clusters evolution..."
gnuplot plot_nclusters

echo "Plotting correlations..."
python3 correlations.py
gnuplot plot_correlations

echo "Plotting system snapshots..."
gnuplot plot_snapshot_beg
gnuplot plot_snapshot_ss

echo "Done!"
