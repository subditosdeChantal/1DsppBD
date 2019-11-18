#!/bin/bash

echo "Plotting cluster size distributions..."
gnuplot plot_size_distr -c 'sim*'
gnuplot plot_size_distr_norm -c 'sim*'

echo "Plotting direction distributions..."
gnuplot plot_dir_distr -c 'sim*'

echo "Plotting total number of clusters evolution..."
gnuplot plot_nclusters -c 'sim*'

echo "Plotting correlations..."
python3 correlations.py 'sim*'
gnuplot plot_correlations -c 'sim*'

echo "Plotting system snapshots..."
gnuplot plot_snapshot_beg -c 'sim*'
gnuplot plot_snapshot_ss -c 'sim*'

echo "Done!"
