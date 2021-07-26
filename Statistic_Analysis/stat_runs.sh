#!/bin/bash

# First generate waveform data of n trials (No ampersand)
python3 Waveform_Generator.py --N_A --N_g --N_f --t0_tf --T --B --trials --seedn --N

for i in {1..10};
do
	python3 statudio.py --trialn --D --N_A --N_g --N_f --t0_tf --T --trials, --run1 --seedn --N_t &

# Wait for a empty directory to fill up before merging its jsons
python3 monitor.py --directory --dir_length

python3 json_stack_keys.py --json_path --merge_path_name &

python3 json_stack_keys.py --json_path --merge_path_name &

python3 json_stack_keys.py --json_path --merge_path_name &

python3 json_update_components.py --json_path --merge_path_name &

python3 json_list_append.py --json_path --merge_path_name &

# Wait for jsons to be merged before using them for the plots
python3 monitor.py --directory --dir_length

python3 test_plotter.py --T --N --var1 --var2 --stat &

python3 ROC_Curve.py --N &

python3 Scatter_plotter.py --thrshld --xvar --yvar --stat --plot &

python3 heatmap.py --var1 --var2 --T --stat --plot1 --plot2 &
