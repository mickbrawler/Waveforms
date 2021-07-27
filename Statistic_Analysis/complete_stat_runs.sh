#!/bin/bash

# Clear the (per trialn) jsons ahead of time

rm -rf Max_BG_TEMP_folder/*
rm -rf output_folder/*
rm -rf Peaks_folder/*
rm -rf Max_OS_folder/*
rm -rf thresholds_folder/*

rm -rf Merged_jsons/*

rm -rf plots/*

echo "correct func"

# First generate waveform data of n trials (No ampersand)
python3 Waveform_Generator.py --N_A 2 --N_g 2 --N_f 2 --t0_tf 4 --T 10 --B 0 --trials 4 --seedn 1 --N 250 --inputfile "input"

count=0
for i in {0..3};
do
	if [ $count -eq 0 ]
	then
		count=$((count+1))
		python3 statudio.py --trialn $i --D .02 --N_A 2 --N_g 2 --N_f 2 --t0_tf 4 --T 10 --trials 4 --run1 True --seedn 1 --N_t 250 --inputfile "input" &
	else
		python3 statudio.py --trialn $i --D .02 --N_A 2 --N_g 2 --N_f 2 --t0_tf 4 --T 10 --trials 4 --run1 False --seedn 1 --N_t 250 --inputfile "input" &
	fi
done

echo "here"

# Wait for a empty directory to fill up before merging its jsons
python3 monitor.py --directory "Max_BG_TEMP_folder/" --dir_length 4

echo "past monitor 1"

python3 json_stack_keys.py --jsons_path "output_folder/" --merge_path_name "Merged_jsons/Merged_output" &
python3 json_stack_keys.py --jsons_path "Peaks_folder/" --merge_path_name "Merged_jsons/Merged_Peaks" &
python3 json_stack_keys.py --jsons_path "Max_OS_folder/" --merge_path_name "Merged_jsons/Merged_Max_OS" &

python3 json_list_append.py --jsons_path "thresholds_folder/" --merge_path_name "Merged_jsons/Merged_thresholds" &


echo "at monitor 2..."

# Checks if all previous merges created
python3 monitor.py --directory "Merged_jsons/" --dir_length 4

python3 json_update_components.py --jsons_path "Max_BG_TEMP_folder/" --merge_path_name "Merged_jsons/Merged_Max_BG_TEMP" &

echo "at monitor 3..."

# Wait for jsons to be merged before using them for the plots
python3 monitor.py --directory "Merged_jsons/" --dir_length 5

echo "past monitor 3!!"

python3 test_plotter.py --T 100 --N 300 --var1 0 --var2 0 --stat 2 --bg_test True --plot1 "test_plot1" --plot2 "test_plot2" & 

python3 ROC_Curve.py --N 100 --outputfile "ROC_test" &

python3 Scatter_plotter.py --thrshld 4 --xvar 2 --yvar 1 --stat 2 --plot "Scatter_Avg_plot" &
python3 Scatter_plotter.py --thrshld 4 --xvar 2 --yvar 0 --stat 2 --plot "Scatter_AvF_plot" &
python3 Scatter_plotter.py --thrshld 3 --xvar 1 --yvar 0 --stat 2 --plot "Scatter_FvG_plot" &

python3 heatmap.py --var1 2 --var2 1 --T 0 --stat 2 --plot1 "heat_AvG_plot1" --plot2 "heat_AvG_plot2" &
python3 heatmap.py --var1 2 --var2 0 --T 0 --stat 2 --plot1 "heat_AvF_plot1" --plot2 "heat_AvF_plot2" &
python3 heatmap.py --var1 1 --var2 0 --T 0 --stat 2 --plot1 "heat_FvG_plot1" --plot2 "heat_FvG_plot2" &
