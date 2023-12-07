# Run HIBEAM.instr with given parameter set 
import math
import numpy as np
from run_hibeam import run_hibeam, run_backprop, optimal_source_positions, analyze_image, output_to_image_data, show_blade_x_scan, colors
from plot_vb import plot_results, show_instr, count_results
from generate_VB import generate_VB_point_focused, generate_VB_bladeFocused, write_VB
from focused_Venbla_xy_to_off import get_focused_blade_array

import csv

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Create Venetian Blinds geometry with given parameters")
	
	parser.add_argument("output_dir", help="Output directory for calculations and image figures")
	parser.add_argument('--vb_v_filename', metavar='vb_v_filename', help='if provided, use pre-generated vertically reflecting venetian blinds filename')
	parser.add_argument('--vb_h_filename', metavar='vb_h_filename', help='if provided, use pre-generated horizontally reflecting venetian blinds filename')
	parser.add_argument('--source_pos_interpolate', metavar='source_pos_file', help='if provided, interpolate blade source position from provided datapoints. *only works for a single input file, => only a single VB_pos')
	parser.add_argument('--detdim', nargs=2, type=float, metavar=('detwidth', 'detheight'),
					help='width and height of detector')
	parser.add_argument('--noShow', action='store_true', help='Show generated figures during calculations')

	# --------------------------------
	# Step 0: Extract parameters, define output files
	# --------------------------------
	args = parser.parse_args()

	output_dir = args.output_dir
	image_dir = output_dir

	noShow = False # Default to showing figures
	if args.noShow:
		noShow = True # Dont show figures, only generate them and store them in output directory

	# Step 0.b: Write output header
	summary_file = "{}/output.csv".format(output_dir)
	with open(summary_file, 'w', newline='') as file:
		writer = csv.writer(file)

		# write header
		header = ['VB_pos', 'VB_length', 'VB_m', 'VB_thickness', 'det_pos x', 'det_pos y', 'VB_filename v reflecting', 'VB_filename h reflecting'] # input variables
		header.extend(['target sum', 'target err', 'ratio', 'ratio_err'])	# tested outputs
		writer.writerow(header)

	# --------------------------------
	# Step 1: Define parameter set 
	# --------------------------------
	n = 1e6
	VB_filenames = [None, None]

	# Step 1.b: Create config_list
	# Baseline configuration
	det_pos = [-0.3, -0.1]
	baseline_config = [10, 0.5, 4, det_pos, 0.0005] # VB_pos, VB_length, VB_m, Det_pos, VB_thickness, VB_filenames
	print(baseline_config)

	# Parameter options
	VB_pos_vals = [8, 10, 15]		# VB pos is tied to geometry file though focusing
	VB_length_vals = [0.3, 0.5]		# VB length is tied to geometry file 
	VB_m_vals = [3, 4]
	VB_thickness_vals = [0.0005]	# VB thickness is tied to geometry file
	Det_pos_x_vals = [det_pos[0]]			# Det pos value is tied to geometry file through focusing
	Det_pos_y_vals = [det_pos[1]]

	# all combinations of configs above 
	grid_array1, grid_array2, grid_array3, grid_array4, grid_array5, grid_array6 = np.meshgrid(VB_pos_vals, VB_length_vals, VB_m_vals, VB_thickness_vals, Det_pos_x_vals, Det_pos_y_vals)
	config_list = np.vstack((grid_array1.ravel(), grid_array2.ravel(), grid_array3.ravel(), grid_array4.ravel(), grid_array5.ravel(), grid_array6.ravel())).T
	print(config_list)

	# --------------------------------
	# Step 2: Run configurations and save output 
	# --------------------------------

	# Step 2.a: Run baseline (without VB)
	no_VB_outDir = run_hibeam(n, baseline_config[0], baseline_config[1], baseline_config[2], baseline_config[3], ("-", "-"), output_dir, with_VB=False, noShow=noShow)

	vertical_image_data = output_to_image_data("{}/v_reflecting_VB_pos_image_midpoint.dat".format(no_VB_outDir)) 
	horizontal_image_data = output_to_image_data("{}/h_reflecting_VB_pos_image_midpoint.dat".format(no_VB_outDir)) 
	target_image_no_vb_data = output_to_image_data("{}/psdt2_large.dat".format(no_VB_outDir)) 

	# Plot results
	plot_results(vertical_image_data, plot_type='full', save_image=f'{image_dir}00_vertically_reflecting_blades_image.pdf', noShow=noShow)
	plot_results(horizontal_image_data, plot_type='full', save_image=f'{image_dir}01_horizontally_reflecting_blades_image.pdf', noShow=noShow)
	plot_results(target_image_no_vb_data, plot_type='full', save_image=f'{image_dir}02_target_image_no_vb.pdf', noShow=noShow)

	# Calculate baseline FoM
	no_vb_sum, no_vb_sum_err = count_results(target_image_no_vb_data, circle=[det_pos[0], det_pos[1], 0.20], save_image=f'{image_dir}03_target_no_vb.pdf', noShow=noShow)
	print(colors.GREEN + f'\nBaseline Calculation: {no_vb_sum} ± {no_vb_sum_err} nT^2/pulse\n' + colors.ENDC)
	print(colors.GREEN + f'\nRatio: {1.00} ± {0.00}\n' + colors.ENDC)

	# Write output to file
	line = ['no vb', 'no vb', 'no vb', 'no vb', det_pos[0], det_pos[1], 'no vb', 'no vb'] 
	line.extend([no_vb_sum, no_vb_sum_err, '1.00', '0.00'])
	with open(summary_file, 'a', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(line)

	# Step 2.b: Run configurations
	n = 1e6
	for i, config in enumerate(config_list):
		# Unpack config
		VB_pos, VB_length, VB_m, VB_thickness, det_pos_x, det_pos_y = config

		# Generate VB Given configuration
		VB_filenames = get_focused_blade_array(VB_pos, VB_length, VB_m, VB_thickness, (det_pos_x, det_pos_y), 55, output_dir, image_dir, noShow=noShow) 

		# Run simulation with config
		yes_VB_outDir = run_hibeam(n, VB_pos, VB_length, VB_m, (det_pos_x, det_pos_y), VB_filenames, output_dir, with_VB=True, noShow=noShow)
	
		# Capture image from output
		target_image_data = output_to_image_data("{}/psdt2_large.dat".format(yes_VB_outDir))
	
		# Show output images/ save to file  
		vb_sum, vb_sum_err = count_results(target_image_data, circle=[det_pos_x, det_pos_y, 0.20], save_image=f'{image_dir}{(4+i):02d}_target_with_vb.pdf', noShow=noShow)
		ratio = vb_sum/no_vb_sum
		ratio_err = ratio*np.sqrt(np.square(no_vb_sum_err/no_vb_sum) + np.square(vb_sum_err/vb_sum))
		print(colors.GREEN + f'\nEstimated improvement: {ratio} ± {ratio_err}\n' + colors.ENDC)

		# Write output to file
		line = [VB_pos, VB_length, VB_m, VB_thickness, det_pos_x, det_pos_y, VB_filenames[0], VB_filenames[1]]
		line.extend([vb_sum, vb_sum_err, ratio, ratio_err])
		with open(summary_file, 'a', newline='') as file:
			writer = csv.writer(file)
			writer.writerow(line)
	
#	# --------------------------------
#	# Step 3: Write output to file 
#	# --------------------------------
#
#	# Step 3.a: Write summary header
#	summary_file = "{}/output.csv".format(output_dir)
#	with open(summary_file, 'w', newline='') as file:
#		writer = csv.writer(file)
#
#		# write header
#		header = ['VB_pos', 'VB_length', 'VB_m', 'det_pos x', 'det_pos y', 'VB_filename v reflecting', 'VB_filename h reflecting'] # input variables
#		header.extend(['ratio', 'ratio_err'])	# tested outputs
#		writer.writerow(header)
#
#	# Step 3.b: Save output for given parameter combinations 
#	for config in config_list:
#		VB_pos, VB_length, VB_m, VB_thickness, det_pos_x, det_pos_y, VB_vr_filename , VB_hr_filename, ratio, ratio_err = config
#
#		# write output to file
#		line = [VB_pos, VB_length, VB_m, det_pos_x, det_pos_y, VB_vr_filename, VB_hr_filename]
#		line.extend([ratio, ratio_err])
#		with open(summary_file, 'a', newline='') as file:
#			writer = csv.writer(file)
#			writer.writerow(line)
