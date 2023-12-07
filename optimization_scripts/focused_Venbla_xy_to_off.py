# Generate focused VB design given .instr file
import math
import os
import numpy as np
from run_hibeam import run_hibeam, run_backprop, optimal_source_positions, analyze_image, output_to_image_data, show_blade_x_scan, source_interpolation, colors
from plot_vb import plot_results, show_instr, count_results
from generate_VB import generate_VB_point_focused, generate_VB_bladeFocused, write_VB

def get_focused_blade_array(VB_pos, VB_length, VB_m, VB_thickness, det_pos, zdet, output_dir, image_dir, noShow, source_pos_interpolate=None):
	# --------------------------------
	# Step 1: Run .instr file with given parameters and extract outputs
	# --------------------------------
	n = 1e6
	VB_filenames = [None, None]
	no_VB_outDir = run_hibeam(n, VB_pos, VB_length, VB_m, det_pos, VB_filenames, output_dir, with_VB=False, dead_monolith=False, noShow=noShow) # Specify reflecting monolith to avoid touching monolith focused
	vertical_image_data = output_to_image_data("{}/v_reflecting_VB_pos_image_midpoint.dat".format(no_VB_outDir)) 
	horizontal_image_data = output_to_image_data("{}/h_reflecting_VB_pos_image_midpoint.dat".format(no_VB_outDir)) 
	target_image_no_vb_data = output_to_image_data("{}/psdt2_large.dat".format(no_VB_outDir)) 

	# --------------------------------
	# Step 2: Find x and y bounds from simulation output
	# --------------------------------
	
	## inner x bounds vertical reflectors: (vx0,vx3)
	## outer x bounds vertical reflectors: (vx1,vx2)
	## inner y bounds vertical reflectors: (vy0,vy3)
	## outer y bounds vertical reflectors: (vy1,vy2)

	## inner x bounds horizontal reflectors: (hx0,hx3)
	## outer x bounds horizontal reflectors: (hx1,hx2)
	## inner y bounds horizontal reflectors: (hy0,hy3)
	## outer y bounds horizontal reflectors: (hy1,hy2)

	# Step 2.a: Find outer bounds by analyzing VB_pos image_data
	vx0, vx3, vy0, vy3 = analyze_image(vertical_image_data)
	hx0, hx3, hy0, hy3 = analyze_image(horizontal_image_data)
		
	# Step 2.b: Find inner bounds by analyzing backpropagated VB_pos image_data (from filtered target)
	# Backpropagate from target
	back_v = run_backprop(no_VB_outDir, VB_pos)
	back_h = run_backprop(no_VB_outDir, VB_pos+VB_length+0.01)

	# Analyze backpropagated image_data to get inner bounds
	backprop_vertical_image_data = output_to_image_data(back_v)
	backprop_horizontal_image_data = output_to_image_data(back_h)

	# Find inner bounds by analyzing backpropagated images
	vx1, vx2, vy1, vy2 = analyze_image(backprop_vertical_image_data)
	hx1, hx2, hy1, hy2 = analyze_image(backprop_horizontal_image_data)

	# Sanity check on target
	back_target = run_backprop(no_VB_outDir, zdet-0.001)
	target_image_data = output_to_image_data(back_target)
	plot_results(target_image_data, plot_type='full', save_image=f'{image_dir}01_filtered_target_image.pdf', noShow=noShow)

	# Sanity check on bounds found
	print('vertical_image (vx0, vx3, vy0, vy3)')
	print(vx0, vx3, vy0, vy3)
	print('horizontal_image (hx0, hx3, hy0, hy3)')
	print(hx0, hx3, hy0, hy3)
	print('backprop_vertical (vx1, vx2, vy1, vy2)')
	print(vx1, vx2, vy1, vy2)
	print('backprop_horizontal (hx1, hx2, hy1, hy2)')
	print(hx1, hx2, hy1, hy2)
	vx = np.array([vx0, vx1, vx2, vx3])
	vy = np.array([vy0, vy1, vy2, vy3])
	hx = np.array([hx0, hx1, hx2, hx3])
	hy = np.array([hy0, hy1, hy2, hy3])
	plot_results(backprop_vertical_image_data, vertical_image_data, plot_type='y', xlims=vx, ylims=vy, save_image=f'{image_dir}02_vr_vb_yproj.pdf', noShow=noShow)
	plot_results(backprop_horizontal_image_data, horizontal_image_data, plot_type='x', xlims=hx, ylims=hy, save_image=f'{image_dir}03_hr_vb_xproj.pdf', noShow=noShow)

	# --------------------------------
	# Step 3: Generate focused Venetian Blinds within found bounds with focusing geometry
	# --------------------------------

	# Step 3.a: Find optimal focus point for y=vy, x=hx blades using backpropagation 
	#show_blade_x_scan(VB_pos, hx, hy, VB_length, no_VB_outDir, image_dir, dx=0.2)
	#vyz0, hxz0 = optimal_source_positions(VB_pos, vx, vy, hx, hy, VB_length, no_VB_outDir, image_dir, noShow) # returns [vy0z0, vy1z0, vy2z0, vy3z0], [hx0z0, hx1z0, hx2z0, hx3z0]
	#print(f'vy: {vy}\nvyz0: {vyz0}')
	#print(f'hx: {hx}\nhxz0: {hxz0}')

	# Step 3.b: Generate Venetian Blinds *assuming linear interpolation of optimal focus point 
	if source_pos_interpolate is None:
		num_gridpoints = 20
		source_pos_file = f"{output_dir}hxz0_{VB_pos}_{VB_length}_{hx[0]}_{hx[3]}_{hy[0]}_{hy[3]}_{num_gridpoints}.dat"
		if os.path.exists(source_pos_file):
			print('source pos file exists!')
			print(source_pos_file)
		else:
			# Generate source pos file for hxz0 interpolation
			source_interpolation(VB_pos, hx, hy, VB_length, num_gridpoints, no_VB_outDir, source_pos_file)
	else:
		source_pos_file = source_pos_interpolate

	#VB_filenames = generate_VB_focused_blades(VB_pos, zdet, det_pos, VB_length, VB_thickness, vx, vy, hx, hy, vyz0, hxz0) # VB with focus linearly interpolated between top and bottom of each blade array 
	#VB_filenames = generate_VB_point_focused(VB_pos, zdet, det_pos, VB_length, VB_thickness, vy, hx) # Venetian Blinds with pointlike source 
	#VB_filenames = ['Venbla_vertically_reflecting_geometry.off', 'Venbla_horizontally_reflecting_geometry.off']
	print(source_pos_file)
	VB_filenames = generate_VB_bladeFocused(VB_pos, zdet, det_pos, VB_length, VB_thickness, vx, vy, hx, hy, source_pos_file) # VB with focus interpolated from file with many source positions 

	return VB_filenames

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Create Venetian Blinds geometry with given parameters")
	
	parser.add_argument("VB_z", type=float, help="Z distance between source and center of vb")
	parser.add_argument("VB_length", type=float, help="Z length of blades")
	parser.add_argument("VB_thickness", type=float, help="Thickness of blades")
	parser.add_argument("VB_m", type=float, default=4, help="VB reflectivity m value, default=4")
	parser.add_argument("det_z", type=float, help="Z distance between center of vb and target")
	parser.add_argument("output_dir", help="Output directory for calculations and image figures")
	parser.add_argument('--source_pos_interpolate', metavar='source_pos_file', help='if provided, interpolate blade source position from provided datapoints.')
	parser.add_argument('--detdim', nargs=2, type=float, metavar=('detwidth', 'detheight'),
					help='width and height of detector')
	parser.add_argument('--detpos', nargs=2, type=float, metavar=('detx', 'dety'),
					help='displacement of detector position')
	parser.add_argument('--noShow', action='store_true', help='Show generated figures during calculations')

	# --------------------------------
	# Step 0: Extract parameters
	# --------------------------------
	args = parser.parse_args()

	output_dir = args.output_dir
	image_dir = output_dir

	VB_length = args.VB_length
	VB_thickness = args.VB_thickness
	VB_pos = args.VB_z
	VB_m = args.VB_m
	zdet = args.det_z

	det_dims = args.detdim
	if args.detpos:
		det_pos = args.detpos
	else:
		det_pos = [0, 0]

	noShow = False # Default to showing figures
	if args.noShow:
		noShow = True # Dont show figures, only generate them and store them in output directory

	# --------------------------------
	# Step 1: Generate VB Configuration 
	# --------------------------------
	VB_filenames = get_focused_blade_array(VB_pos, VB_length, VB_m, VB_thickness, det_pos, zdet, output_dir, image_dir, noShow, source_pos_interpolate=args.source_pos_interpolate)

	# --------------------------------
	# Step 2: Run with generated VB and show output 
	# --------------------------------
	baseline_config = [10, 0.5, 4, (-0.3, -0.1), 0.0005] # VB_pos, VB_length, VB_m, Det_pos, VB_thickness, VB_filenames
	n = 1e6

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

	# Step 2.b: Run with focused VB
	n = 1e6
	yes_VB_outDir = run_hibeam(n, VB_pos, VB_length, VB_m, det_pos, VB_filenames, output_dir, with_VB=True, noShow=noShow)

	# Capture image from output
	target_image_data = output_to_image_data("{}/psdt2_large.dat".format(yes_VB_outDir))

	# Step 2.b: Show output images/ save to file  
	n = 1e6
	VB_filenames = [None, None]

	vb_sum, vb_sum_err = count_results(target_image_data, circle=[det_pos[0], det_pos[1], 0.20], save_image=f'{image_dir}05_target_with_vb.pdf', noShow=noShow)
	ratio = vb_sum/no_vb_sum
	ratio_err = ratio*np.sqrt(np.square(no_vb_sum_err/no_vb_sum) + np.square(vb_sum_err/vb_sum))
	print(f'Estimated improvement: {ratio} ± {ratio_err}')

