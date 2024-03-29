# run hibeam instrument with given parameters and analyze images for greatest slope
import hashlib
import csv
import os
import re
import numpy as np
import mcstasHelper as mc
from plot_vb import plot_results

class colors: 
	RED = '\033[31m'
	ENDC = '\033[m'
	GREEN = '\033[32m'
	YELLOW = '\033[33m'
	BLUE = '\033[34m'

# for unique naming of output files
def generate_hash(*parameters):
	combined_parameters = ''.join(str(param) for param in parameters)
	hash_object = hashlib.sha256(combined_parameters.encode())
	hash_value = hash_object.hexdigest()
	return hash_value

# run simulation with params
def run_hibeam(n, VB_pos, VB_length, VB_m, det_pos, VB_filenames, output_dir, with_VB, dead_monolith=False):
	instr = "hibeam/hibeam_VB.instr"
	instr_no_vb = "hibeam/hibeam_noSecondFocusing.instr"
	instr_no_vb_dead_monolith = "hibeam/hibeam_noSecondFocusing_deadMonolith.instr"
	mpi = 1

	det_x, det_y = det_pos
	v_reflecting_VB_geometry, h_reflecting_VB_geometry = VB_filenames

	# if with_VB = True, then run with generated VB
	if with_VB:
		# const params: VB_pos, VB_length, VB_m 
		print(colors.BLUE + f'==== Running with VB ====\n==== VB_pos:{VB_pos}, VB_length:{VB_length}, VB_m:{VB_m}, det_pos:({det_x}, {det_y}), VB_filenames:{VB_filenames} ====' + colors.ENDC)
		
		hash_value = generate_hash(n, VB_pos, VB_length, VB_m, det_x, det_y, VB_filenames)
		
		filename = "{}run_{}".format(output_dir, hash_value)
		command = "mcrun --mpi={} {} -d {} -s 1 -n {} VB_pos={} VB_length={} VB_m={} target_x={} target_y={} V_r_VB_filename={} H_r_VB_filename={}".format(mpi, instr, filename, n, VB_pos, VB_length, VB_m, det_x, det_y, v_reflecting_VB_geometry, h_reflecting_VB_geometry)
		if dead_monolith:
			print('with VB, dead monolith Option currently not supported')

	# if with_VB = False, then run without VB
	else:
		# const params: VB_pos, VB_length, VB_m 
		print(colors.BLUE + f'==== Running with no VB ====\n==== VB_pos:{VB_pos}, VB_length:{VB_length}, VB_m:{VB_m}, det_pos:({det_x}, {det_y}) ====' + colors.ENDC)
		
		hash_value = generate_hash(n, VB_pos, VB_length, det_x, det_y)
		
		filename = "{}run_{}".format(output_dir, hash_value)
		command = "mcrun --mpi={} {} -d {} -s 1 -n {} VB_pos={} VB_length={} VB_m={} target_x={} target_y={}".format(mpi, instr_no_vb, filename, n, VB_pos, VB_length, VB_m, det_x, det_y)
		if dead_monolith:
			command = "mcrun --mpi={} {} -d {} -s 1 -n {} VB_pos={} VB_length={} VB_m={} target_x={} target_y={}".format(mpi, instr_no_vb_dead_monolith, filename, n, VB_pos, VB_length, VB_m, det_x, det_y)

	# Redirect output to /dev/null (Unix-based) or NUL (Windows)
	print('redirecting output')
	if os.name == "posix":  # Unix-based system (e.g., Linux, macOS)
		command += " > /dev/null 2>&1"
	elif os.name == "nt":   # Windows system
		command += " > NUL 2>&1"
	print(colors.YELLOW + "\nrunning command:\n" + command + colors.ENDC + '\n')
	os.system(command)

	return filename

# run backpropagation cpp file on filtered target image
def run_backprop(dirName, zvb):
	outFile = "{}/histogram_{}.dat".format(dirName, zvb)  

	# Generate mcpl text file from filtered target image for backpropagation
	mcpltool_command = "mcpltool -t {}/target_output.mcpl.gz {}/target_output.dat".format(dirName, dirName)
	#mcpltool_command = "mcpltool -t {}/target_output.mcpl {}/target_output.dat".format(dirName, dirName)
	print(colors.YELLOW + "\nrunning command:\n{}".format(mcpltool_command) + colors.ENDC + '\n')
	os.system(mcpltool_command)

	# Run back_propagate.cpp with given zvb
	backprop_command = "./optimization_scripts/backprop/back_propagate {} {}/target_output.dat {}".format(zvb, dirName, outFile)
	print(colors.YELLOW + "\nrunning command:\n{}".format(backprop_command) + colors.ENDC + '\n')
	os.system(backprop_command)

	# Remove 'target_output.dat' to reduce storage requirements
	os.remove("{}/target_output.dat".format(dirName))

	return outFile

# Backpropagate from blade position and return optimal source position for that point
def optimal_source_positions(zvb, vx, vy, hx, hy, length, outDir, imageDir, noShow):
	# For given blade, find optimal source position
	def blade_optimal_pos(pos, direction): # written with vertical convention
		if direction == 'vertical':
			y = pos

			# Find y extent for blade selection
			ymin = (zvb * y)/ (zvb + length/2) 
			ymax = (zvb * y)/ (zvb - length/2) 
			if ymax < ymin: # make sure ymin < ymax
				tmp = ymax
				ymax = ymin
				ymin = tmp
	
			## In back_propagate_selected.cpp, choose to weight neutrons by ToF^2 or not

			# Generate image of what is being backpropagated
			v_r_VB_blade_file = f'{outDir}/v_r_VB_histogram_{y}_{zvb-0.001}.dat'
			os.system(f"./optimization_scripts/backprop/back_propagate_selected {zvb} {zvb-0.001} {outDir}/v_r_VB_output.dat {v_r_VB_blade_file} --rectangle {vx[0]} {ymin} {vx[3]} {ymax}")
			# Plot what is being backpropagated
			v_r_VB_blade_image_data = output_to_image_data(v_r_VB_blade_file)
			plot_results(v_r_VB_blade_image_data, plot_type='full', save_image=f'{imageDir}05_v_r_VB_blade_{i}_image.pdf', noShow=noShow)
	
			# Backpropagate from selected blade
			v_r_VB_back_file = f'{outDir}/v_r_VB_histogram_{y}_0.dat'
			backprop_command = f"./optimization_scripts/backprop/back_propagate_selected {zvb} 0.01 {outDir}/v_r_VB_output.dat {v_r_VB_back_file} --rectangle {vx[0]} {ymin} {vx[3]} {ymax}"
			print(colors.YELLOW + "\nrunning command:\n{}".format(backprop_command) + colors.ENDC + '\n')
			os.system(backprop_command)
			# Plot result of backpropagation
			v_r_VB_back_image_data = output_to_image_data(v_r_VB_back_file)
			plot_results(v_r_VB_back_image_data, plot_type='full', save_image=f'{imageDir}06_v_r_VB_back_{i}_image.pdf', noShow=noShow)
			plot_results(v_r_VB_back_image_data, plot_type='y', save_image=f'{imageDir}07_v_r_VB_back_y_{i}_image.pdf', noShow=noShow)
	
			# Find weighted average of backpropagated distribution, store in vyz0
			ysrc = find_image_weighted_average(v_r_VB_back_image_data, 'y')
	
			return ysrc
		elif direction == 'horizontal':
			x = pos

			# Find y extent for blade selection
			xmin = (zvb * x)/ (zvb + length/2) 
			xmax = (zvb * x)/ (zvb - length/2) 
			if xmax < xmin: # make sure xmin < xmax
				tmp = xmax
				xmax = xmin
				xmin = tmp
	
			# Generate image of what is being backpropagated
			h_r_VB_blade_file = f'{outDir}/h_r_VB_histogram_{x}_{zvb-0.001}.dat'
			os.system(f"./optimization_scripts/backprop/back_propagate_selected {zvb} {zvb-0.001} {outDir}/h_r_VB_output.dat {h_r_VB_blade_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}")
			# Plot what is being backpropagated
			h_r_VB_blade_image_data = output_to_image_data(h_r_VB_blade_file)
			plot_results(h_r_VB_blade_image_data, plot_type='full', save_image=f'{imageDir}08_h_r_VB_blade_{i}_image.pdf', noShow=noShow)
	
			# Backpropagate from selected blade
			h_r_VB_back_file = f'{outDir}/h_r_VB_histogram_{x}_0.dat'
			backprop_command = f"./optimization_scripts/backprop/back_propagate_selected {zvb} 0.01 {outDir}/h_r_VB_output.dat {h_r_VB_back_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}"
			print(colors.YELLOW + "\nrunning command:\n{}".format(backprop_command) + colors.ENDC + '\n')
			os.system(backprop_command)
			# Plot result of backpropagation
			h_r_VB_back_image_data = output_to_image_data(h_r_VB_back_file)
			plot_results(h_r_VB_back_image_data, plot_type='full', save_image=f'{imageDir}09_h_r_VB_back_{i}_image.pdf', noShow=noShow)
			plot_results(h_r_VB_back_image_data, plot_type='x', save_image=f'{imageDir}10_h_r_VB_back_x_{i}_image.pdf', noShow=noShow)

			# Find weighted average of backpropagated distribution, store in hxz0
			xsrc = find_image_weighted_avg(h_r_VB_back_image_data, 'x')
			
			return xsrc
		else:
			print('incorrect direction specified. specify "vertical" or "horizontal"')
		
	# Non major extent vx, hy used for determining bounds on backprop 
	
	# --------------------------------
	# find vyz0

	# Convert .mcpl to .txt
	mcpltool_command = f'mcpltool -t {outDir}/v_r_VB_output.mcpl.gz {outDir}/v_r_VB_output.dat'
	#mcpltool_command = f'mcpltool -t {outDir}/v_r_VB_output.mcpl {outDir}/v_r_VB_output.dat'
	print(colors.YELLOW + f"\nrunning command:\n{mcpltool_command}" + colors.ENDC + '\n')
	os.system(mcpltool_command)

	# For each blade position, backpropagate only neutrons at that position 
	vyz0 = []
	for y in vy: 
		yz0 = blade_optimal_pos(y, direction='vertical')
		vyz0.append(yz0)
	vyz0 = np.array(vyz0)

	# Remove 'v_r_VB_output.dat' to reduce storage requirements
	os.remove(f"{outDir}/v_r_VB_output.dat")

	# --------------------------------
	# find hxz0

	# Convert .mcpl to .txt
	mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl.gz {outDir}/h_r_VB_output.dat"
	#mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl {outDir}/h_r_VB_output.dat"
	print(colors.YELLOW + f"\nrunning command:\n{mcpltool_command}" + colors.ENDC + '\n')
	os.system(mcpltool_command)

	# For each blade position, backpropagate only neutrons at that position 
	hxz0 = []
	for i, x in enumerate(hx): 
		xz0 = blade_optimal_pos(x, direction='horizontal')
		hxz0.append(xz0)
	hxz0 = np.array(hxz0)

	# Remove 'h_r_VB_output.dat' to reduce storage requirements
	os.remove(f"{outDir}/h_r_VB_output.dat")

	return vyz0, hxz0 

# Shows how source position varies depending on blade x position
def show_blade_x_scan(zvb, hx, hy, length, outDir, histDir, dx=0.2):	# dx [cm]
	# Convert .mcpl to .txt
	mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl.gz {outDir}/h_r_VB_output.dat"
	#mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl {outDir}/h_r_VB_output.dat"
	print(colors.YELLOW + f"\nrunning command:\n{mcpltool_command}" + colors.ENDC + '\n')
	os.system(mcpltool_command)

	# -----------
	# Show image of full backpropagation
	print('Generating full horizontal blade scan:')

	# Perform imaging for first set of blades:
	for i, x in enumerate(np.concatenate((np.arange(hx[0], hx[1], dx), np.arange(hx[2], hx[3], dx)))):
		# Generate image of what is being backpropagated
		blade_file = f'{histDir}/histogram_x_scan_blade_{i:03d}.dat'

		# Find y extent for blade selection
		xmin = (zvb * x)/ (zvb + length/2) 
		xmax = (zvb * x)/ (zvb - length/2) 
		if xmax < xmin: # make sure xmin < xmax
			tmp = xmax
			xmax = xmin
			xmin = tmp

		# Plot what is being backpropagated
		os.system(f"./optimization_scripts/backprop/back_propagate_selected {zvb} {zvb-0.001} {outDir}/h_r_VB_output.dat {blade_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}")
		VB_blade_image_data = output_to_image_data(blade_file)

		# Backpropagate from selected blade
		blade_back_file = f'{histDir}/histogram_x_scan_{i:03d}.dat'
		backprop_command = f"./optimization_scripts/backprop/back_propagate_selected {zvb} 0.001 {outDir}/h_r_VB_output.dat {blade_back_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}"
		print(colors.YELLOW + "\nrunning command:\n{}".format(backprop_command) + colors.ENDC + '\n')
		os.system(backprop_command)

		# For each file, write a line to the top of the output histogram describing the blade xpos:
		with open(blade_back_file, 'r') as file:
			original_content = file.read()
		with open(blade_back_file, 'w') as file:
			file.write(f'# blade_xpos: {x}\n')
			file.write(original_content)

		# These backpropagation images are to be analyzed separately by display_sequence.py, see `optimization_scripts/README.md`
	# -----------
	# Remove 'h_r_VB_output.dat' to reduce storage requirements
	os.remove(f"{outDir}/h_r_VB_output.dat")

# analyze McStas intensity distributions from file 
# return extent, image np array, image err np array as tuple 
def output_to_image_data(inFile):
	I, sigI, N, dataHeader, L = mc.extractMcStasData(inFile)
	# From image, extract numpy array 
	if dataHeader['type'][:8]=="array_2d":
		extent = np.array(dataHeader['xylimits'].split(),dtype=float)
	
		# Call the mcstas2TIFF function to generate .tif files
		imN, imI, sigI = mc.mcstas2TIFF(inFile, save=False)
	else:
		print('inFile incorrect type')
	
	return dataHeader, extent, np.array(imI), np.array(sigI)

# analyze McStas intensity distribution images
# return value corresponding to weighted average of image data summed along plot axis 
def find_image_weighted_average(image_data, plot_type):
	dataHeader, extent, image, imageErr = image_data

	# Find the value of x corresponding to maximum y
	def find_max(x, y):
		max_index = np.argmax(y)
		xmax = x[max_index]
		print(xmax)
		return xmax

	# Find the average of x weighted by y 
	def find_weighted_average(x, y):
		weighted_sum_x = np.sum(x * y)
		sum_y = np.sum(y)
		xavg = weighted_sum_x / sum_y
		
		print(xavg)
		return xavg

	# --------------
	# Find x weighted average:
	if plot_type == 'x':
		dx = (extent[1] - extent[0]) / image.shape[1]
	
		# Generate X cross-section 
		cross_section = np.sum(image, axis=0)  # Sum along the horizontal axis
		err_cross_section = np.sqrt(np.sum(np.square(imageErr), axis=0))
	
		x = np.linspace(extent[0], extent[1], np.size(cross_section))
	
		xmax = find_weighted_average(x, cross_section)
		return xmax

	# --------------
	# Find y weighted average: 
	elif plot_type == 'y':
		dy = (extent[3] - extent[2]) / image.shape[0]
	
		# Generate Y cross-section 
		cross_section = np.sum(image, axis=1)  # Sum along the vertical axis
		err_cross_section = np.sqrt(np.sum(np.square(imageErr), axis=1))
	
		y = np.linspace(extent[3], extent[2], np.size(cross_section))
	
		ymax = find_weighted_average(y, cross_section)
		return ymax
		
	else:
		print('invalid plot type, use plot_type="x" or "y"')
	

# analyze McStas intensity distribution images
# return x and y bounds of max slope
def analyze_image(image_data):
	dataHeader, extent, image, imageErr = image_data

	# Find FWHM values
	def calculate_fwhm(x, y): 
		max_index = np.argmax(y)
		x_max = x[max_index]
		y_max = y[max_index]
	
		half_max_y = y_max/2
		closest_index = np.argmin(np.abs(y - half_max_y))
	
		# Find left_idx where x is less than x_max
		filtered_y = y[:max_index]
		left_idx = np.argmin(np.abs(filtered_y - half_max_y))
		x_left = x[left_idx]
		#print(f'left: {x_left}')
		
		# Find right_idx where x is greater than x_max
		filtered_y = y[max_index:]
		right_idx = max_index + np.argmin(np.abs(filtered_y - half_max_y))
		x_right = x[right_idx]
		#print(f'right: {x_right}')
	
		fwhm = x_left - x_right
		return x_right, x_left
	
	# Find extrema that are at least t_percent% of the maximum value
	t_percent = 1
	def find_extrema(x, y): 
		threshold = 0.01 * t_percent * max(y)  # 1% of the maximum y value
		valid_indices = np.where(y >= threshold)  # Get indices where y values are at least 1% of the maximum

		if len(valid_indices) == 0:
			return None, None  # No valid indices found

		x_valid = x[valid_indices]
		y_valid = y[valid_indices]

		x_min = np.min(x_valid)
		x_max = np.max(x_valid)

		return x_min, x_max

	# --------------
	# Find x bounds:
	dx = (extent[1] - extent[0]) / image.shape[1]

	# Generate X cross-section 
	cross_section = np.sum(image, axis=0)  # Sum along the horizontal axis
	err_cross_section = np.sqrt(np.sum(np.square(imageErr), axis=0))

	x = np.linspace(extent[0], extent[1], np.size(cross_section))

	xmin, xmax = find_extrema(x, cross_section)
	#xmin, xmax = calculate_fwhm(x, cross_section)

	# --------------
	# Find y bounds:
	dy = (extent[3] - extent[2]) / image.shape[0]

	# Generate Y cross-section
	cross_section = np.flip(np.sum(image, axis=1))  # Sum along the vertical axis
	err_cross_section = np.flip(np.sqrt(np.sum(np.square(imageErr), axis=1)))

	y = np.linspace(extent[3], extent[2], np.size(cross_section))

	ymin, ymax = find_extrema(y, cross_section)
	#ymin, ymax = calculate_fwhm(y, cross_section)

	return xmin, xmax, ymin, ymax 

# Scan x blades to find optimal source positions, save to file
# Backpropagate from blade position and return optimal source position for that point
def source_interpolation(zvb, hx, hy, length, num_gridpoints, outDir, vb_summary_file):
	# Write header to outFile
	with open(vb_summary_file, 'a', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(['# blade x pos [cm]','blade optimal source x pos [cm]'])

	# For given blade, find optimal source position
	def blade_optimal_pos(x): # written with vertical convention
		# Find y extent for blade selection
		xmin = (zvb * x)/ (zvb + length/2)
		xmax = (zvb * x)/ (zvb - length/2)
		if xmax < xmin: # make sure xmin < xmax
			tmp = xmax
			xmax = xmin
			xmin = tmp

		# Generate image of what is being backpropagated
		h_r_VB_blade_file = f'{outDir}/h_r_VB_histogram_{x}_{zvb-0.001}.dat'
		os.system(f"./optimization_scripts/backprop/back_propagate_selected {zvb} {zvb-0.001} {outDir}/h_r_VB_output.dat {h_r_VB_blade_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}")
		## Plot what is being backpropagated
		h_r_VB_blade_image_data = output_to_image_data(h_r_VB_blade_file)
		#plot_results(h_r_VB_blade_image_data, plot_type='full', save_image=f'{imageDir}08_h_r_VB_blade_{i}_image.pdf', noShow=noShow)

		# Backpropagate from selected blade
		h_r_VB_back_file = f'{outDir}/histogram_x_scan_{i:03d}.dat'
		backprop_command = f"./optimization_scripts/backprop/back_propagate_selected {zvb} 0.01 {outDir}/h_r_VB_output.dat {h_r_VB_back_file} --rectangle {xmin} {hy[0]} {xmax} {hy[3]}"
		print(colors.YELLOW + "\nrunning command:\n{}".format(backprop_command) + colors.ENDC + '\n')
		os.system(backprop_command)

		## Plot result of backpropagation
		h_r_VB_back_image_data = output_to_image_data(h_r_VB_back_file)
		#plot_results(h_r_VB_back_image_data, plot_type='full', save_image=f'{imageDir}09_h_r_VB_back_{i}_image.pdf', noShow=noShow)
		#plot_results(h_r_VB_back_image_data, plot_type='x', save_image=f'{imageDir}10_h_r_VB_back_x_{i}_image.pdf', noShow=noShow)

		# Write to top of h_r_VB_back_file header info containing x position (for plotting later)
		x_pos_info = f'# blade_xpos: {x}'
		# Read the original content of the file
		with open(h_r_VB_back_file, 'r') as file:
			original_content = file.read()

		# Write the given string to the file
		with open(h_r_VB_back_file, 'w') as file:
			file.write(x_pos_info + '\n')  # Write the string to the top line
			file.write(original_content)  # Write the original content below the string

		# Find peak of backpropagated distribution, store in hxz0
		xsrc = find_image_weighted_average(h_r_VB_back_image_data, 'x')

		return xsrc

	# --------------------------------
	# find horizontal blade optimal source positions

	# Convert .mcpl to .txt
	mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl.gz {outDir}/h_r_VB_output.dat"
	#mcpltool_command = f"mcpltool -t {outDir}/h_r_VB_output.mcpl {outDir}/h_r_VB_output.dat"
	print(colors.YELLOW + f"\nrunning command:\n{mcpltool_command}" + colors.ENDC + '\n')
	os.system(mcpltool_command)

	x_vals = np.linspace(hx[0], hx[3], num_gridpoints)

	# For each blade position, backpropagate only neutrons at that position
	for i, x in enumerate(x_vals):
		print(colors.BLUE + f'interpolating blade {i:03d} at x={x}' + colors.ENDC)
		xz0 = blade_optimal_pos(x)
		#print(f'blade {i:03d} at x={x} has hxz0={xz0}')

		# Write output to outFile
		with open(vb_summary_file, 'a', newline='') as file:
			writer = csv.writer(file, delimiter=' ')
			writer.writerow([x, xz0])


	# Remove 'h_r_VB_output.dat' to reduce storage requirements
	os.remove(f"{outDir}/h_r_VB_output.dat")

