# Functions for generation of focused VB design 
import math
import numpy as np

class colors: 
	RED = '\033[31m'
	ENDC = '\033[m'
	GREEN = '\033[32m'
	YELLOW = '\033[33m'
	BLUE = '\033[34m'

# VB blades individually focused with interpolation from file 
def generate_VB_array_bladeFocused(zvb, zdet, ydet, length, thickness, vy, source_pos_file):
	# Arrays to store final VB characteristics (vestigial implementation from when it was written in C)
	y_vals = np.zeros(10000) 
	angles = np.zeros(10000) 

	# Calculate blade angle exactly
	def calc_blade_angle(ysrc, y_blade, z_blade, y_target, z_target):
		theta_i = math.atan((y_blade - ysrc)/ z_blade)					# angle of incident ray (assumed z_src = 0)
		theta_r = math.atan((y_blade - y_target)/ (z_target - z_blade))	# angle of reflected ray hitting target
		theta_b = 0.5 * (theta_i - theta_r)								# angle of blade for given reflection
		return theta_b

	# Read source_pos_file and store data
	data = np.loadtxt(source_pos_file)
	y_blade, y_src = data[:,0], data[:,1]
	def get_ysrc_from_y(y):
		# y corresponds to y_blade
		# Linearly interpolate ysrc using y in the y_blade array
		y_src_interpolated = 0.01*np.interp(y*100, y_blade, y_src) # file is in [cm], convert to [m]
		return y_src_interpolated 

	# --------------------------------
	# Generate upper array

	# inner angle, outer angle for vb blades
	inner_angle = math.atan(vy[2] / zvb)
	outer_angle = math.atan(vy[3] / zvb)

	# 1. Generate top blade
	i = 0
	y_vals[i] = vy[3]  # center of blade = (zvb, vy[3])
	ysrc = get_ysrc_from_y(y_vals[i])	# optimal source position for top blade
	angles[i] = calc_blade_angle(ysrc, y_vals[i], zvb, ydet, zdet)  # mirror angle required for reflection from (0, ysrc) to (zdet, 0)

	# 2. Generate successive blades with following: 
	while (inner_angle <= math.atan(y_vals[i] / zvb) <= outer_angle):
		i += 1
		
		# Use ysrc of previous blade for solid angle calculation
		ysrc = get_ysrc_from_y(y_vals[i - 1])

		# 2.a: Use previous blade and previous ysrc to find angle for top left of current blade (accounting for thickness of blades)
		y_bottom_right_previous = y_vals[i - 1] + ((length / 2.) * math.tan(angles[i - 1])) - ((thickness / 2.) * math.cos(angles[i - 1]))
		theta = math.atan((y_bottom_right_previous - ysrc) / (zvb + length / 2.))	# angle from (z=0, y=ysrc) to bottom right of previous blade = angle to top left of current blade

		# Top left of new blade (accounting for thickness)
		y_top_left = (zvb - (length / 2.)) * math.tan(theta) + ysrc
		y0 = y_top_left - thickness		# current estimate for bottom left of blade

		yrc = y0			# yrc is y center of reflection on blade, initially defined for angle=0
		yrctmp = yrc + 1	# yrctmp is previous yrc, used to determine when to end loop. Define initial as arbitrarily greater than precision limit

		# 2.b: Interpolate ysrc from data file for blade focusing calculation 
		ysrc = get_ysrc_from_y(yrc)		# ysrc optimal is approximated by using ysrc for theta=0 

		# 2.c: Generate blade focused from ysrc, to y_blade, to y_target 
		j=0
		max_iter = 10
		#TODO make sure this calculation holds for ysrc!=0
		#TODO maybe try to implement ydet!=0
		while ((0.000001 < abs(yrctmp - yrc)) and (j < max_iter)):
			j+=1
			yrctmp = yrc

			delta = 0.5 * (math.atan(zdet / yrc) - math.atan(zvb / (yrc - ysrc)))	# delta is small angle change
			yrc = y0 + (math.tan(delta) * (length / 2.))
		if j == max_iter:
			print("exceeded max iteration count!")

		angles[i] = math.atan((yrc - y0) / (length / 2.))			# Blade angle
		y_vals[i] = yrc + (thickness / 2.) * math.cos(angles[i])	# Blade center position, !not center of reflective surface

	# --------------------------------
	# Generate lower array

	#TODO make this work for the lower array
	# inner angle, outer angle for vb blades
	inner_angle = math.atan(vy[1] / zvb)
	outer_angle = math.atan(vy[0] / zvb)
	print('lower blade array:')
	print(inner_angle, outer_angle)

	# 1. Generate bottom blade
	# write over blade i, it is lower than <inner_angle>
	y_vals[i] = vy[0]  # center of blade = (zvb, vy[0])
	ysrc = get_ysrc_from_y(y_vals[i])	# optimal source position for top blade
	angles[i] = calc_blade_angle(ysrc, y_vals[i], zvb, ydet, zdet)  # mirror angle required for reflection from (0, ysrc) to (zdet, 0)

	# 2. Generate successive blades with following: 
	while (outer_angle <= math.atan(y_vals[i] / zvb) <= inner_angle):
		print('i value lower blades:')
		print(i)
		i += 1
		
		# Use ysrc of previous blade for solid angle calculation
		ysrc = get_ysrc_from_y(y_vals[i - 1])
		print(f'yval: {y_vals[i-1]}, ysrc: {ysrc}')

		# 2.a: Use previous blade and previous ysrc to find angle for bottom left of current blade (accounting for thickness of blades)
		y_top_right_previous = y_vals[i - 1] + ((length / 2.) * math.tan(angles[i - 1])) + ((thickness / 2.) * math.cos(angles[i - 1]))
		theta = math.atan((y_top_right_previous - ysrc) / (zvb + length / 2.))	# angle from (z=0, y=ysrc) to top right of previous blade = angle to bottom left of current blade
		print('y top right previous')
		print(y_top_right_previous)

		# Bottom left of new blade (accounting for thickness)
		y_bot_left = (zvb - (length / 2.)) * math.tan(theta) + ysrc
		y0 = y_bot_left + thickness		# current estimate for top left of blade

		yrc = y0			# yrc is y center of reflection on blade, initially defined for angle=0
		yrctmp = yrc + 1	# yrctmp is previous yrc, used to determine when to end loop. Define initial as arbitrarily greater than precision limit

		# 2.b: Interpolate ysrc from data file for blade focusing calculation 
		ysrc = get_ysrc_from_y(yrc)		# ysrc optimal is approximated by using ysrc for theta=0 
		print(f'yrc: {yrc}')
		print(f'ysrc: {ysrc}')

		# 2.c: Generate blade focused from ysrc, to y_blade, to y_target 
		j=0
		max_iter = 10
		#TODO see 2.c top array
		while ((0.000001 < abs(yrctmp - yrc)) and (j < max_iter)):
			j+=1
			yrctmp = yrc

			delta = 0.5 * (math.atan(zdet / yrc) - math.atan(zvb / (yrc - ysrc)))
			yrc = y0 + math.tan(delta) * (length / 2.)
			print(f'==========\n{j} iteration')
			print(f'delta {delta}, yrc {yrc}, y0 {y0}')
		if j == max_iter:
			print("exceeded max iteration count!")

		angles[i] = math.atan((yrc - y0) / (length / 2.))			# Blade angle
		y_vals[i] = yrc + (thickness / 2.) * math.cos(angles[i])	# Blade center position, !not center of reflective surface
		print(y_vals[i])
		print(outer_angle,math.atan(y_vals[i] / zvb),inner_angle)

	print('top and lower blades')
	print(y_vals)

	# 3. Add small change in blade angles to achieve targeted reflection to ydet
#	for j in range(numBlades):
#		# if j < i, reflecting surface on bottom of blade
#		# if j >= i, reflecting surface on top of blade
#		if j < i:
#			yrc = y_vals[j] - (thickness / 2.) * math.cos(angles[j])
#		else:
#			yrc = y_vals[j] + (thickness / 2.) * math.cos(angles[j])
#
#		# zvb != zrc because blades are rotated about the center of the blade,
#		# zrc = zvb + (thickness/2.) * abs(sin(angles[j])).
#		# approximation of zrc=zvb is valid for small angles[], small thickness
#
#		# shift all blades by d_angle to target reflection to (zdet, ydet) relative (zvb, 0)
#		d_angle = 0.5 * (math.atan(yrc / zdet) - math.atan((yrc - ydet) / zdet))
#		angles[j] += d_angle

	#TODO perform checks to make sure blades dont overlap
	if y_vals[i] > y_vals[i - 1] - thickness:
		print("ERROR: blade overlap! max overlap @center=", y_vals[i] - (y_vals[i - 1] - thickness))
	#if numBlades > 9999:
	#	print("ERROR: numBlades greater than blade limit!")
	
	# remove trailing zeros from y_vals, angles
	return np.trim_zeros(y_vals, 'b'), np.trim_zeros(angles, 'b')

def generate_VB_bladeFocused(VB_pos, zdet, det_pos, VB_length, VB_thickness, vx, vy, hx, hy, source_pos_file):
	print(colors.BLUE + f'VB parameters:\nVB_pos={VB_pos}, zdet={zdet}, det_pos={det_pos}, VB_length={VB_length}, VB_thickness={VB_thickness}, vy={vy}, hx={hx}, source_pos_filename={source_pos_file}' + colors.ENDC)

	# vy, hx [cm]
	vy *= 0.01
	hx *= 0.01
	# xdet, ydet [m]
	xdet, ydet = det_pos

	hvb = vy[3]-vy[0]
	wvb = hx[3]-hx[0]
	y_center = 0.5*(vy[3] + vy[0])
	x_center = -0.5*(hx[3] + hx[0])

	# --------------------------------
	# Create vertically reflecting blade geometry with point focus!!!!!!!!
	y_vals_v, angles_v = generate_VB_array_pointFocused(VB_length, VB_thickness, VB_pos, zdet, 100, 0.01, ydet) # hvb = 100, hdet = 0.01 

	# Apply restrictions to blade geometry:
	mask1_v = np.logical_and(y_vals_v >= vy[0], y_vals_v <= vy[1])
	mask2_v = np.logical_and(y_vals_v >= vy[2], y_vals_v <= vy[3])
	final_mask_v = np.logical_or(mask1_v, mask2_v)

	# apply the mask to both x_coords and y_coords
	y_vals_v = y_vals_v[final_mask_v]
	angles_v = angles_v[final_mask_v]

	# Write generated geometry to file
	v_blades = []

	print(f'wvb: {wvb}, x_center: {x_center}')
	for i in range(len(y_vals_v)):
		# Append each combination of parameters to v_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		v_blades.append((wvb, VB_thickness, VB_length, x_center, y_vals_v[i], 0.0, angles_v[i])) # zvb := 0.0 in relation to center of component

	print(y_vals_v)
	print(angles_v)
	print(f'number of v_blades: {y_vals_v.size}')
	v_reflecting_venbla_geometry = f"Venbla_vertically_reflecting_geometry_pointFocused_{VB_pos}_{zdet}_{det_pos[0]}_{det_pos[1]}_{VB_length}_{VB_thickness}.off"
	write_VB(v_blades, v_reflecting_venbla_geometry)

	# --------------------------------
	# Create horizontally reflecting blade geometry with interpolated individual blade focusing!!!!!!!!!
	x_vals_h, angles_h = generate_VB_array_bladeFocused(VB_pos, zdet, xdet, VB_length, VB_thickness, hx, source_pos_file)

	# Write the blade geometry to file:
	h_blades = []

	print(f'hvb: {hvb}, y_center: {y_center}')
	for i in range(len(x_vals_h)):
		# Append each combination of parameters to h_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		h_blades.append((hvb, VB_thickness, VB_length, y_center, x_vals_h[i], 0.0, angles_h[i])) # zvb := 0.0 in relation to center of component

	print(x_vals_h)
	print(angles_h)
	print(f'number of h_blades: {x_vals_h.size}')
	h_reflecting_venbla_geometry = f"Venbla_horizontally_reflecting_geometry_bladeFocused_{VB_pos}_{zdet}_{det_pos[0]}_{det_pos[1]}_{VB_length}_{VB_thickness}.off"
	write_VB(h_blades, h_reflecting_venbla_geometry)

	return v_reflecting_venbla_geometry, h_reflecting_venbla_geometry 

# VB blades with pointlike source
def generate_VB_point_focused(zvb, zdet, det_pos, length, thickness, vy, hx):
	print(colors.BLUE + f'VB parameters:\nVB_pos={zvb}, zdet={zdet}, det_pos={det_pos}, VB_length={length}, VB_thickness={thickness}, vy={vy}, hx={hx}' + colors.ENDC)

	# vy, hx [cm]
	vy0, vy1, vy2, vy3 = 0.01*vy
	hx0, hx1, hx2, hx3 = 0.01*hx
	# det_pos [m]
	xdet, ydet = det_pos

	hvb = vy3-vy0
	wvb = hx3-hx0
	y_center = 0.5*(vy3 + vy0)
	x_center = -0.5*(hx3 + hx0)

	# --------------------------------
	# Create vertically reflecting blade geometry
	y_vals_v, angles_v = generate_VB_array_pointFocused(length, thickness, zvb, zdet, 100, 0.01, ydet) # hvb = 100, hdet = 0.01 

	# Apply restrictions to blade geometry:
	mask1_v = np.logical_and(y_vals_v >= vy0, y_vals_v <= vy1)
	mask2_v = np.logical_and(y_vals_v >= vy2, y_vals_v <= vy3)
	final_mask_v = np.logical_or(mask1_v, mask2_v)

	# apply the mask to both x_coords and y_coords
	y_vals_v = y_vals_v[final_mask_v]
	angles_v = angles_v[final_mask_v]

	# Write the blade geometry to file:
	v_blades = []

	print(f'wvb: {wvb}, x_center: {x_center}')
	for i in range(len(y_vals_v)):
		# Append each combination of parameters to v_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		v_blades.append((wvb, thickness, length, x_center, y_vals_v[i], 0.0, angles_v[i])) # zvb := 0.0 in relation to center of component

	print(y_vals_v)
	print(angles_v)
	print(f'number of v_blades: {y_vals_v.size}')
	v_reflecting_venbla_geometry = f"Venbla_vertically_reflecting_geometry_pointFocused_{zvb}_{zdet}_{det_pos[0]}_{det_pos[1]}_{length}_{thickness}.off"
	write_VB(v_blades, v_reflecting_venbla_geometry)

	# --------------------------------
	# Create horizontally reflecting blade geometry
	x_vals_h, angles_h = generate_VB_array_pointFocused(length, thickness, zvb, zdet, 100, 0.01, xdet) # wvb = 100, wdet = 0.01

	# Apply restrictions to blade geometry:
	mask1_h = np.logical_and(x_vals_h >= hx0, x_vals_h <= hx1)
	mask2_h = np.logical_and(x_vals_h >= hx2, x_vals_h <= hx3)
	final_mask_h = np.logical_or(mask1_h, mask2_h)

	# apply the mask to both x_coords and y_coords
	x_vals_h = x_vals_h[final_mask_h]
	angles_h = angles_h[final_mask_h]

	# Write the blade geometry to file:
	h_blades = []

	print(f'hvb: {hvb}, y_center: {y_center}')
	for i in range(len(x_vals_h)):
		# Append each combination of parameters to h_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		h_blades.append((hvb, thickness, length, y_center, x_vals_h[i], 0.0, angles_h[i])) # zvb := 0.0 in relation to center of component
		#h_blades.append((hvb, thickness, length, x_vals_h[i], 0.0, angles_h[i])) # zvb := 0.0 in relation to center of component

	print(x_vals_h)
	print(angles_h)
	print(f'number of h_blades: {x_vals_h.size}')
	h_reflecting_venbla_geometry = f"Venbla_horizontally_reflecting_geometry_bladeFocused_{zvb}_{zdet}_{det_pos[0]}_{det_pos[1]}_{length}_{thickness}.off"
	write_VB(h_blades, h_reflecting_venbla_geometry)

	return v_reflecting_venbla_geometry, h_reflecting_venbla_geometry 

# VB blades individually focused 
def generate_VB_focused_blades(VB_pos, zdet, det_pos, VB_length, VB_thickness, vx, vy, hx, hy, vyz0, hxz0):
	print(colors.BLUE + f'VB parameters:\nVB_pos={VB_pos}, zdet={zdet}, det_pos={det_pos}, VB_length={VB_length}, VB_thickness={VB_thickness}, vy={vy}, hx={hx}, vyz0={vyz0}, hxz0={hxz0}' + colors.ENDC)

	# vy, hx [cm]
	vy *= 0.01
	hx *= 0.01
	# xdet, ydet [m]
	xdet, ydet = det_pos

	hvb = vy3-vy0
	wvb = hx3-hx0
	y_center = 0.5*(vy3 + vy0)
	x_center = -0.5*(hx3 + hx0)

	# --------------------------------
	# Create vertically reflecting blade geometry
	y_vals_v, angles_v = generate_VB_array_pointFocused(length, thickness, zvb, zdet, ydet, vy, vyz0) 

	# Write generated geometry to file
	v_blades = []

	print(f'wvb: {wvb}, x_center: {x_center}')
	for i in range(len(y_vals_v)):
		# Append each combination of parameters to v_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		v_blades.append((wvb, thickness, length, x_center, y_vals_v[i], 0.0, angles_v[i])) # zvb := 0.0 in relation to center of component

	print(y_vals_v)
	print(angles_v)
	print(f'number of v_blades: {y_vals_v.size}')
	v_reflecting_venbla_geometry = f"Venbla_vertically_reflecting_geometry_bladeFocused_{zvb}_{zdet}_{det_pos[0]}_{det_pos[1]}_{length}_{thickness}.off"
	write_VB(v_blades, v_reflecting_venbla_geometry)

	# --------------------------------
	# Create horizontally reflecting blade geometry
	x_vals_h, angles_h = generate_VB_array_bladeFocused(length, thickness, zvb, zdet, xdet, hx, hxz0)

	# Write the blade geometry to file:
	h_blades = []

	print(f'hvb: {hvb}, y_center: {y_center}')
	for i in range(len(x_vals_h)):
		# Append each combination of parameters to h_blades list
		# x extent [m], y extent [m], z extent [m], xpos[m], ypos [m], zpos [m], angle [deg]
		h_blades.append((hvb, thickness, length, y_center, x_vals_h[i], 0.0, angles_h[i])) # zvb := 0.0 in relation to center of component

	print(x_vals_h)
	print(angles_h)
	print(f'number of h_blades: {x_vals_h.size}')
	h_reflecting_venbla_geometry = f"Venbla_horizontally_reflecting_geometry_bladeFocused_{zvb}_{zdet}_{det_pos[0]}_{det_pos[1]}_{length}_{thickness}.off"
	write_VB(h_blades, h_reflecting_venbla_geometry)

	return v_reflecting_venbla_geometry, h_reflecting_venbla_geometry 

# VB focused to y=0
# assumed ysrc = 0
def generate_VB_array_pointFocused(length, thickness, zvb, zdet, hvb, hdet, ydet):
	y_vals = np.zeros(10000) 
	angles = np.zeros(10000) 

	# inner angle, outer angle for vb blades
	inner_angle = math.atan((hdet / 2.) / (zvb + zdet))
	outer_angle = math.atan((hvb / 2.) / zvb)

	# define top blade
	i = 0
	y_vals[i] = hvb / 2.  # center of blade = (zvb, hvb/2.)
	angles[i] = 0.5 * (math.atan((hvb / 2.) / zvb) - math.atan((hvb / 2.) / zdet))  # mirror angle required for reflection to (zdet, 0)

	# define upper blade array
	while (inner_angle <= math.atan(y_vals[i] / zvb) <= outer_angle):
		i += 1

		# angle for ray to max z of previous blade (bottom right accounting for thickness)
		theta = math.atan((y_vals[i - 1] + (length / 2.) * math.tan(angles[i - 1]) - (thickness / 2.) * math.cos(angles[i - 1])) / (zvb + length / 2.))

		# top left of new blade (accounting for thickness)
		y0 = (zvb - (length / 2.)) * math.tan(theta)
		y0 -= thickness

		yrc = y0			# yrc is y center of reflection on blade, initially defined for angle=0
		yrctmp = yrc + 1	# yrctmp is previous yrc, used to determine when to end loop. Define initial as arbitrarily greater
		ysrc = 0			# placeholder for y pos of source

		j=0
		max_iter = 10
		while ((0.000001 < abs(yrctmp - yrc)) and (j < max_iter)):
			j+=1
			yrctmp = yrc

			delta = 0.5 * (math.atan(zdet / yrc) - math.atan(zvb / (yrc - ysrc)))
			yrc = y0 + math.tan(delta) * (length / 2.)
		if j == max_iter:
			print("exceeded max iteration count!")

		angles[i] = math.atan((yrc - y0) / (length / 2.))
		y_vals[i] = yrc + (thickness / 2.) * math.cos(angles[i])  # y_vals stores center of blade, not center of reflective surface

	numBlades = 2 * i

	# write over blade i, it is lower than <inner_angle>
	for j in range(i):
		# define lower blades
		angles[j + i] = -angles[j]
		y_vals[j + i] = -y_vals[j]

	# targeted reflection
	for j in range(numBlades):
		# if j < i, reflecting surface on bottom of blade
		# if j >= i, reflecting surface on top of blade
		if j < i:
			yrc = y_vals[j] - (thickness / 2.) * math.cos(angles[j])
		else:
			yrc = y_vals[j] + (thickness / 2.) * math.cos(angles[j])

		# zvb != zrc because blades are rotated about the center of the blade,
		# zrc = zvb + (thickness/2.) * abs(sin(angles[j])).
		# approximation of zrc=zvb is valid for small angles[], small thickness

		# shift all blades by d_angle to target reflection to (zdet, ydet) relative (zvb, 0)
		d_angle = 0.5 * (math.atan(yrc / zdet) - math.atan((yrc - ydet) / zdet))
		angles[j] += d_angle

	if y_vals[i] > y_vals[i - 1] - thickness:
		print("ERROR: blade overlap! max overlap @center=", y_vals[i] - (y_vals[i - 1] - thickness))
	if numBlades > 9999:
		print("ERROR: numBlades greater than blade limit!")
	
	# remove trailing zeros from y_vals, angles
	return np.trim_zeros(y_vals, 'b'), np.trim_zeros(angles, 'b')

# Functions for writing geometry to OFF file
def rotate_vertex(vertex, rotation_angles):
	a = rotation_angles[0] #math.radians(rotation_angles[0])
	cos_a = math.cos(a)
	sin_a = math.sin(a)

	# Rotation matrix for x-axis
	rotation_matrix = [
		[1, 0, 0],
		[0, cos_a, -sin_a],
		[0, sin_a, cos_a]
	]

	return [sum(rotation_matrix[i][j] * vertex[j] for j in range(3)) for i in range(3)]

def write_rectangular_prism(vertices, faces, filename):
	# Write the OFF file
	with open(filename, "w") as file:
		file.write("OFF\n")
		file.write(f"{len(vertices)} {len(faces)} 0\n")

		# Write vertices
		for vertex in vertices:
			file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

		# Write faces
		for face in faces:
			file.write(f"{len(face)} {' '.join(map(str, face))}\n")

def generate_rectangular_prism(x_extent, y_extent, z_extent, xpos, ypos, zpos, a):
	# Rectangular prism vertices
	vertices = [
		[-x_extent / 2, -y_extent / 2, -z_extent / 2],  # Vertex 0
		[x_extent / 2, -y_extent / 2, -z_extent / 2],   # Vertex 1
		[x_extent / 2, -y_extent / 2, z_extent / 2],	# Vertex 2
		[-x_extent / 2, -y_extent / 2, z_extent / 2],   # Vertex 3
		[-x_extent / 2, y_extent / 2, -z_extent / 2],   # Vertex 4
		[x_extent / 2, y_extent / 2, -z_extent / 2],	# Vertex 5
		[x_extent / 2, y_extent / 2, z_extent / 2],	 # Vertex 6
		[-x_extent / 2, y_extent / 2, z_extent / 2]	 # Vertex 7
	]

	# Rotate and translate the rectangular prism vertices
	rotated_vertices = [rotate_vertex(vertex, [-a, 0.0, 0.0]) for vertex in vertices]
	translated_vertices = [[vertex[0] + xpos, vertex[1] + ypos, vertex[2] + zpos] for vertex in rotated_vertices]

	# Rectangular prism faces (vertices are 0-indexed)
	faces = [
		[0, 1, 2, 3],  # Bottom face
		[0, 4, 5, 1],  # Side faces
		[1, 5, 6, 2],
		[2, 6, 7, 3],
		[3, 7, 4, 0],
		[4, 7, 6, 5],  # Top face
	]

	return translated_vertices, faces

def write_VB(rectangular_prisms, filename):
	all_vertices = []
	all_faces = []

	for prism_info in rectangular_prisms:
		x_extent, y_extent, z_extent, xpos, ypos, zpos, a = prism_info
		vertices, faces = generate_rectangular_prism(x_extent, y_extent, z_extent, xpos, ypos, zpos, a)
		num_vertices = len(all_vertices)

		# Add current prism vertices to the combined list
		all_vertices.extend(vertices)

		# Add current prism faces to the combined list
		all_faces.extend([face_index + num_vertices for face_index in face] for face in faces)

	write_rectangular_prism(all_vertices, all_faces, filename)
