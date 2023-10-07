import math
import numpy as np

# Define geometry
def generate_venbla(length, thickness, zvb, zdet, hvb, hdet, ydet):
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
	if numBlades > 999:
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

def generate_rectangular_prism(x_extent, y_extent, z_extent, ypos, zpos, a):
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
	translated_vertices = [[vertex[0], vertex[1] + ypos, vertex[2] + zpos] for vertex in rotated_vertices]

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

def write_off(rectangular_prisms, filename):
	all_vertices = []
	all_faces = []

	for prism_info in rectangular_prisms:
		x_extent, y_extent, z_extent, ypos, zpos, a = prism_info
		vertices, faces = generate_rectangular_prism(x_extent, y_extent, z_extent, ypos, zpos, a)
		num_vertices = len(all_vertices)

		# Add current prism vertices to the combined list
		all_vertices.extend(vertices)

		# Add current prism faces to the combined list
		all_faces.extend([face_index + num_vertices for face_index in face] for face in faces)

	write_rectangular_prism(all_vertices, all_faces, filename)

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Create Venetian Blinds geometry with given parameters")

	parser.add_argument("length", type=float, help="Z length of blades")
	parser.add_argument("thickness", type=float, help="Thickness of blades")
	parser.add_argument("zvb", type=float, help="Z distance between source and center of vb")
	parser.add_argument("zdet", type=float, help="Z distance between center of vb and target")
	parser.add_argument('--detdim', nargs=2, type=float, metavar=('detwidth', 'detheight'),
                    help='width and height of detector')
	parser.add_argument('--detpos', nargs=2, type=float, metavar=('detx', 'dety'),
                    help='displacement of detector position')
	parser.add_argument('--xbounds', nargs=4, type=float, metavar=('x0', 'x1', 'x2', 'x3'),
                    help='X bounds for creation of 2 sets of venetian blinds')
	parser.add_argument('--ybounds', nargs=4, type=float, metavar=('y0', 'y1', 'y2', 'y3'),
                    help='Y bounds for creation of 2 sets of venetian blinds')

	# Parse arguments
	args = parser.parse_args()

	length = args.length
	thickness = args.thickness
	zvb = args.zvb
	zdet = args.zdet

	det_dims = args.detdim
	det_pos = args.detpos
	x_bounds = args.xbounds
	y_bounds = args.ybounds

	# --------------------------------
	# Extract parameters
	wvb = 2*max(np.abs(x_bounds[0]), np.abs(x_bounds[3]))
	hvb = 2*max(np.abs(y_bounds[0]), np.abs(y_bounds[3]))

	wdet, hdet = det_dims
	xdet, ydet = det_pos

	# --------------------------------
	# Create vertically reflecting blade geometry
	y_vals_v, angles_v = generate_venbla(length, thickness, zvb, zdet, 100, 0.01, ydet)
	#y_vals_v, angles_v = generate_venbla(length, thickness, zvb, zdet, hvb, hdet, ydet)

	# Apply restrictions to blade geometry:
	mask1_v = np.logical_and(y_vals_v >= y_bounds[0], y_vals_v <= y_bounds[1])
	mask2_v = np.logical_and(y_vals_v >= y_bounds[2], y_vals_v <= y_bounds[3])
	final_mask_v = np.logical_or(mask1_v, mask2_v)

	# apply the mask to both x_coords and y_coords
	y_vals_v = y_vals_v[final_mask_v]
	angles_v = angles_v[final_mask_v]

	# Write the blade geometry to file:
	# x extent [m], y extent [m], z extent [m], ypos [m], zpos [m], angle [deg]
	v_blades = []
	
	for i in range(len(y_vals_v)):
		# Append each combination of parameters to v_blades list
		v_blades.append((wvb, thickness, length, y_vals_v[i], 0.0, angles_v[i])) # zvb := 0.0 in relation to center of component
	
	print(y_vals_v)
	print(angles_v)
	print(f'number of v_blades: {y_vals_v.size}')
	filename = "Venbla_vertically_reflecting_geometry.off"
	write_off(v_blades, filename)

	# --------------------------------
	# Create horizontally reflecting blade geometry
	x_vals_h, angles_h = generate_venbla(length, thickness, zvb+(length+0.01), zdet-(length+0.01), 100, 0.01, xdet)
	#x_vals_h, angles_h = generate_venbla(length, thickness, zvb+(length+0.01), zdet-(length+0.01), wvb, wdet, xdet)

	# Apply restrictions to blade geometry:
	mask1_h = np.logical_and(x_vals_h >= x_bounds[0], x_vals_h <= x_bounds[1])
	mask2_h = np.logical_and(x_vals_h >= x_bounds[2], x_vals_h <= x_bounds[3])
	final_mask_h = np.logical_or(mask1_h, mask2_h)
	
	# apply the mask to both x_coords and y_coords
	x_vals_h = x_vals_h[final_mask_h]
	angles_h = angles_h[final_mask_h]
	
	# Write the blade geometry to file:
	# x extent [m], y extent [m], z extent [m], ypos [m], zpos [m], angle [deg]
	h_blades = []
	
	for i in range(len(x_vals_h)):
		# Append each combination of parameters to h_blades list
		h_blades.append((hvb, thickness, length, x_vals_h[i], 0.0, angles_h[i])) # zvb := 0.0 in relation to center of component
	
	print(x_vals_h)
	print(angles_h)
	print(f'number of h_blades: {x_vals_h.size}')
	filename = "Venbla_horizontally_reflecting_geometry.off"
	write_off(h_blades, filename)
