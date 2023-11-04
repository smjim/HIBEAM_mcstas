import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import os
import re

def show_instr():
	os.system('mcdisplay target40cm_off.instr')
	#os.system('mcdisplay-webgl target40cm_off.instr')

# Plot any number of distributions on 1d projection
def plot_results(*image_data_list, plot_type='full', save_image=None, xlims=None, ylims=None, noShow=False):
	ax = None

	for i, image_data in enumerate(image_data_list):
		if ax is None:
			fig, ax = plt.subplots(figsize=(8, 6))
		else:
			fig = ax.figure

		dataHeader, extent, image, imageErr = image_data
		position = dataHeader['position']
		component = dataHeader['component']
		title = f"{component}; ({position})m"
		#print(f'position: {position}, i: {i}')

		if plot_type == "x":
			unit = re.findall(r"\[(.*?)\]", dataHeader['xlabel'])
			dx = (extent[1] - extent[0]) / image.shape[1]

			# Generate X cross-section cross_section
			cross_section = np.sum(image, axis=0)  # Sum along the horizontal axis
			err_cross_section = np.sqrt(np.sum(np.square(imageErr), axis=0))

			# Plot X bounds
			if xlims is not None:
				for x in xlims:
					plt.axvline(x, color='black', linestyle='--')

			x = np.linspace(extent[0], extent[1], np.size(cross_section))
			plt.errorbar(x, cross_section, err_cross_section, capsize=2, label=f'{component}')
			plt.xlabel(dataHeader['xlabel'])
			plt.ylabel('Intensity [n/s]/ '+"{:.2e}".format(dx)+' ['+unit[0]+']')
			plt.title("X Cross-Section: "+title)

		elif plot_type == "y":
			unit = re.findall(r"\[(.*?)\]", dataHeader['ylabel'])
			dx = (extent[3] - extent[2]) / image.shape[0]

			# Generate Y cross-section
			cross_section = np.flip(np.sum(image, axis=1))  # Sum along the horizontal axis
			err_cross_section = np.flip(np.sqrt(np.sum(np.square(imageErr), axis=1)))

			# Plot Y bounds
			if ylims is not None:
				for y in ylims:
					plt.axvline(y, color='black', linestyle='--')

			y = np.linspace(extent[3], extent[2], np.size(cross_section))
			plt.errorbar(y, cross_section, err_cross_section, capsize=2, label=f'{component}')
			plt.xlim(extent[2], extent[3])
			plt.xlabel(dataHeader['ylabel'])
			plt.ylabel('Intensity [n/s]/ '+"{:.2e}".format(dx)+' ['+unit[0]+']')
			plt.title("Y Cross-Section: "+title)

		elif plot_type == "full":
			unit1 = re.findall(r"\[(.*?)\]", dataHeader['xlabel'])
			unit2 = re.findall(r"\[(.*?)\]", dataHeader['ylabel'])
			# Calculate bin size
			dx = (extent[1] - extent[0]) / (image.shape[1] - 1)
			dy = (extent[3] - extent[2]) / (image.shape[0] - 1)

			#plt.imshow(image, extent=extent, cmap='plasma', norm='log')
			plt.imshow(np.flipud(image), extent=extent, cmap='plasma', norm='log')
			plt.colorbar().set_label('Intensity [n/s]/ '+"{:.2e}".format(dx*dy)+' ['+unit1[0]+'*'+unit2[0]+']')
			plt.xlabel(dataHeader['xlabel'])
			plt.ylabel(dataHeader['ylabel'])
			plt.title(title, pad=10)

			plt.grid()
			plt.tight_layout()
			plt.savefig(save_image, format='pdf')
			if not noShow:
				plt.show()
			return

		else:
			print("Invalid plot type. Please specify 'x', 'y', or 'full'.")


	plt.grid()
	plt.legend()
	plt.tight_layout()
	plt.savefig(save_image, format='pdf')
	if not noShow:
		plt.show()

def count_results(image_data, square=None, circle=None, noShow=False, save_image=None):
	dataHeader, extent, image, imageErr = image_data

	dx = (extent[1] - extent[0]) / (image.shape[1] - 1)
	dy = (extent[3] - extent[2]) / (image.shape[0] - 1)
	mask = np.zeros_like(image, dtype=bool)
	roi_area = 0

	if square:
		x0, x1, y0, y1 = square
		x_indices = np.where((x0 <= extent[0] + np.arange(image.shape[1]) * dx) & (extent[0] + np.arange(image.shape[1]) * dx < x1))
		y_indices = np.where((-y1 <= extent[2] + np.arange(image.shape[0]) * dy) & (extent[2] + np.arange(image.shape[0]) * dy < -y0))
		mask[y_indices[0][:, np.newaxis], x_indices[0]] = True
		roi_area = (x1 - x0) * (y1 - y0)

	if circle:
		x0, y0, radius = circle
		x_indices, y_indices = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
		mask = ((extent[0] + x_indices * dx - x0) ** 2 + (extent[2] + y_indices * dy - y0) ** 2 <= radius ** 2)
		roi_area = np.pi * radius ** 2

	roi_sum = np.sum(image[mask])
	sum_err = np.sqrt(np.sum(np.square(imageErr[mask])))
	unit1 = re.findall(r"\[(.*?)\]", dataHeader['xlabel'])
	unit2 = re.findall(r"\[(.*?)\]", dataHeader['ylabel'])

	print("Sum within ROI: ", "{:.2e}".format(roi_sum), " ± ", "{:.2e}".format(sum_err))
	print("Area within ROI: ", "{:.2e}".format(roi_area) + ' [' + unit1[0] + '*' + unit2[0] + ']\n')

	# Show plot
	fig, ax = plt.subplots()
	img = ax.imshow(np.flipud(image), extent=extent, cmap='plasma')
	ax.set_title(f"{dataHeader['component']}; ({dataHeader['position']})m")
	ax.set_xlabel(dataHeader['xlabel'])
	ax.set_ylabel(dataHeader['ylabel'])
	cbar = fig.colorbar(img, ax=ax)
	cbar.set_label(dataHeader['zvar'] + '/ ' + "{:.2e}".format(dx * dy) + ' [' + unit1[0] + '*' + unit2[0] + ']')

	if square:
		square = Rectangle((x0, y0), (x1 - x0), (y1 - y0), fill=False, color='red', linewidth=2)
		ax.add_patch(square)
	if circle:
		circle = Circle((x0, y0), radius, fill=False, color='red', linewidth=2)
		ax.add_patch(circle)

	# Display ROI sum and error as text on the plot
	roi_label = f"ROI Sum: {roi_sum:.2e} ± {sum_err:.2e}\nROI Area: {roi_area:.2e} [{unit1[0]}*{unit2[0]}]"
	ax.text(0.05, 0.95, roi_label, transform=ax.transAxes, fontsize=12, va='top', ha='left', backgroundcolor='white')

	#plt.grid()
	plt.tight_layout()
	plt.savefig(save_image, format='pdf')

	if not noShow:
		plt.show()
	
	return roi_sum, sum_err
