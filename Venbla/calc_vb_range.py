import matplotlib.pyplot as plt
import numpy as np

# constants
hdet = 0.8 # height of detector [m]

# points for calculation: (z, y) [m]
top_source = np.array([0.01, 0.015])
mid_source = np.array([0.01, 0.000])
bot_source = np.array([0.01, -0.015])

top_corner_ellipse = np.array([5.374, 0.035])
bot_corner_ellipse = np.array([5.374, -0.100])

top_corner_target = np.array([65, hdet/2.])
bot_corner_target = np.array([65, -hdet/2.])

# Calculate equations of lines
def calculate_line_equation(point1, point2):
	slope = (point2[1] - point1[1]) / (point2[0] - point1[0])
	intercept = point1[1] - slope * point1[0]
	return slope, intercept

# Calculate intersections of lines and plot
def plot_intersections(z_value):
	plt.figure(figsize=(10, 6))

	plt.plot(*top_source, 'ro', label='Top Source')
	plt.plot(*mid_source, 'go', label='Mid Source')
	plt.plot(*bot_source, 'bo', label='Bot Source')
	plt.plot(*top_corner_ellipse, 'mx', label='Top Corner Ellipse')
	plt.plot(*bot_corner_ellipse, 'cx', label='Bot Corner Ellipse')
	plt.plot(*top_corner_target, 'ms', label='Top Corner Target')
	plt.plot(*bot_corner_target, 'cs', label='Bot Corner Target')

	plt.vlines(x=z_value, ymin=-1.0, ymax=1.0, label=f'zvb={z_value}')

	for eq in equations:
		x_vals = np.linspace(0, 70, 100)
		y_vals = eq[0] * x_vals + eq[1]
		plt.plot(x_vals, y_vals)

		y_intersection = eq[0] * z_value + eq[1]
		plt.plot(z_value, y_intersection, 'kx')  # Plot intersection point
		print(y_intersection)

	plt.xlabel('Z [m]')
	plt.ylabel('Y [m]')
	plt.title('Lines, Points, and Intersections')
	plt.legend()
	plt.grid(True)
	plt.show()

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Find Venbla Optimal Geometry")
	parser.add_argument("zvb", type=float, help="Z displacement of venetian blinds from source")
	args = parser.parse_args()
	zvb = args.zvb

	# Calculate equations for the specified lines
	lines = [
		(top_source, bot_corner_ellipse),
		(mid_source, bot_corner_ellipse),
		(bot_source, top_corner_ellipse),
		(mid_source, top_corner_ellipse),
		(top_corner_ellipse, top_corner_target),
		(bot_corner_ellipse, bot_corner_target)
	]

	equations = [calculate_line_equation(p1, p2) for p1, p2 in lines]

	print(plot_intersections(zvb))
