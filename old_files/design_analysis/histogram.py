# histograms detections for x, y
import matplotlib.pylab as plt
import numpy as np

data = np.loadtxt('data/TOF_target_radial.dat')
#data = np.loadtxt('data/TOF_l50cm.dat')
X = data[:,1]
Y = data[:,2]

hist, xedges, yedges = np.histogram2d(Y, X, bins=100)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

fig = plt.figure()
hist_plot = fig.add_subplot(111)
plt.xlabel('Horizontal Position (m)')
plt.ylabel('Vertical Position (m)')

pos = hist_plot.imshow(hist, extent=extent, origin='lower', aspect='auto', cmap='plasma')
fig.colorbar(pos,label='Total Counts')

plt.show()
