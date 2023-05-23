# arguments: <filename> <y medium>
import numpy as np
import matplotlib.pyplot as plt
import sys

# p, x, y, r_xy, tof 
inFile = sys.argv[1]
data = np.loadtxt(inFile, dtype='float')

y_m = 0.001 #0.95
i = 100
size = 0.01 #0.5
ymin = y_m-size
ymax = y_m+size
dy = (ymax - ymin)/ (2*i)
print('resolution [m]: ', str(dy))
y_vals = np.linspace(ymin, ymax, i)
fomDat = np.zeros((i, 3))
fomDat[:,0] += y_vals
# fomDat = [[y_vals, fom, weight],...]

for line in data: 
	weight = line[0]	# neutron weight
	x_val = line[1]		# neutron x detection
	y_val = line[2]		# neutron y detection
#	radius = line[3]	# detection radius (xy)
	time_f = line[4]	# time of flight [s]
	for j in range(i):
		y = y_vals[j]
		if ((y_val - dy < y) & (y < y_val + dy)):
			fomDat[j,1] += (time_f**2)*weight
			fomDat[j,2] += weight

def find_fwhm(y_values, weights):
    max_y = np.max(y_values)
    hm_y = max_y / 2
    left_idx = np.argmin(np.abs(y_values[:np.argmax(y_values)] - hm_y))
    right_idx = np.argmin(np.abs(y_values[np.argmax(y_values):] - hm_y)) + np.argmax(y_values)
    fwhm = weights[right_idx] - weights[left_idx]

    return fwhm, left_idx, right_idx 

fwhm, left_idx, right_idx = find_fwhm(fomDat[:,2], y_vals)
# find total neutron weight within center +- fwhm
fom = 0
weight = 0
for i in range(left_idx, right_idx):
	fom += fomDat[i, 1]
	weight += fomDat[i, 2]
center = fomDat[np.argmax(fomDat[:,2]),0]

print('FWHM, center y, weight within center +-fwhm, fom within center +-fwhm:')
print(fwhm, center, weight, fom)
print('total of all weights:')
print(np.sum(fomDat[:,2]))
print('total fom:') 
print(np.sum(fomDat[:,1]))

# total neutrons at target
plt.plot(fomDat[:,0],fomDat[:,2])
#plt.semilogy()

plt.vlines(x=center, ymin=0, ymax=np.max(fomDat[:,2]), colors='black', ls=':', lw=1, label='center')
#plt.vlines(x=[center-fwhm/2., center+fwhm/2.], ymin=0, ymax=np.max(fomDat[:,2]), colors='purple', ls='--', lw=2, label='center +- fwhm')
#plt.hlines(xmin=center-fwhm/2., xmax=center+fwhm/2., y=[np.max(fomDat[:,2]),np.max(fomDat[:,2])/2.], colors='purple', ls='--', lw=2, label='max and half max')
plt.vlines(x=[y_vals[right_idx], y_vals[left_idx]], ymin=0, ymax=np.max(fomDat[:,2]), colors='purple', ls='--', lw=2)
plt.hlines(xmin=y_vals[right_idx], xmax=y_vals[left_idx], y=[np.max(fomDat[:,2]),np.max(fomDat[:,2])/2.], colors='purple', ls='--', lw=2, label='max and half max')

plt.xlabel('y [m]')
plt.ylabel('weight')
plt.show()
