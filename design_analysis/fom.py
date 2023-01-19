import numpy as np
import matplotlib.pyplot as plt

# p, x, y, r_xy, tof 
data = np.loadtxt('data/TOF_target_radial.dat')

radii = np.linspace(0, 1, 100)
fom_dat = []
sum_dat = []
for r in radii:
	fom_tmp = 0
	sum_tmp = 0
	for line in data: 
		weight = line[0]	# neutron weight
		radius = line[3]	# detection radius (xy)
		time_f = line[4]	# time of flight [s]
		if (radius < r):
			fom_tmp += (time_f**2)*weight
			sum_tmp += weight
	fom_dat.append(fom_tmp)
	sum_dat.append(sum_tmp)
	print(str(r)+' done')
fom_dat = np.array(fom_dat)
sum_dat = np.array(sum_dat)

plt.plot(radii,fom_dat) 

plt.xlabel('radii [m]')
plt.ylabel('FOM sum')
plt.show()

plt.plot(radii,sum_dat) 
plt.xlabel('radii [m]')
plt.ylabel('weight sum')
plt.show()

plt.plot(radii,np.sqrt(fom_dat/sum_dat)) 
plt.xlabel('radii [m]')
plt.ylabel('avg TOF [s]')
plt.show()
