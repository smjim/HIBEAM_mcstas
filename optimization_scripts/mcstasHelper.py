#!/bin/python3
#McStas Helpers

from PIL import Image
from numpy import *#genfromtxt
import matplotlib.pyplot as plt
import tifffile

def extractMcStasData(filename):
	data = genfromtxt(filename)
	raw = []
	with open(filename, 'r') as f:
		while True:
			r = f.readline()
			if (r[0]=="#"): raw.append(r[2:-1])
			else: break
	dataHeader = {}
	for i in raw:
		t = i.split(": ")
		#print(t[0])
		try:
			dataHeader[t[0]]=t[1]	
		except:
			break
	if dataHeader['type'][:8]=="array_2d":
		s = shape(data)
		I = data[:s[0]//3]
		sigI = data[s[0]//3:s[0]//3*2]
		N = data[2*s[0]//3:]
		L = []
	elif dataHeader['type'][:8]=="array_1d":
		L = data[:,0]
		I = data[:,1]
		sigI = data[:,2]
		N = data[:,3]
	else:
		print("Unknown Data Type.")
		return -1
	return I, sigI, N, dataHeader, L

def rebin(a,I):
	newa = zeros(len(a)//I)
	for i in arange(I): newa+=a[i::I]
	return newa	

# For 1d array data
def mcstas2np(filename, statsOnly=True):
	if statsOnly:
		I, sigI, N, dataHeader, L = extractMcStasData(filename)
		arr_I, arr_sigI, arr_N, arr_L = array(I), array(sigI), array(N), array(L) 
		return arr_I, arr_sigI, arr_N, arr_L 
		
# For 2d array data
def mcstas2TIFF(filename, statsOnly=True, save=False):
	if statsOnly:
		I, sigI, N, dataHeader, L = extractMcStasData(filename)
		tiffFilename = filename+".tif"
		xmin,xmax,ymin,ymax = array(dataHeader['xylimits'].split(),dtype=float)
		s = shape(N)
		xres = s[0]/(xmax-xmin)/100
		yres = s[1]/(ymax-ymin)/100
		#print(xres,yres)
		imN = Image.fromarray(array(N,dtype=uint16),mode='I;16')
		imI = Image.fromarray(array(I,dtype=float32),mode='F')
		if save:
			imN.save(filename+"_N.tif", resolution_unit=3, x_resolution=xres, y_resolution=yres)
			imI.save(filename+"_I.tif", resolution_unit=3, x_resolution=xres, y_resolution=yres)
		else:
			return imN, imI, sigI

def show_tiff(filename, extent, xlabel, ylabel):
	# Load tiff file data
	data = tifffile.imread(filename, is_ome=False)

	# Display the image
	plt.imshow(data, extent=extent, cmap='plasma')
	plt.colorbar()
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(filename, pad=10)
	plt.show()

