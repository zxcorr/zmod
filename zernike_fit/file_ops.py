# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2022-01

Author: João Alberto
'''
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.File Handling.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#


import numpy as np
from astropy.io import fits as pyfits


def remove_spaces(line):
	'''
	Reads a string and returns a split list without the spaces.
	'''
	
	line = line.split(" ")
	while True:
		try: line.remove("")
		except: break
		
	return line



def grasp_spher_grd_file(filename, shift_center=True, verbose=False):
	'''
	This functions opens a .grd file containing information about
	the spherical grid and returns the grid limits, the number of
	points in which the grid is divided, and the columns which
	contain the values from the field in different points.
	Ref.: GRASP manual, p. 1060-1066.
	
	"shif_center" should always be True, unless you want grid_center
	to be returned as the center provided in the .grd file (and not
	at the beam maximum).
	'''

			
	###
	# dataset strcuture:
	# 1st index: frequency selection
	# 2nd index: position selection
	# 3rd index: field component selection
	#
	# 1st comp: co real part
	# 2nd comp: co imag part
	# 3rd comp: cx real part
	# 4th comp: cx imag part
	##
	

	if verbose: print("\nCollecting data from {}".format(filename))

	grd_file = open(filename,"r")
	
	
	# Going over the headers
	
	for line in grd_file:
		# Only considers the last registered frequency.
		if "FREQUENCY:" in line: freqs = remove_spaces(line)
		if line=="++++\n": break
	grd_file.readline()
	
	freqs = np.array([float(freqs[i]) for i in range(len(freqs)) if i%2==1 and i!=0])
	
	# NSET:  number of field sets or beams.
	# ICOMP: denotes field components. If NCOMP = 3, linear co and cx components are given.
	# NCOMP: number of components. If NCOMP = 2, two field components for each point are given.
	# IGRID: control parameter of field grid type. If IGRID = 1, grid is uv type.
	config = remove_spaces(grd_file.readline())
	nset = int(config[0])
	if config[1:4] != ["3","2","1\n"]:
		if input("Different configurations found (NSET, ICOMP, NCOMP and/or IGRID). Are you sure you want to continue? (y/n)") != "y":
			print("Stopping execution")
			return False # arrumar esse return lá embaixo
	
	
	grid_center = [[]]*nset
	
	if not shift_center:
	
		for f in range(nset):
			grid_center[f] = remove_spaces(grd_file.readline())
			for i in range(2): grid_center[f][i] = float(grid_center[f][i])
		
	else: 
	
		for f in range(nset): grd_file.readline()
	
		
	grid_lim = []
	Npoints = []
	data = []
	
	for f in range(nset):
		
		line = remove_spaces(grd_file.readline())
		grid_lim.append([])
		for i in range(4): grid_lim[f].append(float(line[i]))
		
		# Getting the number of points
		
		line = remove_spaces(grd_file.readline())
		Npoints.append([int(line[0]), int(line[1])])
			
		
		# Getting data
		
		data.append([])
		
		for i in range(Npoints[f][0]*Npoints[f][1]):
		
			line = remove_spaces(grd_file.readline())
			data[f].append([float(comp) for comp in line])
			
			if i==Npoints[f][0]*Npoints[f][1]: print(data[f][i])
				
			#cols[0].append(float(data[0])) # Real part co
			#cols[1].append(float(data[1])) # Imaginary part co
			#cols[2].append(float(data[2])) # Real part cx
			#cols[3].append(float(data[3])) # Imaginary part cx
			
			
		if shift_center:
		
			co_real = np.array(data)[f][:,0] # Real co
			co_imag = np.array(data)[f][:,1] # Imag co
			co_2 = co_real**2 + co_imag**2
			
			max_point = np.where(co_2==max(co_2))[0][0]    # Selects the first maximum found
			#max_point = np.unravel_index(np.argmax(co_2, axis=None), co_2.shape)
			grid_center[f] = [0,0]
			
			grid_center[f][0] = grid_lim[f][0] + (max_point%Npoints[f][0])*(grid_lim[f][2]-grid_lim[f][0])/(Npoints[f][0]-1)
			grid_center[f][1] = grid_lim[f][1] + (max_point//Npoints[f][1])*(grid_lim[f][3]-grid_lim[f][1])/(Npoints[f][1]-1)
	
		
	return np.array(data), np.array(grid_lim), np.array(grid_center), np.array(Npoints), np.array(freqs)



def dualcut_file(file_name):
	'''
	Takes a .cut file and returns the data containing in it.
	Only works properly for (co,cx) polarisations, constant
	phi, and two field components. 
	'''
	
	data_set = 0    # There may be multiple datasets containing different
					# values of fixed phi in one single file.
					
	theta_i, delta_theta, Npoints, phi, data = [], [], [], [], []
	
	
	arq = open(file_name,"r")
	
	text_line = arq.readline()
	line_count = 1
	
	while text_line!="":
			
		# The second line contains the following information:
		#
		# V_INI - initial value (usually theta_ini)
		# V_INC - increment (usually delta_theta)
		# V_NUM - number of values in cut
		# C     - constant (usually fixed value of phi)
		# ICOMP - polarisation (if ==3, co and cx)
		# ICUT  - polar (=1, phi=const.) or conical (=2, theta=const.)
		# NCOMP - number of field components
		
		line = arq.readline().split()
		line_count+=1
		
		if line[4:7]!=["3","1","2"]:
			print("Cut configurationis different than  '3, 1, 2'. Are you sure you want to continue? (y/n)")
			if input()!="y": return False
			
		theta_i.append(float(line[0]))
		delta_theta.append(float(line[1]))
		Npoints.append(int(line[2]))
		phi.append(float(line[3]))
		
		
		data.append([])
			
		for line_count in range(line_count,line_count+Npoints[data_set]):
			
			line = arq.readline().split()
			data[data_set].append([float(comp) for comp in line])

		data_set+=1
		text_line = arq.readline()
	
	arq.close()
	
	return np.array(data), np.array(theta_i), np.array(delta_theta), np.array(Npoints), np.array(phi)



def fits_writing(filename, coefficients, freqs, R):
	'''
	This function receives a filename, a (Nfreq x Ncoeff x 3)
	coefficients array,	a frequency list, and writes the data
	into a FITS file. Each line of the coefficients array must
	correspond to a frequency-wise set of coefficients,
	containing the value of each coefficient in the first column,
	and its respective beta and alpha indices in the second and
	third columns. Frequencies must be in GHz.
	
	E.g.:
	
	coefficients = [[[5.2,0,0],   [4.8,2,0] ],
					 [-13.,1,-1], [0.5,3,3]]]
					 
	freqs = [0.98, 0.99] 
	'''

	hdu = pyfits.PrimaryHDU(coefficients)
	hdu.header["ttype1"] = "coefficients"
	hdu.header["ttype2"] = "beta"
	hdu.header["ttype3"] = "alpha"
	hdu.header.comments["ttype2"] = "radial index"
	hdu.header.comments["ttype3"] = "azimuthal index"
	hdu.header["radius"] = str(R)
	
	hdu_f = pyfits.BinTableHDU.from_columns([pyfits.Column(name="frequencies", format="D", array=freqs)])
	
	hdul = pyfits.HDUList([hdu, hdu_f])
	hdul.writeto(filename,output_verify="warn")



def fits_writing_old(filename, coefficients, freqs):
	'''
	WARNING: deprecated.
	
	This function receives a filename, a (Nfreq x Ncoeff x 3)
	coefficients array,	a frequency list, and writes the data
	into a FITS file. Each line of the coefficients array must
	correspond to a frequency-wise set of coefficients,
	containing the value of each coefficient in the first column,
	and its respective beta and alpha indices in the second and
	third columns. Frequencies must be in GHz.
	
	E.g.:
	
	coefficients = [[[5.2,0,0],   [4.8,2,0] ],
					 [-13.,1,-1], [0.5,3,3]]]
					 
	freqs = [0.98, 0.99] 
	'''

	hdul = []

	#hdu = pyfits.PrimaryHDU(freqs)
	#hdu.header["comment"] = "frequency list (GHz)"
	#hdul.append(hdu)
	
	for i in range(len(freqs)):
	
		col1 = pyfits.Column(name="coefficients", \
							 format="D", \
							 array=coefficients[i,:,0])
		col2 = pyfits.Column(name="beta", \
							 format="J", \
							 array=coefficients[i,:,1])
		col3 = pyfits.Column(name="alpha", \
							 format="J", \
							 array=coefficients[i,:,2])

		hdu = pyfits.BinTableHDU.from_columns([col1, col2, col3])
		hdu.header.comments["TTYPE2"] = "radial index"
		hdu.header.comments["TTYPE3"] = "azimuthal index"
		hdu.header.set("FREQUENCY",freqs[i])
		hdul.append(hdu)
	
	hdul = pyfits.HDUList(hdul)
	hdul.writeto(filename,output_verify="warn")
	
	return True




