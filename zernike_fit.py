# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 20/02/2021
Last update 10/05/2021

Author: João Alberto
'''

import numpy as np
from math import factorial

c = 299792458


#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Zernike Polynomials.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def normalisation_factor(beta,alpha):

	if alpha==0: return np.sqrt(beta+1)
	
	else: return np.sqrt(2*beta+2)



def radial_polynomial(r,beta,alpha):

	if r<0 or r>1: raise ValueError("The radius value must be betweeen 0 and 1.")

	if alpha<0: alpha *= -1
	
	R = 0
	
	for s in range(int((beta-alpha)/2)+1):
	
		R += (-1)**s * factorial(beta-s) * r**(beta-2*s)/ \
			 (factorial(s) * factorial((beta+alpha)/2-s) * factorial((beta-alpha)/2-s))
	
	return R
	


def zernike_polynomial(r,theta,beta,alpha):
	'''
	Returns the numerical value of the Zernike polynomials in the 
	given radius and azimuth (polar coordinates), to chosen values
	of beta, which is the degree of the radial polynomial, and
	alpha, which is the azimutal frequency.
	'''
	
	if beta!=int(beta) or alpha!=int(alpha):
		raise TypeError("Alpha and beta must be integers.")
		
	if beta<0: return 0
	if (beta-alpha)%2!=0: return 0
	
	
	if alpha<0:
		return -normalisation_factor(beta,alpha) * \
			   radial_polynomial(r,beta,alpha)  * \
			   np.sin(alpha*theta)

	elif alpha>=0:
		return normalisation_factor(beta,alpha)* \
			   radial_polynomial(r,beta,alpha)  * \
			   np.cos(alpha*theta)
			   

			   
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.File Handling.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def remove_spaces(line):
	'''
	Reads a string and returns a split list without the spaces.
	'''
	
	line = line.split(" ")
	while True:
		try: line.remove("")
		except: break
		
	return line



def grasp_spher_grd_file(filename, shift_center=False, verbose=False):
	'''
	This functions opens a .grd file containing information about
	the spherical grid and returns the grid limits, the number of
	points in which the grid is divided, and the columns which
	contain the values from the field in different points.
	Ref.: GRASP manual, p. 1060-1066.
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
	

	if verbose: print("\nCollecting data from {}".format(file_path))

	grd_file = open(filename,"r")
	
	
	# Going over the headers
	
	for line in grd_file:
		# Only considers the last registered frequency.
		if "FREQUENCY:" in line: freqs = remove_spaces(line)
		if line=="++++\n": break
	grd_file.readline()
	
	freqs = np.array([float(freqs[i]) for i in range(len(freqs)) if i%2==1 and i!=0])
	
	config = remove_spaces(grd_file.readline())
	nset = int(config[0])
	if config[1:4] != ["3","2","1\n"]:
	
		# NSET:  number of field sets or beams.
		# ICOMP: denotes field components. If NCOMP = 3, linear co and cx components are given.
		# NCOMP: number of components. If NCOMP = 2, two field components for each point are given.
		# IGRID: control parameter of field grid type. If IGRID = 1, grid is uv type.
		
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
			
			max_point = np.where(co_2==max(co_2))[0][0]
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



def fits_writing(filename, coefficients, freqs):
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
					 
	freqs = [980, 990] 
	'''

	from astropy.io import fits as pyfits

	hdu = pyfits.PrimaryHDU(coefficients)
	hdu.header["ttype1"] = "coefficients"
	hdu.header["ttype2"] = "beta"
	hdu.header["ttype3"] = "alpha"
	hdu.header.comments["ttype2"] = "radial index"
	hdu.header.comments["ttype3"] = "azimuthal index"
	
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
					 
	freqs = [980, 990] 
	'''

	from astropy.io import fits as pyfits
	
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



#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Grid Operations.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def rect_to_polar(rec_coord, center):
	'''
	This function takes rectangular coordinates (x,y) and transforms
	them into polar coordinates (r,theta).
	0 <= theta < 2*pi
	'''
	
	x = rec_coord[0] - center[0]
	y = rec_coord[1] - center[1]
	r = np.sqrt( x**2 + y**2 )
	
	
	if x!=0 and y!=0:
	
		theta = np.arctan(y/x)
		if x<0: theta = theta-np.sign(theta)*(np.pi)
	
	else:
	
		if y==0:
			if x<0: theta = np.pi
			else: theta = 0
			
		else: # x==0
			if y>0: theta = np.pi/2
			else: theta = -np.pi/2 
	
	if theta<0: theta = 2*np.pi + theta
	
	return [r,theta]



def polar_grid(grid_lim, Npoints, center, verbose=False):
	'''
	This function takes the limits from a grid and calculates the
	coordinates of all of the points, and then returns the angle and 
	radius respective to each point in a one-dimensional structure.
	'''
	
	if verbose: print("\nGenerating grid...")
	
	X = np.linspace(grid_lim[0],grid_lim[2],Npoints[0])
	Y = np.linspace(grid_lim[1],grid_lim[3],Npoints[1])
	
	XX, YY = np.meshgrid(X,Y)
	
	coordinates = []
	
	for i in range(len(XX)):
		for j in range(len(XX[0])):
			coordinates.append( rect_to_polar([XX[i][j], YY[i][j]], center) )
			
	return np.array(coordinates)



def uv_to_polar(rec_coord, center):
	'''
	This function takes uv coordinates (u,v) and transforms
	them into polar coordinates (r,theta).
	0 <= theta < 2*pi
	'''
	
	u = rec_coord[0] - center[0]
	v = rec_coord[1] - center[1]
	r = np.arctan(np.sqrt( (u**2 + v**2)/(1-u**2-v**2) ))
	
	
	if u!=0 and v!=0:
	
		theta = np.arctan(v/u)
		if u<0: theta = theta-np.sign(theta)*(np.pi)
	
	else:
	
		if v==0:
			if u<0: theta = np.pi
			else: theta = 0
			
		else: # x==0
			if v>0: theta = np.pi/2
			else: theta = -np.pi/2 
	
	if theta<0: theta = 2*np.pi + theta
	
	return [r,theta]



def polar_uv_grid(grid_lim, Npoints, center, verbose=False):
	'''
	This function takes the limits from a uv-grid and calculates the
	coordinates of all of the points, and then returns the angle and 
	radius respective to each point in a one-dimensional structure.
	'''

	if verbose: print("\nGenerating grid...")
	
	U = np.linspace(grid_lim[0],grid_lim[2],Npoints[0])
	V = np.linspace(grid_lim[1],grid_lim[3],Npoints[1])
	
	UU, VV = np.meshgrid(U,V)
	
	coordinates = []
	
	for i in range(len(UU)):
		for j in range(len(VV[0])):
			coordinates.append( uv_to_polar([UU[i][j], VV[i][j]], center) )
			
	return np.array(coordinates)
	


def circular_mask(data,coordinates,radius,normalisation=False, verbose=False):
	'''
	Takes a dataset and its respective coordinates, and masks all of
	the points outside a given radius.
	Return: new dataset.
	'''
	
	if verbose: print("\nSelecting circular area around beam maximum...")
	
	# Making new arrays with only the data inside min_radius

	data = np.ma.masked_where(coordinates[:,0]>radius, data)
	
	if normalisation: coordinates[:,0] = coordinates[:,0]/radius
				
	return data, coordinates
				
		

			 
			   
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Beam Functions'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def zernike_indices(ncoeff):

	indices = []
	
	n = int(np.sqrt(2*ncoeff+0.25)-0.5)
	
	for beta in range(n+1):
		for alpha in range(-beta, beta+1,2):
			
			if beta==n and (beta*(beta+1)/2 + (beta+alpha)/2)==ncoeff: break
			indices.append([beta, alpha])
				
	return np.array(indices)



def zernike_values(coordinates, indices, verbose=False):
	'''
	This function takes a one-dimensional dataset containing polar 
	coordinates and the number of coefficients to be accounted for, and
	calculates the numerical values of the Zernike polynomials in the
	respective coordinates.
	
	Return: (ncoeff X len(coordinates)) matrix.
	
	ncoeff recommended: n(n+1)/2, in which n = beta_max-1
	'''
	
	if verbose: print("\nCalculating zernike polynomials in each point of the grid...")
	
	indices = np.array(indices)
	p=0

	Z = []
	npix = len(coordinates)
	
	for i in range(npix):
		
		r = coordinates[i][0]
		theta = coordinates[i][1]
		Z.append([])

			
		if verbose:
			if int(100*i/npix)!=int(p):
				print(int(100*i/npix),"%",end="\r")
			p = 100*i/npix

		
		for index in indices:
			
			Z[i].append(zernike_polynomial(r, theta, index[0], index[1]))

		
	return np.array(Z)



def zernike_coeff(data, Z, verbose=False):
	'''
	C = inv( Z.T * Z ) * Z.T * data
	'''

	if verbose: print("\nEstimating coefficients...")
	
	Z = np.array(Z)
	data = np.array(data)	
	C = np.dot(np.dot( np.linalg.inv(np.dot(Z.T, Z)), Z.T), data)
	
	C_info = []
	
	beta = 0
	alpha = 0
	
	for i in range(len(C)):
	
		C_info.append([C[i],beta,alpha])
		
		if alpha==beta:
			beta+=1
			alpha = -beta
			
		else: alpha+=2
	
	return np.array(C_info)
	


def beam_reconstruction(C,Npoints,verbose):
	'''
	This function takes Zernike coefficients and reconstructs
	its correspondent beam on a unit circle inside a 2x2 square
	with the given number of pixels Npoints=[Nx,Ny].
	Returns a one-dimensional array containing the beam data.
	'''
	
	coordinates = polar_grid([-1,-1,1,1], (Npoints[1],Npoints[0]), [0,0])
	
	coord_mask = [0 if coordinates[i,0]>1 else 1 for i in range(Npoints[0]*Npoints[1])]
	coordinates = np.array([coordinates[i] if coord_mask[i]==1 else [np.nan, np.nan] for i in range(len(coordinates))])
	
	indices = C[:,1:3]
	Z = np.array(zernike_values(coordinates, indices, verbose=verbose))
	beam = np.dot(Z,C[:,0])
	
	return beam
	
	
	
def hpbeam(C, nside, center, beam_radius):
	'''
	This function receives the coefficients that describe a beam
	and reconstructs it in a healpix format, centered at center and
	with the given radius. 
	
	center = (theta, phi) (rad)
	beam_radius: angle    (rad)
	'''

	import healpy as hp
	
	mapa = np.zeros(12*nside**2)
	
	center_vec = hp.ang2vec(center)
	disk_pixels = hp.query_disc(nside, center_vec, radius=beam_radius)
	
	
	# TO BE FINISHED


def dualcut_plot(data, theta_i, delta_theta, Npoints, phi, polarisation="co"):
	'''
	
	'''

	if polarisation=="co": amplitude = 10*np.log10(abs(data[:,0]**2 + data[:,1]**2))
	if polarisation=="cx": amplitude = 10*np.log10(abs(data[:,2]**2 + data[:,3]**2))

	
	theta_f = theta_i + delta_theta*(Npoints-1)
	thetas = np.arange(theta_i, theta_f, delta_theta)
	thetas = np.linspace(theta_i,theta_f,Npoints)
	# linspace introduces a tiny little error in delta_theta

	import matplotlib.pyplot as plt
		
	plt.figure(0)
	plt.plot(thetas, amplitude, marker="o")
	plt.title(r"$\phi$ = {}{}".format(phi,chr(176))) # chr(176) = degree sign
	plt.xlabel(r"$\theta$ ({})".format(chr(176)))
	plt.ylabel("Amplitude (dB)")
	plt.show()



#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Not-so-general Functions'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def beam_zernike_coefficients(cols, grid_lim, grid_center, Npoints, freq, radius, beta_max, verbose=False, show_plots=False):
	'''
	This functions takes a GRASP file and makes a fit using
	Zernike Polynomials within a given radius around the beam 
	maximum. It does so using only those coefficient whose
	absolute value is greeater than msc, considering a max
	value for beta, which is the degree of the polynomial.
	
	Return: estimated coefficients.
	'''
	
	#===
	# Grid Operations
	#===
		
	coordinates = polar_uv_grid(grid_lim, Npoints,grid_center,verbose)
	
	data = np.sqrt(cols[:,0]**2 + cols[:,1]**2)
	data, coordinates = circular_mask(data, coordinates, radius, normalisation=True, verbose=verbose)
	
	# Turned off
	if show_plots and False:
	
		print("\nPlotting grid...")
		
		from matplotlib import pyplot as plt
	
		fig = plt.figure(0)
		
		X = np.linspace(grid_lim[0],grid_lim[2],Npoints[0],endpoint=True)
		Y = np.linspace(grid_lim[1],grid_lim[3],Npoints[1],endpoint=True)
		XX, YY = np.meshgrid(X,Y)
		ax1 = fig.add_subplot(211)
		c1 = ax1.scatter(XX,YY)		
		
		ax2 = fig.add_subplot(212, projection="polar")
		c2 = ax2.scatter(coordinates[:,1],coordinates[:,0])
		
		plt.show()	
	
	valid_data = data[~data.mask]
	valid_coordinates = coordinates[~data.mask]
	
	
	#===
	# Zernike Coefficients
	#===

	coeff_number = (beta_max+1)*(beta_max+2)/2
	# This is flexible so that you can input the desired indices list
	indices = zernike_indices(coeff_number)

	Z = np.array(zernike_values(valid_coordinates, indices, verbose))
	C = zernike_coeff(valid_data,Z,verbose)
	
	if verbose:
	
		print("Estimated coefficients ({}):".format(len(C)))
		print(C)
	
	reconstructed_beam = np.dot(Z,C[:,0])
	
	masked_beam = np.ma.array(np.zeros(data.shape),mask=data.mask)
	k=0
	for i in range(len(masked_beam)):
		if masked_beam[i] is not np.ma.masked:
			masked_beam[i], k = reconstructed_beam[k], k+1


	#===
	# Error Analysis
	#===

	residuals = data - masked_beam

	NRME = np.sqrt(np.sum(residuals**2)/len(valid_coordinates))
	if verbose: print("\nNormalised Root-Mean-Squared Error (NRME):",NRME)

	# Reconstructed Power Fraction
	# According to GRASP, the power in far-field cases is P = (k|E|)^2
	
	data_matrix = data.reshape((Npoints[0],Npoints[1]))
	beam_matrix = masked_beam.reshape((Npoints[0],Npoints[1]))
	residuals_matrix = residuals.reshape((Npoints[0],Npoints[1]))
	
	k = 2*np.pi*10**6*float(freq)/c    # wavenumber
	
	from scipy.integrate import simps
	
	# Here, we integrate over x and then over y
	
	P1_x = [simps(data_matrix[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
	P_original =  k**2*simps(P1_x, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
	
	P2_x = [simps(beam_matrix[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
	P_reconstructed =  k**2*simps(P2_x, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
	
	P0_x = [simps(residuals_matrix[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
	P_res =  k**2*simps(P0_x, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
	
	rec_power = P_reconstructed/P_original
	res_power = P_res/P_original
	
	if verbose: print("\nRelative reconstructed power (reconstructed/original): {}%".format(100*rec_power))
	if verbose: print("Relative residual power (residuals/original): {}%".format(100*res_power))
		
	
	#===
	# Beam Plotting
	#===
	
	if show_plots:
	
		print("\nPlotting original and reconstructed beams...")
	
		from matplotlib import pyplot as plt
		fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
		#plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Helvetica"]})
				
		
		print("\nConverting polar to rectangular coordinates...")
		
		r = coordinates[:,0]
		theta = coordinates[:,1]
		
		XX = (r*np.cos(theta)).reshape(Npoints)
		YY = (r*np.sin(theta)).reshape(Npoints)
		
		
		vmin = np.min(data_matrix)
		vmax = np.max(data_matrix)
		
		# Original Beam
		
		axs[0].set_title("Original Beam")
		c1 = axs[0].pcolormesh(XX,YY,20*np.log10(abs(data_matrix)), shading="auto")
		cbar = fig.colorbar(c1, ax=axs[0])
		cbar.set_label("Amplitude (dB)")
		
		
		# Reconstructed Beam
		
		axs[1].set_title(r"Reconstructed Beam with $\beta_{max}$=" + "{}".format(beta_max))
		c2 = axs[1].pcolormesh(XX,YY,20*np.log10(abs(beam_matrix)), vmin=20*np.log10(vmin), vmax=20*np.log10(vmax), shading="auto")
		cbar = fig.colorbar(c2, ax=axs[1])
		cbar.set_label("Amplitude (dB)")
		
		
		# Residuals

		axs[2].set_title("Residuals")
		c3 = axs[2].pcolormesh(XX,YY,20*np.log10(abs(residuals_matrix)), vmin=20*np.log10(vmin), vmax=20*np.log10(vmax), shading="auto")
		cbar = fig.colorbar(c3,ax=axs[2])
		cbar.set_label("Amplitude (dB)")
		
		fig.suptitle("{} GHz Beam Resconstruction".format(freq))
		plt.show()

	return C, NRME, rec_power, res_power



def multifreq_zernike(filename, radius, beta_max, msc, verbose=False, show_plots=False):
	'''
	This function takes a .grd file (which may contain datasets from
	one or more frequencies), a radius, and returns Cs, freqs, NRMEs,
	rec_powers and res_powers.
	'''

	cols, grid_lims, grid_centers, Nps, freqs = grasp_spher_grd_file(filename, shift_center=True, verbose=False)
	
	# Sorting in frequency-increasing order
	indices = np.argsort(freqs)
	freqs, cols, grid_lims, grid_centers, Nps = freqs[indices], cols[indices], grid_lims[indices], grid_centers[indices], Nps[indices]
	
	
	Cs = []
	NRMEs = []
	rec_powers = []
	res_powers = []

	if verbose: print("\nFrequencies encountered in current file (GHz):\n {}".format(freqs))

	for f in range(len(freqs)):
	
		C, NRME, rec_power, res_power = beam_zernike_coefficients(cols[f], grid_lims[f], grid_centers[f], Nps[f], freqs[f], radius, beta_max, verbose, show_plots)
		
		Cs.append(C)
		NRMEs.append(NRME)
		rec_powers.append(rec_power)
		res_powers.append(res_power)
		
	return np.array(Cs), freqs, np.array(NRMEs), np.array(rec_powers), np.array(res_powers)
	
	
	
def multifile_zernike(file_list, radius, beta_max, msc, verbose=False, show_plots=False, final_txt=None, final_fits=None):
	'''
	Receives a list of .grd files and performs the zernike fit for
	each of them.
	
	final_txt: if None, does nothing; if str, writes a txt with the
	final results, saving it with the given name.
	final_fits: if None, does nothing; if str, writes a fits file
	with the calculated coefficients and frequencies.
	'''

	file_name = file_list[0]
	if verbose:
		print("\n\n========================================")
		print("****************************************")
		print("<> Extracting data from file ({}/{}):\n{}\n".format(1,len(file_list),file_name))
	
	Cs, freqs, NRMEs, rec_powers, res_powers = multifreq_zernike(file_name, radius, beta_max, msc, verbose, show_plots)

	for i in range(1,len(file_list)):
	
		file_name = file_list[i]
		
		if verbose:
			print("\n\n==============================================")
			print("**********************************************")
			print("<> Extracting data from file({}/{}):\n{}\n".format(i+1,len(file_list),file_name))
	
		Cs_i, freqs_i, NRMEs_i, rec_powers_i, res_powers_i = multifreq_zernike(file_name, radius, beta_max, msc, verbose, show_plots)
		
		Cs         = np.concatenate((Cs,Cs_i))
		freqs      = np.concatenate((freqs,freqs_i))
		NRMEs      = np.concatenate((NRMEs,NRMEs_i))
		rec_powers = np.concatenate((rec_powers,rec_powers_i))
		res_powers = np.concatenate((res_powers,res_powers_i))

	
	if verbose:
		print("\n\n FINAL RESULTS:")
		print(Cs)
		print("Radius:",radius)
		print("Frequencies (GHz):\n",freqs)
		print("NRME:\n",NRMEs)
		print("Reconstructed Power Fractions:\n",100*rec_powers)
		print("Residual Power Fraction:\n",100*res_powers)
		
	if final_txt!=None:
		f_txt = open(final_txt,"w")
		f_txt.write("Radius: {}\n".format(radius))
		f_txt.write("Frequencies (GHz):\n{}\n".format(freqs))
		f_txt.write("NRME:\n{}\n".format(NRMEs))
		f_txt.write("Reconstructed Power Fractions:\n{}\n".format(100*rec_powers))
		f_txt.write("Residual Power Fraction:\n{}\n".format(100*res_powers))
		f_txt.write("Coefficients:\n{}".format(Cs))
		
	if final_fits!=None:
		zk.fits_writing(final_fits, Cs, freqs)


	return Cs, freqs, NRMEs, rec_powers, res_powers



def spectral_plot(Cs, freqs, msc=0, verbose=False):
	'''
	Takes a multifrequency list of Zernike coefficients and
	plots those which are above MSC value (Minimum Signifficant
	Coefficient) *at least* in one frequency channel.
	'''
	
	C = []
	for i in range(len(Cs)):
		C.append([Cs[i,k] for k in range(len(Cs[i])) if (True in np.array(abs(Cs[:,k,0])>=msc))])
	C = np.array(C)
	
	if verbose: print(C)


	import matplotlib.pyplot as plt
	
	for i in range(C.shape[1]):
		plt.plot(freqs, abs(C[:,i,0]), label=str(C[0,i,1:3]))
		plt.legend(bbox_to_anchor=(1,1) ,loc="upper left")
		
	plt.title("Spectral Representation of Zernike Coefficients")
	plt.xlabel("Frequency (GHz)")
	plt.ylabel(r"|C($\nu$)|")
		
	plt.show()
	
	return True
	
	
	
def singlecut_interp(file, theta_max, phi, N_interp,show_plots=False, verbose=False):
	'''
	This function takes a .grd file containing 2D data from a beam and
	performs an interpolation at the given phi for a N_interp number of
	points between -theta_max < theta < theta_max.
	This was originally created only to compare the results with the single
	cut data produced by GRASP, to assert if the correct geometric convertions
	(uv -> theta-phi) were being applied.
	
	theta_max, phi -> radians
	0 < phi < pi
	'''
	
	cols, grid_lims, grid_centers, Nps, freqs = grasp_spher_grd_file(file_path, shift_center=True, verbose=verbose)
		
	indices = np.argsort(freqs)
	freqs, cols, grid_lims, grid_centers, Nps = freqs[indices], cols[indices], grid_lims[indices], grid_centers[indices], Nps[indices]
		
	grid_lim = grid_lims[0]
	grid_center = grid_centers[0]
	Npoints = Nps[0]
	
	if verbose: print("Grid center = {}, Npoints = {}".format(grid_center, Npoints))
		
		
	import matplotlib.pyplot as plt
		
	coordinates = polar_uv_grid(grid_lim, Npoints,grid_center,verbose=verbose)
	extra_points = np.array([ [theta,-np.pi+phi] for theta in np.linspace(0, theta_max, N_interp)][::-1] +\
							[ [theta,       phi] for theta in np.linspace(0, theta_max, N_interp)])
	
	X = coordinates[:,0]*np.cos(coordinates[:,1])
	Y = coordinates[:,0]*np.sin(coordinates[:,1])
	X_new = extra_points[:,0]*np.cos(extra_points[:,1])
	Y_new = extra_points[:,0]*np.sin(extra_points[:,1])
	
	if show_plots:
		plt.figure(0)	
		plt.scatter(np.degrees(X),np.degrees(Y))
		plt.scatter(np.degrees(X_new),np.degrees(Y_new),marker="x",color="red")
			
	data = np.sqrt(cols[0,:,0]**2 + cols[0,:,1]**2)
	
	if show_plots and False:
		plt.figure(1)
		plt.pcolormesh(X.reshape(Npoints),Y.reshape(Npoints),20*np.log10(abs(data.reshape(Npoints))))
		plt.colorbar()
		
		
	if verbose: print("Performing interpolation...")
		
	import scipy.interpolate
	#scipy.interpolate.interp2d(X,Y,data)(X_new,Y_new)
	data = 20*np.log10(data)
	tck = scipy.interpolate.bisplrep(X,Y,data)  # m-sqrt(2m) < s < m+sqrt(2m)
	data_interp = scipy.interpolate.bisplev(X_new,Y_new,tck)
		
	d = np.array([data_interp[i][i] for i in range(data_interp.shape[0])])
		
	if show_plots:
		theta = np.array([extra_points[i][0] if extra_points[i][1]>0 else -extra_points[i][0] for i in range(len(extra_points))])
		plt.figure(2)
		plt.plot(np.degrees(theta),abs(d),marker="o")
		plt.show()	
	
	return d
	
	

#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Executing Tests.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



#if __name__ == "__main__":

