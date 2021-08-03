# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2021-06

Author: João Alberto
'''

import numpy as np
from math import factorial
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})

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
	alpha, which is the azimuthal frequency.
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

	from astropy.io import fits as pyfits

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
	them into polar coordinates (r,phi).
	0 <= phi < 2*pi
	'''
	
	x = rec_coord[0] - center[0]
	y = rec_coord[1] - center[1]
	r = np.sqrt( x**2 + y**2 )
	
	
	if x!=0 and y!=0:
	
		phi = np.arctan(y/x)
		if x<0: phi = phi-np.sign(phi)*(np.pi)
	
	else:
	
		if y==0:
			if x<0: phi = np.pi
			else: phi = 0
			
		else: # x==0
			if y>0: phi = np.pi/2
			else: phi = -np.pi/2 
	
	if phi<0: phi = 2*np.pi + phi
	
	return [r,phi]



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



def uv_to_polar(rec_coord, center=[0,0]):
	'''
	This function takes uv coordinates (u,v) and transforms
	them into polar coordinates (r,phi).
	0 <= phi < 2*pi
	'''
	
	u = rec_coord[0] - center[0]
	v = rec_coord[1] - center[1]
	r = np.arctan(np.sqrt( (u**2 + v**2)/(1-u**2-v**2) ))
	
	
	if u!=0 and v!=0:
	
		phi = np.arctan(v/u)
		if u<0: phi = phi-np.sign(phi)*(np.pi)
	
	else:
	
		if v==0:
			if u<0: phi = np.pi
			else: phi = 0
			
		else: # x==0
			if v>0: phi = np.pi/2
			else: phi = -np.pi/2 
	
	if phi<0: phi = 2*np.pi + phi
	
	return [r,phi]



def polar_uv_grid(UV, center=[0,0], verbose=False):
	'''
	This function takes the limits from a uv-grid and calculates the
	coordinates of all of the points, and then returns the angle and 
	radius respective to each point in a one-dimensional structure.
	'''

	if verbose: print("\nConverting from uv to theta-phi...")
	
	coordinates = []
	
	for uv in UV:
		coordinates.append( uv_to_polar(uv, center) )
			
	return np.array(coordinates)
	


def circular_mask(data,coordinates,radius,normalisation=False, verbose=False):
	'''
	Takes a dataset and its respective coordinates, and masks all of
	the points outside a given radius.
	Return: masked dataset.
	'''
	
	if verbose: print("\nSelecting circular area around beam maximum...")
	
	# Making new arrays with only the data inside min_radius

	data = np.ma.masked_where(coordinates[:,0]>radius, data)
	
	if normalisation: coordinates[:,0] = coordinates[:,0]/radius
				
	return data, coordinates
				


def approx_zero(x, error=10**-14):
	'''
	Verifies if a given value is near zero inside
	an error range. This is mostly useful to prevent
	floating point arithmetic errors.
	'''
	return abs(x)<error
	
	

def polar_rot(coords, phi, theta, verbose=False): 
	'''
	Takes a specific point on the surface of the unit sphere
	described by two angles theta (coords[0]) and phi (coords[1])
	and performs a rotation which is equivalent to:
	
	1) a phi   (rad) rotation about the z-axis;
	2) a theta (rad) rotation about the new y-axis;
	3) a -phi  (rad) rotation about the new z-axis.
	
	Returns: final coords np.array([theta_new, phi_new])
	'''

	from scipy.spatial.transform import Rotation as R

	if verbose: print("Performing rotations to centralise polar coordinates...")
	x = np.sin(coords[0])*np.cos(coords[1])
	y = np.sin(coords[0])*np.sin(coords[1])
	z = np.cos(coords[0])
	
	vec_0 = np.array([x,y,z])
	
	R1 = R.from_rotvec(phi*np.array([0,0,1])).as_matrix() # Rotating about the z-axis
	y_axis = R1[:,1] # New y-axis
	R2 = R.from_rotvec(-theta*y_axis).as_matrix() # Rotating about the new y-awis
	# Rotating phi back
	z_axis = np.linalg.inv(np.dot(R1,R2))[:,2]
	R3 = R.from_rotvec(-phi*z_axis).as_matrix()
	
	#vec = np.dot(np.dot(R1,R2),vec_0)
	vec = np.dot(np.dot(np.dot(R1,R2),R3),vec_0)
	
	assert approx_zero(vec[0]**2+vec[1]**2+vec[2]**2-1), "Coordinates do not point to unit sphere."

	if approx_zero(vec[0]): # vec[0]==0
		if approx_zero(vec[1]): # vec[1]==0
			if   approx_zero(vec[2]-1): return np.array([0,0])     # vec[2]==1
			elif approx_zero(vec[2]+1): return np.array([np.pi,0]) # vec[2]==-1
			else: print("Error: coordinates do not point to unit sphere")
		elif vec[1]<0: phi_new = 3*np.pi/2
		else:          phi_new = np.pi/2
	
	else:
		if approx_zero(vec[1]):
			if vec[0]>0: phi_new = 0
			else:        phi_new = np.pi 	
		elif vec[1]>0: phi_new = np.arctan(vec[1]/vec[0])%np.pi
		elif vec[1]<0: phi_new = np.pi + (np.arctan(vec[1]/vec[0]))%np.pi
			
	theta_new = np.arccos(vec[2])
	
	return np.array([theta_new, phi_new])



def rotate_coords(coords, center, verbose=False):
	'''
	This function applies polar_rot to a whole array of coordinates.
	
	Returns: rotated coordinates.
	'''
	
	if verbose: print("\nRotating coordinates to centralize on beam maximum...")
	
	p=0
	N=len(coords)
	
	rot_coords = []
	
	for i in range(N):
		
		coord = coords[i]
	
		if verbose:
			if int(100*i/N)!=int(p):
				print(int(100*i/N),"%",end="\r")
			p = 100*i/N
			
		rot_coords.append(polar_rot(coord, center[1], center[0]))
		
	return np.array(rot_coords)



def sphere_separation(d1, a1, d2, a2):
    """
	great circle distance http://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas

	:param d1: dec 1
	:param a1: ra 1
	:param d2: dec 2
	:param a2:ra 2

    Taken from HIDE.
    """
    
    sin_d1, cos_d1 = np.sin(d1), np.cos(d1) #sin_cos(d1) #np.sin(d1), np.cos(d1)
    sin_d2, cos_d2 = np.sin(d2), np.cos(d2) #sin_cos(d2) #np.sin(d2), np.cos(d2)
    sin_a2a1, cosa2a1 = np.sin(a2-a1), np.cos(a2-a1) #sin_cos(a2-a1) #np.sin(a2-a1), np.cos(a2-a1)
    return (np.arctan2((np.sqrt(cos_d2**2 * sin_a2a1**2 + (cos_d1 * sin_d2 - sin_d1 * cos_d2 * cosa2a1)**2)), (sin_d1 * sin_d2 + cos_d1 * cos_d2 * cosa2a1)))




def spherical_arccos(a,b,c):
	'''
	Spherical law of cosines:
	https://en.wikipedia.org/wiki/Spherical_trigonometry#Derivation_of_the_cosine_rule
	'''
	
	if b!=0 and c!=0:
		cosA = (np.cos(a)-np.cos(b)*np.cos(c))/(np.sin(b)*np.sin(c))
	else: return 0
	
	print(a,b,c,cosA)
		
	return np.arccos(cosA)
	
	


def uv_integrate(data, thetas_matrix, Npoints, grid_lim):
	'''
	This function integrates a dataset over u and then over v.
	The data must be given in a matrix format.
	
	The jacobian in terms of theta-phi is 
	sin(theta)*d(theta)*d(phi) = cos(theta)*du*dv
	'''
	
	try: data[data.mask==1] = 0 # taking out masked elements
	except: pass
	
	from scipy.integrate import simps

	jacob = 1/np.cos(thetas_matrix)
	uv_area = (grid_lim[3]-grid_lim[1])*(grid_lim[2]-grid_lim[0])
	data = data*jacob
	data_max = np.max(data)
	
	# This method seems to have an intrinsic error
	#I_u = [simps(data[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
	#I_uv =  simps(I_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
	
	
	# 2D Monte Carlo integration
	
	N = 100 # number of tries per pixel
	total = N*Npoints[0]*Npoints[1]
	counts = 0
	samples = data_max*np.random.rand(Npoints[0],Npoints[1],N)

	for i in range(Npoints[0]):
		for j in range(Npoints[1]):
			y = data[i][j]
			counts += np.sum(samples[i,j,:]<=y)

	fraction = counts/total
	I_uv = data_max * uv_area * fraction
	
	return I_uv

			 
			   
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
	
	#indices = np.array(indices)
	p=0

	Z = []
	npix = len(coordinates)
	
	for i in range(npix):
		
		r = coordinates[i][0]
		phi = coordinates[i][1]
		Z.append([])

			
		if verbose:
			if int(100*i/npix)!=int(p):
				print(int(100*i/npix),"%",end="\r")
			p = 100*i/npix

		
		for index in indices:
			
			Z[i].append(zernike_polynomial(r, phi, index[0], index[1]))

		
	return np.array(Z)



def zernike_coeff(data, Z, indices, verbose=False):
	'''
	C = inv( Z.T * Z ) * Z.T * data
	'''

	if verbose: print("\nEstimating coefficients...")
	
	Z = np.array(Z)
	data = np.array(data)	
	Coeff = np.dot(np.dot( np.linalg.inv(np.dot(Z.T, Z)), Z.T), data)
	
	Coeff_info = []
	
	for i in range(len(Coeff)):
	
		Coeff_info.append([Coeff[i],indices[i][0],indices[i][1]])
	
	return np.array(Coeff_info)
	


def beam_reconstruction(Coeff,Npoints,verbose=False):
	'''
	This function takes Zernike coefficients ((ncoeffs,3) matrix)
	and reconstructs its correspondent beam on a unit circle inside
	a 2x2 square with the given number of pixels Npoints=[N_x,N_y].
	
	Returns a one-dimensional array containing the beam data.
	'''
	
	coordinates = polar_grid([-1,-1,1,1], (Npoints[1],Npoints[0]), [0,0])
	
	coord_mask = [0 if coordinates[i,0]>1 else 1 for i in range(Npoints[0]*Npoints[1])]
	coordinates = np.array([coordinates[i] if coord_mask[i]==1 else [np.nan, np.nan] for i in range(len(coordinates))])
	
	indices = Coeff[:,1:3]
	Z = np.array(zernike_values(coordinates, indices, verbose=verbose))
	beam = np.dot(Z,Coeff[:,0])
	
	return beam
	
	
	
def hp_beam(Coeff, nside, center, beam_radius):
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
	disk_angs = hp.pix2ang(disk_pixels)
	
	
	# rotate to get polar
	
	# reconstruct with Z
	
	# put signal in pixels
	
	
	
	# TO BE FINISHED



def radec_beam(Coeff, RA, DEC, R):
	'''
	This function receives the coefficients that describe a beam
	and two meshgrids containing angular coordinates (RA, dec),
	and	reconstructs the beam for the given directions.
	RA and DEC must be centered at 0 (as how is done in HIDE).
	'''

	coordinates = []
	for i in range(len(RA)):
		for j in range(len(DEC[0])):
			# Transforming to polar coordinates (dec, ra) -> (r, phi)
			r = np.sqrt(DEC[i][j]**2 + RA[i][j]**2)
			phi = -uv_to_polar([DEC[i][j], RA[i][j]])[1] #-np.arctan(RA[i][j]/DEC[i][j]) if DEC[i][j]!=0 else 0
			coordinates.append([r,phi])  
	
	coordinates = np.array(coordinates)
	min_dist = min(np.max(RA),np.max(DEC))
	assert min_dist>=R, "The parameters beam_azimut ({}) or beam_elevation ({}) should be greater than the beam radius ({}).".format(RA[0][-1], DEC[-1][0], R)
	coordinates[:,0] = coordinates[:,0]/R   # Normalization


	# 1) with 0's and 1's
	#coord_mask = [0 if coordinates[i,0]>1 else 1 for i in range(RA.shape[0]*RA.shape[1])]
	#coordinates = np.array([coordinates[i] if coord_mask[i]==1 else [0, 0] for i in range(len(coordinates))])
	
	# or 2) with np.masked elements
	coord_mask = np.ones(len(coordinates))
	coord_mask = np.ma.masked_where(coordinates[:,0]>1, coord_mask)
	valid_coordinates = coordinates[~coord_mask.mask]
		
	indices = Coeff[:,1:3]

	Z = zernike_values(valid_coordinates, indices, verbose=False)
	reconstructed_beam = np.dot(Z,Coeff[:,0])
	
	beam = reconstruct_masked(reconstructed_beam, coord_mask)

	return beam.reshape(RA.shape)



def reconstruct_masked(data, masked, keep_shape="False"):
	'''
	Returns data according to the mask given in "masked".
	'''
	
	if keep_shape: shape_0 = data.shape
	data.flatten()
	masked.flatten()
	
	masked_data = np.ma.array(np.zeros(masked.shape), mask=masked.mask)
	k=0
	for i in range(len(masked_data)):
		if masked_data[i] is not np.ma.masked:
			masked_data[i], k = data[k], k+1
			
	if not keep_shape: return masked_data.reshape(shape_0)
	return masked_data



def gauss_2d_wrapper(radius=-1):

	def gauss_2d(coords, A, mu_x, mu_y, sigma):
	
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( -r2 /sigma**2 /2 )
		else: return A*np.exp(-r2 /sigma**2 /2) * (r2<=r02) #np.piecewise(r2, [r2<=r02, r2>r02], [A*np.exp(-r2 /sigma**2 /2), 0])
		
	return gauss_2d



######



def dualcut_plot(data, theta_i, delta_theta, Npoints, phi, polarisation="co"):
	'''
	
	'''

	if polarisation=="co": amplitude = 10*np.log10(abs(data[:,0]**2 + data[:,1]**2))
	if polarisation=="cx": amplitude = 10*np.log10(abs(data[:,2]**2 + data[:,3]**2))

	
	theta_f = theta_i + delta_theta*(Npoints-1)
	thetas = np.arange(theta_i, theta_f, delta_theta)
	thetas = np.linspace(theta_i,theta_f,Npoints)
	# linspace introduces a tiny little error in delta_theta

	plt.figure(0)
	plt.plot(thetas, amplitude, marker="o")
	plt.title(r"$\phi$ = {}{}".format(phi,chr(176))) # chr(176) = degree sign
	plt.xlabel(r"$\theta$ ({})".format(chr(176)))
	plt.ylabel("Amplitude (dB)")
	plt.show()



def beam_plot(filename, grid="theta-phi", pol="co", title=None, verbose=False, fig_path=None, fig_filename=None,f=0):
	'''
	Takes .grd file and plots the 'co' polarisation data in it.
	If there are multiple freqs in one file, the parameter f will
	determine which frequency should be plotted.
	
	If figpath_name=None, it just shows the figure instead of saving it.
	'''

	cols, grid_lims, grid_centers, Nps,     freqs = grasp_spher_grd_file(filename, shift_center=True, verbose=verbose)
	cols, grid_lim,  grid_center,  Npoints, freq  = cols[f], grid_lims[f], grid_centers[f], Nps[f], freqs[f]
	
	if pol=="co":   data = np.sqrt(cols[:,0]**2 + cols[:,1]**2).reshape(Npoints)
	elif pol=="cx": data = np.sqrt(cols[:,2]**2 + cols[:,3]**2).reshape(Npoints)
	
	if   grid=="uv": coordinates = polar_grid(grid_lim, Npoints, grid_center, verbose)
	elif grid=="theta-phi": coordinates = polar_uv_grid(grid_lim, Npoints, grid_center, verbose)
	else: print("grid parameter must be either 'uv' or 'theta-phi'.")

	r = coordinates[:,0]
	phi = coordinates[:,1]
	XX = (r*np.cos(phi)).reshape(Npoints)
	YY = (r*np.sin(phi)).reshape(Npoints)
	
	if title==None: title="{} GHz - {} points".format(freq,Npoints)
	
	plt.figure()
	plt.pcolormesh(XX, YY, 20*np.log(data), shading="auto")
	plt.title(title)
	plt.colorbar()
	if fig_filename==None: fig_filename="freq{}_npts{}.png".format(int(1000*freq),Npoints)
	if fig_path!=None: plt.savefig(fig_path+fig_filename)
	else: plt.show()
	plt.close()
	
	return True



#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Not-so-general Functions'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



def spectral_plot(Coeffs, freqs, msc=0, verbose=False, fig_path=None):
	'''
	Takes a multifrequency list of Zernike coefficients and
	plots those which are above MSC value (Minimum Signifficant
	Coefficient) *at least* in one frequency channel.
	'''
	
	Coeff = []
	for i in range(len(Coeffs)):
		Coeff.append([Coeffs[i,k] for k in range(len(Coeffs[i])) if (True in np.array(abs(Coeffs[:,k,0])>=msc))])
	Coeff = np.array(Coeff)
	
	if verbose:
		print("\nCoefficients greater than MSC={} ({}/{})):\n".format(msc,len(Coeff[0]),len(Coeffs[0])),Coeff[0,:,1:3])

	
	plt.figure(0,figsize=(13,7.2))
	for i in range(Coeff.shape[1]):
		plt.plot(freqs, abs(Coeff[:,i,0]), label=str(Coeff[0,i,1:3]))
		plt.legend(bbox_to_anchor=(1,1) ,loc="upper left")
		
	plt.title("Spectral Representation of Zernike Coefficients - MSC={} ({}/{})".format(msc,len(Coeff[0]),len(Coeffs[0])))
	plt.xlabel("Frequency (GHz)")
	plt.ylabel(r"|C($\nu$)|")
		
	if fig_path==None: plt.show()
	else: plt.savefig(fig_path + "spectral_coefficients_{:.2f}_to_{:.2f}.png".format(freqs[0],freqs[-1]))
	
	return Coeff
	
	
	
def singlecut_interp(file_path, theta_max, phi, N_interp, title=None, show_plots=False, verbose=False):
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
		plt.title("Interpolated Region")
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
		if title==None: title = r"Dual cut Interpolation $phi$={}{}".format(np.degrees(phi),chr(176))
		plt.title(title)
		plt.show()	
	
	return d
	
	

#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Executing Tests.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#



#if __name__ == "__main__":

