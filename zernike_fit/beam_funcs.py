# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2022-01

Author: João Alberto
'''
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Beam Functions'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#


import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})

from . import grid_ops
from . import zernike_pols


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
			Z[i].append(zernike_pols.zernike_polynomial(r, phi, index[0], index[1]))
		
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
	
	coordinates = grid_ops.polar_grid([-1,-1,1,1], (Npoints[1],Npoints[0]), [0,0])
	
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
	and two arrays containing angular coordinates (RA, dec),
	and	reconstructs the beam for the given directions.
	RA and DEC must be centered at 0 (as how it is done in HIDE).
	'''

	coordinates = []
	assert len(RA)==len(DEC)
	for i in range(len(RA)):
		# Transforming to spherical coordinates (dec, ra) -> (theta, phi)
		theta = np.sqrt(DEC[i]**2 + RA[i]**2)
		phi = -np.arctan2(DEC[i],RA[i])
		coordinates.append([theta,phi])  
	
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



def gauss_circular_wrapper(radius=-1):
	def gauss_circular(coords, A, mu_x, mu_y, sigma):
	
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( -r2 /sigma**2 /2 )
		else: return A*np.exp(-r2 /sigma**2 /2) * (r2<=r02)
		
	return gauss_circular
	
	
	
def gauss_elliptical_wrapper(radius=-1):
	def gauss_elliptical(coords, A, mu_x, mu_y, sigma_x, sigma_y, alpha):
	
		a = (np.cos(alpha)/sigma_x)**2/2  + (np.sin(alpha)/sigma_y)**2/2
		b = -np.sin(2*alpha)/sigma_x**2/4 + np.sin(2*alpha)/sigma_y**2/4
		c = (np.sin(alpha)/sigma_x)**2/2  + (np.cos(alpha)/sigma_y)**2/2
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( a*x_0**2 + 2*b*x_0*y_0 + c*y_0**2 )
		else: return A*np.exp(-( a*x_0**2 + 2*b*x_0*y_0 + c*y_0**2 )) * (r2<=r02)
		
	return gauss_elliptical



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


