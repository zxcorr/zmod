# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2022-01

Author: João Alberto
'''
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Miscellaneous Functions.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})


def spectral_plot(Coeffs, freqs, msc=0, verbose=False, fig_path=None):
	'''
	Takes a multifrequency list of Zernike coefficients and
	plots those which are above MSC value (Minimum Signifficant
	Coefficient) *at least* in one frequency channel.
	'''
	
	assert len(Coeffs)==len(freqs), "Freqs and Coeffs lengths must be equal."
	
	# Sorting arrays
	f_idxs = np.argsort(freqs)
	freqs = freqs[f_idxs]
	Coeffs = Coeffs[f_idxs]
	
	Coeff = []
	for i in range(len(Coeffs)):
		Coeff.append([Coeffs[i,k] for k in range(len(Coeffs[i]))
					  if (True in np.array(abs(Coeffs[:,k,0])>=msc))])
	Coeff = np.array(Coeff)
	
	if verbose:
		print("\nCoefficients greater than MSC={} ({}/{})):"
			  "\n".format(msc,len(Coeff[0]),len(Coeffs[0])),Coeff[0,:,1:3])
	
	plt.figure(0,figsize=(13,7.2))
	for i in range(Coeff.shape[1]):
		plt.plot(freqs, abs(Coeff[:,i,0]), label=str(Coeff[0,i,1:3]))
		plt.legend(bbox_to_anchor=(1,1) ,loc="upper left")
	plt.title("Spectral Representation of Zernike Coefficients - "
			  "MSC={} ({}/{})".format(msc,len(Coeff[0]),len(Coeffs[0])))
	plt.xlabel("Frequency (GHz)")
	plt.ylabel(r"|C($\nu$)|")
		
	if fig_path==None: plt.show()
	else: plt.savefig(fig_path)
	
	return Coeff
	
	
	
# NOT REVISED
"""
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
	
	cols, grid_lims, grid_centers, Nps, freqs = file_ops.grasp_spher_grd_file(file_path, shift_center=True, verbose=verbose)
		
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
	
"""
