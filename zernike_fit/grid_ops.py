# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2022-01

Author: João Alberto
'''
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Grid Operations.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#


import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.integrate import simps



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
	


def uv2thetaphi(uv):
	'''
	Receives the coordinates (u,v) and transforms it to the equivalent theta-phi coords
	with theta<pi/2 (upper sphere).
	'''
	u,v = uv[0],uv[1]
	theta = np.arcsin(np.sqrt(u**2+v**2))
	phi = np.arctan2(v,u)%(2*np.pi)
	return [theta,phi]



def uv2thetaphigrid(UV, verbose=False):
	'''
	This function takes a uv-grid and calculates the equivalent angles in spherical
	coordinates (theta-phi).
	
	UV: [[u1,v1], ..., [uN,vN]] - (N,2) array
	Returns: (N,2) array
	'''

	if verbose: print("\nConverting from uv to theta-phi...")
	coordinates = []
	for uv in UV:
		coordinates.append(uv2thetaphi(uv))
			
	return np.array(coordinates)
	
	

def radec2thetaphi(radec):
	'''
	This function takes a radec array and calculates the equivalent thetaphi coordinates.
	
	radec: [[ra1,dec1], ..., [raN,decN]] - (N,2) array
	Returns: (N,2) array
	'''
	
	ras = radec.T[0]
	decs = radec.T[1]
	
	thetas = np.sqrt(ras**2+decs**2)
	phis = np.arctan2(ras,decs)
	
	return np.vstack((thetas,phis)).T
	
	
	
def thetaphi2radec(thetaphi): pass
	


def circular_mask(data,coordinates,radius,normalisation=False, verbose=False):
	'''
	Takes a dataset and its respective coordinates, and masks all of
	the points outside a given radius.
	Return: masked dataset.
	'''
	
	if verbose: 
		print("\nSelecting circular area around beam maximum...")
	# Making new arrays with only the data inside min_radius
	data_masked = np.ma.masked_where(coordinates[:,0]>radius, data)
	if normalisation: 
		coordinates[:,0] = coordinates[:,0]/radius
	
	return data_masked, coordinates
				


def approx_zero(x, error=10**-14):
	'''
	Verifies if a given value is near zero inside
	an error range. This is mostly useful to prevent
	floating point arithmetic errors.
	'''
	return abs(x)<error

	

def sphere_rot(coords, phi, theta, verbose=False): 
	'''
	Takes a specific point on the surface of the unit sphere
	described by two angles theta (coords[0]) and phi (coords[1])
	and performs a rotation which is equivalent to:
	
	1) a phi   (rad) rotation about the z-axis;
	2) a theta (rad) rotation about the new y-axis;
	3) a -phi  (rad) rotation about the new z-axis.
	
	Returns: final coords np.array([theta_new, phi_new])
	'''

	if verbose: print("Performing rotations to centralise spherical coordinates...")
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
	vec = np.dot(np.dot(np.dot(R1,R2),R3),vec_0)
	
	assert approx_zero(vec[0]**2+vec[1]**2+vec[2]**2-1), "Coordinates do not point to unit sphere."

	if approx_zero(vec[0]): # vec[0]==0
		if approx_zero(vec[1]): # vec[1]==0
			if   approx_zero(vec[2]-1): return np.array([0,0])	 # vec[2]==1
			elif approx_zero(vec[2]+1): return np.array([np.pi,0]) # vec[2]==-1
			else: print("Error: coordinates do not point to unit sphere")
		elif vec[1]<0: phi_new = 3*np.pi/2
		else:		  phi_new = np.pi/2
	else:
		if approx_zero(vec[1]):
			if vec[0]>0: phi_new = 0
			else:		phi_new = np.pi 	
		elif vec[1]>0: phi_new = np.arctan(vec[1]/vec[0])%np.pi
		elif vec[1]<0: phi_new = np.pi + (np.arctan(vec[1]/vec[0]))%np.pi
			
	theta_new = np.arccos(vec[2])
	
	return np.array([theta_new, phi_new])


def sphere_rot_vec(coords, phi, theta, verbose=False):
	"""
	Vectorized: Applies spherical rotation to multiple coordinates at once.
	"""
	if verbose:
		print("Performing rotations to centralize spherical coordinates...")
	
	# Convertendo coordenadas esféricas para Cartesianas
	x = np.sin(coords[:, 0]) * np.cos(coords[:, 1])
	y = np.sin(coords[:, 0]) * np.sin(coords[:, 1])
	z = np.cos(coords[:, 0])
	vec_0 = np.stack((x, y, z), axis=-1)
	
	# Matrizes de rotação compostas
	R_phi = R.from_euler('z', phi).as_matrix()
	R_theta = R.from_euler('y', -theta).as_matrix()
	R_phi_neg = R.from_euler('z', -phi).as_matrix()
	R_total = R_phi @ R_theta @ R_phi_neg
	
	# Aplicar a matriz de rotação em massa
	vec_rot = vec_0 @ R_total.T  # Rotação de vetores (multiplicação matriz-vetor)
	
	# Converter de volta para coordenadas esféricas
	theta_new = np.arccos(vec_rot[:, 2])
	phi_new = np.arctan2(vec_rot[:, 1], vec_rot[:, 0]) % (2 * np.pi)

	return np.stack((theta_new, phi_new), axis=-1)


def rotate_coords(coords, center, verbose=False, optimize=False):
	'''
	This function applies sphere_rot to a whole array of theta-phi coordinates.
	
	Returns: rotated coordinates.
	'''
	
	if verbose: print("\nRotating coordinates to centralize on beam maximum...")

	if optimize:
		return sphere_rot_vec(coords, center[1], center[0], verbose=verbose)

	p=0
	N=len(coords)
	rot_coords = []
	for i in range(N):
		coord = coords[i]
		if verbose:
			if int(100*i/N)!=int(p):
				print(int(100*i/N),"%",end="\r")
			p = 100*i/N
		rot_coords.append(sphere_rot(coord, center[1], center[0]))
		
	return np.array(rot_coords)
	

def uv_integrate_optimized(data, thetas_matrix, grid_lim, eps=1e-6):
	'''
	Integrates a dataset over u and v using the adaptive quadrature method.
	'''
	import numpy as np
	from scipy.integrate import quad
	# Ensure masked elements are set to zero
	try:
		data[data.mask] = 0
	except AttributeError:
		pass

	# Compute Jacobian
	jacob = 1 / np.cos(thetas_matrix)
	data *= jacob

	# Define grid spacings
	u_min, v_min, u_max, v_max = grid_lim
	u = np.linspace(u_min, u_max, data.shape[0])
	v = np.linspace(v_min, v_max, data.shape[1])

	# Function to integrate over v for each u (i.e., for each row)
	def integrand_v(v, u_idx):
		v_grid = np.linspace(v_min, v_max, data.shape[1])
		return np.interp(v, v_grid, data[u_idx, :])  # Interpolate data for the given v

	# Integrate over v first (for each u), using adaptive quadrature
	def integrate_v(u_idx):
		result, _ = quad(integrand_v, v_min, v_max, args=(u_idx,), epsabs=eps, epsrel=eps)
		return result

	# Integrate over u (columns) using adaptive quadrature
	I_v = np.array([integrate_v(u_idx) for u_idx in range(data.shape[0])])

	# Now integrate over u (columns)
	def integrand_u(u, I_v):
		u_grid = np.linspace(u_min, u_max, len(I_v))
		return np.interp(u, u_grid, I_v)  # Interpolate I_v for the given u

	# Integrate along u using adaptive quadrature
	I_uv, _ = quad(integrand_u, u_min, u_max, args=(I_v,), epsabs=eps, epsrel=eps)

	return I_uv


def uv_integrate(data, thetas_matrix, grid_lim, optimize=False):
	'''
	This function integrates a dataset over u and then over v.
	The data must be given in a matrix format.
	
	The jacobian in terms of theta-phi is 
	sin(theta)*d(theta)*d(phi) = cos(theta)*du*dv
	'''

	if optimize:
		return uv_integrate_optimized(data, thetas_matrix, grid_lim)

	Npoints = data.shape
	
	try: data[data.mask==1] = 0 # removing masked elements
	except: pass

	jacob = 1/np.cos(thetas_matrix)
	uv_area = (grid_lim[3]-grid_lim[1])*(grid_lim[2]-grid_lim[0])
	data = data*jacob
	data_max = np.max(data)
	
	# Rectangular Integration over uv
	# This method seems to have an intrinsic error.
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

