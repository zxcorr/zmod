import numpy as np
import os
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})

from zmod import zernike_fit as zk


class Grasp_File(object):

	def __init__(self, filename, path="./"):
	
		assert len(filename.split(".")) == 2
		#assert not ("." in filename.replace())
		self.path = path
		self.filename = filename
		self.filepath = os.path.join(path,filename)
	


class Grd_File(Grasp_File):

	def __init__(self, path="./"):
		
		#assert filename.split(".")[1] == "grd"
		filename = "spherical_grid.grd"		
		Grasp_File.__init__(self,filename, path)
		self.beams        = None
		self.frequencies  = []
		self.grid_lims    = [] 
		self.grid_centers = [] 
		self.Nps          = []
		
		
	def extract_data(self, verbose=False):
		'''
		Returns: cols, grid_lims, grid_centers, Nps, freqs
		
		Since one .grd file may have data from more than one frequency,
		the returned arrays first index will refer to the selection of
		the desired frequency. 
		
		E.g. for .grd file containing only one frequency, one can get
		the data simply by selecting
		cols[0], grid_lims[0], grid_centers[0], Nps[0], freqs[0].
		'''
		
		beams = []
				
		cols, grid_lims, grid_centers, Nps, freqs = zk.grasp_spher_grd_file(self.filepath,
																		shift_center=True,
																		verbose=verbose)
		
		idx = np.argsort(freqs) # Sorting in frequency-increasing order
		cols = cols[idx]
		self.frequencies  = freqs[idx]
		self.grid_lims    = grid_lims[idx]
		self.grid_centers = grid_centers[idx]
		self.Nps          = Nps[idx]
		
		for f in range(len(freqs)):
			beams.append(Beam_Data(cols[f], grid_lims[f], grid_centers[f],
								   Nps[f], freqs[f]))
								   
		self.frequencies  = np.array(self.frequencies)
		self.grid_lims    = np.array(self.grid_lims) 
		self.grid_centers = np.array(self.grid_centers)
		self.Nps          = np.array(self.Nps)
		self.beams = beams




class Beam_Data(object):

	def __init__(self, cols, grid_lims, grid_center, Npoints, frequency):
	
		#assert ...
	
		self.field_components = cols.T
		self.co = np.sqrt(cols[:,0:2].T[0]**2 + cols[:,0:2].T[1]**2)
		self.cx = np.sqrt(cols[:,2:4].T[0]**2 + cols[:,2:4].T[1]**2)

		self.grid_lims   = grid_lims
		self.grid_center = grid_center  # center in UV coordinates
		self.Npoints     = Npoints
		self.frequency   = frequency
		
		# Grids must be generated through self.generate_grid().
		self.uv_grid       = None
		self.thetaphi_grid = None
		self.polar_grid    = None
		
		# valid_data refers to the data inside a specified radius.
		# These must be generated through self.circular_mask(radius).
		# Warning: this attribute will only hold information about
		# only one polarization.
		self.valid_data       = None
		self.valid_polar_grid = None
	
	
	def generate_grid(self, fast=False, verbose=False):
	
		grid_lim    = self.grid_lims
		grid_center = self.grid_center
		Npoints     = self.Npoints
	
		# v must be iterated before u
		UV = np.array([ [u,v] for v in np.linspace(grid_lim[1],grid_lim[3],Npoints[1])
							  for u in np.linspace(grid_lim[0],grid_lim[2],Npoints[0]) ])
		thetaphi_coordinates = zk.grid_ops.uv2thetaphigrid(UV, verbose=verbose)

		# Quicker option but without proper rotations
		if fast:
			recentered_UV = np.array(UV.T[0]-grid_center[0],UV.T[1]-grid_center).T
			coordinates = zk.grid_ops.uv2thetaphigrid(recentered_UV, verbose=verbose)
		else:
			thetaphi_center = zk.grid_ops.uv2thetaphi(grid_center)
			coordinates = zk.grid_ops.rotate_coords(thetaphi_coordinates, thetaphi_center,
													verbose=verbose)
		self.uv_grid       = UV
		self.thetaphi_grid = thetaphi_coordinates
		self.polar_grid    = coordinates
		
		
	def circular_mask(self, radius, pol="co", verbose=False):
	
		if pol == "co": data = self.co
		elif pol == "cx": data = self.cx
		coordinates = self.polar_grid
		data, coordinates = zk.grid_ops.circular_mask(data, coordinates, radius,
													  normalisation=True, verbose=verbose)
		self.valid_data = data[~data.mask]
		self.valid_polar_grid = coordinates[~data.mask]
		
		return data
		
		
	def interpolate(self, Npoints_new):
	
		pass
		
	
	def plot_beam(self, pol="co", title=None, fig_path=None, verbose=False):
	
		if verbose: print("\nPlotting beam...")
		fig, axs = plt.subplots(ncols=1, sharex=True, sharey=True)#, figsize=(16,4))
				
		if verbose: print("\nConverting theta-phi to Alt-Az coordinates...")
		theta = self.thetaphi_grid[:,0]
		phi = self.thetaphi_grid[:,1]
		XX = (-theta*np.cos(phi)).reshape(self.Npoints)
		YY = (theta*np.sin(phi)).reshape(self.Npoints)
		data_matrix = self.co.reshape(self.Npoints)
	
		if title==None: title="Beam Data"
		axs.set_title(title)
		c1 = axs.pcolormesh(XX,YY,20*np.log10(abs(data_matrix)), shading="auto")
		cbar = fig.colorbar(c1, ax=axs)
		cbar.set_label("Amplitude (dB)")
		axs.set_xlabel("Azimuth (rad)")
		axs.set_ylabel("Elevation (rad)")

		if fig_path==None: plt.show()
		else: plt.savefig(fig_path)
		
		
	def plot_grid(self):
	
		pass
	

