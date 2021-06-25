import numpy as np
import os
import zernike_fit as zk
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})



class Grasp_File(object):



	def __init__(self, filename, path="./"):
	
		assert len(filename.split(".")) == 2
		#assert not ("." in filename.replace())
	
		self.path = path
		self.filename = filename
		self.filepath = os.path.join(path,filename)
		
		
		
	def __str__(self):
		
		return str(str.filepath)
		



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
				
		cols, grid_lims, grid_centers, Nps, freqs = zk.grasp_spher_grd_file(self.filepath, shift_center=True, verbose=verbose)
		
		idx = np.argsort(freqs) # Sorting in frequency-increasing order
		freqs, cols, grid_lims, grid_centers, Nps = freqs[idx], cols[idx], grid_lims[idx], grid_centers[idx], Nps[idx]
		
		self.frequencies  = freqs
		self.grid_lims    = grid_lims 
		self.grid_centers = grid_centers 
		self.Nps          = Nps
		
		for f in range(len(freqs)):
		
			beams.append(Beam_Data(cols[f], grid_lims[f], grid_centers[f], Nps[f], freqs[f]))
		
		self.frequencies  = np.array(self.frequencies)
		self.grid_lims    = np.array(self.grid_lims) 
		self.grid_centers = np.array(self.grid_centers)
		self.Nps          = np.array(self.Nps)
		self.beams = beams
		
		return True




#class Beam_Data(object): pass




class Beam_Data(object):



	def __init__(self, cols, grid_lims, grid_center, Npoints, frequency):
	
		#assert ...
	
		self.field_components = cols.T
		self.co = np.sqrt(cols[:,0:2].T[0]**2 + cols[:,0:2].T[1]**2)
		self.cx = np.sqrt(cols[:,2:4].T[0]**2 + cols[:,2:4].T[1]**2)
		
		self.grid_lims   = grid_lims
		self.grid_center = grid_center
		self.Npoints     = Npoints
		self.frequency   = frequency
		
		# warning
		self.uv_grid       = None
		self.thetaphi_grid = None
		self.polar_grid    = None
		
		self.valid_data       = None
		self.valid_polar_grid = None
		
	
	
	def generate_grid(self, fast=False, verbose=False):
	
		grid_lim   = self.grid_lims
		grid_center = self.grid_center
		Npoints     = self.Npoints
	
		UV = np.array([ [u,v] for v in np.linspace(grid_lim[1],grid_lim[3],Npoints[1]) for u in np.linspace(grid_lim[0],grid_lim[2],Npoints[0]) ]) # v must be iterated before u

		thetaphi_coordinates = zk.polar_uv_grid(UV, verbose=verbose)

		# Quicker option but without proper rotations
		if fast:
			coordinates = zk.polar_uv_grid(UV, center=grid_center, verbose=verbose)
		else:
			grid_center = zk.uv_to_polar(grid_center)
			coordinates = zk.rotate_coords(thetaphi_coordinates, grid_center, verbose)
			
		self.uv_grid       = UV
		self.thetaphi_grid = thetaphi_coordinates
		self.polar_grid    = coordinates
		
		return True
		
		
		
	def circular_mask(self, radius, pol="co", verbose=False):
	
		if pol == "co": data = self.co
		elif pol == "cx": data = self.cx
		
		coordinates = self.polar_grid
	
		data, coordinates = zk.circular_mask(data, coordinates, radius, normalisation=True, verbose=verbose)
		
		self.valid_data = data[~data.mask]
		self.valid_polar_grid = coordinates[~data.mask]
		
		return data
		
		
		
	def interpolate(self, Npoints_new):
	
		pass
		
	
	
	def beam_plot(self):
	
		pass
		
		
		
	def grid_plot(self):
	
		pass
	

		
		
class Beam_Zernike_Coefficients(object):

	def __init__(self, beam_data, radius, beta_max, verbose=False, pol="co", indices=None, show_plot=False, fig_path=None):
	
		# assert ...
	
		self.pol        = pol
		self.beta_max   = beta_max
		self.radius     = radius
		self.indices    = indices
		self.data       = None
		self.reconstructed_beam = None
		
		self.coeffs     = None
		self.NRMS       = None
		self.P_original = None
		self.P_reconstructed = None
		self.P_res      = None
		self.P_gauss    = None
		self.rec_power  = None
		self.res_power  = None
		
	
		# colocar em um objeto de uma nova classe Beam_Data
		# diferente da Beam_Data que ja existe
		# Beam_Data != Grd_Data
		
		if pol=="co":   self.data = np.sqrt(beam_data.co[0]**2 + beam_data.co[1]**2)
		elif pol=="cx": self.data = np.sqrt(beam_data.cx[0]**2 + beam_data.cx[1]**2)
		
		self.coeffs, self.reconstructed_beam = self.zernike_run(beam_data, verbose)
		
		self.NRMS, self.P_original, self.P_reconstructed, self.P_res= self.error_analysis(beam_data, verbose)
		self.rec_power, self.res_power = self.P_reconstructed/self.P_original, self.P_res/self.P_original
		
		
		
		if show_plot: self.plot_results(beam_data, fig_path, verbose)
		
		
		
	def zernike_run(self, beam_data, verbose=False):
	
	
		# interpolate data
	
		beam_data.generate_grid(fast=False, verbose=verbose)
		self.data = beam_data.circular_mask(self.radius, self.pol, verbose=verbose)
		
		
		if type(self.indices)==type(None):
			coeff_number = (self.beta_max+1)*(self.beta_max+2)/2
			self.indices = zk.zernike_indices(coeff_number)

		Z     = zk.zernike_values(beam_data.valid_polar_grid, self.indices, verbose)
		Coeff = zk.zernike_coeff(beam_data.valid_data, Z, self.indices, verbose) # Fitting with E field
		
		if verbose:
	
			print("Estimated coefficients ({}):".format(len(Coeff)))
			print(Coeff)
		
		masked_beam        = np.dot(Z,Coeff[:,0])
		reconstructed_beam = zk.reconstruct_masked(masked_beam, self.data)
		
		return Coeff, reconstructed_beam
		
		
		
	def error_analysis(self, beam_data, verbose=False):
	
		Npoints = beam_data.Npoints		
		grid_lim = beam_data.grid_lims
		thetaphi_coordinates = beam_data.thetaphi_grid
		
		data = self.data
		reconstructed_beam = self.reconstructed_beam

	
		# According to GRASP, the power in far-field cases is P = (k|E|)^2
		
		residuals = data - reconstructed_beam
		
		data_matrix = data.reshape((Npoints[0],Npoints[1]))
		beam_matrix = reconstructed_beam.reshape((Npoints[0],Npoints[1]))
		residuals_matrix = residuals.reshape((Npoints[0],Npoints[1]))
		

		k = 1#2*np.pi*10**9*float(freq)/c    # wavenumber


		NRMS = np.sqrt(np.sum(residuals**2)/reconstructed_beam.count())   # ~ average residual value per pixel
		if verbose: print("\nNormalised Root-Mean-Square Error (NRMS):",NRMS)


		# Reconstructed Power Fraction
		
		# Here, we integrate over u and then over v.
		# The jacobian in terms of theta-phi is 
		# sin(theta)*d(theta)*d(phi) = cos(theta)*du*dv
		
		from scipy.integrate import simps
		
		# Multiplying data for cos(theta) before integrating
		jacobian = np.cos(thetaphi_coordinates[:,0]).reshape((Npoints[0],Npoints[1]))
		data_jacob = data_matrix*jacobian
		beam_jacob = beam_matrix*jacobian
		residuals_jacob = residuals_matrix*jacobian
		
		
		#if interp_factor<=1:
		P1_u = [simps(data_jacob[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
		P_original =  k**2*simps(P1_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
		
		P2_u = [simps(beam_jacob[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
		P_reconstructed =  k**2*simps(P2_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
		
		P0_u = [simps(residuals_jacob[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
		P_res =  k**2*simps(P0_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))
		
		if verbose: print("\nOriginal Power: {}, Reconstructed Power: {}, Resisdual Power: {}".format(P_original,P_reconstructed,P_res))
		
		
		rec_power = P_reconstructed/P_original
		res_power = P_res/P_original
		
		if verbose: print("Relative reconstructed power (reconstructed/original): {}%".format(100*rec_power))
		if verbose: print("Relative residual power (residuals/original): {}%".format(100*res_power))
		
		
		return NRMS, P_original, P_reconstructed, P_res
		
	
	
	def gaussian_fit(self, coordinates, verbose=False, show_plots=False):
	
		popt, pcov = zk.gaussian_fit(self.data, coordinates, self.radius, verbose, show_plots)
		
		z_fit = gauss_2d_wrapper(radius)((x,y), popt[0],popt[1],popt[2],popt[3])
		z_res = z_fit-data


		# PLOTS
		
		if show_plots:


			# cut
			# x fixo
			
			#x0 = Npoints[1]//2-1
			#yi = x0*Npoints[0]
			#dy = int(0.5*Npoints[0]/2)
			#yf = x0+dy
			
			#print(list(y[yi:yf]).sort())
			
			#fig = plt.figure(0)
			#ax1 = fig.add_subplot(211)
			#ax1.scatter(y[yi:yf], data[yi:yf]  , c="black")
			#ax1.scatter(y[yi:yf], z_fit[yi:yf], c="b")
			#ax1.scatter(y[yi:yf,1], z_res[yi:yf], c="g")
			
			#ax2 = fig.add_subplot(212)
			#ax2.scatter(y[yi:yf], z_res[yi:yf], c="r")
			#plt.show()
			
			

			X = x.reshape(Npoints)
			Y = y.reshape(Npoints)
			data  = data.reshape(Npoints)
			z_fit = z_fit.reshape(Npoints)
			z_res = z_res.reshape(Npoints)
			
			
			plot_2d = True
			if plot_2d:
			
				fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
				
				#fig.suptitle(r"A = {:.02f}, $\mu_x$ = {:.02f}, $\mu_y$ = {:.02f}, $\sigma$ = {:.02f}".format(popt[0], popt[1], popt[2], popt[3]))
				
				c1 = axs[0].pcolormesh(X, Y, data, shading="auto")
				axs[0].set_title("Electric Field Data")
				cbar = fig.colorbar(c1, ax=axs[0])

				axs[0].set_xlim(-radius,radius)
				axs[0].set_ylim(-radius,radius)
				
				c3 = axs[1].pcolormesh(X, Y, z_fit, shading="auto")
				axs[1].set_title("Fit")
				cbar = fig.colorbar(c3, ax=axs[1])
				
				c4 = axs[2].pcolormesh(X, Y, z_res, shading="auto")
				axs[2].set_title("Residuals")
				cbar = fig.colorbar(c4, ax=axs[2])
				
				plt.show()

			else:
			
				fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
				fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
				fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})
				
				ax1.plot_surface(X, Y, data)
				ax2.plot_surface(X, Y, z_fit)
				ax3.plot_surface(X, Y, z_res)
				
				plt.show()

			
		
		
	def plot_results(self, beam_data, fig_path=None, verbose=False):
	
		coordinates = beam_data.polar_grid
		Npoints     = beam_data.Npoints
		frequency   = beam_data.frequency
		radius      = self.radius
		rec_power   = self.rec_power
		res_power   = self.res_power
		beta_max    = self.beta_max
		
		data_matrix      = self.data.reshape(Npoints)
		beam_matrix      = self.reconstructed_beam.reshape(Npoints)
		residuals_matrix = (self.data-self.reconstructed_beam).reshape(Npoints)
	
		if verbose: print("\nPlotting original and reconstructed beams...")
		fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
				
		
		if verbose: print("\nConverting polar to rectangular coordinates...")
		
		r = coordinates[:,0]
		phi = coordinates[:,1]
		
		XX = (radius*r*np.cos(phi)).reshape(Npoints)
		YY = (radius*r*np.sin(phi)).reshape(Npoints)
		
		vmin = np.min(data_matrix)
		vmax = np.max(data_matrix)
		
		
		# Original Beam
		
		axs[0].set_title("Original Beam")
		c1 = axs[0].pcolormesh(XX,YY,20*np.log10(abs(data_matrix)), shading="auto") #data_matrix
		cbar = fig.colorbar(c1, ax=axs[0])
		cbar.set_label("Amplitude (dB)")
		axs[0].set_xlim(-radius,radius)
		axs[0].set_ylim(-radius,radius)
		#axs[0].set_ylim(-1.3*radius,1.3*radius)
		
		
		# Reconstructed Beam
		
		axs[1].set_title("Reconstructed Beam (" + r"$P_{rec}$=" + "{:.2f}%)".format(100*rec_power))
		c2 = axs[1].pcolormesh(XX,YY,20*np.log10(abs(beam_matrix)), shading="auto", vmin=20*np.log10(vmin), vmax=20*np.log10(vmax)) #beam_matrix
		cbar = fig.colorbar(c2, ax=axs[1])
		cbar.set_label("Amplitude (dB)")
		
		
		# Residuals

		axs[2].set_title("Residuals (" + r"$P_{res}$=" + "{:.2f}%)".format(100*res_power))
		c3 = axs[2].pcolormesh(XX,YY,20*np.log10(abs(residuals_matrix)), shading="auto") # vmin=20*np.log10(vmin), vmax=20*np.log10(vmax)) #residuals_matrix
		cbar = fig.colorbar(c3,ax=axs[2])
		cbar.set_label("Amplitude (dB)")
		
		
		fig.suptitle("{} GHz Beam Resconstruction with R={} and  ".format(frequency,radius) + r"$\beta_{max}$=" + "{}".format(beta_max))
		
		if fig_path==None: plt.show()
		else:
			if verbose: print("\nSaving figure as:\n{}".format(fig_path))
			plt.savefig(fig_path)
			
			
			
			
class Zernike_Analysis(object):

	def __init__(self, radius, beta_max, filepath = "./", verbose=False, pol="co", indices=None, show_plot=False, fig_path=None, final_txt=None, final_fits=None):
		
		#filename        = "spherical_grid.grd"
		#self.filename   = filename
		
		self.radius     = radius
		self.beta_max   = beta_max
		self.pol        = pol
		self.beam_zernike_coefficients = []
		self.Coeffs     = []
		self.NRMSs      = []
		self.rec_powers = []
		self.res_powers = []
		
		grd_file = Grd_File(filepath)
		grd_file.extract_data(verbose)
		
		for i in range(len(grd_file.beams)):
			beam = grd_file.beams[i]
			
			self.beam_zernike_coefficients.append(Beam_Zernike_Coefficients(beam, radius, beta_max, verbose, pol, indices, show_plot, fig_path))

			self.Coeffs.append(self.beam_zernike_coefficients[i].coeffs)
			self.NRMSs.append(self.beam_zernike_coefficients[i].NRMS)
			self.rec_powers.append(self.beam_zernike_coefficients[i].rec_power)
			self.res_powers.append(self.beam_zernike_coefficients[i].res_power)

		self.Coeffs     = np.array(self.Coeffs) 
		self.NRMSs      = np.array(self.NRMSs)
		self.rec_powers = np.array(self.rec_powers)
		self.res_powers = np.array(self.res_powers)
		
		if verbose:
			print("\n\n FINAL RESULTS:")
			print(self.Coeffs)
			print("Radius:",self.radius)
			print("Frequencies (GHz):\n",grd_file.frequencies)
			print("NRMS:\n",self.NRMSs)
			print("Reconstructed Power Fractions:\n",100*self.rec_powers)
			print("Residual Power Fraction:\n",100*self.res_powers)
		
		if final_txt!=None:
			if verbose: print("Writing .txt file...")
			f_txt = open(final_txt,"w")
			f_txt.write("Input File Info:\n{}\n".format(filepath))
			f_txt.write("No. points = {}\nGrid limits = {}\nMaximum at = {}\n\n".format(grd_file.Nps, grd_file.grid_lims, grd_file.grid_centers))
			f_txt.write("Radius: {}\n".format(self.radius))
			f_txt.write("Frequencies (GHz):\n{}\n".format(grd_file.frequencies))
			f_txt.write("NRMS:\n{}\n".format(self.NRMSs))
			f_txt.write("Reconstructed Power Fractions:\n{}\n".format(self.rec_powers))
			f_txt.write("Residual Power Fraction:\n{}\n".format(self.res_powers))
			f_txt.write("Coefficients:\n{}".format(self.Coeffs))
			
		if final_fits!=None:
			if verbose: print("Writing .fits file...")
			zk.fits_writing(final_fits, self.Coeffs, grd_file.frequencies, radius)
			
		
		
		
		
	
