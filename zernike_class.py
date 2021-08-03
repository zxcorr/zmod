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
	

		
		
class Beam_Zernike_Fit(object):

	def __init__(self, beam_data, radius, beta_max, verbose=False, pol="co", indices=None, show_plot=False, fig_path=None):
	
		# assert ...
	
		self.pol       = pol
		self.beta_max  = beta_max
		self.radius    = radius
		self.indices   = indices
		self.data = None
		self.reconstructed_beam = None
		self.coeffs    = None
		self.NRMS      = None
		self.rec_power = None
		self.res_power = None
		
	
		# ideia: colocar em um objeto de uma nova classe
		# Beam_Data diferente da Beam_Data que ja existe.
		# Beam_Data != Grd_Data
		
		if pol=="co":   self.data = np.sqrt(beam_data.co[0]**2 + beam_data.co[1]**2)
		elif pol=="cx": self.data = np.sqrt(beam_data.cx[0]**2 + beam_data.cx[1]**2)
		
		self.coeffs, self.reconstructed_beam = self.zernike_run(beam_data, verbose)
		
		self.NRMS, self.res_power = self.error_analysis(beam_data, verbose)
		self.rec_power = 1-self.res_power
		
		if show_plot or not isinstance(fig_path, type(None)): self.plot_results(beam_data, fig_path, verbose)
		
		
		
	def zernike_run(self, beam_data, verbose=False):
	
		if isinstance(beam_data.thetaphi_grid, type(None)):
			beam_data.generate_grid(fast=False, verbose=verbose)
		self.data = beam_data.circular_mask(self.radius, self.pol, verbose=verbose)
		
		
		if type(self.indices)==type(None):
			coeff_number = (self.beta_max+1)*(self.beta_max+2)/2
			self.indices = zk.zernike_indices(coeff_number)

		Z     = zk.zernike_values(beam_data.valid_polar_grid, self.indices, verbose)
		Coeff = zk.zernike_coeff(beam_data.valid_data, Z, self.indices, verbose) # Fitting with |E| field
		
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

		residuals = data - reconstructed_beam
		
		data_matrix = data.reshape((Npoints[0],Npoints[1]))
		beam_matrix = reconstructed_beam.reshape((Npoints[0],Npoints[1]))
		residuals_matrix = residuals.reshape((Npoints[0],Npoints[1]))
		

		NRMS = np.sqrt(np.sum(residuals**2)/reconstructed_beam.count())   # ~ average residual value per pixel
		if verbose: print("\nNormalised Root-Mean-Square Error (NRMS):",NRMS)


		# Reconstructed Power Fraction
		
		thetas = thetaphi_coordinates[:,0].reshape(Npoints)
		P_original      = zk.uv_integrate(data_matrix**2, thetas, Npoints, grid_lim)
		#P_reconstructed = zk.uv_integrate(beam_matrix**2, thetas, Npoints, grid_lim)
		P_res           = zk.uv_integrate(residuals_matrix**2, thetas, Npoints, grid_lim)
		
		if verbose: print("\nOriginal Power: {}, Resisdual Power: {}".format(P_original,P_res))
		
		res_power = P_res/P_original
		
		#if verbose: print("Relative reconstructed power (reconstructed/original): {}%".format(100*rec_power))
		if verbose: print("Relative residual power (residuals/original): {}%".format(100*res_power))
		
		
		return NRMS, res_power
			
		
		
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
		
		self.radius     = radius
		self.beta_max   = beta_max
		self.pol        = pol
		self.beam_zernike_fit = []
		self.Coeffs     = []
		self.NRMSs      = []
		self.rec_powers = []
		self.res_powers = []
		
		grd_file = Grd_File(filepath)
		grd_file.extract_data(verbose)
		
		for i in range(len(grd_file.beams)):
			beam = grd_file.beams[i]
			
			self.beam_zernike_fit.append(Beam_Zernike_Fit(beam, radius, beta_max, verbose, pol, indices, show_plot, fig_path))

			self.Coeffs.append(self.beam_zernike_fit[i].coeffs)
			self.NRMSs.append(self.beam_zernike_fit[i].NRMS)
			self.rec_powers.append(self.beam_zernike_fit[i].rec_power)
			self.res_powers.append(self.beam_zernike_fit[i].res_power)

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
			f_txt.write("No. points = \n{}\nGrid limits = \n{}\nMaximum at = \n{}\n\n".format(grd_file.Nps, grd_file.grid_lims, grd_file.grid_centers))
			f_txt.write("Radius: \n{}\n".format(self.radius))
			f_txt.write("Frequencies (GHz):\n{}\n".format(grd_file.frequencies))
			f_txt.write("NRMS:\n{}\n".format(self.NRMSs))
			f_txt.write("Reconstructed Power Fractions:\n{}\n".format(self.rec_powers))
			f_txt.write("Residual Power Fraction:\n{}\n".format(self.res_powers))
			f_txt.write("Coefficients:\n{}".format(self.Coeffs))
			
		if final_fits!=None:
			if verbose: print("Writing .fits file...")
			#zk.fits_writing(final_fits, self.Coeffs, grd_file.frequencies, radius)
		
			from astropy.io import fits as pyfits

			hdu = pyfits.PrimaryHDU(self.Coeffs)
			hdu.header["ttype1"] = "coefficients"
			hdu.header["ttype2"] = "beta"
			hdu.header["ttype3"] = "alpha"
			hdu.header.comments["ttype2"] = "radial index"
			hdu.header.comments["ttype3"] = "azimuthal index"
			hdu.header["radius"] = str(self.radius)
			hdu.header.comments["radius"] = "angular radius (rad)"
			hdu.header["rec"] = str(self.rec_powers)
			
			import datetime as dt
			hdu.header["date"] = str(dt.date.today())

			hdu_f = pyfits.BinTableHDU.from_columns([pyfits.Column(name="frequencies", format="D", array=grd_file.frequencies)])

			hdul = pyfits.HDUList([hdu, hdu_f])
			hdul.writeto(final_fits,output_verify="warn")
			
		
		
		
class Beam_Gaussian_Fit(object):

	def __init__(self, beam_data, radius, verbose=False, show_plot=False, fig_path=None, final_txt=None, pol="co"):
	
		self.beam_data = beam_data
		self.radius = radius
		
		self.params = None
		self.cov    = None
		self.fit_data = None
		
		self.NRMS   = None
		self.rec_power = None
		self.res_power = None
		self.P_data = None
		self.P_fit  = None
		self.P_res  = None
		
		self.gauss_run(radius, verbose, show_plot, fig_path, final_txt, pol)
		
		
		
	def gauss_run(self, radius, verbose=False, show_plot=False, fig_path=None, final_txt=None, pol="co"):

		if   pol=="co": data_0 = self.beam_data.co
		elif pol=="cx": data_0 = self.beam_data.cx
		
		grid_lims = self.beam_data.grid_lims
		Npoints = self.beam_data.Npoints
		
		if isinstance(self.beam_data.thetaphi_grid, type(None)):
			self.beam_data.generate_grid(fast=False, verbose=verbose)
		coordinates = self.beam_data.thetaphi_grid
		
		polar_center = zk.uv_to_polar(self.beam_data.grid_center)
		x_center = polar_center[0]*np.cos(polar_center[1])
		y_center = polar_center[0]*np.sin(polar_center[1])
		x = coordinates[:,0]*np.cos(coordinates[:,1])
		y = coordinates[:,0]*np.sin(coordinates[:,1])
		X,Y = x.reshape(Npoints), y.reshape(Npoints)
	

		data = np.where(self.beam_data.polar_grid[:,0]**2<=radius**2, data_0, 0) 
		A_est, mu_x_est, mu_y_est, sigma_est = np.max(data_0), x_center, y_center, 0.01
		Data = data.reshape(Npoints)


		if verbose: print("\nStarting analysis...")

		from scipy.optimize import curve_fit
		popt, pcov = curve_fit(zk.gauss_2d_wrapper(radius), (x,y), data, p0=[A_est, mu_x_est, mu_y_est, sigma_est])
		self.params, self.cov = popt, pcov

		if verbose:
			print("\nParameters (A, mu_x, mu_y, sigma):\n",popt)
			print("Cov:\n",pcov)

		self.fit_data = zk.gauss_2d_wrapper(radius)((x,y), popt[0],popt[1],popt[2],popt[3])
		z_res = self.fit_data-data
		
		

		#NRMS = np.sqrt(np.sum(z_res**2)/self.fit_data.count())   # ~ average residual value per pixel
		#if verbose: print("\nNormalised Root-Mean-Square Error (NRMS):",NRMS)
		#self.NRMS = NRMS
		
		
		z_fit = self.fit_data.reshape(Npoints)
		z_res = z_res.reshape(Npoints)
		
				
		thetas = self.beam_data.thetaphi_grid[:,0].reshape(Npoints)
				
		P_data = zk.uv_integrate(Data**2,  thetas, self.beam_data.Npoints, self.beam_data.grid_lims)
		P_fit  = zk.uv_integrate(z_fit**2, thetas, self.beam_data.Npoints, self.beam_data.grid_lims)
		P_res  = zk.uv_integrate(z_res**2, thetas, self.beam_data.Npoints, self.beam_data.grid_lims)
		
		#self.P_data = P_data
		#self.P_fit  = P_fit
		#self.P_res  = P_res
		self.rec_power = P_fit/P_data
		self.res_power = P_res/P_data
		
		if verbose: print("Total power:\nP_data = {}; P_fit = {}; P_res = {}".format(P_data, P_fit, P_res))
		if verbose: print("Relative powers:\nr_fit = {}; r_res = {}".format(self.rec_power, self.res_power))
		
		
		if not isinstance(final_txt, type(None)):
			if verbose: print("Writing .txt file...")
			f_txt = open(final_txt,"w")
			#f_txt.write("Input File Info:\n{}\n".format(filepath))
			f_txt.write("No. points = \n{}\nGrid limits = \n{}\nMaximum at = \n{}\n".format(self.beam_data.Npoints, self.beam_data.grid_lims, self.beam_data.grid_center))
			f_txt.write("Radius: \n{}\n".format(self.radius))
			f_txt.write("Frequency (GHz):\n{}\n".format(self.beam_data.frequency))
			#f_txt.write("NRMS:\n{}\n".format(self.NRMS))
			f_txt.write("Reconstructed Power Fractions:\n{}\n".format(self.rec_power))
			f_txt.write("Residual Power Fraction:\n{}\n".format(self.res_power))
			f_txt.write("Params:\n{}\n".format(self.params))
			f_txt.write("Cov:\n{}\n".format(self.cov))
			
		
		
		if show_plot or not isinstance(fig_path, type(None)):
		
			fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
		
			#fig.suptitle(r"A = {:.02f}, $\mu_x$ = {:.02f}, $\mu_y$ = {:.02f}, $\sigma$ = {:.02f}".format(popt[0], popt[1], popt[2], popt[3]))
			
			c1 = axs[0].pcolormesh(X, Y, Data, shading="auto")
			axs[0].set_title("Electric Field Data")
			cbar = fig.colorbar(c1, ax=axs[0])

			axs[0].set_xlim(popt[1]-radius,popt[1]+radius)
			axs[0].set_ylim(popt[2]-radius,popt[2]+radius)
			
			c3 = axs[1].pcolormesh(X, Y, z_fit, shading="auto")
			axs[1].set_title("Fit")
			cbar = fig.colorbar(c3, ax=axs[1])
			
			c4 = axs[2].pcolormesh(X, Y, z_res, shading="auto")
			axs[2].set_title("Residuals")
			cbar = fig.colorbar(c4, ax=axs[2])
			
			if isinstance(fig_path, type(None)):
				plt.show()
			else:
				plt.savefig(fig_path)
	
	
	
	
class Gaussian_Zernike_Analysis(object):

	def __init__(self, filepath, radius, beta_max, indices=None, verbose=False, show_plot=False, gauss_fig_path=None, zernike_fig_path=None, final_txt=None, final_fits=None, pol="co"):
		
		self.radius = radius
		self.beta_max = beta_max
		self.pol = pol
		self.indices = indices
		
		self.gauss_fits   = [] # Beam_Gaussian_Fit objects
		self.zernike_fits = [] # Beam_Zernike_Fit
		
		grd_file = Grd_File(filepath)
		grd_file.extract_data()
		
		for i in range(len(grd_file.beams)):
		
			beam_data = grd_file.beams[i]
			
			self.gauss_fits.append(Beam_Gaussian_Fit(beam_data, radius, verbose, show_plot, gauss_fig_path))
			gauss_fit = self.gauss_fits[i]
			
			if   pol=="co": beam_data.co = beam_data.co - gauss_fit.fit_data
			elif pol=="cx": beam_data.cx = beam_data.cx - gauss_fit.fit_data
			
			self.zernike_fits.append(Beam_Zernike_Fit(beam_data, radius, beta_max, verbose, pol, indices, show_plot, zernike_fig_path))
		
		
		self.gauss_rec_powers = np.array([ gauss.rec_power for gauss in self.gauss_fits ])
		self.gauss_res_powers = np.array([ gauss.res_power for gauss in self.gauss_fits ])
		self.gauss_params = np.array([ gauss.params for gauss in self.gauss_fits ])
		self.gauss_covs = np.array([ gauss.cov for gauss in self.gauss_fits ])
		
		self.zernike_rec_powers = np.array([ zernike.rec_power for zernike in self.zernike_fits ])
		self.zernike_res_powers = np.array([ zernike.rec_power for zernike in self.zernike_fits ])
		self.zernike_NRMSs  = np.array([ zernike.NRMS for zernike in self.zernike_fits ])
		self.zernike_coeffs = np.array([ zernike.coeffs for zernike in self.zernike_fits ])
			
			
		if verbose and False:
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
			f_txt.write("No. points = \n{}\nGrid limits = \n{}\nMaximum at = \n{}\n".format(grd_file.Nps, grd_file.grid_lims, grd_file.grid_centers))
			f_txt.write("Radius: \n{}\n".format(self.radius))
			f_txt.write("Frequencies (GHz):\n{}\n".format(grd_file.frequencies))
			f_txt.write("Gaussian Fit:\n")
			f_txt.write("Reconstructed Power Fractions:\n{}\n".format(self.gauss_rec_powers))
			f_txt.write("Residual Power Fraction:\n{}\n".format(self.gauss_res_powers))
			f_txt.write("Params:\n{}\n".format(self.gauss_params))
			f_txt.write("Cov:\n{}\n".format(self.gauss_covs))
			f_txt.write("Zernike Residuals Fit:\n")
			f_txt.write("NRMS:\n{}\n".format(self.zernike_NRMSs))
			f_txt.write("Reconstructed Power Fractions:\n{}\n".format(self.zernike_rec_powers))
			f_txt.write("Residual Power Fraction:\n{}\n".format(self.zernike_res_powers))
			f_txt.write("Coefficients:\n{}".format(self.zernike_coeffs))
		
