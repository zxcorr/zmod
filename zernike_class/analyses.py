import numpy as np
import os
import matplotlib.pyplot as plt
import datetime as dt
import json
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit

from zmod import zernike_fit as zk
from zmod.zernike_class import base_classes as bclass


plt.rcParams.update({"font.family": "serif", "font.serif": "Times New Roman"})

		
class Beam_Zernike_Fit(object):

	def __init__(self, beam_data, radius, beta_max, verbose=False, pol="co",
				 indices=None, show_plot=False, fig_path=None, record_file=None):
	
		assert isinstance(beam_data, bclass.Beam_Data)
		assert radius>0
		assert isinstance(beta_max, int)
		assert isinstance(verbose,bool)
		assert pol in ["co","cx"]
		#assert indices
		assert isinstance(show_plot,bool)
		assert isinstance(fig_path,str)
	
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
		
		if pol=="co":   self.data = np.sqrt(beam_data.co[0]**2 + beam_data.co[1]**2)
		elif pol=="cx": self.data = np.sqrt(beam_data.cx[0]**2 + beam_data.cx[1]**2)
		
		self.coeffs, self.reconstructed_beam = self.zernike_run(beam_data, verbose)
		self.NRMS, self.res_power = self.error_analysis(beam_data, verbose)
		self.rec_power = 1-self.res_power
		
		if show_plot or not isinstance(fig_path, type(None)):
			self.plot_results(beam_data, fig_path, verbose)
			
		if not isinstance(record_file, type(None)):
			self.fits_record(beam_data, record_file, verbose)
		
		
		
	def zernike_run(self, beam_data, verbose=False):
	
		if isinstance(beam_data.thetaphi_grid, type(None)):
			beam_data.generate_grid(fast=False, verbose=verbose)
		self.data = beam_data.circular_mask(self.radius, self.pol, verbose=verbose)

		if type(self.indices)==type(None):
			coeff_number = (self.beta_max+1)*(self.beta_max+2)/2
			self.indices = zk.beam_funcs.zernike_indices(coeff_number)

		Z = zk.beam_funcs.zernike_values(beam_data.valid_polar_grid,
							 				  self.indices, verbose)
		Coeff = zk.beam_funcs.zernike_coeff(beam_data.valid_data, Z,
											  self.indices, verbose)
		
		if verbose:
			print("Estimated coefficients ({}):".format(len(Coeff)))
			print(Coeff)
		
		masked_beam        = np.dot(Z,Coeff[:,0])
		reconstructed_beam = zk.beam_funcs.reconstruct_masked(masked_beam, self.data)
		
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
		
		# NRMS ~ average residual value per pixel
		NRMS = np.sqrt(np.sum(residuals**2)/reconstructed_beam.count())
		if verbose: print("\nNormalised Root-Mean-Square Error (NRMS):",NRMS)

		# Reconstructed Power Fraction
		thetas = thetaphi_coordinates[:,0].reshape(Npoints)
		P_original = zk.grid_ops.uv_integrate(data_matrix**2, thetas, grid_lim)
		P_res = zk.grid_ops.uv_integrate(residuals_matrix**2, thetas, grid_lim)
		#P_reconstructed = zk.grid_ops.uv_integrate(beam_matrix**2, thetas, grid_lim)
		if verbose:
			print("\nOriginal Power: {}, Resisdual Power: {}".format(P_original, P_res))
		res_power = P_res/P_original
		if verbose:
			print("Relative residual power (residuals/original): {}%".format(
																	  100*res_power))
		
		return NRMS, res_power
			
		
		
	def plot_results(self, beam_data, fig_path=None, verbose=False):
	
		coordinates = beam_data.thetaphi_grid
		grid_center = beam_data.grid_center
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
				
		if verbose: print("\nConverting theta-phi to Alt-Az coordinates...")
		theta = coordinates[:,0]
		phi = coordinates[:,1]
		thetaphi_center = zk.grid_ops.uv2thetaphi(grid_center)
		XX = (-theta*np.cos(phi)).reshape(Npoints) # azimuth
		YY = (theta*np.sin(phi)).reshape(Npoints)  # altitude
		
		altaz_center = [-thetaphi_center[0]*np.cos(thetaphi_center[1]),
						 thetaphi_center[0]*np.sin(thetaphi_center[1])]
		azmin, azmax = altaz_center[0]-radius, altaz_center[0]+radius
		altmin, altmax = altaz_center[1]-radius, altaz_center[1]+radius
		vmin = np.min(data_matrix)
		vmax = np.max(data_matrix)
		
		axs[0].set_title("Original Beam")
		axs[1].set_title("Reconstructed Beam (" + r"$P_{rec}$=" + 
						 "{:.2f}%)".format(100*rec_power))
		axs[2].set_title("Residuals (" + r"$P_{res}$=" +
						 "{:.2f}%)".format(100*res_power))
		c1 = axs[0].pcolormesh(XX,YY,20*np.log10(abs(data_matrix)), shading="auto")
		c2 = axs[1].pcolormesh(XX,YY,20*np.log10(abs(beam_matrix)), shading="auto", 
							   vmin=20*np.log10(vmin), vmax=20*np.log10(vmax))
		c3 = axs[2].pcolormesh(XX,YY,20*np.log10(abs(residuals_matrix)), shading="auto",
							   vmin=20*np.log10(vmin), vmax=20*np.log10(vmax))
		cbar1 = fig.colorbar(c1, ax=axs[0])
		cbar1.set_label("Amplitude (dB)")
		cbar2 = fig.colorbar(c2, ax=axs[1])
		cbar2.set_label("Amplitude (dB)")
		cbar3 = fig.colorbar(c3,ax=axs[2])
		cbar3.set_label("Amplitude (dB)")
		axs[0].set_xlim(azmin, azmax)
		axs[0].set_ylim(altmin, altmax)
		axs[0].set_xlabel("Azimuth (rad)")
		axs[0].set_ylabel("Elevation (rad)")
		fig.suptitle("{} GHz Beam Resconstruction with R={} and  ".format(frequency,
					 radius) + r"$\beta_{max}$=" + "{}".format(beta_max))
		
		if fig_path==None: plt.show()
		else:
			if verbose: print("\nSaving figure as:\n{}".format(fig_path))
			plt.savefig(fig_path)
			
			
		def fits_record(self, beam_data, record_file, verbose=False):
		
			if verbose: print("Writing .json file...")
			json_data = {}
			json_data["input_file"] = "File path not given."
			json_data["Npoints"] = beam_data.Npoints.tolist()
			json_data["grid_lims"] = beam_data.grid_lims.tolist()
			json_data["maximum_pos"] = beam_data.grid_center.tolist()
			json_data["radius"] = self.radius
			json_data["frequencies"] = beam_data.frequency.tolist()
			json_data["NRMS"] = self.NRMS.tolist()
			json_data["rec_powers"] = self.rec_power.tolist()
			json_data["res_powers"] = self.res_power.tolist()
			with open(record_file+".json","w+") as f_json:
				json.dump(json_data,f_json)
				
			hdu = pyfits.PrimaryHDU(self.coeffs)
			hdu.header["ttype1"] = "coefficients"
			hdu.header["ttype2"] = "beta"
			hdu.header["ttype3"] = "alpha"
			hdu.header.comments["ttype2"] = "radial index"
			hdu.header.comments["ttype3"] = "azimuthal index"
			hdu.header["radius"] = str(self.radius)
			hdu.header.comments["radius"] = "angular radius (rad)"
			hdu.header["rec"] = str(self.rec_powers)
			hdu.header["date"] = str(dt.date.today())
			hdu_f = pyfits.BinTableHDU.from_columns([pyfits.Column(name="frequencies",
													 format="D",
													 array=np.array(beam_data.frequency))])
			hdul = pyfits.HDUList([hdu, hdu_f])
			hdul.writeto(record_file+".fits",output_verify="warn")
			
			
			
			
class Zernike_Analysis(object):

	def __init__(self, radius, beta_max, filepath = "./", verbose=False, pol="co",
				 indices=None, show_plot=False, record_file=None):
		
		self.radius     = radius
		self.beta_max   = beta_max
		self.pol        = pol
		self.beam_zernike_fit = []
		self.Coeffs     = []
		self.NRMSs      = []
		self.rec_powers = []
		self.res_powers = []
		
		grd_file = bclass.Grd_File(filepath)
		grd_file.extract_data(verbose)
		
		self.grd_file = grd_file
		
		for i in range(len(grd_file.beams)):
			beam = grd_file.beams[i]
			fig_path = record_file + ".png" # TODO: problem for multifreq grd files
			self.beam_zernike_fit.append(Beam_Zernike_Fit(beam, radius, beta_max, verbose,
										 pol, indices, show_plot, fig_path))

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
		
		if record_file!=None:
			if verbose: print("Writing .json file...")
			json_data = {}
			json_data["input_file"] = filepath
			json_data["Npoints"] = grd_file.Nps.tolist()
			json_data["grid_lims"] = grd_file.grid_lims.tolist()
			json_data["maximum_pos"] = grd_file.grid_centers.tolist()
			json_data["radius"] = self.radius
			json_data["frequencies"] = grd_file.frequencies.tolist()
			json_data["NRMS"] = self.NRMSs.tolist()
			json_data["rec_powers"] = self.rec_powers.tolist()
			json_data["res_powers"] = self.res_powers.tolist()
			with open(record_file+".json","w+") as f_json:
				json.dump(json_data,f_json)
				
			hdu = pyfits.PrimaryHDU(self.Coeffs)
			hdu.header["ttype1"] = "coefficients"
			hdu.header["ttype2"] = "beta"
			hdu.header["ttype3"] = "alpha"
			hdu.header.comments["ttype2"] = "radial index"
			hdu.header.comments["ttype3"] = "azimuthal index"
			hdu.header["radius"] = str(self.radius)
			hdu.header.comments["radius"] = "angular radius (rad)"
			hdu.header["rec"] = str(self.rec_powers)
			hdu.header["date"] = str(dt.date.today())
			hdu_f = pyfits.BinTableHDU.from_columns([pyfits.Column(name="frequencies",
													 format="D",
													 array=grd_file.frequencies)])
			hdul = pyfits.HDUList([hdu, hdu_f])
			hdul.writeto(record_file+".fits",output_verify="warn")
			
		
		
		
class Beam_Gaussian_Fit(object):

	def __init__(self, beam_data, radius, elliptical=False, verbose=False, show_plot=False,
				 fig_title=None, record_file=None, pol="co"):
	
		self.beam_data = beam_data
		self.radius = radius
		self.params = None
		self.cov    = None
		self.ellipticity = None
		self.fit_data = None
		self.NRMS   = None
		self.rec_power = None
		self.res_power = None
		self.P_data = None
		self.P_fit  = None
		self.P_res  = None
		
		self.gauss_run(beam_data, radius, elliptical=elliptical, verbose=verbose, 
					   show_plot=show_plot, fig_title=fig_title, record_file=record_file,
					   pol=pol)
		
		
	def gauss_run(self, beam_data, radius, elliptical=False, verbose=False, 
				  show_plot=False, fig_title=None, record_file=None, pol="co"):

		if   pol=="co": data_0 = self.beam_data.co
		elif pol=="cx": data_0 = self.beam_data.cx
		
		grid_lims = self.beam_data.grid_lims
		Npoints = self.beam_data.Npoints
		
		if isinstance(self.beam_data.thetaphi_grid, type(None)):
			self.beam_data.generate_grid(fast=False, verbose=verbose)
		coordinates = self.beam_data.thetaphi_grid
		
		polar_center = zk.grid_ops.uv2thetaphi(self.beam_data.grid_center)
		x_center = -polar_center[0]*np.cos(polar_center[1])
		y_center = polar_center[0]*np.sin(polar_center[1])
		x = -coordinates[:,0]*np.cos(coordinates[:,1]) # Az
		y = coordinates[:,0]*np.sin(coordinates[:,1]) # Alt
		X,Y = x.reshape(Npoints), y.reshape(Npoints)
	
		data = np.where(self.beam_data.polar_grid[:,0]**2<=radius**2, data_0, 0)
		Data = data.reshape(Npoints)

		if verbose: print("\nStarting analysis...")
		if elliptical:
			A_est, mu_x_est, mu_y_est, sigma_est, alpha_est = (np.max(data_0), x_center,
															   y_center, 0.01, 0.)
			popt, pcov = curve_fit(zk.beam_funcs.gauss_elliptical_wrapper(radius),
								   (x,y), data, 
								   p0=[A_est, mu_x_est, mu_y_est,
								   	   sigma_est, sigma_est, alpha_est],
								   bounds=([0,-np.inf,-np.inf,0,0,0],
								   		   [np.inf,np.inf,np.inf,np.inf,np.inf,2*np.pi]))
			self.params, self.cov = popt, pcov
			self.fit_data = zk.beam_funcs.gauss_elliptical_wrapper(radius)((x,y), popt[0],
																				  popt[1],
																				  popt[2],
																				  popt[3],
																				  popt[4],
																				  popt[5])
			z_res = self.fit_data-data
			if popt[4]>popt[3]:
				self.ellipticity = np.sqrt(1-(popt[3]/popt[4])**2)
			else:
				self.ellipticity = np.sqrt(1-(popt[4]/popt[3])**2)
			if verbose:
				print("\nParameters (A, mu_x, mu_y, sigma_x, sigma_x, alpha):\n",popt)
				print("Cov:\n",pcov)
				print("Ellipticity (sqrt(1-(a/b)^2)): {}".format(self.ellipticity))
				
		else:
			A_est, mu_x_est, mu_y_est, sigma_est = np.max(data_0), x_center, y_center, 0.01
			popt, pcov = curve_fit(zk.beam_funcs.gauss_circular_wrapper(radius),
								   (x,y), data, p0=[A_est, mu_x_est, mu_y_est, sigma_est])
			self.params, self.cov = popt, pcov
			if verbose:
				print("\nParameters (A, mu_x, mu_y, sigma):\n",popt)
				print("Cov:\n",pcov)

			self.fit_data = zk.beam_funcs.gauss_circular_wrapper(radius)((x,y),
																		 popt[0], 
																		 popt[1], 
																		 popt[2], 
																		 popt[3])
		z_res = self.fit_data-data
		
		# ~ average residual value per pixel
		NRMS = np.sqrt(np.sum(z_res**2)/np.count_nonzero(self.fit_data))
		if verbose: print("\nNormalised Root-Mean-Square Error (NRMS):",NRMS)
		self.NRMS = NRMS
		
		z_fit = self.fit_data.reshape(Npoints)
		z_res = z_res.reshape(Npoints)
				
		thetas = self.beam_data.thetaphi_grid[:,0].reshape(Npoints)
		P_data = zk.grid_ops.uv_integrate(Data**2, thetas, self.beam_data.grid_lims)
		P_fit  = zk.grid_ops.uv_integrate(z_fit**2, thetas, self.beam_data.grid_lims)
		P_res  = zk.grid_ops.uv_integrate(z_res**2, thetas, self.beam_data.grid_lims)
		
		self.res_power = P_res/P_data
		self.rec_power = 1-self.res_power
		
		if verbose and False: print("Total power:\nP_data = {}; "
									"P_fit = {}; P_res = {}".format(P_data, P_fit, P_res))
		if verbose: print("Relative powers:\nr_fit = {}; "
						  "r_res = {}".format(self.rec_power, self.res_power))
		

		if show_plot:
			if fig_title==None:
				popt = self.beam_gaussian_fit.params
				if elliptical:
					fig_title = (r"A = {:.02f} dB, "
								 r"($\mu_x$, $\mu_y$) = ({:.02f}, {:.02f}), "
								 r"($\sigma_x$, $\sigma_y$) = ({:.02f}, {:.02f}) arcmin, "
								 r"$\alpha$ = {:.02f} °").format(20*np.log10(popt[0]),
									 						   popt[1], popt[2],
									 						   60*np.degrees(popt[3]),
									 						   60*np.degrees(popt[4]),
									 						   np.degrees(popt[5]))
				else:
					fig_title = (r"A = {:.02f} dB, "
								 r"($\mu_x$, $\mu_y$) = ({:.02f}, {:.02f}) rad, "
								 r"$\sigma$ = {:.02f} arcmin").format(20*np.log10(popt[0]),
									 						   popt[1], popt[2],
									 						   60*np.degrees(popt[3]))
									 						   
			fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
			fig.suptitle(fig_title)
			
			c1 = axs[0].pcolormesh(X, Y, 20*np.log10(Data), shading="auto")
			c3 = axs[1].pcolormesh(X, Y, 20*np.log10(z_fit), shading="auto")
			c4 = axs[2].pcolormesh(X, Y, 20*np.log10(abs(z_res)), shading="auto")
			axs[0].set_title("Beam Data")
			axs[0].set_xlim(popt[1]-radius,popt[1]+radius)
			axs[0].set_ylim(popt[2]-radius,popt[2]+radius)
			axs[0].set_xlabel("Azimuth (rad)")
			axs[0].set_ylabel("Elevation (rad)")
			axs[1].set_title("Fit")
			axs[2].set_title("Residuals")
			cbar = fig.colorbar(c1, ax=axs[0])
			cbar = fig.colorbar(c3, ax=axs[1])
			cbar = fig.colorbar(c4, ax=axs[2])
			
			
			if isinstance(record_file, type(None)):
				plt.show()
			else:
				plt.savefig(record_file+".png")
		
	
	
	
class Gaussian_Analysis(object):

	def __init__(self, filepath, radius, elliptical=False, verbose=False, show_plot=False,
				 fig_title=None, record_file=None, pol="co"):
				 
		# TODO: problem for multifreq grd files (e.g. titles, record_files etc)
		grd_file = bclass.Grd_File(filepath)
		grd_file.extract_data(verbose)
		self.beam_data = grd_file.beams[0]
		#self.beam_gaussian_fit = []
		self.beam_gaussian_fit = Beam_Gaussian_Fit(beam_data=self.beam_data, radius=radius, elliptical=elliptical, verbose=verbose, show_plot=show_plot,
				 fig_title=fig_title, record_file=record_file, pol=pol)
										  
		if record_file!=None:
			if verbose: print("Writing .json file...")
			json_data = {}
			json_data["input_file"] = filepath
			json_data["Npoints"] = self.beam_data.Npoints.tolist()
			json_data["grid_lims"] = self.beam_data.grid_lims.tolist()
			json_data["maximum_pos"] = self.beam_data.grid_center.tolist()
			json_data["radius"] = self.beam_gaussian_fit.radius
			json_data["frequencies"] = self.beam_data.frequency.tolist()
			json_data["params"] = self.beam_gaussian_fit.params.tolist()
			if elliptical:
				json_data["params_order"] = ["A", "mu_x", "mu_y",
											 "sigma_x", "sigma_y", "alpha"]
				json_data["ellipticity"] = self.beam_gaussian_fit.ellipticity
			else:
				json_data["params_order"] = ["A", "mu_x", "mu_y", "sigma"]
			json_data["cov"] = self.beam_gaussian_fit.cov.tolist()
			json_data["NRMS"] = self.beam_gaussian_fit.NRMS.tolist()
			json_data["rec_powers"] = self.beam_gaussian_fit.rec_power.tolist()
			json_data["res_powers"] = self.beam_gaussian_fit.res_power.tolist()
			with open(record_file+".json","w+") as f_json:
				json.dump(json_data,f_json)
	
	
	
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
			
			
			
			coordinates = beam_data.polar_grid
			r = coordinates[:,0]
			phi = coordinates[:,1]
			XX = (radius*r*np.cos(phi)).reshape(beam_data.Npoints)
			YY = (radius*r*np.sin(phi)).reshape(beam_data.Npoints)
			
			if   pol=="co": beam_data.co = beam_data.co - gauss_fit.fit_data
			elif pol=="cx": beam_data.cx = beam_data.cx - gauss_fit.fit_data
			else: print("Wrong polarisation option")
			
			data_test = beam_data.co - gauss_fit.fit_data
			'''
			plt.figure(1)
			plt.pcolormesh(XX,YY, gauss_fit.fit_data.reshape(beam_data.Npoints), shading="auto")
			plt.colorbar()
			plt.figure(2)
			plt.pcolormesh(XX,YY, beam_data.co.reshape(beam_data.Npoints), shading="auto")
			plt.colorbar()
			plt.figure(3)
			plt.pcolormesh(XX,YY, (beam_data.co - gauss_fit.fit_data).reshape(beam_data.Npoints), shading="auto")
			plt.colorbar()
			plt.show()
			'''
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
		
