import numpy as np
import sys
import os
import gc
import astropy.io.fits as pyfits
from zmod import zernike_class as zclass
from zmod import zernike_fit as zfit


############################################
############# INPUT DATA ###################
############################################

# Fitting options
verbose = True
radius   = 0.02  #[rad]
beta_max = 22
pol = "mueller" # "co" for co-polarized, "cx" for cross-polarized, "mueller" for full mueller matrix
indices     = None # Choose which indices should be fit
show_plot   = True
overwrite = True
shift_center = False
optimize = True # Whether to optimize the rotation of the grid and the integration method for faster execution.

mueller_names = ["m11","m12","m13","m14",
                "m21","m22","m23","m24",
                "m31","m32","m33","m34",
                "m41","m42","m43","m44"]

for mueller_name in ["m12"]: #mueller_names:
        # Path to directory containing .grd file
        grasp_dir = f"/home/gcosta/fisica/bingo_freq/mueller/{mueller_name}/" 
        # Directory to contain final results
        results_path = f"/home/gcosta/fisica/repos/zernike/zmod/teste/{mueller_name}/" 
        # File name for the final fits containing the coefficients
        final_fits = f"mueller_{mueller_name}.fits"
        title = f'mueller_{mueller_name}'
        record_file = results_path + title



        filepath = grasp_dir

        print("\n\n1) Executing analysis.\n")
        Coeffs = []
        freqs = []
        rec_powers = []

        try:
                print("\n\n\nInput File: {}\nResults File: {}\n"
                        "Title: {}\n\n".format(grasp_dir, record_file, title))
                
                zernike_analysis = zclass.Zernike_Analysis(radius, beta_max, filepath=filepath,
                                                        verbose=verbose, pol=pol, indices=indices,
                                                        show_plot=show_plot,record_file=record_file, optimize=optimize)

                print("\n\n2) Gathering data into one final fits file.\n")

                Coeffs = []
                rec_powers = []
                freqs = []

                for indice in range(len(zernike_analysis.Coeffs)):
                        print(f"Coefficients for frequency {indice}: {zernike_analysis.Coeffs[indice]}")
                        print(f"Reconstructed power for frequency {indice}: {zernike_analysis.rec_powers[indice]}")
                        Coeffs.append(zernike_analysis.Coeffs[indice])
                        freqs.append(zernike_analysis.grd_file.frequencies[indice])
                        rec_powers.append(zernike_analysis.rec_powers[indice])
                # Deleting data to avoid overflow
                del zernike_analysis
                gc.collect()


                ############################################
                ############ GATHERING DATA ################
                ############################################


                print("\n\n2) Gathering data into one final fits file.\n")

                Coeffs = np.array(Coeffs)
                rec_powers = np.array(rec_powers)
                freqs = np.array(freqs)

                # Creating FITS file

                hdu = pyfits.PrimaryHDU(Coeffs)
                hdu.header["ttype1"] = "coefficients"
                hdu.header["ttype2"] = "beta"
                hdu.header[" ttype3"] = "alpha"
                hdu.header.comments["ttype2"] = "radial index"
                hdu.header.comments["ttype3"] = "azimuthal index"
                hdu.header["radius"] = radius
                hdu.header.comments["radius"] = "angular radius (rad)"
                hdu.header["date"] = None

                hdu_f = pyfits.BinTableHDU.from_columns([pyfits.Column(name="frequencies",
                                                                                                format="D",
                                                                                                array=freqs)])
                hdu_r = pyfits.BinTableHDU.from_columns([pyfits.Column(name="rec_powers",
                                                                                                format="D",
                                                                                                array=rec_powers)])
                hdul = pyfits.HDUList([hdu, hdu_f, hdu_r])
                hdul.writeto(os.path.join(results_path+"{}/".format(mueller_name),final_fits),output_verify="warn",overwrite=overwrite)


        except Exception as e:
                print(f"Erro inesperado: {e}")
                print("\nPróximo...")
