# Zernike Module (zmod)

This program contains functions that are suited for performing
a fit of the beam pattern data produced by GRASP, more specifically, 
it takes data from a .grd file and finds the Zernike coefficients
for this set of values, although it might be flexible enough to
be used as a general beam pattern data analyser.

The main objects to be used are the classes contained in the
submodule zmod.zernike_class.

Example:

1) Plot Zernike Polynomials

```
import numpy as np
import matplotlib.pyplot as plt
from zmod import zernike_fit as zfit

Z_beta_alpha = np.array([[1, 5, 1]]) # (beta, alpha) = (5, 1)
grid_size = [101, 101]
zernike_polynomial = zfit.beam_reconstruction(Z_beta_alpha, grid_size,
											  verbose=True).reshape(grid_size)

plt.figure(0)
plt.title(r"$Z_{beta}^{alpha}$ Zernike Polynomial".format(beta=Z_beta_alpha[0][1],
														  alpha=Z_beta_alpha[0][2]))
plt.pcolormesh(zernike_polynomial)
plt.show()
```


2) Perform Zernike Fit

```
from zmod import zclass

# Fitting options
verbose = True
radius   = 0.03 #radians
beta_max = 22
pol = "co"
show_plot   = True
record_file = "my_first_run"
filepath = "./spherical_grid.grd"

zernike_analysis = zclass.Zernike_Analysis(radius, beta_max, filepath=filepath,
										   verbose=verbose, pol=pol, indices=indices,
										   show_plot=show_plot,record_file=record_file)
										   
Coeffs = zernike_analysis.Coeffs[0]
freqs = zernike_analysis.grd_file.frequencies[0]
rec_powers = zernike_analysis.rec_powers[0]
```
