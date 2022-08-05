from setuptools import setup, find_packages

import numpy as np
import os
import matplotlib.pyplot as plt
import datetime as dt
import json
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit

requires = ["datetime", "numpy", "scipy", "matplotlib", "astropy" ]

setup(
    name='zmod',
    version='0.2.1',
    description="A package desgined to perform operations over data produced by the GRASP.",
    author_email="joao.alb@usp.br",
    url="https://github.com/zxcorr/zmod",
    license='GPLv3',
    keywords="zernike",
    package_dir={"zmod":"zmod"},
    classifiers = [
		"Development Status :: 2 - Pre-Alpha",
		"Intended Audience :: Science/Research",
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
		"Natural Language :: English",
		"Operating System :: OS Independent",
		"Topic :: Scientific/Engineering :: Image Processing",
	],
    install_requires=requires,
    packages=find_packages(
    	where=".",
    	include="*"
    ),
)
