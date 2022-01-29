# -*- coding: cp1252 -*-
'''
This program aims to make a numerical fitting of the beam
calculated by GRASP through the Zernike polynomials.

Created on 2021-02
Last update 2022-01

Author: João Alberto
'''
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.Zernike Polynomials.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#
#.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'.'#


import numpy as np
from math import factorial


def normalisation_factor(beta,alpha):

	if alpha==0: return np.sqrt(beta+1)
	
	else: return np.sqrt(2*beta+2)



def radial_polynomial(r,beta,alpha):

	if r<0 or r>1: raise ValueError("The radius value must be betweeen 0 and 1.")

	if alpha<0: alpha *= -1
	
	R = 0
	
	for s in range(int((beta-alpha)/2)+1):
	
		R += (-1)**s * factorial(beta-s) * r**(beta-2*s)/ \
			 (factorial(s) * factorial((beta+alpha)/2-s) * factorial((beta-alpha)/2-s))
	
	return R
	


def zernike_polynomial(r,theta,beta,alpha):
	'''
	Returns the numerical value of the Zernike polynomials in the 
	given radius and azimuth (polar coordinates), to chosen values
	of beta, which is the degree of the radial polynomial, and
	alpha, which is the azimutal frequency.
	'''
	
	if beta!=int(beta) or alpha!=int(alpha):
		raise TypeError("Alpha and beta must be integers.")
		
	if beta<0: return 0
	if (beta-alpha)%2!=0: return 0
	
	
	if alpha<0:
		return -normalisation_factor(beta,alpha) * \
			   radial_polynomial(r,beta,alpha)  * \
			   np.sin(alpha*theta)

	elif alpha>=0:
		return normalisation_factor(beta,alpha)* \
			   radial_polynomial(r,beta,alpha)  * \
			   np.cos(alpha*theta)
