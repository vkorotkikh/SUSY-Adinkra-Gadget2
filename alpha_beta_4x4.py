# ******************************************************************************
# Name:    Hardcoded Alpha & Beta matrices - times 2i
# Author:  Vadim Korotkikh
# Email:   Vadim.Korotkikh
# Date:    November 2016
# Version: 1.3
#
# Description:
#
# ******************************************************************************


# ******************************************************************************
# Begin Imports
import sys, math, itertools
import numpy as np
import numpy.matlib
from numpy import array

def main():
	pass

# ******************************************************************************
# Main() function.
# def illuminator_of_elles():
def illuminator_of_elfes():

	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are results of outerproducts inbetween
		Pauli matrices and Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.matrix([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.matrix([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.matrix([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])
	
	return [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]


# ******************************************************************************
# Calculate all the possible +1 and -1 l coefficient combinations
def sign_combimutations():
	combo_l = []

	lsigns = [1,1,1,1,1,1]
	for signs in itertools.product([-1,1], repeat=len(lsigns)):
		temp = np.array([a*sign for a,sign in zip(lsigns,signs)])
		combo_l.append(temp)
		# print(temp)
	# print(len(combo_l))
	return combo_l

# ******************************************************************************
# Execute main
def verify_billionspaper():

	"""Starting with {P4} ~ {VM1}
		According to the draft paper the tetrad order is the following:

		{ (1432), (24), (1234), (13) }  - which is equivalent to the
		following in "bracket_overbar" notation:

		{ <4123>, <1432>, <2341>, <3214> }

		Using the following tetrad as 1st test case

		<1432> <23-4-1> <3-2-14> <4-12-3> but in the above order so
		<4-12-3>, <1432>, <23-4-1>, <3-2-14>
	"""



	l1 = np.matrix([[0,0,0,1], [-1,0,0,0], [0,1,0,0], [0,0,-1,0]])
	l2 = np.matrix([[1,0,0,0], [0,0,0,1], [0,0,1,0], [0,1,0,0]])
	l3 = np.matrix([[0,1,0,0], [0,0,1,0], [0,0,0,-1], [-1,0,0,0]])
	l4 = np.matrix([[0,0,1,0], [0,-1,0,0], [-1,0,0,0], [0,0,0,1]])

	tl1 = np.transpose(l1)
	tl2 = np.transpose(l2)
	tl3 = np.transpose(l3)
	tl4 = np.transpose(l4)

	v12 = np.matmul(tl1,l2) - np.matmul(tl2,l1)
	v13 = np.matmul(tl1,l3) - np.matmul(tl3,l1)
	v14 = np.matmul(tl1,l4) - np.matmul(tl4,l1)
	v23 = np.matmul(tl2,l3) - np.matmul(tl3,l2)
	v24 = np.matmul(tl2,l4) - np.matmul(tl4,l2)
	v34 = np.matmul(tl3,l4) - np.matmul(tl4,l3)

	print("Vij - 12 result mat:")
	print(v12)
	print("")
	print("Vij - 13 result mat:")
	print(v13)
	print("")
	print("Vij - 14 result mat")
	print(v14)
	print("")
	print("Vij - 23 result mat:")
	print(v23)
	print("")
	print("Vij - 24 result mat:")
	print(v24)
	print("")
	print("Vij - 34 result mat:")
	print(v34)

	l1_13 = np.matrix([[0,0,1,0], [0,1,0,0], [1,0,0,0], [0,0,0,1]])
	l1_24 = np.matrix([[1,0,0,0], [0,0,0,1], [0,0,1,0], [0,1,0,0]])
	print(np.dot(l1_13, l1_24))

# main()
