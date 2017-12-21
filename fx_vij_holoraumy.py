# ******************************************************************************
# Name:    Functions for computing Vij matrices
# Author:  Vadim Korotkikh
# Date:    September 2017
# Description:  Scripts for calculating Vij matrices for each one of 36864
# unique Adinkra tetrads.
#
# ******************************************************************************

# Library Imports
import itertools
import numpy as np

#>******************************************************************************
def alphas_betas():
	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are tensor products of Pauli spin matrices
 		Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""
	# Redoing the matrices with np.array
	# Alphas
	# alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha1i	= np.array([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.array([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.array([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.array([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.array([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.array([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])

	tlist = [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]
	return tlist


#>******************************************************************************
def calc_fermi_vij(adinkra_list):
	""" Calculates the (6) Fermionic Holoraumy matrices for each Adinkra in the
		adinkra_list.
		adinkra_list - list containing n lists - each nested list contains 6 tuples
		( int value - L matrix index, 4x4 np.array - L matrix)
	"""

	loc_adinkras = []
	if not isinstance(adinkra_list, list):
		loc_adinkras = [adinkra_list]
	else:
		print("Adinkra-list length: ", len(adinkra_list))
		loc_adinkras = adinkra_list
		# pass

	''' Store Fermionic Holoraumy matrices '''
	vij_fermi	= []
	ij_inds		= list(itertools.combinations([0,1,2,3], 2))
	for idink in loc_adinkras:
		temp_six = []
		for ijtup in ij_inds:
			limat 		= idink[ijtup[0]][1]
			ljmat 		= idink[ijtup[1]][1]
			rimat		= np.transpose(limat)
			rjmat		= np.transpose(ljmat)
			""" Vij eq from 1601.00 (3.2) """
			""" Probably needs 1/2	"""
			holo_mat	= np.divide((np.dot(rimat, ljmat) - np.dot(rjmat, limat)),2)
			# holo_mat	= np.divide((np.dot(rjmat, limat) - np.dot(rimat, ljmat)),2)
			temp_six.append(holo_mat)
		vij_fermi.append(temp_six)
	return vij_fermi

#>******************************************************************************
def calc_ultrafermi_vij(adinkra_list):
	""" Calculates the 24 Fermionic Holoraumy matrices for each Adinkra in the
		adinkra_list. Not coded to handle more than that.
		adinkra_list - list containing n lists - each nested list contains 6 tuples
		( int value - L matrix index, 4x4 np.array - L matrix)
	"""

	loc_adinkras = []
	if not isinstance(adinkra_list, list):
		loc_adinkras = [adinkra_list]
	else:
		print("Adinkra-list length: ", len(adinkra_list))
		loc_adinkras = adinkra_list

	vij_ultrafermi	= []
	ij_inds = list(itertools.permutations([0,1,2,3], 2)) #  24 possible permutations
	for idink in loc_adinkras:
		temp_dozen = []
		for ijtup in ij_inds:
			limat 		= idink[ijtup[0]][1]
			ljmat 		= idink[ijtup[1]][1]
			rimat		= np.transpose(limat)
			rjmat		= np.transpose(ljmat)
			""" Vij eq from 1601.00 (3.2) - has the 1/2 factor"""
			holo_mat	= np.divide((np.dot(rimat, ljmat) - np.dot(rjmat, limat)),2)
			temp_dozen.append(holo_mat)
		vij_ultrafermi.append(temp_dozen)
	return vij_ultrafermi


#>******************************************************************************
def fermionic_holomats(adinkra):
	return calc_fermi_vij(adinkra)


# ******************************************************************************
def bosonic_holomats(adinkra_list):
	""" Calculates the (6) Bosonic Holoraumy matrices for each Adinkra in the
		adinkra_list.
		adinkra_list - list containing n lists - each nested list contains 6 tuples
		( int value - L matrix index, 4x4 np.array - L matrix)
	"""
	loc_adinkras = []
	if not isinstance(tetrad_list, list):
		loc_adinkras = [adinkra_list]
	else:
		pass
	''' Store n Vij bosonic matrices in vij_bosonic	'''
	vij_bosonic	= []
	ij_indices = list(itertools.combinations([0,1,2,3], 2))

	for idink in loc_adinkras:
		temp_six	= []
		for ijtup in ij_indices:

			limat	= idink[ijtup[0]]
			ljmat	= idink[ijtup[1]]
			rimat	= np.transpose(limat)
			rjmat	= np.transpose(ljmat)

			""" Vij bosonic eq	1501.00101 3.5	"""
			holo_mat = np.divide((np.dot(limat, rjmat) - np.dot(ljmat, rimat)),2)
			temp_six.append(holo_mat)
		vij_bosonic.append(temp_six)

	return vij_bosonic


#>******************************************************************************
def calc_vij_alphabeta(main_tetrad_list):
	"""	Note: main_tetrad_list is a list of lists,
		with each list containing four tuples, with tuples being an integer
		index (0 to 383) of the L matrix and the 4x4 np.array representing
		the L sign permutation matrix
	"""

	vij_possibilities = alphas_betas()
	ij_indices	= list(itertools.combinations([0,1,2,3], 2))
	vij_mats	= []
	debug		= 0
	print("							")
	print("Calculating Vij matrices")
	for ti, teti in enumerate(main_tetrad_list):
		vij_tempset = []
		for ijtup in ij_indices:
			limat 		= teti[ijtup[0]][1]
			ljmat 		= teti[ijtup[1]][1]
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= ij_temp
			tr_limat, tr_ljmat	= np.transpose(limat), np.transpose(ljmat)
			""" Vij eq from 1601.00 (3.2) """
			temp_mat	= np.dot(tr_limat, ljmat) - np.dot(tr_ljmat, limat)
			""" Compare against the 6 possible matrix solutions """
			tf_bool = 0
			# Compare against the 6 possible matrix solutions
			for xi, ijx in enumerate(vij_possibilities):
				ijx_neg = np.multiply(ijx, -1)
				if np.array_equal(temp_mat, ijx):
					tf_bool = 1
					tmint = np.int(1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						# vij_tempset.append([tmp_str, ijstr, tmint])
						vij_tempset.append((tmp_str, ijstr, tmint))
						continue
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						# vij_tempset.append([tmp_str, ijstr, tmint])
						vij_tempset.append((tmp_str, ijstr, tmint))
						continue
				elif np.array_equal(temp_mat, ijx_neg):
					tf_bool = 1
					tmint = np.int(-1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						# vij_tempset.append([tmp_str, ijstr, tmint])
						vij_tempset.append((tmp_str, ijstr, tmint))
						continue
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						# vij_tempset.append([tmp_str, ijstr, tmint])
						vij_tempset.append((tmp_str, ijstr, tmint))
						continue
				else:
					if tf_bool == 0 and xi >= 5:
						if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
							print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
							print("Anomaly found:",xi)
							print(temp_mat)
			tf_bool = 0
		vij_mats.append(vij_tempset)

	return vij_mats

# ******************************************************************************
# Do the final Vij calculation
def gadget_one_abcalc(vij_fermi_list):
	""" Remember that the vij_fermi_list is a list of lists,
		with each list containing four tuples, with tuples being an integer
		index (0 to 383) of the L matrix and the 4x4 np.array representing
		the L sign permutation matrix
	"""

	vij_possibil = alphas_betas()
	ij_indices	= list(itertools.combinations([0,1,2,3], 2))
	coef_vals	= []

	# anomaly_switch  = 0
	debug			= 0

	""" Here each Adinkra in vij_fermi_list is a tetrad, ie 4 L matrices
		for ind, adinkra in enumerate(vij_fermi_list):	"""
	for ind, vijset in enumerate(vij_fermi_list):

		vij_tempset 	= []
		vij_tempcoefs 	= []
		tf_bool = 0
		for vijtemp in vijset:
			# vijtemp = np.multiply(vijmat,2)
			for xi, ijx in enumerate(vij_possibil):
				ijx_neg = np.multiply(ijx, -1)
				# print(xi)
				if np.array_equal(vijtemp, ijx):
					tf_bool = 1
					# if debug:
						# print("*************$$$$$$$$$$$$$$$$$$ ")
						# print("l-solution found:", ijx)
						# print(ijx)
					tmint = np.int(1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						vij_tempset.append((tmint, tmp_str))
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						vij_tempset.append((tmint, tmp_str))
				elif np.array_equal(vijtemp, ijx_neg):
					tf_bool = 1
					# if debug:
						# print("*************$$$$$$$$$$$$$$$$$$ ")
						# print("l-solution found: ", ijx_neg)
						# print(ijx_neg)
					tmint = np.int(-1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						vij_tempset.append((tmint, tmp_str))
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						vij_tempset.append((tmint, tmp_str))
				else:
					if tf_bool == 0 and xi >= 5:
						if not(np.array_equal(vijtemp, ijx)) or not np.array_equal(vijtemp, ijx_neg):
							print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
							print("Anomaly found:", ind)
							print(vijtemp)
							# anomaly_switch = 1
			tf_bool = 0

		coef_vals.append(vij_tempset)

	print("*************$$$$$$$$$$$$$$$$$$ ")
	print("Vij Matrix Coefficients Results:")
	print("Num#: ", len(coef_vals))
	for mvals in coef_vals:
		if any(x for x in mvals if x[1].startswith('alpha')) and any(y for y in mvals if y[1].startswith('beta')):
			print("MIXED ALPHA_BETA ERROR")
			print(mvals)
		else:
			 continue
			# print(mvals)
	return coef_vals
