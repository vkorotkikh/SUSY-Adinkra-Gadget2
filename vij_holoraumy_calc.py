# ******************************************************************************
# Name:    Calculate Vij matrices and Adinkra Gadget values
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    November 2016
# Version: 1.3
#
# Description: Scripts for calculating Vij matrices for each one of 36864
# unique Adinkra tetrads and scripts for calculating the Gadget values from the
# Vij matrices
#
# ******************************************************************************


# ******************************************************************************
# Begin Imports
import numpy as np
import numpy.matlib
import itertools
# from numpy import array
# from numpy.linalg import inv

# import matrix_outerprod_calc
def illuminator_of_elfes():

	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are results of outerproducts inbetween
		Pauli matrices and Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	alpha1i = np.array([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.array([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.array([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.array([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.array([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.array([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])

	return [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]

# ******************************************************************************
# Do the final Vij calculation
def calc_vij_matrices(main_tetrad_list):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """

	vij_possibilities = []
	''' vij_possibilities should be a function input argument instead of it
		being hardcoded in inside the function definition/scope	'''
	vij_possibilities = illuminator_of_elfes()

	print("							")
	print("Calculating Vij matrices")
	print("							")
	vij_alphas 		= []
	vij_betas  		= []

	vij_mats		= []

	vij_matrices	= []

	anomaly_switch  = 0
	debug			= 0

	for ti, teti in enumerate(main_tetrad_list):
		if debug:
			print("# ********************************")
			print("								     ")
			print("Tetrad i: ", ti)
			# print(teti[0][1][0,:], teti[1][1][0,:], teti[2][1][0,:], teti[3][1][0,:])
			# print(teti[0][1][1,:], teti[1][1][1,:], teti[2][1][1,:], teti[3][1][1,:])
			# print(teti[0][1][2,:], teti[1][1][2,:], teti[2][1][2,:], teti[3][1][2,:])
			# print(teti[0][1][3,:], teti[1][1][3,:], teti[2][1][3,:], teti[3][1][3,:])
			print("								     ")

		vij_tempset = []

		""" Store 6 Vij matrices in temp_vijmat"""
		temp_vijmat		= []

		""" This section does a double loop over the same tetrad to calculate
		the set of 6 Vij matrices for the tetrad.
		So for each matrix in the tetrad its checked against all the possible others,
	 	bypassing the duplicate calculations
		"""
		for i, li in enumerate(teti):
			# print(li[1])
			bigli = li[1]
			tr_bigli = np.transpose(bigli)
			temp_combos = []

			for j, lj in enumerate(teti):
				biglj = lj[1]
				ij_temp = [i, j]
				ij_temp.sort()
				ir = i + 1
				jr = j + 1
				ijstr = str(ir) + str(jr)
				if ij_temp not in temp_combos and i != j:
					# print("Vij matrix i-j vals:", ij_temp)
					# print("Vij matrix i-j vals:", ijstr)
					temp_combos.append(ij_temp)
					tr_biglj = np.transpose(biglj)
					# temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
					""" Vij eq from 1601.00 (3.2) """
					# temp_mat = np.matmul(tr_biglj, bigli) - np.matmul(tr_bigli, biglj)
					temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
					tf_bool = 0
					# Compare against the 6 possible matrix solutions
					for xi, ijx in enumerate(vij_possibilities):
						ijx_neg = np.multiply(ijx, -1)
						# print(xi)
						if np.array_equal(temp_mat, ijx):
							tf_bool = 1
							temp_vijmat.append(temp_mat)
							if debug:
								print("*************$$$$$$$$$$$$$$$$$$ ")
								print("l-solution found:")
								print(ijx)
							tmint = np.int(1)
							if xi < 3:
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								vij_tempset.append([tmp_str, ijstr, tmint])
							elif xi >= 3:
								tmp_str = "beta" + str((xi - 2))
								vij_tempset.append([tmp_str, ijstr, tmint])
						elif np.array_equal(temp_mat, ijx_neg):
							tf_bool = 1
							temp_vijmat.append(temp_mat)
							if debug:
								print("*************$$$$$$$$$$$$$$$$$$ ")
								print("l-solution found:")
								print(ijx_neg)
							# xint = (xi + 1) * ( -1)
							tmint = np.int(-1)
							if xi < 3:
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								vij_tempset.append([tmp_str, ijstr, tmint])
							elif xi >= 3:
								tmp_str = "beta" + str((xi - 2))
								vij_tempset.append([tmp_str, ijstr, tmint])
						else:
							if i != j and tf_bool == 0 and xi >= 5:
								if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
									print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
									print("Anomaly found:",i,j)
									print(temp_mat)
									anomaly_switch = 1
					tf_bool = 0

		# vij_matrices.append(temp_vijmat)
		vij_mats.append(vij_tempset)

	print("*************$$$$$$$$$$$$$$$$$$ ")
	print("Vij Matrix Coefficients Results:")
	print("")
	for mvals in vij_mats:
		if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
			print("MIXED ALPHA_BETA ERROR")
			print(mvals)
		else:
			continue
			# print(mvals)

	return vij_mats


#>**************************************************************************
def gadget_two_calc(vijmat_list):

	gadget_vals		= []
	one_count 		= 0
	ptre_count		= 0
	ntre_count		= 0
	zero_count		= 0
	if not anomaly_switch:
		for fi, ijf in enumerate(vijmat_list):

			for xj, ijx in enumerate(vijmat_list):
				# ind_temp = [fi, xj]
				# ind_temp.sort()
				# x = [val]
				if ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2] and ijf[2][0:2] == ijx[2][0:2]:
					# als = ijf[0][3] * ijx[0][3]
					gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in range(0, len(ijf))])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print(ijf)
						print(ijx)
						print("Gadget ERROR 1:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("****** Gadget calculation ******")
					# print("Calc #:", calc_count)
					# print(div_const)
					# print("G values:", gadget_vals)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2]:
					gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [0, 5]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 2:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2]:
					# print(ijf, ijx)
					gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [1, 4]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 3:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[2][0:2] == ijx[2][0:2]:
					gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in [2, 3]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 4:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2] and ijf[2][0:2] != ijx[2][0:2]:
					gadget_sum = 0
					zero_count += 1

					div_const = gadget_sum / 6
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				else:
					print("ERROR**********")
					print(ijf)
					print(ijx)
			print("zero count", zero_count)
			print(" 1/3 count", ptre_count)
			print("-1/3 count", ntre_count)
			print("  1  count", one_count)
			print(gadget_vals)
	else:
		pass

	print("################################################")
	print(" Printing final Gadget values and counts        ")
	print("							")
	print("zero count", zero_count)
	print(" 1/3 count", ptre_count)
	print("-1/3 count", ntre_count)
	print("  1  count", one_count)
	print(gadget_vals)
