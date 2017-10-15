# ******************************************************************************
# Name:    Functions for computing Gadget values
# Author:  Vadim Korotkikh
# Date:    October 2017
# Description:  Contains functions for calculating different Gadget types
#
# ******************************************************************************

# Library Imports
import itertools
from itertools import permutations
import numpy as np

#>******************************************************************************
# def main():
# 	nx = 4
# 	perms_list = list(itertools.permutations([0,1,2,3],nx))
# 	# parity_list = perm_parity(perms_list)
# 	for p in perms_list:
# 		l = list(p)
# 		# print("%2i %r" % (perm_parity(l), p))
# 		print("%r %2i" % (p, perm_parity(l)))

#>******************************************************************************
def perm_parity(lst):
	'''
	*** Fails if list doesn't start w/ 0 ***
	Given a permutation of the digits 0..N in order as a list,
	returns its parity (or sign): +1 for even parity; -1 for odd.
	'''
	parity = 1
	for i in range(0,len(lst)-1):
		if lst[i] != i:
			parity *= -1
			mn = min(range(i,len(lst)), key=lst.__getitem__)
			lst[i],lst[mn] = lst[mn],lst[i]
	return parity

#>******************************************************************************
def org_gadgetcalc(vij_holomat_list, abcoefs)

	lenvijlt = len(vij_holomat_list)
	print("Len vijlist: ", lenvijlt)
	vijlist = vij_holomat_list
	exelen = int(lenvijlt/32)
	totcount = 0
	probcount = 0
	if lenvijlt > 1 and isinstance(vij_holomat_list, list):
		for xind in range(0, exelen):
			# for zind in range(xind+1, exelen+1):
			for zind in range(xind+1, lenvijlt):
				# print("A:",xind,zind)
				# g1val = gadget_one_trace(vijlist[xind], vijlist[zind])
				# gorig = original_coef_gadget(abcoefs[xind], abcoefs[zind])
				gtrace = newgadget_trace(vijlist[xind], vijlist[zind])
				gnewab = newgadget_abcoef(abcoefs[xind], abcoefs[zind])
				# print("Adinkras:",xind,zind,"-> Gnew", gtrace, "GnewAB", gnewab)
				if gtrace == gnewab:
					totcount+=1
					# print("Adinkras:",xind,zind,"-> Gnew", gtrace,"GnewAB", gnewab)
					print("Adinkras",xind,zind,"-> GnewTrace", gtrace, "GnewAB", gnewab)
				elif gtrace != gnewab:
					probcount+=1
					print("ERROR: ",xind,zind,"-> GnewTrace", gtrace, "GnewAB", gnewab)

	print("################################################")
	print("Gadget #:", totcount)
	print("Error  #:", probcount)


#>******************************************************************************
def newgadget_trace(vij_holomats1, vij_holomats2):
	""" Calculate the new Gadget equation
		G(a,a') = (1/N)eps^(IJKL)*Tr[V(a)IJ, V(a')KL]
		vij_holomats1 = List of 6 vij matrices, each 4x4 np.array
		vij_holomats2 = List of 6 vij matrices, each 4x4 np.array

		vij6 = [v12, v13, v14, v23, v24, v34]
		calc = Tr[ v12v34 - v13v24 + v14v23 + v23v14 - v24v13 + v34v12 ]
	"""
	ev	= [1, -1, 1, 1, -1, 1]
	# ev = [-1, 1, -1, -1, 1, -1]
	vh1 = vij_holomats1
	vh2 = vij_holomats2
	somev = [np.multiply(np.dot(vh1[i],vh2[-i-1]),ev[i]) for i in range(0,len(ev))]
	trvals = [np.trace(vx) for vx in somev]
	gadgetval = int(sum(trvals)/(8))
	# print("NewGadg trace", trvals, " Trace Sum", gadgetval)
	# matc = np.dot(vh1[0],vh2[5]) - np.dot(vh1[1],vh2[4]) + \
	# + np.dot(vh1[2],vh2[3]) + np.dot(vh1[3],vh2[2]) + \
	# - np.dot(vh1[4], vh2[1]) + np.dot(vh1[5],vh2[0])
	# rawtrace = np.trace(matc)
	# print("Trace value:", rawtrace)
	return gadgetval

#>******************************************************************************
def gadget_deux_calc(holomat_list):
	pass

#>******************************************************************************
def gadget_one_trace(vij_holomats1, vij_holomats2):
	""" G(a,a') = (-1/N)Tr[V(a)IJ, V(a')IJ]

		vij_holomats1 = List of 6 vij matrices, each 4x4 np.array
		vij_holomats2 = List of 6 vij matrices, each 4x4 np.array

		vij6 = [v12, v13, v14, v23, v24, v34]
	"""
	ev	= [1, 1, 1, 1, 1, 1]
	vh1 = vij_holomats1
	vh2 = vij_holomats2
	''' For some reason Trace[classic gadget] must be multiplied by -1 '''
	somev = [np.multiply(np.matmul(vh1[i],vh2[i]),1) for i in range(0,6)]
	mkadd1 = [np.add(somev[ix*2],somev[(ix*2)+1]) for ix in range(0,3)]
	finsum = np.add(np.add(mkadd1[0], mkadd1[1]),mkadd1[2])
	trvals = [np.trace(vx) for vx in somev]

	gadgetval = int(np.trace(finsum)/(-8))
	# print("Gadget1 trace", trvals, " Trace Sum", gadgetval)
	return gadgetval

#>******************************************************************************
def original_coef_gadget(coef_l1, coef_l2):
	gadget_vals		= 0
	div_factor		= 2
	# one_count, ptre_count, ntre_count, zero_count = 0, 0, 0, 0

	ijf = coef_l1
	ijx = coef_l2

	if ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2] and ijf[2][0:2] == ijx[2][0:2]:
		gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in range(0, len(ijf))])
		# if gadget_sum == 2:
		# 	ptre_count += 1
		# elif gadget_sum == -2:
		# 	ntre_count += 1
		# elif gadget_sum == 6:
		# 	one_count += 1
		# elif gadget_sum == 0:
		# 	zero_count += 1
		# else:
		# 	print("Gadget ERROR 1:",gadget_sum)
		gadget_vals = gadget_sum / div_factor

	elif ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2]:
		gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [0, 5]])
		gadget_vals = gadget_sum / div_factor

	elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2]:
		gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [1, 4]])
		gadget_vals = gadget_sum / div_factor

	elif ijf[0][0:2] != ijx[0][0:2] and ijf[2][0:2] == ijx[2][0:2]:
		gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in [2, 3]])
		gadget_vals = gadget_sum / div_factor

	elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2] and ijf[2][0:2] != ijx[2][0:2]:
		gadget_sum = 0
		gadget_vals = gadget_sum / div_factor

	else:
		print("ERROR**********")
		print(ijf)
		print(ijx)
	return gadget_vals

#>******************************************************************************
def newgadget_abcoef(coef_l1, coef_l2):

	gadgetval	= 0
	# div_factor	= 2

	ijf = coef_l1.copy()
	ijx = coef_l2.copy()
	ev	= [1, -1, 1, 1, -1, 1]
	rng6 = range(0,6)
	rng3 = range(0,3)
	if all(ijf[i][0] == ijx[-i-1][0] for i in rng6):
		# print("ijf", coef_l1)
		# print("ijx", coef_l2[::-1])
		coefv = sum([(ijf[ix][2] * ijx[-ix-1][2])*(ev[ix]) for ix in rng6])
		# gadgetval = coefv/div_factor
		# print("NGadget AB:", int(coefv/(-2)))
		return int(coefv/(-2))
	elif any(ijf[i][0] == ijx[-i-1][0] for i in rng3):
		for ix in rng3:
			if (ijf[ix][0] == ijx[-ix-1][0]):
				coefv = sum([ijf[ix][2] * ijx[5-ix][2] * ev[ix], ijf[5-ix][2] * ijx[ix][2] * ev[5-ix]])
				# print("NGadget AB:", int(coefv/(-2)))
				return int(coefv/(-2))
		# print("ijf", coef_l1[0:3])
		# print("ijx", coef_l2[:2:-1])
	else:
		return gadgetval
