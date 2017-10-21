# ******************************************************************************
# Name:    MP some Calculations
# Author:  Vadim Korotkikh
# Date:    October 2017
# Description:  BLANK
#
# ******************************************************************************

# Library Imports
import random, string
import os, time, itertools
import numpy as np
# Specific imports -< FOR DEV >-
import logging
import multiprocessing as mp

import adinkra_nxn_constructor

#>******************************************************************************
def info(title):
	print(title)
	# print('module name:', __name__)
	print('parent process:', os.getppid())
	print('process id:', os.getpid())
	print("")

#>******************************************************************************
# def org_gadgetcalc(vij_holomat_list, abcoefs):
def mp_gadgetcalc_abonly(abcoef_list, numpaks=8):
	""" Do the new Gadget calc only for Alpha-Beta coefficients	"""
	print("Executing ", __name__)
	start_time = time.time()

	lenablist = len(abcoef_list)
	ablist = abcoef_list[:]

	numpaks	= 8	# Number of calculation packs
	nx		= int(lenablist/64)
	pool 	= mp.Pool(processes=64)
	paksize = int(lenablist/numpaks)	# Splitting 36k A into n sets of len paksize
	indpaks	= [paksize*n for n in range(0,numpaks)]

	for ix, pak in enumerate(indpaks):
		islice = pak + paksize
		abcoefpak 	= ablist[pak:islice]
		print("")
		# print("Length ab ", len(abcoefpak))
		print("Pack/Slice -", pak, ":", islice)
		for ind in range(0,paksize):
			adjadinknum = ind + pak
			print("Adinkra:", adjadinknum)
			icof	= abcoefpak[ind]
			complt	= []
			cpacked = [ablist[i:i+nx] for i in range(ind + pak, lenablist, nx)]
			cpklen 	= len(cpacked)
			print(cpklen)
			abcalc = [pool.apply_async(newgadget_abcoefs, args=(icof, cpacked[xi],xi)) for xi in range(0,cpklen)]
			for px in range(0,cpklen):
				abcofpak = abcalc[px].get()
				for gtval in abcofpak:
					gval = "Gadget val: " + str(gtval)
					complt.append(gval)

			# acalc = "GadgetVal/Adinkra_" + str(ind+islice) + "_GnewVal.txt"
			acalc = "GadgetVal/Adinkra_" + str(adjadinknum) + ".txt"
			with open(acalc, 'w') as wfile:
				wfile.write("Adinkra: %s \n" % ind)
				for cval in complt:
					wfile.write("%s \n" % cval)

	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))

#>******************************************************************************
# def org_gadgetcalc(vij_holomat_list, abcoefs):
def mp_gadgetcalc_abtr(vij_holomat_list, abcoefs):
	print("Executing ", __name__)
	start_time = time.time()

	lenvijlt = len(vij_holomat_list)
	vijlist = vij_holomat_list

	numpaks	= 8	# Number of calculation packs
	nx		= int(lenvijlt/64)
	pool 	= mp.Pool(processes=64)
	paksize = int(lenvijlt/numpaks)	# Splitting 36k A into n sets of len paksize
	indpaks	= [paksize*n for n in range(0,numpaks)]

	for ix, pak in enumerate(indpaks):

		islice = pak + paksize
		tracespak 	= vijlist[pak:islice]
		abcoefpak 	= abcoefs[pak:islice]
		# print("Length ab ", len(abcoefpak))
		# print("Pack/Slice -", pak, ":", islice)
		for ind in range(0,paksize):
			adjadinknum = ind + pak
			print("Adinkra:", adjadinknum)

			imat	= tracespak[ind]
			icof	= abcoefpak[ind]
			complt	= []
			npacked = [vijlist[i:i+nx] for i in range(islice+ind, lenvijlt, nx)]
			cpacked = [abcoefs[i:i+nx] for i in range(islice+ind, lenvijlt, nx)]
			npklen 	= len(npacked)
			mtracs = [pool.apply_async(newgadget_mtraces, args=(imat, npacked[xi], xi)) for xi in range(0,npklen)]
			abcalc = [pool.apply_async(newgadget_abcoefs, args=(icof, cpacked[xi],xi)) for xi in range(0,npklen)]
			pr_sw  = 0
			for px in range(0,npklen):
				mtracpak = mtracs[px].get()
				abcofpak = abcalc[px].get()
				# pr_sw = (mtracpak == abcofpak)
				# cstr = "Gadget val # " + str(px) + " Result: " + str(mtracpak == abcofpak)
				for gtval in abcofpak:
					gval = "Gadget val: " + str(gtval)
					complt.append(gval)

			# for px in chkval:
			# 	tl = px.get()
			# 	complt.append(tl)
			# acalc = "GadgetVal/Adinkra_" + str(ind+islice) + "_GnewVal.txt"
			acalc = "GadgetVal/Adinkra_" + str(adjadinknum) + ".txt"
			with open(acalc, 'w') as wfile:
				wfile.write("Adinkra: %s \n" % ind)
				for cval in complt:
					wfile.write("%s \n" % cval)

	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))


#>******************************************************************************
def newgadget_mtraces(vij_holomat, vij_holomats, adind):
	""" Calculate the new Gadget equation
		G(a,a') = (1/N)eps^(IJKL)*Tr[V(a)IJ, V(a')KL]
		vij_holomat = One of the vij matrices, each 4x4 np.array
		vij_holomats2 = List of 6 vij matrices, each 4x4 np.array

		vij6 = [v12, v13, v14, v23, v24, v34]
		calc = Tr[ v12v34 - v13v24 + v14v23 + v23v14 - v24v13 + v34v12 ]
	"""
	gadgetvals = []
	ev	= [1, -1, 1, 1, -1, 1]
	vh = vij_holomat
	vhs = vij_holomats
	for iv in vij_holomats:
		somev = [np.multiply(np.dot(vh[i],iv[-i-1]),ev[i]) for i in range(0,len(ev))]
		trvals = [np.trace(vx) for vx in somev]
		gadgetval = int(sum(trvals)/(8))
		gadgetvals.append(gadgetval)
	return gadgetvals

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
def newgadget_abcoefs(coef_l1, abcoefs, xind):

	gadgetval	= 0
	gadgetvals 	= []
	ijf = coef_l1.copy()
	abind	= xind*len(abcoefs)
	# div_factor	= 2

	ev	= [1, -1, 1, 1, -1, 1]
	rng6 = range(0,6)
	rng3 = range(0,3)
	for ind, abcoef in enumerate(abcoefs):
		ijx = abcoef
		if all(ijf[i][0] == ijx[-i-1][0] for i in rng6):
			# print("ijf", coef_l1)
			# print("ijx", coef_l2[::-1])
			coefv = sum([(ijf[ix][2] * ijx[-ix-1][2])*(ev[ix]) for ix in rng6])
			# gadgetval = coefv/div_factor
			# print("NGadget AB:", int(coefv/(-2)))
			gadgetvals.append((int(coefv/(-2)), abind))
		elif any(ijf[i][0] == ijx[-i-1][0] for i in rng3):
			for ix in rng3:
				if (ijf[ix][0] == ijx[-ix-1][0]):
					coefv = sum([ijf[ix][2] * ijx[5-ix][2] * ev[ix], ijf[5-ix][2] * ijx[ix][2] * ev[5-ix]])
					# print("NGadget AB:", int(coefv/(-2)))
					gadgetvals.append((int(coefv/(-2)), abind))
			# print("ijf", coef_l1[0:3])
			# print("ijx", coef_l2[:2:-1])
		else:
			gadgetvals.append((0, abind))
	return gadgetvals

#>******************************************************************************
def mpp_org_gadgetcalc(vij_holomat_list, abcoefs):
	print("Executing ", __name__)
	start_time = time.time()
	lenvijlt = len(vij_holomat_list)
	vijlist = vij_holomat_list
	exelen = int(lenvijlt/32)

	pack	= 32
	nx		= int(lenvijlt/16)
	pool = mp.Pool(processes=8)
	fonesix = vijlist[0:pack]
	# for ind, xad in enumerate(fonesix):
	rescalcs = [pool.apply_async(newgadget_mtraces, args=(fonesix[xs], vijlist[xs:], xs)) for xs in range(0,pack)]
	for p in rescalcs:
		allgs = p.get()

	# for ind, xad in enumerate(fonesix):
	# 	print("Adinkra:", ind)
	# 	npacked = [vijlist[i:i+nx] for i in range(ind, lenvijlt, nx)]
	# 	# acalc 	= [pool.apply(newgadget_mtraces, args=(xad, npacked[xi], xi)) for xi in range(0,16)]
	# 	acalc = [pool.apply_async(newgadget_mtraces, args=(xad, npacked[xi], xi)) for xi in range(0,len(npacked))]
	# 	for px in acalc:
	# 		somex = px.get()
			# print(len(somex))
		# acalc = "Adinkra_" + str(ind+1) + "_Gadgets.txt"
		# with open(acalc, 'w') as rf:


#>******************************************************************************
def rand_string(length, output):
	info('function rand_string')
	rand_str = ''.join(random.choice(string.ascii_lowercase + string.ascii_uppercase
									+ string.digits) for i in range(length))
	output.put(rand_str)

#>******************************************************************************
def test_mpp():
	random.seed(123)
	# output Queue
	output = mp.Queue()
	processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]

	# Run the processes, starting sequentially
	for pt in processes:
		pt.start()
	# Exit the completed processes
	for pt in processes:
		pt.join()

	resutstr = [output.get() for p in processes]
	print(resutstr)

if __name__ == "__main__":
	# test_mpp()
	start_time = time.time()
	test_mpp()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
