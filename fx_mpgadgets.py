# ******************************************************************************
# Name:    MP some Calculations
# Author:  Vadim Korotkikh
# Date:    October 2017
# Description:  BLANK
#
# ******************************************************************************

# Library Imports
import os, time, itertools, shutil
import numpy as np
# Specific imports -< FOR DEV >-
import logging
import multiprocessing as mp

logger = logging.getLogger(__name__)

#>******************************************************************************
def info(xtitle):
	print('parent process:', os.getppid(), 'title', xtitle)
	print('process id:', os.getpid())
	print("")

#>******************************************************************************
def mp_gcalcab_valcount(abcoeflist, mpcount=64, numpaks=8):
	return mp_gadgetcalc_abonly(abcoeflist, mpcount, numpaks)

#>******************************************************************************
# def org_gadgetcalc(vij_holomat_list, abcoefs):
def mp_gadgetcalc_abonly(abcoef_list, mpcount=64, numpaks=8):
	""" This function executs the New Gadget calculation aka Gadget mk II using
	passed in Alpha-Beta elle and ~elle coefficients. Uses the multiprocessing
	module to speed up the Gadget calculation.
	abcoef_list - list of lists, w/ each list containg 6 int coefficients
	mpcount - number of threads for multiprocessing.Pool
	numpaks - default = 8. num of calculation sets or packs. The 36864 are
	divided into this amount of equally sized sets, if possible.
	"""
	# print("Executing ", __name__)
	logger.info("Running %s " % (__name__))
	start_time = time.time()
	lenablist = len(abcoef_list)
	ablist = abcoef_list[:]

	numpaks	= 8	# Number of calculation packs
	nx		= int(lenablist/mpcount)
	pool 	= mp.Pool(processes=mpcount)
	paksize = int(lenablist/numpaks)	# Splitting 36k A into n sets of len paksize
	indpaks	= [paksize*n for n in range(0,numpaks)]
	pgvals	= {}

	''' Splitting Gadget calculation into 8 equally sized sets of length 4608.
	This is done to minimize the number of Adinkra files per results directory
	'''
	for ix, pak in enumerate(indpaks): # for ix 0-7, pak 0-32256, steps 4608
		islice = pak + paksize
		abcoefpak 	= ablist[pak:islice]
		start_pakt = time.time()
		print("Pack/Slice: ", pak, ":", islice)
		for ind in range(0,paksize): # paksize = 4608
			adjnum = ind + pak
			adjanumstr	= str(adjnum)
			cprintstr1	= "(" + adjanumstr + ", "
			icof	= abcoefpak[ind] # i adinkra a/b coef integers
			complt	= []
			rngval 	= list(range(adjnum, lenablist, nx)) # lenablist = 36k, nx = 576
			cpacked = [ablist[i:i+nx] for i in rngval if (i+nx) < lenablist]
			cpacked.append(ablist[rngval[-1]:])
			cpklen 	= len(cpacked)
			abcalc 	= [pool.apply_async(ng_abc, args=(icof, cpacked[xi],xi, adjnum)) for xi in range(0,cpklen)]
			for px in range(0,cpklen):
				abcofpak = abcalc[px].get()
				# abcofpak = abcalc[px][0]
				for gtval in abcofpak:
					if gtval[0] not in pgvals:
						pgvals[gtval[0]] = 1
					elif gtval[0] in pgvals:
						pgvals[gtval[0]]+=1

			# 		indstr = cprintstr1 + str(gtval[1]) + ") -> "
			# 		gval = indstr + gtval[0]
			# 		complt.append(gval)
			''' Gadget vals written to text files '''
			# acalc = "GadgetVal/Adinkra_" + adjanumstr + ".txt"
			# with open(acalc, 'w') as wfile:
			# 	for cval in complt:
			# 		wfile.write("%s \n" % cval)
		calcrestxt 	= "GadgetVal/gvalstats" + str(ix+1) + "of8.txt"
		pakctime	= time.time() - start_pakt
		adinpermin 	= int(paksize/(pakctime/60))
		with open(calcrestxt, 'w') as wres:
			wres.write(("Adinkra Slice: " + str(pak) + " : " + str(islice) + "\n"))
			wres.write("\n")
			wres.write(" Gadget - Counts \n")
			for key, gval in pgvals.items():
				gvalmul = int(gval*2)
				wres.write(" %s	= %s \n" % (str(key), str(gvalmul)))
			wres.write("---- Execution time ----\n")
			wres.write(("-- " +str(pakctime) + " seconds --\n"))
			wres.write(("-- " + str(adinpermin) + " Adinkras / minute --\n"))
		# ixdir_name	= "GadgetVal" + str(ix+1) + "of8"
		# resfiles_p	= os.path.join(mkpath, 'GadgetVal')
		# fpath_ixdir	= os.path.join(mkpath, ixdir_name)
		# if not os.path.isdir(fpath_ixdir):
		# 	shutil.move(resfiles_p, fpath_ixdir)
		print("New Gadget Values", pgvals)
	tval = time.time() - start_time
	print("New Gadget Values", pgvals)
	print("-- Execution time --")
	print("---- %s seconds ----" % tval)
	print("---- %s Adinkras / minute ----" % int(paksize/(tval/60)))

#>******************************************************************************
# def org_gadgetcalc(vij_holomat_list, abcoefs):
def mp_gadgetcalc_abtrace(vij_holomat_list, abcoefs, numpaks=8):
	""" Performs new Gadget calculation using both the Trace method and
		Alpha-Beta coefficient method
	"""
	print("Executing ", __name__)
	start_time = time.time()

	lenvijlt = len(vij_holomat_list)
	vijlist = vij_holomat_list[:]

	numpaks	= 8	# Number of calculation packs
	divcntr = 64
	nx		= int(lenvijlt/divcntr)
	pool 	= mp.Pool(processes=divcntr)
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
def check_or_makedir(dirname, dirpath=''):
	print("Dirpath: ", dirpath)
	runad_fullpath = os.path.realpath(__file__)
	runad_dirpath  = os.path.dirname(runad_fullpath)
	if dirpath != '' and dirpath == runad_dirpath:
		pass
	elif dirpath != '':
		runad_dirpath = dirpath
	# print("dirname :", runad_dirpath)
	chkdirpath = runad_dirpath + "/" + dirname
	if not os.path.isdir(chkdirpath):
		os.makedirs(chkdirpath)
		return runad_dirpath
	else:
		return runad_dirpath

#>******************************************************************************
def newgadget_mtraces(vij_holomat, vij_holomats_list, adind):
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
	# vhs = vij_holomats_list
	for iv in vij_holomats_list:
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
def ng_abc(coef, locabc, xi, sti):
	return newgadget_abcoefs(coef, locabc, xi, sti)

#>******************************************************************************
def newgadget_abcoefs(coef_l1, locabcoefs, xind,  stind):
	""" Performs the Gadget Mk II aka the New Gadget calculation by using the
	Alpha-Beta coefficients to calculate the final Gadget value for each Adinkra
	"""

	gadgetvals 	= []
	ijf = coef_l1.copy()
	# startind	= xind*len(locabcoefs) + stind
	startind 	= xind * 576 + stind
	div_factor	= -6

	ev	= [1, -1, 1, 1, -1, 1]
	rng6 = range(0,6)
	rng3 = range(0,3)
	for xab, abcoef in enumerate(locabcoefs):
		reali = startind + xab
		ijx = abcoef
		gadgetval	= "0"
		if all(ijf[i][0] == ijx[-i-1][0] for i in rng6):
			# print("ijf", coef_l1)
			# print("ijx", coef_l2[::-1])
			coefv = sum([(ijf[ix][2] * ijx[-ix-1][2])*(ev[ix]) for ix in rng6])/div_factor
			if coefv == (-1/3):
				gadgetval = "-1/3"
			elif coefv == (1/3):
				gadgetval = "1/3"
			elif coefv == 1:
				gadgetval = "1"
			elif coefv == -1:
				gadgetval = "-1"
			# gadgetval = coefv/div_factor
			# print("NGadget AB:", int(coefv/(-2)))
			gadgetvals.append((gadgetval, reali))
		elif any(ijf[i][0] == ijx[-i-1][0] for i in rng3):
			for ix in rng3:
				if (ijf[ix][0] == ijx[-ix-1][0]):
					coefv = sum([ijf[ix][2] * ijx[5-ix][2] * ev[ix], ijf[5-ix][2] * ijx[ix][2] * ev[5-ix]])/div_factor
					# print("NGadget AB:", int(coefv/(-2)))
					if coefv == (-1/3):
						gadgetval = "-1/3"
					elif coefv == (1/3):
						gadgetval = "1/3"
					elif coefv == 1:
						gadgetval = "1"
					elif coefv == -1:
						gadgetval = "-1"
					gadgetvals.append((gadgetval, reali))
			# print("ijf", coef_l1[0:3])
			# print("ijx", coef_l2[:2:-1])
		else:
			gadgetvals.append((gadgetval, reali))
	# return gadgetvals, (startind, len(locabcoefs)) # for index value testing
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
