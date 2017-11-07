# ******************************************************************************
# Name:    Compare Gadget Values
# Author:  Vadim K
# Date:    November 2017
# Description:  Parse text files and compare Gadget values
#
# ******************************************************************************

import os, sys, glob
import multiprocessing as mp

#>******************************************************************************
def main(gvalpath1="", gvalpath2=""):

	filecnt	= 36864
	numpaks	= 8
	paksize	= int(filecnt/numpaks)
	indpaks	= [paksize*n for n in range(0,numpaks)]

	gvalpath1 = "/Users/vkorotki/Movies/CodArk-Adinkra/VMatrixNewGadget_Parallel/values"
	gvalpath2 = "/Users/vkorotki/Movies/CodArk-Adinkra/SUSY-Adinkra-NewGadget"

	gval_dirstr	= [("GadgetVal" + str(ix) + "of8") for ix in range(1,9)]
	gvaldirpath2 = gvalpath2
	gvalallpaths2 = [(gvalpath2 + "/" + x) for x in gval_dirstr]

	for ig in (gvalpath1, gvalpath2):
		if check_path(ig):
			pass
		else:
			sys.exit("DIRECTORY DNE")
			# pass

	gval_list1 = []
	gval_list2 = []

	# gval_list1 = parse_gadgetvals(gvalpath1)
	# print("Size", sys.getsizeof(gval_list1)/1024)

	for ix, pak in enumerate(indpaks):
		if ix >= 6:
			continue
		islice	= pak + paksize

		temp_vkpath = gvalpath2 + "/" + gval_dirstr[ix]
		print(temp_vkpath)
		lk_gval	= parse_lgadgetvals(gvalpath1, pak, islice)
		vk_gval	= parse_lgadgetvals(temp_vkpath, pak, islice, "Adinkra_")
		gvals_compare(vk_gval, lk_gval, pak, islice)
		# print("Size", sys.getsizeof(lk_gval)/1024)



#>******************************************************************************
def parse_gadgetvals(gvalpath, minnum=0, maxnum=4608, fst = ""):

	temp_list	= []
	frange 	= list(range(minnum,maxnum))
	maxline	= 26
	act_filecnt	= 0

	for ix in frange:
		valfile	= gvalpath + "/" + fst + str(ix) + ".txt"
		gval_dict	= {}
		gval_dict[ix] = []
		if os.path.isfile(valfile):
			act_filecnt+=1
		else:
			continue
		with open(valfile, 'r') as rf:
			print("Parsing file:", valfile)
			for inum, line in enumerate(rf):
				if inum == maxline:
					break
				gadgetval = line.split(">")[1].strip()
				gadgetnum = line.split("-")[0].strip()
				gadgetnum = gadgetnum.strip("()").split(",")
				# anum1, anum2 = gadgetnum[0], gadgetnum[1]
				anum2	= gadgetnum[1]
				# gval_dict[ix].append((anum1, anum2, gadgetval))
				gval_dict[ix].append((int(anum2), gadgetval))
		temp_list.append(gval_dict)
	return temp_list


#>******************************************************************************
def parse_lgadgetvals(gvalpath, minnum=0, maxnum=4608, fst = ""):

	# maxnum	= 250
	temp_list	= []
	frange 	= list(range(minnum,maxnum))
	gfilend	= ".txt"
	maxline	= 26
	act_filecnt	= 0
	
	for ix in frange:
		valfile	= gvalpath + "/" + fst + str(ix) + ".txt"
		gval_dict	= {}
		gval_dict[ix] = []
		if os.path.isfile(valfile):
			act_filecnt+=1
		else:
			continue
		with open(valfile, 'r') as rf:
			for inum, line in enumerate(rf):
				if inum == maxline:
					break
				gadgetval = line.split(">")[1].strip()
				gadgetnum = line.split("-")[0].strip()
				gadgetnum = gadgetnum.strip("()").split(",")
				anum1, anum2 = gadgetnum[0], gadgetnum[1]
				# gval_dict[ix].append((anum1, anum2, gadgetval))
				gval_dict[ix].append((int(anum2), gadgetval))
		# print(gval_dict)
		temp_list.append(gval_dict)
	return temp_list

#>******************************************************************************
def gvals_compare(val_list1, val_list2, minr, maxr):

	print("Types", type(val_list1), type(val_list2))
	print("Lengths", len(val_list1), len(val_list2))
	val_len = 0
	if len(val_list1) == len(val_list2):
		val_len = len(val_list1)

	dif_toggle	= 0
	dif_list	= []
	for ix in range(0, val_len):
		val_dict1 = val_list1[ix]
		val_dict2 = val_list2[ix]
		keyx = minr + ix
		# print("dict1", type(val_dict1[ix]))
		# print("dict2", val_dict2[ix])
		val_vals1	= set(val_dict1[keyx])
		val_vals2	= set(val_dict2[keyx])
		isect12	= val_vals1 ^ val_vals2
		# print("Number-diff: ",  len(isect12))
		if len(isect12) > 0:
			dif_toggle = 1
			print("Value difference")
			print("Adinkra %s", keyx)
			dif_list.append(keyx)
	if dif_toggle:
		for adinkra in dif_list:
			print(adinkra)


#>******************************************************************************
def check_path(direcpath):

	if os.path.isdir(direcpath):
		return direcpath
	else:
		return 0
		# sys.exit("DIRECTORY DNE")

#>******************************************************************************
if __name__ == "__main__":

	try:
		main(sys.argv[1], sys.argv[2])
	except IndexError:
		main()
