# ******************************************************************************
# Name:    Testing
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:	   February 2017
#
# Description: Algorithm for generating n-colour sized Adinkras of n-node
# matrix representations

import os, time, itertools
from cython.parallel import *

def main():
	nx = 4
	perms_list = list(itertools.permutations([0,1,2,3],nx))
	# parity_list = perm_parity(perms_list)
	for p in perms_list:
		l = list(p)
		# print("%2i %r" % (perm_parity(l), p))
		print("%r %2i" % (p, perm_parity(l)))
	print("")
	gadgetVV()

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

# def gadgetVV(allVtilde1, allVtilde2):
def gadgetVV():
	value = 0
	LC = [[[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  1],[ 0,  0, -1,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  0, -1],[ 0,  0,  0,  0],[ 0,  1,  0,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  1,  0],[ 0, -1,  0,  0],[ 0,  0,  0,  0]]],

	[[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0, -1],[ 0,  0,  1,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]],
	[[ 0,  0,  0,  1],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[-1,  0,  0,  0]],
	[[ 0,  0, -1,  0],[ 0,  0,  0,  0],[ 1,  0,  0,  0],[ 0,  0,  0,  0]]],

	[[[ 0,  0,  0,  0],[ 0,  0,  0,  1],[ 0,  0,  0,  0],[ 0, -1,  0,  0]],
	[[ 0,  0,  0, -1],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 1,  0,  0,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]],
	[[ 0,  1,  0,  0],[-1,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]]],

	[[[ 0,  0,  0,  0],[ 0,  0, -1,  0],[ 0,  1,  0,  0],[ 0,  0,  0,  0]],
	[[ 0,  0,  1,  0],[ 0,  0,  0,  0],[-1,  0,  0,  0],[ 0,  0,  0,  0]],
	[[ 0, -1,  0,  0],[ 1,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]],
	[[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0],[ 0,  0,  0,  0]]]]
	perms_list = list(itertools.permutations([0,1,2,3],4))
	ps = perms_list[:]

	print("ϵ^(ijkl) values")
	for ix in ps:
		LCtemp = LC[ix[0]][ix[1]][ix[2]][ix[3]]
		ixmod = (ix[0]+1, ix[1]+1, ix[2]+1, ix[3]+1)
		ixmodstr = "".join([str(x) for x in ixmod])
		eijkl = "ϵ^(" + ixmodstr + ") = "
		print(eijkl, LCtemp)
	print("")
	for i in prange(4):
		print("prange i:", i)
		for j in prange(4):
			for k in prange(4):
				for l in prange(4):
					LCtemp = LC[i-1][j-1][k-1][l-1]
					if LCtemp != 0:
						print("e ijkl:",i+1,j+1,k+1,l+1," - ", LCtemp)


if __name__ == "__main__":
	start_time = time.time()
	main()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
