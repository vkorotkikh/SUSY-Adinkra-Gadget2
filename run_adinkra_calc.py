# ******************************************************************************
# Name:    Caculate All Possible BC4-based Adinkras
# Author:  Vadim Korotkikh
# Date:    December 2016
# Description: The code calculates all all unique 36,864 ordered BC4-based
# adinkras with four colors, four open-nodes and four closed nodes and then
# calculates the Vij holoraumy matrices and the Gadget values
#
# ******************************************************************************

# Library Imports
import os, sys, time
import numpy as np

# Function Imports
import fx_gadgets
import vij_holoraumy_calc
import adinkra_nxn_constructor
import cls_adinkra_set
import fx_mpgadgets

#>******************************************************************************
def main():
	print("#>***********************************************************************")
	print("# Name:    Caculate all BC4-based Adinkras")
	print("# Author:  Vadim Korotkikh	")
	# print("# Email:   va.korotki@gmail.com")
	print("# Date:    December 2016		")
	print("# Version: 2.0 Update in Progress Sept 2017")
	print("#							")
	print("# Description: Calculates all unique 36,864 ordered BC4-based adinkras")
	print("# with four colors, four open-nodes and four closed nodes.             ")
	print("#	")
	print("#>***********************************************************************")
	print("		")
	adinkra_list	= []
	adinkra_list	= adinkra_nxn_constructor.create_adinkras(4,4)

	# check_or_makedir("GadgetVal", '')
	if len(adinkra_list) == 36864:
		NewAdink = cls_adinkra_set.AdinkraSet(4,4,adinkra_list)
		NewAdink.exe_fermiorder()
		print("Executing New Gadget calc")
		vijmat_list = NewAdink.ret_fermiholo()
		abcoef_list = NewAdink.ret_fermi_abcoef()
		# ultra_fermi = NewAdink.ret_ultrafermi()

		# fx_mpgadgets.mporg_gadgetcalc(vijmat_list, abcoef_list)
		fx_mpgadgets.mp_gadgetcalc_abonly(abcoef_list)
		# fx_gadgets.org_gadgetcalc(vijmat_list, abcoef_list,ultra_fermi)
		# NewAdink.exe_fermiorder()

#>******************************************************************************
def check_or_makedir(dirname, dirpath=''):

	runad_fullpath = os.path.realpath(__file__)
	runad_dirpath  = os.path.dirname(runad_fullpath)
	if dirpath != '' and dirpath == runad_dirpath:
		pass
	elif dirpath != '':
		runad_dirpath = dirpath
	elif dirpath == '':
		pass
	# print("dirname :", runad_dirpath)
	chkdirpath = runad_dirpath + "/" + dirname
	print("Does this directory exist: ", chkdirpath)
	if not os.path.isdir(chkdirpath):
		os.makedirs(chkdirpath)
	else:
		pass

#>******************************************************************************
def user_options():
	""" Terminal user options defined/coded here w/ help from inner functions
	"""

	# **************************************************************************
	def core_options():
		print("Choose one of the following calculation options:")
		print("")
		print(" < 1 >  -  Calculate all BC4 Adinkras")
		print(" < 2 >  -  Calculate BC4 Adinkras Gadget I values")
		print(" < 3 >  -  Calculate BC4 Adinkras Gadget II values")
		print(" < 4 >  -  Calculate BC4 Adinkras Gadget III values")
		# print(" < 5 >  -  Set output file string")
		print(" < 5 >  =  Nevermind. Get me outa here! Exit")
		print("")
		return input(": ")

	# Set loopcount = 0 of no arg is supplied for first time
	def option_one(loopcount=0):
		counter	= 0
		print("")
		print(" < 1 >  -  Display entire BC4 CG Small Library")
		print(" < 2 >  -  Display select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			# Lets make the P-Set options a function
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_one(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	def option_two(loopcount=0):
		print("")
		print(" < 1 >  -  Verify entire BC4 CG Library")
		print(" < 2 >  -  Verify select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			# for now
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			print(usr_pset)
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_two(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	# Set loopcount = 0 of no arg is supplied for first time
	def option_three(loopcount=0):
		print("")
		print(" < 1 >  -  Calculate all P-sets")
		print(" < 2 >  -  Calculate select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			# bc4_validation_organizer('PALL', 'Vmats', 'fermi')
			bc4_validation_organizer('PALL', 'mats', 'fermi')
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			print(usr_pset)
			bc4_validation_organizer(usr_pset, 'mats', 'fermi')
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_three(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')


	def option_five():
		print("NOT ACTIVATED (code not finished)")
		print("Going back to main menu")
		option_activator('core')

	def option_six():
		print("")
		print("Quiting script. Are you sure (yes/no)?")
		ninput = input(": ")
		if ninput.lower() == 'yes':
			sys.exit("EXITING BC4 CG Library Utility")
		elif ninput.lower() == 'no':
			print("Going back to main menu")
			option_activator('core')
		else:
			pass

	# **************************************************************************
	def option_activator(input_str, mcounter=0):
		""" option_activator - Inner function of user_options that depending on
			input_str either executes other inner functions or rexecutes itself
			with core_options() func. providing input string.
		"""
		if input_str.strip() == '1':
			option_one()
		elif input_str.strip() == '2':
			option_two()
		elif input_str.strip() == '3':
			option_three()
		# elif input_str.strip() == '4':
		# 	option_four()
		# elif input_str.strip() == '5':
			option_five()
		elif input_str.strip() == '4' or input_str.lower() == 'exit':
			option_six()
		elif input_str.strip() == 'core' or input_str.lower() == 'core':
			option_activator(core_options())
		else:
			print("UNRECOGNIZED SELECTION: ", input_str)
			print("Try again...")
			if mcounter <= 6:
				mcounter+=1
				option_activator(core_options(), mcounter)
			elif mcounter >= 7:
				print("Too many attempts. Try again later.")
				# print("EXITING BC4 CG Library Utility")

				sys.exit("EXITING BC4 CG Library Utility")

	uinput = core_options()
	option_activator(uinput)

#>******************************************************************************
if __name__ == "__main__":
	start_time = time.time()
	main()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
