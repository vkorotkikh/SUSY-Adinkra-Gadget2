# ******************************************************************************
# Name:    Adinkra sets class
# Author:  Vadim Korotkikh
# Date:    September 2017
# Tab delimited = 4 spaces.
# Description: Storage class for Adinkras and their various Vij matrix lists
#
# ******************************************************************************

# Library Imports
import itertools
import numpy as np

# Function Imports
import fx_vij_holoraumy

#>******************************************************************************
class AdinkraSet():

	@classmethod
	def aset_classmethod(cls):
		print("Class Name: %s" % cls.__name__)

	@staticmethod
	def aset_method(tbd_var):
		print(tbd_var, type(tbd_var))

		# pass
	#>**************************************************************************
	def __init__(self, dim, nodes, adinkra_list):

		self.dim 		= dim
		self.nodes		= nodes

		self.adinkra_list = adinkra_list
		self.fermi_mats		= []
		self.ultrafermi_mt	= []
		self.fermi_abcoefs 	= []
		# self.abcoef_list	= []

	#>**************************************************************************
	def exe_fermiorder(self):		# class instance method
		self.get_fermiholo()
		self.get_fermi_abcoef()
		# self.get_ultrafermi()

	#>**************************************************************************
	def get_fermiholo(self):
		templist = fx_vij_holoraumy.calc_fermi_vij(self.adinkra_list)
		self.fermi_mats = templist

	#>**************************************************************************
	def get_fermi_abcoef(self):
		templist = fx_vij_holoraumy.calc_vij_alphabeta(self.adinkra_list)
		self.fermi_abcoefs = templist

	#>**************************************************************************
	def ret_fermiholo(self):
		return self.fermi_mats

	#>**************************************************************************
	def ret_ultrafermi(self):
		return self.ultrafermi_mt

	#>**************************************************************************
	def ret_fermi_abcoef(self):
		return self.fermi_abcoefs
