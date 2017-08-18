#!/usr/bin/python

#Equations for Pharm regression
#Author: Evan Su
#Version: 1.0
#Date: 6/19/2017 
from lmfit  import Parameters #optimize parameters
import math #math functions


##class implementation#############################################
#class, NVDataSet#######################################
#purpose: container for all variables
#note: This class is specifically created for lmfit libraries.
#contents:
##variables: 
'''
kAtt
d
L
kAd
kT
ATP
ton
kDet
r
Vlittle
'''
##functions:
'''
___init___(self)
__str__(self)
recalc(self)
set_*(self,value) 
		#value: desired value to set
'''
class BioEquClass:

	def __init__(self):


		self.param = Parameters()
		self.param.add('N', value = 1, 
					min = 0.0, vary = False)

		self.param.add('kAtt', value = 1, 
					min = 0.0, vary = False)

		self.param.add('d', value = 1, 
					min = 0.0, vary = False)

		self.param.add('L', value = 1, 
					min = 0.0, vary = False)

		self.param.add('kAd', value = 1, 
					min = 0.0, vary = False)

		self.param.add('kT', value = 1, 
					min = 0.0, vary = False)

		self.param.add('ATP', value = 1, 
					min = 0.0, vary = False)
		
		self.ton = 1
		self.kDet = 1
		self.r = 1
		self.Vlittle = 1
		self.VResult = 1

	def __str__(self):
		string = ''

		string = string +"N (s^-1): %.3f" % (self.param['N'].value)
		if(self.param['N'].vary) and not (self.param['N'].stderr is None):
			string = string +" +/- %.3f" % (self.param['N'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"kAtt (s^-1): %.3f" % (self.param['kAtt'].value)
		if(self.param['kAtt'].vary) and not (self.param['kAtt'].stderr is None):
			string = string +" +/- %.3f" % (self.param['kAtt'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nd (nm): %.3f" % (self.param['d'].value)
		if(self.param['d'].vary) and not (self.param['d'].stderr is None):
			string = string +" +/- %.3f" % (self.param['d'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nL(nm): %.3f" % (self.param['L'].value)
		if(self.param['L'].vary) and not (self.param['L'].stderr is None):
			string = string +" +/- %.3f" % (self.param['L'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nk-AD(s^-1): %.3f" % (self.param['kAd'].value)
		if(self.param['kAd'].vary) and not (self.param['kAd'].stderr is None):
			string = string +" +/- %.3f" % (self.param['kAd'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nkT (uM^-1 * s^-1): %.3f" % (self.param['kT'].value)
		if(self.param['kT'].vary) and not (self.param['kT'].stderr is None):
			string = string +" +/- %.3f" % (self.param['kT'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nATP (uM): %.3f" % (self.param['ATP'].value)
		if(self.param['ATP'].vary) and not (self.param['ATP'].stderr is None):
			string = string +" +/- %.3f" % (self.param['ATP'].stderr)
		else:
			string = string +" (fixed)"


		string = string +"\nton (s): %.3f" % (self.ton)
		string = string +"\nkDet (s^-1): %.3f" % (self.kDet)
		string = string +"\nr: %.3f" % (self.r)
		string = string +"\nVlittle (s^-1): %.3f" % (self.Vlittle)
		string = string +"\nVResult: %.3f" % self.VResult
		string = string +"\n"

		return string


	def recalc(self):
		#recalc dependend variables
		self.ton = calc_ton(self.param['kAd'], self.param['kT'], self.param['ATP'])
		self.kDet = calc_kDet(self.ton)
		self.r = calc_r(self.ton, self.param['kAtt'])
		self.Vlittle = calc_Vlittle(  self.param['kAtt'],self.param['kAd'],
											 	self.param['kT'], self.param['ATP'])
		Vprev = 0.0001
		self.VResult = calc_Vn(self.r, self.param['L'], self.ton, Vprev, 
								self.param['N'], self.Vlittle, self.param['d'])

	#set value functions in no particular order
	def set_N(self, value):
		self.param['N'].value = value
	def set_kAtt(self, value):
		self.param['kAtt'].value = value
	def set_d(self, value):
		self.param['d'].value = value
	def set_L(self, value):
		self.param['L'].value = value
	def set_kAd(self, value):
		self.param['kAd'].value = value
	def set_kT(self, value):
		self.param['kT'].value = value
	def set_ATP(self, value):
		self.param['ATP'].value = value


	def set_NVary(self, value):
		self.param['N'].vary = value
	def set_kAttVary(self, value):
		self.param['kAtt'].vary = value
	def set_dVary(self, value):
		self.param['d'].vary = value
	def set_LVary(self, value):
		self.param['L'].vary = value
	def set_kAdVary(self, value):
		self.param['kAd'].vary = value
	def set_kTVary(self, value):
		self.param['kT'].vary = value
	def set_ATPVary(self, value):
		self.param['ATP'].vary = value


##function implementation#############################################

def calc_Vn(r, L, ton, Vprev, N, Vlittle, d):
	'''
	#for debugging
	print r
	print L
	print ton
	print Vprev
	print N
	print Vlittle
	print d
	'''
	#calls the actual recursive function
	result = calc_Vn_recursive(r, L, ton, Vprev, N, Vlittle, d, 0)
	return result

def calc_Vn_recursive(r, L, ton, Vprev, N, Vlittle, d, iterations):

	p = calc_p(r,L,ton,Vprev)
	#calculates the next iteration
	result = float( 
					math.pow( (1 - p), N) * N * Vlittle * d 
					+ (1- (math.pow((1 - p), N)) ) * 
					1/( 1/(N*Vlittle*d)  + ton/L  )
					)

	#check if recurision is done
	if (math.fabs(Vprev - result) < 0.00001) or iterations > 50:
		return result
	else:
		return calc_Vn_recursive(r, L, ton, result, N, 
								Vlittle, d, iterations + 1)

def calc_r(ton, kAtt):
	result = ton/ ( ton + math.pow(kAtt, -1) )
	return result
	

def calc_ton(kAd, kT, ATP):
	result = ( (1/float(kAd)) + ( 1/(float(kT) *float(ATP)) ) )
	return result
	

def calc_kDet(ton):
	result = 1/ton
	return result
	

def calc_Vlittle(kAtt, kAd, kT, ATP):
	result = math.pow(  (1/float(kAtt)) + (1/float(kAd)) + (1/float((kT*ATP))) 
						, -1 )
	return result
	

def calc_p(r, L, ton, V):
	result = r*math.pow(  math.e,( -1*L/float(ton*V) )  )
	return result




