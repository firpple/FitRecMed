#!/usr/bin/python

#Main driver for Pharm regression
#Author: Evan Su
#Version: 2.0
#Date: 6/19/2017
#tab space: 4


from lmfit  import Parameters, minimize #optimize parameters, levenberg
import matplotlib.pyplot as plt #grapher
import numpy as np #array handler
import math #math functions
import csv #excel reading format
import Tkinter #Menu
import time #Debugging
import os.path #check for valid files
#library within the folder
#contains equations to calculate recursion
import PharmEquations as equations#BioEquClass


##main function####################################################
#created for modularity
def main():
	#menu system
	master = Tkinter.Tk()
	menuhandler = NVMenuHandler(master)
	master.mainloop()
	return	


##class implementation#############################################
#class: NVDataSetClass
#purpose: container for all variables
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
residual(self,params, x_array, y_array)
readFrom(self,name) #name: name of file
set_*(self,value) #value: desired value to set
'''
class NVRsquaredDataSetClass(equations.BioEquClass):

	def __init__(self):
		equations.BioEquClass.__init__(self)

		#RSquared Value(coefficient of determination)
		self.cod = 0

		#calc curve
		self.x_fit =  np.empty(shape=(0)) 
		self.y_fit = np.empty(shape=(0))

		#creates empty arrays to read observed values
		self.N_array = np.empty(shape=(0))
		self.V_array = np.empty(shape=(0))



	def __str__(self):
		string = ''

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
		string = string +"\nR^2: %.3f" % (self.cod)
		string = string +"\n"

		return string


	def calcMin(self):
		#calc correct parameters if possible
		if (self.N_array.size > 0) and (self.V_array.size > 0):

			#calcs mim parameters
			out = minimize(self.residual,self.param, args=(self.N_array, self.V_array))

			#reassign variables

			self.param = out.params
			self.cod = 1- out.redchi/np.var(self.V_array)

			#print 1- out.redchi/np.var(self.V_array, ddof=2)
			#print out.redchi
			#print out.chisqr
		self.recalc()

		return

	#function that works with lmfit's minimize function
	def residual(self, params, x_array, y_array):
	
		#pulldata
		kAtt = params['kAtt']
		d = params['d']
		L = params['L']
		kAd = params['kAd']
		kT = params['kT']
		ATP = params['ATP']

		#recalc dependend variables
		ton = equations.calc_ton(kAd, kT, ATP)
		kDet = equations.calc_kDet(ton)
		r = equations.calc_r(ton, kAtt)
		Vlittle = equations.calc_Vlittle(kAtt,kAd, kT, ATP)

		#calcs model numbers
		Vprev = 0.001
		model = np.empty(shape =(0))
		for i in x_array:
			model = np.append(model,equations.calc_Vn(  r, L, ton, 
														Vprev, i, Vlittle, d))

		return np.subtract(model,y_array)




	def readFrom(self, name):
		#creates empty arrays to read observed values
		self.N_array = np.empty(shape=(0))
		self.V_array = np.empty(shape=(0))

		#reads values from file
		ifile = open(name, "rb")
		reader = csv.reader(ifile)
		reader.next()
		for row in reader:
			self.N_array = np.append(self.N_array, float(row[0]))
			self.V_array = np.append(self.V_array, float(row[1]))
			#print self.V_array
		
		ifile.close()
		
		if(self.V_array[0] < 50):
			self.V_array = np.multiply(self.V_array, 1000)

	
	def get_L(self):
		return self.param['L'].value
	def get_kAd(self):
		return self.param['kAd'].value
	def get_RSquared(self):
		return self.cod

#class:
#purpose: container for all variables
#contents:
##variables: 
'''
self.kAttVary
self.dVary
self.LVary
self.kAdVary
self.kTVary
self.ATPVary
self.master #menu 
self.set #dataset

#menu components
self.l_variables
self.l_value
self.l_fixed

#kAtt
self.l_kAtt
self.e_kAtt
self.c_kAtt

#d
self.l_d
self.e_d
self.c_d

#L
self.l_L
self.e_L
self.c_L

#kAd
self.l_kAd
self.e_kAd
self.c_kAd


#kT
self.l_kT
self.e_kT
self.c_kT

#ATP
self.l_ATP
self.e_ATP
self.c_ATP

#filename
self.l_fileName
self.e_fileName

#calculate button
self.c_kAtt
self.b_calc
'''
##functions:
'''
def __init__(self, master) #master: Tkinter.Tk() #menu
def ignore(self):
def calc_button_press(self):
def dataTransferTo(self):
def makeGraph(self):
'''
class NVMenuHandler:
	#local menu system function
	def __init__(self, master):

		self.LMin = 0
		self.LMax = 10

		self.kAdMin = 0
		self.kAdMax = 10
		self.kAdStep = .1

		self.RMatrix = np.zeros(shape=(0,0))

		self.set = NVRsquaredDataSetClass()
		############creating menu system##############
	
		self.master = master
		self.master.wm_title("Rsquared Contour Maker")

		#title labels

		self.l_variables = Tkinter.Label(master, text = "Variables")
		self.l_variables.grid(row = 0, column = 0)
		self.l_value = Tkinter.Label(master, text = "Values")
		self.l_value.grid(row =0, column = 1)

		#kAtt
		self.l_kAtt = Tkinter.Label(master, text = "kAtt (s^-1)", justify = 'left')
		self.l_kAtt.grid(row = 1, column = 0, sticky = Tkinter.W)

		self.e_kAtt = Tkinter.Entry(master)
		self.e_kAtt.grid(row = 1, column = 1)
		self.e_kAtt.insert(0,10)


		#d
		self.l_d = Tkinter.Label(master, text = "d (nm)", justify = 'left')
		self.l_d.grid(row = 2, column = 0, sticky = Tkinter.W)

		self.e_d = Tkinter.Entry(master)
		self.e_d.grid(row = 2, column = 1)
		self.e_d.insert(0,8)



		#L
		self.l_L = Tkinter.Label(master, text = "L Range (nm)", justify = 'left')
		self.l_L.grid(row = 3, column = 0, sticky = Tkinter.W)

		self.e_min_L = Tkinter.Entry(master)
		self.e_min_L.grid(row = 3, column = 1)
		self.e_min_L.insert(0,0)

		self.e_max_L = Tkinter.Entry(master)
		self.e_max_L.grid(row = 3, column = 2)
		self.e_max_L.insert(0,100)

		#kAd
		self.l_kAd = Tkinter.Label(master, text = "k-Ad Range(s^-1)", justify = 'left')
		self.l_kAd.grid(row = 5, column = 0, sticky = Tkinter.W)

		self.e_min_kAd = Tkinter.Entry(master)
		self.e_min_kAd.grid(row = 5, column = 1)
		self.e_min_kAd.insert(0,50)

		self.e_max_kAd = Tkinter.Entry(master)
		self.e_max_kAd.grid(row = 5, column = 2)
		self.e_max_kAd.insert(0,100)

		self.l_Step_kAd = Tkinter.Label(master, text = "k-Ad Step Size", justify = 'left')
		self.l_Step_kAd.grid(row = 6, column = 0, sticky = Tkinter.W)

		self.e_Step_kAd = Tkinter.Entry(master)
		self.e_Step_kAd.grid(row = 6, column = 1)
		self.e_Step_kAd.insert(0,1)

	
		#kT
		self.l_kT = Tkinter.Label(master, text = "kT(uM^-1 * s^-1)", justify = 'left')
		self.l_kT.grid(row = 7, column = 0, sticky = Tkinter.W)

		self.e_kT = Tkinter.Entry(master)
		self.e_kT.grid(row = 7, column = 1)
		self.e_kT.insert(0,2)


		#ATP
		self.l_ATP = Tkinter.Label(master, text = "ATP (uM)", justify = 'left')
		self.l_ATP.grid(row = 8, column = 0, sticky = Tkinter.W)

		self.e_ATP = Tkinter.Entry(master)
		self.e_ATP.grid(row = 8, column = 1)
		self.e_ATP.insert(0,1000)


		#filename
		self.l_fileName = Tkinter.Label(master, 
										text = "Data File Name",
										justify = 'left')
		self.l_fileName.grid(row = 9, column = 0, sticky = Tkinter.W)

		self.e_fileName = Tkinter.Entry(master)
		self.e_fileName.grid(row = 9, column = 1)
		self.e_fileName.insert(0, "sk v vs n.csv")

		#save filename
		self.l_saveName = Tkinter.Label(master, 
										text = "Save File Name",
										justify = 'left')
		self.l_saveName.grid(row = 10, column = 0, sticky = Tkinter.W)

		self.e_saveName = Tkinter.Entry(master)
		self.e_saveName.grid(row = 10, column = 1)
		self.e_saveName.insert(0, "Output.csv")

		#calculate button

		self.b_calc = Tkinter.Button(   master, 
										text = "Calculate", 
										command = self.calc_button_press)
		self.b_calc.grid(row = 11, column = 0, columnspan =1)

		#save button

		self.b_save = Tkinter.Button(   master, 
										text = "Save File", 
										command = self.save_button_press)
		self.b_save.grid(row = 11, column = 1, columnspan =1)


		return


	

	def ignore(self):
		print "ignore"
		return



	def calc_button_press(self):
		#validates input
		if not self.validateInput():
			return
		#disable button
		##kind of works
		self.b_calc.configure(text = "Please Wait", command = self.ignore)
		self.b_save.configure(text = "Please Wait", command = self.ignore)

		self.master.update()

		self.primaryCalc()

		self.printSetTerminalOutput()

		#display graph
		self.makeGraph()

		#reenable button
		self.b_calc.configure(	text = "Calculate",
							    command = self.calc_button_press)
		self.b_save.configure(	text = "Save File",
							    command = self.save_button_press)

		#updates GUI
		self.master.update()
		return


	def save_button_press(self):
		#validates input
		if not self.validateInput():
			return
		#disable button
		##kind of works
		self.b_calc.configure(text = "Please Wait", command = self.ignore)
		self.b_save.configure(text = "Please Wait", command = self.ignore)
		self.master.update()

		self.primaryCalc()
		#Formats output
		WriteMatrix = self.FormatSaveOutput()

		#writing is performed here
		fout = open('OUTPUT/' + self.e_saveName.get(), 'w')
		fout.truncate()
		
		for row in WriteMatrix:
			for col in row:
				if(type(col) is str):
					fout.write(col)
				fout.write(',')
			fout.write('\n')
		
		fout.close()
		print "file saved"

		self.printSetTerminalOutput()

		self.makeGraph()

		#reenable button
		self.b_calc.configure(	text = "Calculate",
							    command = self.calc_button_press)
		self.b_save.configure(	text = "Save File",
							    command = self.save_button_press)

		#updates GUI
		self.master.update()
		return

	def FormatSaveOutput(self): #change this
		#not 100% memory effiecent, extra 2 blank columns exist 
		#for the purpose of formating

		#sets column index for variables
		NCol = 0
		VCol = 1

		KCol = 3
		LCol = 4
		RCol = 5

		VariableCol = 7
		ValueCol = 8
		SEECol = 9

		#indexes for navigating through output Matrix
		heightIndex = 0
		widthIndex = 0

		arrayHeight = 20 #default height of output matrix
		if(self.set.N_array.size > arrayHeight):
			arrayHeight = self.set.N_array.size
		if(len(self.RMatrix)*3 > arrayHeight):
			arrayHeight = len(self.RMatrix)*3

		arrayHeight += 2 #+2 for the header

		arrayWidth = 10 #default width of output matrix

		#creates matrix
		outputMatrix = [[0 for width in range(arrayWidth)]for height in range(arrayHeight)]

		#filling matrix
		####original data###################

		#data labels
		outputMatrix[0][NCol] = 'N (uM)'
		outputMatrix[0][VCol] = 'V (um/s)'

		#Setup
		heightIndex = 1
		NIter = np.nditer(self.set.N_array)
		VIter = np.nditer(self.set.V_array)

		#data points
		while not NIter.finished:
			outputMatrix[heightIndex][NCol] = str(float(NIter[0]))
			outputMatrix[heightIndex][VCol] = str(float(VIter[0])/1000)
			NIter.iternext()
			VIter.iternext()
			heightIndex += 1

		####Rsquared data######################
		#data labels
		outputMatrix[0][KCol] = 'K-Ad (s^-1)'
		outputMatrix[0][LCol] = 'L (nm)'
		outputMatrix[0][RCol] = 'R Squared'

		#data input
		heightIndex = 1
		for index in range(len(self.RMatrix)):
			for jndex in range(len(self.RMatrix[index])/2):
				outputMatrix[heightIndex][KCol] = str(float(index*self.kAdStep + self.kAdMin))
				outputMatrix[heightIndex][LCol] = str(self.RMatrix[index][jndex*2])
				outputMatrix[heightIndex][RCol] = str(self.RMatrix[index][jndex*2 + 1])
				heightIndex += 1
		####Parameter Data##################
		outputMatrix[0][VariableCol] = 'Variable Name'
		outputMatrix[0][ValueCol] = 'Value'
		outputMatrix[0][SEECol] = 'S.E.E.'
		
		#can be automized with for loop?
		#manual entry here as the output for each is very specific
		
		heightIndex = 1
		#automaize attempt
		paramHeight = 4
		paramWidth = 2
		paramArray = [[0 for x in range(paramWidth)]for x in range(paramHeight)]

		#[row][0] = CodeName, [row][1] = outputName
		paramArray[0] = ['kAtt', 'kAtt(s^-1)']
		paramArray[1] = ['d', 'd(nm)']
		paramArray[2] = ['kT', 'kT(uM^-1 * s^-1)']
		paramArray[3] = ['ATP', 'ATP(uM)']

		#print paramArray

		#setup for output
		CODENAME = 0
		OUTPUTNAME = 1

		for row in paramArray:
			outputMatrix[heightIndex][VariableCol] = row[OUTPUTNAME]
			outputMatrix[heightIndex][ValueCol] = str(
					self.set.param[row[CODENAME]].value)

			if self.set.param[str(row[CODENAME])].vary:
				outputMatrix[heightIndex][SEECol] = str(
						self.set.param[row[CODENAME]].stderr)
			else:
				outputMatrix[heightIndex][SEECol] = '(Fixed)'
			heightIndex += 1
		

		
		#outputMatrix
		#print outputMatrix
		return outputMatrix

	def dataTransferTo(self):
		#sets initial value
		self.set.set_kAtt(float(self.e_kAtt.get()))
		self.set.set_d(float(self.e_d.get()))
		self.set.set_kT(float(self.e_kT.get()))
		self.set.set_ATP(float(self.e_ATP.get()))

		#sets fixed or floating
		self.set.set_kAttVary(False)
		self.set.set_dVary(False)
		self.set.set_kTVary(False)
		self.set.set_ATPVary(False)
		return



	def makeGraph(self):#change this
		
		#plotting the numbers

		##plot min max

		#data input
		XPoints = np.zeros(shape=(len(self.RMatrix)*3))
		YPoints = np.zeros(shape=(len(self.RMatrix)*3))
		heightIndex = 0 
		for index in range(len(self.RMatrix)):
			for jndex in range(len(self.RMatrix[index])/2):
				XPoints[heightIndex] = float(index*self.kAdStep + self.kAdMin)
				YPoints[heightIndex] = self.RMatrix[index][jndex*2]
				heightIndex += 1 


		##plot original data
		ax = plt.figure()
		plt.scatter(XPoints, YPoints, color = 'b')
		#sets borders of graph
		plt.xlim(min(XPoints) - 5, max(XPoints) + 5)
		plt.ylim(min(YPoints) - .1, max(YPoints) + .1)

		#enables axis lines
		plt.axhline(y =0, color = 'k')
		plt.axvline(x = 0, color = 'k')

		#sets title
		plt.suptitle(self.e_fileName.get())

		#sets axis labels
		plt.xlabel('k-AD(s^-1)')
		plt.ylabel('L(nm)')


		#sets the info legend texts
		tempString = 'kAtt(s^-1):' + \
					str(self.set.param['kAtt'].value) + \
					'\nd(nm):'+ \
					str(self.set.param['d'].value) +\
					'\nkT(uM^-1 * s^-1):' +\
					str(self.set.param['kT'].value) +\
					'\nATP(uM):'+\
					str(self.set.param['ATP'].value)
					
		plt.annotate(tempString, 
					horizontalalignment = 'right', 
					verticalalignment = 'top', 
					xy = (1,1), 
					xycoords = 'axes fraction',
					fontweight = 'bold')
		plt.show(block = False)


	def printSetTerminalOutput(self):
		print '##############End Set#########################'
		print self.set
		print '##############End Set#########################'

	def primaryCalc(self):		
		self.set.readFrom('INPUT/' + self.e_fileName.get())
		print 'Read file'
		#sends the data to the dataset#sends the data to the dataset
		self.dataTransferTo()

		self.LMin = float(self.e_min_L.get())
		self.LMax = float(self.e_max_L.get())
		#self.LStep = float(self.e_Step_L.get())

		self.kAdMin = float(self.e_min_kAd.get())
		self.kAdMax = float(self.e_max_kAd.get())
		self.kAdStep = float(self.e_Step_kAd.get()) 
		#LElements = (self.LMax - self.LMin)/self.LStep
		kAdElements = (self.kAdMax - self.kAdMin)/self.kAdStep


		self.RMatrix = np.zeros(shape= (int(kAdElements), 6))
		self.set.set_L(float(100 + .000001))
		self.set.set_kAd(float(1000 + .000001))
		
		#put smart algrorithm here
		for index in range(len(self.RMatrix)):
			self.set.set_LVary(True)
			self.set.set_L(1)
			kAdVal = index*self.kAdStep + self.kAdMin
			if kAdVal == 0:
				kAdVal = 0.0001
			print kAdVal
			self.set.set_kAd(kAdVal)
			self.dataTransferTo()
			self.set.calcMin()
			#print index
			#print self.set
			#print 'pizza'
			if self.set.get_RSquared() < 0:
				self.RMatrix[index][0] = (self.LMin+self.LMax) / 2.0
				self.RMatrix[index][1] = 0
				self.RMatrix[index][2] = self.LMin
				self.RMatrix[index][3] = 0
				self.RMatrix[index][4] = self.LMax
				self.RMatrix[index][5] = 0
				continue #skips this iteration as nothing more can be done
			#print (self.set.get_L(),self.set.get_RSquared())
			self.RMatrix[index][0] = self.set.get_L()
			self.RMatrix[index][1] = self.set.get_RSquared()
			##<><><><<><><><>><><><><><<
			#calc lower bound
			lower_L = self.LMin
			upper_L = self.set.get_L()
			mid_L = 0
			result_R = 0 
			for x in range(16):
				mid_L = (lower_L+upper_L) / 2.0
				self.set.set_LVary(False)
				self.set.set_L(mid_L)
				self.set.calcMin()
				result_R = self.set.get_RSquared()
				#print (lower_L, mid_L ,upper_L, result_R)
				if result_R < 0.0001 and result_R > 0.0:
					break
				if result_R > 0:
					upper_L = mid_L
				else:
					lower_L = mid_L


			#save results
			'''
			self.set.set_L(upper_L)
			self.set.recalc()
			result_R = self.set.get_RSquared()
			'''
			#print (upper_L,result_R)
			
			self.RMatrix[index][2] = mid_L
			self.RMatrix[index][3] = result_R
		
			#calc upper bound
			lower_L = self.set.get_L()
			upper_L = self.LMax
			mid_L = 0
			result_R = 0 
			for x in range(16):
				mid_L = (lower_L+upper_L) / 2.0
				self.set.set_LVary(False)
				self.set.set_L(mid_L)
				self.set.calcMin()
				result_R = self.set.get_RSquared()
				#print (lower_L, mid_L ,upper_L, result_R)
				if result_R < 0.0001 and result_R > 0.0:
					break
				if result_R > 0:
					lower_L = mid_L
				else:
					upper_L = mid_L

			#save results
			'''
			self.set.set_L(lower_L)
			self.set.recalc()
			result_R = self.set.get_RSquared()
			'''
			#print (lower_L,result_R)
			
			self.RMatrix[index][4] = mid_L
			self.RMatrix[index][5] = result_R

			
			print (kAdVal,self.RMatrix[index])

		#calcs parameters
		self.set.calcMin()
		#print self.set
		print self.set.get_RSquared()

		return


	def validateInput(self):
		testPassed = True
		try:
			test = float(self.e_kAtt.get())
		except ValueError:
			print 'Invalid kAtt value!'
			testPassed = False

		try:
			test = float(self.e_d.get())
		except ValueError:
			print 'Invalid d value!'
			testPassed = False


		try:
			test = float(self.e_min_L.get())
		except ValueError:
			print 'Invalid min L value!'
			testPassed = False
		try:
			test = float(self.e_max_L.get())
		except ValueError:
			print 'Invalid max L value!'
			testPassed = False
		try:
			test = float(self.e_Step_kAd.get())
		except ValueError:
			print 'Invalid k-Ad step value!'
			testPassed = False
		try:
			test = float(self.e_min_kAd.get())
		except ValueError:
			print 'Invalid min k-Ad value!'
			testPassed = False
		try:
			test = float(self.e_max_kAd.get())
		except ValueError:
			print 'Invalid max k-Ad value!'
			testPassed = False

		try:
			test = float(self.e_kT.get())
		except ValueError:
			print 'Invalid kT value!'
			testPassed = False

		try:
			test = float(self.e_ATP.get())
		except ValueError:
			print 'Invalid ATP value!'
			testPassed = False

		if(not (os.path.isfile('INPUT/' + self.e_fileName.get())) ):
			print 'Invalid Data File Name!'
			testPassed = False
		
		try:
			f = open('OUTPUT/' + self.e_saveName.get(), 'w')
			f.close()
		except IOError:
			print 'Invalid Save File Name!'
			testPassed = False
		
		if not testPassed:
			print 'Please check your input!\n'

		return testPassed


##Runtime Function Calls###########################################
#this is what actual gets ran
#main()
main()
print("hello world")


