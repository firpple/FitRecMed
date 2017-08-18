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
class NVDataSetClass(equations.BioEquClass):

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


	def calcCurve(self):
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

		#calculates best fit curve
		self.x_fit = np.array(range(1,  int( max(self.N_array)),  1))
		self.y_fit = np.empty(shape =(0))
		for i in self.x_fit:
		#sorry ranout of room
			self.y_fit =np.append(self.y_fit,
									equations.calc_Vn(  self.r, self.param['L'],
													self.ton, 0.001,  
													i, self.Vlittle, 
														self.param['d']))
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
		self.kAttFixed = Tkinter.BooleanVar()
		self.dFixed = Tkinter.BooleanVar()
		self.LFixed = Tkinter.BooleanVar()
		self.kAdFixed = Tkinter.BooleanVar()
		self.kTFixed = Tkinter.BooleanVar()
		self.ATPFixed = Tkinter.BooleanVar()


		self.set = NVDataSetClass()
		############creating menu system##############
	
		self.master = master
		self.master.wm_title("N vs V regression fitter")

		#title labels

		self.l_variables = Tkinter.Label(master, text = "Variables")
		self.l_variables.grid(row = 0, column = 0)
		self.l_value = Tkinter.Label(master, text = "Values")
		self.l_value.grid(row =0, column = 1)
		self.l_fixed = Tkinter.Label(master, text = "Fixed?")
		self.l_fixed.grid(row = 0, column = 2)

		#kAtt
		self.l_kAtt = Tkinter.Label(master, text = "kAtt (s^-1)", justify = 'left')
		self.l_kAtt.grid(row = 1, column = 0, sticky = Tkinter.W)

		self.e_kAtt = Tkinter.Entry(master)
		self.e_kAtt.grid(row = 1, column = 1)
		self.e_kAtt.insert(0,1)

		self.c_kAtt = Tkinter.Checkbutton(master, variable = self.kAttFixed)
		self.c_kAtt.grid(row = 1, column = 2)

		#d
		self.l_d = Tkinter.Label(master, text = "d (nm)", justify = 'left')
		self.l_d.grid(row = 2, column = 0, sticky = Tkinter.W)

		self.e_d = Tkinter.Entry(master)
		self.e_d.grid(row = 2, column = 1)
		self.e_d.insert(0,8)

		self.c_d = Tkinter.Checkbutton(master, variable = self.dFixed)
		self.c_d.grid(row = 2, column = 2)


		#L
		self.l_L = Tkinter.Label(master, text = "L (nm)", justify = 'left')
		self.l_L.grid(row = 3, column = 0, sticky = Tkinter.W)

		self.e_L = Tkinter.Entry(master)
		self.e_L.grid(row = 3, column = 1)
		self.e_L.insert(0,25.8)

		self.c_L = Tkinter.Checkbutton(master, variable = self.LFixed)
		self.c_L.grid(row = 3, column = 2)

		#kAd
		self.l_kAd = Tkinter.Label(master, text = "k-Ad(s^-1)", justify = 'left')
		self.l_kAd.grid(row = 4, column = 0, sticky = Tkinter.W)

		self.e_kAd = Tkinter.Entry(master)
		self.e_kAd.grid(row = 4, column = 1)
		self.e_kAd.insert(0,125)

		self.c_kAd = Tkinter.Checkbutton(master, variable = self.kAdFixed)
		self.c_kAd.grid(row = 4, column = 2)
		
	
		#kT
		self.l_kT = Tkinter.Label(master, text = "kT(uM^-1 * s^-1)", justify = 'left')
		self.l_kT.grid(row = 5, column = 0, sticky = Tkinter.W)

		self.e_kT = Tkinter.Entry(master)
		self.e_kT.grid(row = 5, column = 1)
		self.e_kT.insert(0,2)

		self.c_kT = Tkinter.Checkbutton(master, variable = self.kTFixed)
		self.c_kT.grid(row = 5, column = 2)

		#ATP
		self.l_ATP = Tkinter.Label(master, text = "ATP (uM)", justify = 'left')
		self.l_ATP.grid(row = 6, column = 0, sticky = Tkinter.W)

		self.e_ATP = Tkinter.Entry(master)
		self.e_ATP.grid(row = 6, column = 1)
		self.e_ATP.insert(0,1000)

		self.c_ATP = Tkinter.Checkbutton(master, variable = self.ATPFixed)
		self.c_ATP.grid(row = 6, column = 2)


		#filename
		self.l_fileName = Tkinter.Label(master, 
										text = "Data File Name",
										justify = 'left')
		self.l_fileName.grid(row = 7, column = 0, sticky = Tkinter.W)

		self.e_fileName = Tkinter.Entry(master)
		self.e_fileName.grid(row = 7, column = 1)
		self.e_fileName.insert(0, "sk v vs n.csv")

		#save filename
		self.l_saveName = Tkinter.Label(master, 
										text = "Save File Name",
										justify = 'left')
		self.l_saveName.grid(row = 8, column = 0, sticky = Tkinter.W)

		self.e_saveName = Tkinter.Entry(master)
		self.e_saveName.grid(row = 8, column = 1)
		self.e_saveName.insert(0, "Output.csv")

		#calculate button

		self.b_calc = Tkinter.Button(   master, 
										text = "Calculate", 
										command = self.calc_button_press)
		self.b_calc.grid(row = 9, column = 0, columnspan =1)

		#save button

		self.b_save = Tkinter.Button(   master, 
										text = "Save File", 
										command = self.save_button_press)
		self.b_save.grid(row = 9, column = 1, columnspan =1)


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
		#reenable button
		self.b_calc.configure(	text = "Calculate",
							    command = self.calc_button_press)
		self.b_save.configure(	text = "Save File",
							    command = self.save_button_press)

		#updates GUI
		self.master.update()
		return

	def FormatSaveOutput(self):
		#not 100% memory effiecent, extra 2 blank columns exist 
		#for the purpose of formating

		#sets column index for variables
		NCol = 0
		VCol = 1

		XCol = 3
		YCol = 4

		VariableCol = 6
		ValueCol = 7
		SEECol = 8

		#indexes for navigating through output Matrix
		heightIndex = 0
		widthIndex = 0

		arrayHeight = 20 #default height of output matrix
		if(self.set.N_array.size > arrayHeight):
			arrayHeight = self.set.N_array.size
		if(self.set.x_fit.size > arrayHeight):
			arrayHeight = self.set.x_fit.size 

		arrayHeight += 2 #+2 for the header

		arrayWidth = 9 #default width of output matrix

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

		####Curve data######################
		#data Labels
		outputMatrix[0][XCol] = 'Fit X (uM)'
		outputMatrix[0][YCol] = 'Fit Y (um/s)'

		#Setup
		heightIndex = 1
		XIter = np.nditer(self.set.x_fit)
		YIter = np.nditer(self.set.y_fit)

		#Data Points
		while not XIter.finished:
			outputMatrix[heightIndex][XCol] = str(float(XIter[0]))
			outputMatrix[heightIndex][YCol] = str(float(YIter[0])/1000)
			XIter.iternext()
			YIter.iternext()
			heightIndex += 1

		####Parameter Data##################
		outputMatrix[0][VariableCol] = 'Variable Name'
		outputMatrix[0][ValueCol] = 'Value'
		outputMatrix[0][SEECol] = 'S.E.E.'
		
		#can be automized with for loop?
		#manual entry here as the output for each is very specific
		
		heightIndex = 1
		#automaize attempt
		paramHeight = 6
		paramWidth = 2
		paramArray = [[0 for x in range(paramWidth)]for x in range(paramHeight)]

		#[row][0] = CodeName, [row][1] = outputName
		paramArray[0] = ['kAtt', 'kAtt(s^-1)']
		paramArray[1] = ['d', 'd(nm)']
		paramArray[2] = ['L', 'L(nm)']
		paramArray[3] = ['kAd', 'k-Ad(s^-1)']
		paramArray[4] = ['kT', 'kT(uM^-1 * s^-1)']
		paramArray[5] = ['ATP', 'ATP(uM)']

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
		

		#print extra's variables stuffs#####
		#ton
		outputMatrix[heightIndex][VariableCol] = 'ton (s)'
		outputMatrix[heightIndex][ValueCol] = str(self.set.ton)
		heightIndex += 1

		#kDet
		outputMatrix[heightIndex][VariableCol] = 'kDet (S^-1)'
		outputMatrix[heightIndex][ValueCol] = str(self.set.kDet)
		heightIndex += 1

		#r
		outputMatrix[heightIndex][VariableCol] = 'r'
		outputMatrix[heightIndex][ValueCol] = str(self.set.r)
		heightIndex += 1

		#Vlittle
		outputMatrix[heightIndex][VariableCol] = 'v (ATPase)(S^-1)'
		outputMatrix[heightIndex][ValueCol] = str(self.set.Vlittle)
		heightIndex += 1

		#R^2
		outputMatrix[heightIndex][VariableCol] = 'R^2'
		outputMatrix[heightIndex][ValueCol] = str(self.set.cod)
		heightIndex += 1
		
		#outputMatrix
		return outputMatrix

	def dataTransferTo(self):
		#sets initial value
		self.set.set_kAtt(float(self.e_kAtt.get()))
		self.set.set_d(float(self.e_d.get()))
		self.set.set_L(float(self.e_L.get()))
		self.set.set_kAd(float(self.e_kAd.get()))
		self.set.set_kT(float(self.e_kT.get()))
		self.set.set_ATP(float(self.e_ATP.get()))

		#sets fixed or floating
		self.set.set_kAttVary(not self.kAttFixed.get())
		self.set.set_dVary(not self.dFixed.get())
		self.set.set_LVary(not self.LFixed.get())
		self.set.set_kAdVary(not self.kAdFixed.get())
		self.set.set_kTVary(not self.kTFixed.get())
		self.set.set_ATPVary(not self.ATPFixed.get())
		return



	def makeGraph(self):
		
		#plotting the numbers
		##plot original data
		ax = plt.figure()
		plt.scatter(self.set.N_array, 
					np.divide(self.set.V_array,1000), color = 'r')
		##plot best fit curve
		plt.scatter(self.set.x_fit, np.divide(self.set.y_fit, 1000), color = 'b')

		#sets borders of graph
		plt.xlim(min(self.set.N_array) - 10, max(self.set.N_array) + 5)
		plt.ylim(  min(np.divide(self.set.V_array,1000)) - .5, 
					max(np.divide(self.set.V_array, 1000)) + 1)

		#enables axis lines
		plt.axhline(y =0, color = 'k')
		plt.axvline(x = 0, color = 'k')

		#sets title
		plt.suptitle(self.e_fileName.get())

		#sets axis labels
		plt.xlabel('N(uM)')
		plt.ylabel('V(um/s)')


		#sets the info legend texts
		plt.annotate(self.set, 
					horizontalalignment = 'right', 
					verticalalignment = 'bottom', 
					xy = (1,0), 
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
		#sends the data to the dataset
		self.dataTransferTo()
		self.set.recalc()
		self.set.calcCurve()
		print 'Finish Levenburg Marquard'
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
			test = float(self.e_L.get())
		except ValueError:
			print 'Invalid L value!'
			testPassed = False

		try:
			test = float(self.e_kAd.get())
		except ValueError:
			print 'Invalid kAd value!'
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


