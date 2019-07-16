#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copyright (C) S. Merkel, Universite de Lille, France

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


# System functions, to manipulate command line arguments
import sys
import optparse
import argparse
import os.path
# Time functions
import time
import datetime
# string module contains a number of functions that are useful for manipulating strings
import string
# Mathematical stuff (for data array)
import math
import numpy as np
import scipy as sp
import scipy.optimize
# StringIO used for strings behaves like a file
import linecache
# to search for a certain line in the file
from StringIO import StringIO
# Plotting library
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
font = {'family' : 'sans-serif','weight' : 'normal', 'size'   : 12}
matplotlib.rc('font', **font)
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon


#################################################################
#
# Global variables, some are set later when starting the program
#
#################################################################

# Pe - pressure of point on pv / (pv+ppv) phase boundary, in GPa
# Te - temperature of point on pv / (pv+ppv) phase boundary, in K
# clapeyron - Clapeyron slope in GPa / K
Pe = 0.0
Te = 0.0
clapeyron = 0.0

# Limits for pressure in plot
Pmin = 110.
Pmax = 140.

# Wait time when showing the plots
waittime = 1.

#################################################################
#
# Support functions
#
#################################################################

#####
# Function to convert depth into pressures for D''
#   depth is in km, pressure is in GPa
# 
def pressure(x): 
	p = []
	for i in range(0,len(x)):
		if (x[i] > 2891.):
			p.append(-200.94 + 0.1265 * x[i] - 3.5384e-06 * np.power(x[i],2))
		else:
			p.append(5.1676 + 0.031296 * x[i] + 4.7313e-06 * np.power(x[i],2))
	return p


#####
# Shear Transformation model
#  x[:,0] is temperature in K
#  x[:,1] is pressure in GPa
# A is ln(k2) constant in paper
# E is activation energy
# C is activation volume
# 
def shearModel(x,E,A,C): 
	DeltaP = (x[:,1]-((Pe-clapeyron*Te)+(x[:,1]+clapeyron*x[:,0]))/2)
	return  np.exp(A)*np.exp(E/(x[:,0]*8.314))/np.sinh(C*DeltaP/(x[:,0]*8.314))



#####
# Nucleation and growth model
#  x[:,0] is temperature in K
#  x[:,1] is pressure in GPa
#  x[:,2] is DeltaV is A3 (maybe, there are weird conversions)

# A is ln(k1) constant in paper
# E is activation energy
# B is activation volume
# 
def nuclGrowthModel(x,E,A,B):
	DeltaP = (x[:,1]-((Pe-clapeyron*Te)+(x[:,1]+clapeyron*x[:,0]))/2)
	preC = np.exp(A)
	Arhenius = np.exp((E-B*DeltaP)/(x[:,0]*8.314))
	DistanceTerm = (1- np.exp(-x[:,2]*1E-24*DeltaP/(x[:,0]*1.32E-23)))
	return preC*(1/x[:,0])*Arhenius/DistanceTerm
	
	# Old version from Chris
	# ureturn np.exp(A)*(1/x[:,0])*np.exp((E-B*(x[:,1]-((Pe-clapeyron*Te)+(x[:,1]+clapeyron*x[:,0]))/2))/(x[:,0]*8.314))/(1- np.exp(-x[:,3]*1E-24*((x[:,1]-((Pe-clapeyron*Te)+(x[:,1]+clapeyron*x[:,0]))/2))/(x[:,0]*1.32E-23)))


#################################################################
#
# Plotting subroutines
#
#################################################################


def plotshearModel(E, A, C, Pexp, Texp, cold, warm, hot,factor,noneobservation):
	"""
	Plot results of the Shear transformation model 
	- E: activation energy 
	- A: constant
	- C: activation volume
	- Pexp, Texp: experimental data
	- cold, warm, hot: geotherms
	- factor - multiplication factor for times in plotting of kinetics (in log)
	- noneobservation - File with negative observations (can be set to none)
	"""

	# Build a grid on the explored range of 1/T and P
	inverseT2 = [1./(5000-5*m) for m in range(800)]
	P2 = [Pmin+m*(Pmax-Pmin)/800.for m in range(800)]
	inverseT2, P2 = np.meshgrid(inverseT2, P2)

	# Calculation of kinetics at each grid point
	Z = np.log10(np.exp(A)*np.exp(E/((1/inverseT2)*8.314))/np.sinh(C*(P2-((Pe-clapeyron*Te)+(P2+clapeyron*(1/inverseT2)))/2)/((1/inverseT2)*8.314)))+factor

	# Create a plot with matplotlib
	fig = plt.figure(figsize=(7.5,5))
	ax = fig.add_subplot(111)

	# Definition of phase boundary between pure Pv and Pv + pPv mix
	inverseT = [1./(5000-5*m) for m in range(800)]
	Tclap = []
	phaseboundary = []
	for i in range (0,len(inverseT)):
		Tclap.append(float(1/float(inverseT[i])))
		phaseboundary.append(float(clapeyron*(1/float(inverseT[i])-Te)+Pe))
	
	# Plot experimental data as red dots 
	plt.plot(Texp,Pexp,'bo',color='red')
	
	# Plot negative observations, if any
	if (noneobservation != None):
		f = open(noneobservation, 'r')
		datcontent = [line.strip() for line in f.readlines()]
		f.close()
		Tn = []
		Pn = []
		for line in datcontent:
			if (line[0] != '#'):
				a = line.split(' ')
				Tn.append(float(a[1]))
				Pn.append(float(a[0]))
		plt.plot(Tn,Pn,'bo',color='blue')
	
	# Plot phase boundary
	plt.plot(Tclap[:], phaseboundary[:],color='black')
	
	# Plot geotherms
	plt.plot(hot[:,1],hot[:,0],'--',color='red', zorder=5)
	plt.plot(warm[:,1],warm[:,0],'--',color='green', zorder=5)
	plt.plot(cold[:,1],cold[:,0],'--',color='blue', zorder=5)
	
    # Plot kinetics as contours
	levels = np.arange(-2, 10, 1)
	CS = ax.contourf(1/inverseT2, P2, Z, levels, origin='lower', cmap=cm.jet_r, extend='both')
	CS2 = ax.contour(CS, levels,
					colors='black',
					origin='lower',
					hold='on',linewidths=1)

	# Add label at some location so it looks better
	# Painful to find where to label. Removed it for now
	# manual_locations = [(1750,98),(1750,105), (1750,108),(1750,111),(1750,114),(1750, 116),(1750, 118.5),(1750, 121),(1750, 124.2), (1750, 127), (1750,130), (1750,132.8), (1750,135.75), (1750,138.4),(1750, 140)]
	# plt.clabel(CS2, inline=1, fontsize=12, fmt='%2.0f', color='black', manual=manual_locations)
	test = plt.clabel(CS2, inline=1, fontsize=11, fmt='%2.0f', color='black')
	#for x in test:
	#	print x.get_text(), x.get_position()


	# Adding the CMB to hide the rest later
	#ax.plot([1500, 4500], [135.2, 135.2], 'k-', zorder=4,color='black')
	rect = Rectangle((1500, 135.2), 3000, 4.8, edgecolor=None, facecolor='0.7', zorder=4)
	ax.add_artist(rect)
	
	# Arrow with D''
	verts = [(4500, 136), (4100, 136), (4100, 126), (4300, 124), (4500, 126), (4500, 136)]
	poly = Polygon(verts, facecolor='0.7', edgecolor=None, zorder=4)
	ax.add_patch(poly)

	
	
	#plot colorbar
	CB = plt.colorbar(CS, shrink=0.8, pad=0.2)
	CB.ax.yaxis.major.formatter._useMathText = True
	CB.set_label(r'$\mathrm{log(\tau)}$', rotation=270, fontsize=10)
	
	
	# Secondary labelling with twiny, for inverse temperatures
	ax1 = ax.twiny()
	ax1.plot(hot[:,1],hot[:,0],'--',color='none', zorder=5) # dummy trace
	ax1.set_xlim(1500,4500)
	ax1.invert_xaxis()
	ax1.set_xticks([1000./0.25, 1000./.30, 1000./.40, 1000./.5, 1000./0.6])
	ax1.set_xticklabels(["0.25", "0.30", "0.40", "0.50", "0.60"], fontsize=10)
	ax1.set_xlabel("$\mathrm{1000/T (K^{-1})}$", fontsize=10)
	ax1.xaxis.tick_bottom()
	ax1.xaxis.set_label_position("bottom")
	
    # Flipping the x-axis to have warm temperatures on the left and adding labels
	ax.set_xlim(1500,4500)
	ax.invert_xaxis()
	ax.set_xlabel('$\mathrm{Temperature (K)}$', fontsize=10)
	ax.xaxis.set_ticks_position('top')
	ax.xaxis.set_label_position('top')
	ax.tick_params(axis='both', which='major', labelsize=10)
  
	
	# Secondary y-axis
	ax2 = ax.twinx()
	ax2.set_ylim(Pmin,Pmax)
	ax2.set_yticks([135.2, 130, 124, 118.4, 112.5])
	ax2.set_yticklabels(["2891", "2800", "2700", "2600", "2500"])
	ax2.set_ylabel("Depth (km)", fontsize=10)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	ax2.tick_params(axis='both', which='major', labelsize=10)
	
	#py-axis label
	ax.set_ylim(Pmin,Pmax)
	ax.set_ylabel('Pressure (GPa)', fontsize=10)
   
	#plot title
	label = "Shear model\n" + r"$\mathrm{P_e = %d  GPa, T_e = %d  K, \alpha = %.1f  MPa/K}$" % (Pe, Te, clapeyron*1000)
	plt.text(0.05, 0.05, label, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes, fontsize=9, bbox=dict(edgecolor=(0.682, 0.682, 0.682), facecolor='white', alpha=1.0))
	
	plt.text(3000, 135.5, "Core-Mantle Boundary", horizontalalignment='center', verticalalignment='bottom', fontsize=10)
	plt.text(4300, 130, "D''\nlayer", horizontalalignment='center', verticalalignment='center', fontsize=10)

	#save plot in svg type with the name : "shearModel"
	plt.savefig('shearModel.svg')
	plt.savefig('shearModel.png')

	#plot
	plt.show(block=False)
	plt.pause(waittime)
	plt.close()
    
	return     



def plotnuclGrowthModel(E,A,B,Pexp,Texp,cold,warm,hot,factor,noneobservation):
	"""
	Plot results of the Nucleation and Growth transformation model 
	 - E: activation energy 
	 - A: constant
	 - C: activation volume
	 - Pexp, Texp: experimental data
	 - cold, warm, hot: geotherms
	 - factor - multiplication factor for times in plotting of kinetics (in log)
	- noneobservation - File with negative observations (can be set to none)
	"""
	
	# Build a grid on the explored range of 1/T and P
	inverseT2 = [1./(4500-3.75*m) for m in range(800)]
	P2 = [Pmin+m*(Pmax-Pmin)/800.for m in range(800)]
	inverseT2, P2 = np.meshgrid(inverseT2, P2)
	
	# Calculation of kinetics at each grid point
	DeltaP = (P2-((Pe-clapeyron*Te)+(P2+clapeyron*(1/inverseT2)))/2)
	arrhenius = np.exp((E-B*DeltaP)/((1/inverseT2)*8.314))
	Z = np.log10(np.exp(A)*(1/(1/inverseT2))*arrhenius/(1- np.exp(-1.4*1E-24*(DeltaP)/((1/inverseT2)*1.32E-23)))) + factor
	
	# Preparing a figure with matplotlib
	fig = plt.figure(figsize=(7.5,5))
	ax = fig.add_subplot(111)
	
	# Adding a line for phase boundary
	inverseT = [1./(4500-3.75*m) for m in range(800)]
	Tclap = []
	phaseboundary = []
	for i in range (0,len(inverseT)):
		Tclap.append(float(1/float(inverseT[i])))
		phaseboundary.append(float(clapeyron*(1/float(inverseT[i])-Te)+Pe))
	plt.plot(Tclap[:], phaseboundary[:],color='black')

	# Plotting experimental data
	plt.plot(Texp, Pexp,'bo',color='red')
	#Pexp = np.array(Pexp)
	#Texp = np.array(Texp)
	#DeltaP = (Pexp-((Pe-clapeyron*Te)+(Pexp+clapeyron*(Texp)))/2)
	#print DeltaP
	
	# Plot negative observations, if any
	if (noneobservation != None):
		f = open(noneobservation, 'r')
		datcontent = [line.strip() for line in f.readlines()]
		f.close()
		Tn = []
		Pn = []
		for line in datcontent:
			if (line[0] != '#'):
				a = line.split(' ')
				Tn.append(float(a[1]))
				Pn.append(float(a[0]))
		plt.plot(Tn,Pn,'bo',color='blue')
	
	# Geotherms
	plt.plot(hot[:,1],hot[:,0],'--',color='red', zorder=5)
	plt.plot(warm[:,1],warm[:,0],'--',color='green', zorder=5)
	plt.plot(cold[:,1],cold[:,0],'--',color='blue', zorder=5)

	# Contours for kinetics
	levels = np.arange(-2, 10, 1)
	CS = ax.contourf(1/inverseT2, P2, Z,levels, origin='lower', cmap=cm.jet_r, extend='both')
	CS2 = ax.contour(CS, levels,
					colors='black',
					origin='lower',
					hold='on',linewidths=1)
	
	# Add label at some location so it looks better
	# Painful to find where to label. Removed it for now
	# manual_locations = [(1750,98),(1750,105), (1750,108),(1750,111),(1750,114),(1750, 116),(1750, 118.5),(1750, 121),(1750, 124.2), (1750, 127), (1750,130), (1750,132.8), (1750,135.75), (1750,138.4),(1750, 140)]
	# plt.clabel(CS2, inline=1, fontsize=12, fmt='%2.0f', color='black', manual=manual_locations)
	test = ax.clabel(CS2, inline=1, fontsize=11, fmt='%2.0f', color='black')
	#for x in test:
	#	print x.get_text(), x.get_position()


	# Adding the CMB to hide the rest later
	#ax.plot([1500, 4500], [135.2, 135.2], 'k-', zorder=4,color='black')
	rect = Rectangle((1500, 135.2), 3000, 4.8, edgecolor=None, facecolor='0.7', zorder=4)
	ax.add_artist(rect)
	
	# Arrow for D''
	verts = [(4500, 136), (4100, 136), (4100, 126), (4300, 124), (4500, 126), (4500, 136)]
	poly = Polygon(verts, facecolor='0.7', edgecolor=None, zorder=4)
	ax.add_patch(poly)

	
	
	#plot colorbar
	CB = plt.colorbar(CS, shrink=0.8, pad=0.2)
	CB.ax.yaxis.major.formatter._useMathText = True
	CB.set_label(r'$\mathrm{log(\tau)}$', rotation=270, fontsize=10)
	
	
	# Secondary labelling with twiny, for inverse temperatures
	ax1 = ax.twiny()
	ax1.plot(hot[:,1],hot[:,0],'--',color='none', zorder=5) # dummy trace
	ax1.set_xlim(1500,4500)
	ax1.invert_xaxis()
	ax1.set_xticks([1000./0.25, 1000./.30, 1000./.40, 1000./.5, 1000./0.6])
	ax1.set_xticklabels(["0.25", "0.30", "0.40", "0.50", "0.60"], fontsize=10)
	ax1.set_xlabel("$\mathrm{1000/T (K^{-1})}$", fontsize=10)
	ax1.xaxis.tick_bottom()
	ax1.xaxis.set_label_position("bottom")
	
    # Flipping the x-axis to have warm temperatures on the left and adding labels
	ax.set_xlim(1500,4500)
	ax.invert_xaxis()
	ax.set_xlabel('$\mathrm{Temperature (K)}$', fontsize=10)
	ax.xaxis.set_ticks_position('top')
	ax.xaxis.set_label_position('top')
	ax.tick_params(axis='both', which='major', labelsize=10)
  
	
	# Secondary y-axis
	ax2 = ax.twinx()
	ax2.set_ylim(Pmin,Pmax)
	ax2.set_yticks([135.2, 130, 124, 118.4, 112.5])
	ax2.set_yticklabels(["2891", "2800", "2700", "2600", "2500"])
	ax2.set_ylabel("Depth (km)", fontsize=10)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	ax2.tick_params(axis='both', which='major', labelsize=10)
	
	#py-axis label
	ax.set_ylim(Pmin,Pmax)
	ax.set_ylabel('Pressure (GPa)', fontsize=10)
   
	#plot title
	label = "Nucleation and growth model\n" + r"$\mathrm{P_e = %d  GPa, T_e = %d  K, \alpha = %.1f  MPa/K}$" % (Pe, Te, clapeyron*1000)
	plt.text(0.05, 0.05, label, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes, fontsize=9, bbox=dict(edgecolor=(0.682, 0.682, 0.682), facecolor='white', alpha=1.0))
	
	plt.text(3000, 135.5, "Core-Mantle Boundary", horizontalalignment='center', verticalalignment='bottom', fontsize=10)
	plt.text(4300, 130, "D''\nlayer", horizontalalignment='center', verticalalignment='center', fontsize=10)

	#save plot in svg type with the name : "shearModel"
	plt.savefig('nuclGrowthModel.svg')
	plt.savefig('nuclGrowthModel.png')

	#plot
	plt.draw()
	plt.show(block=False)
	plt.pause(waittime)
	plt.close()
	
	return     



###################################################################################################################


def testmodel(datafile,factor,noneobservation):
	"""
	Fits and plots kinetics models
		
	Parameters:
	  datafile - file with experimental data
	  factor - multiplication factor for times in plotting of kinetics (in log)
	  noneobservation - File with negative observations (can be set to none)
	"""

	# Read data file
	f = open(datafile, 'r')
	datcontent = [line.strip() for line in f.readlines()]
	f.close()
	T = []
	tau = []
	P = []
	DeltaV = []  
	ratio = []
	for line in datcontent:
		if (line[0] != '#'):
			a = line.split(' ')                      #T K tau ratio P 1/T G=k*d ln(G)
			# Old version for file format from Christopher, includes redundant information
			# T.append(float(a[0]))
			# P.append(float(a[4]))
			# DeltaV.append(float(a[8]))
			# tau.append(float(a[2]))
			# New version for minimal data format with no redundant information to avoid confusion
			T.append(float(a[0]))
			P.append(float(a[1]))
			DeltaV.append(float(a[5]))
			tau.append(float(a[2]))

	# Reading geotherms
	g = open('hernlundhot.dat', 'r')
	e = open('hernlundwarm.dat', 'r')
	h = open('hernlundcold.dat', 'r')
	hernlundhot = [line.strip() for line in g.readlines()]
	hernlundwarm = [line.strip() for line in e.readlines()]
	hernlundcold = [line.strip() for line in h.readlines()]
	g.close()
	e.close()
	h.close()
	hot = []
	warm = []
	cold = []
	j = 1
	for i in range(0,len(hernlundhot)-2):
		j += 1  
		hotline = hernlundhot[j]
		hot1 = hotline.split(' ')
		warmline = hernlundwarm[j]
		warm1 = warmline.split(' ')
		coldline = hernlundcold[j]
		cold1 = coldline.split(' ')
		hot.append(map(float,hot1))
		warm.append(map(float,warm1))
		cold.append(map(float,cold1))
	hot = np.array(hot)    
	warm = np.array(warm)
	cold = np.array(cold)
	hot[:,0] = pressure(hot[:,0])    #####convert depth(km) in P
	warm[:,0] = pressure(warm[:,0])  #####convert depth(km) in P
	cold[:,0] = pressure(cold[:,0])  #####convert depth(km) in P
  
	# Fitting transformation models
	
	# Preparing array with input data
	#x[:0] is temperature
	#x[:1] is pressure
	#x[:2] is DeltaV
	x = []
	for i in range(0,len(P)):
		x.append([T[i], P[i], DeltaV[i]])

	# Starting parameters
	E0 = 400.*1000
	V0 = 15.*1000
	A0 = -15.
	paramsShear, covShear = sp.optimize.curve_fit(shearModel, x, tau, p0=(E0,A0,V0)) 
	paramsNG, covNG = sp.optimize.curve_fit(nuclGrowthModel, x, tau, p0=(E0,A0,V0))
	print(covShear)
	print(covNG)

	perrShear = np.sqrt(np.diag(covShear))   #standard deviation
	perrNG = np.sqrt(np.diag(covNG))   #standard deviation


	# Print, save, and plot results

	fi = open('summarymodel.dat', 'w')
	print>>fi, 'model', 'parameters', 'standard dev',"\n"
	print>>fi,"Shear model",paramsShear, perrShear,"\n"
	print>>fi,"Nucleation and growth model",paramsNG, perrNG,"\n"
	fi.close()

	Eshear, Ashear, Vshear = paramsShear  
  
	# V* divided by 1000 (conversion of pressures to Pa and volumes to cm3)
	# Q0 divided by 1000 to express into in kJ/mol
	print "\nResults of shear model"
	print "Q0 = %d (%d) kJ/mol" % (float(Eshear)/1000, float(perrShear[0])/1000)
	print "ln(k) = %.1f (%.1f)" % (Ashear, perrShear[1])
	print "V* = %.1f (%.1f) cm3/mol" % (float(Vshear)/1000,float(perrShear[2])/1000)
	print('\n')

	plotshearModel(Eshear,Ashear,Vshear,P,T,cold,warm,hot,factor,noneobservation) 
  
	ENG, ANG, VNG = paramsNG
	
	print "\nResults of nucleation and growth model"
	print "Q0 = %d (%d) kJ/mol" % (float(ENG)/1000, float(perrNG[0])/1000)
	print "ln(k) = %.1f (%.1f)" % (ANG, perrNG[1])
	print "V* = %.1f (%.1f) cm3/mol" % (float(VNG)/1000,float(perrNG[2])/1000)
	print ('\n')

	plotnuclGrowthModel(ENG,ANG,VNG,P,T,cold,warm,hot,factor,noneobservation)



#################################################################
#
# Main subroutines
#
#################################################################

class MyParser(argparse.ArgumentParser):
	"""
	Extend the regular argument parser to show the full help in case of error
	"""
	def error(self, message):
		
		sys.stderr.write('\nError : %s\n\n' % message)
		self.print_help()
		sys.exit(2)


def main(argv):
	"""
	Main subroutine
	"""
	
	global Pe, Te, clapeyron

	parser = MyParser(usage='%(prog)s [options] datafile', description="Fit of a shear or nucleation and growth model for the perovskite to post-perovskite transformation.\nAdapted from C. Langrand PhD thesis.")
	
	# Required arguments
	parser.add_argument('datafile',  help="Data file (required)")
	
	# Optionnal arguments
	parser.add_argument('-P', '--Pe', required=False, type=float, help="Pressure of point on pv / (pv+ppv) phase boundary in GPa. Default is %(default)s", default=128.)
	parser.add_argument('-T', '--Te', required=False, type=float, help="Temperature of point on pv / (pv+ppv) phase  in K. Default is %(default)s", default=3300.)
	parser.add_argument('-c', '--clapeyron', required=False, type=float, help="Clapeyron's slope in GPa/K. Default is %(default)s", default=6.7e-3)
	parser.add_argument('-f', '--factor', required=False, type=float, help="Multiply times by 10^x in plotting. Default is %(default)s", default=0.)
	parser.add_argument('-n', '--noneobservation', required=False, help="File with negative observations. Default is %(default)s", default=None)

	args = vars(parser.parse_args())

	datafile = args['datafile']
	Pe = args['Pe']
	Te = args['Te']
	clapeyron = args['clapeyron']
	factor = args['factor']
	noneobservation=args['noneobservation']

	testmodel(datafile, factor, noneobservation)


################################################################################################################### 

# Calling main, if necessary
if __name__ == "__main__":
    main(sys.argv[1:])
