# TheisSolutionStepPumpTest.py

# This code similates pumping induced drawdown for designing a pump test

# Scripted by Casey Gierke of Lee Wilson & Associates on 7/2/2013.
# Updated 10/9/2013
# Updated 6/13/2018

# With Notepad++, use F5 then copy this into box
# C:\Python27\python.exe -i "$(FULL_CURRENT_PATH)"

# ------------------------------------------------------
# IMPORTS
# ------------------------------------------------------

import os
import math
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import pylab as p
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# ------------------------------------------------------
# DEFINE FUNCTIONS
# ------------------------------------------------------

# Define last position finder
def find_last(s,t):
	last_pos = -1
	while True:
		pos = s.find(os.sep, last_pos +1)
		if pos == -1:
			return last_pos
		last_pos = pos

# Create finder function
def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

# Define Theis function 
def Theis(Q, S, T, r, t):
	
	# Create a file to write to
	fout = open(path+os.sep+'Outfile- '+str(Q)+' gpm.txt','w')
	# foutU = open(path+os.sep+'OutfileU- '+str(Q)+' gpm.txt','w')
	# foutWu = open(path+os.sep+'OutfileWu- '+str(Q)+' gpm.txt','w')
	# foutS = open(path+os.sep+'OutfileSeq- '+str(Q)+' gpm.txt','w')

	# Convert pumpint to ft^3/day
	Q = Q*192.5
	
	# Build matrices populated with zeros
	u = [[0 for _ in range(len(r))] for _ in range(len(t))]
	Wu = [[0 for _ in range(len(r))] for _ in range(len(t))]
	s = [[0 for _ in range(len(r))] for _ in range(len(t))]
	# Seq = [[0 for _ in range(len(r))] for _ in range(len(t))]

	# Open loop to populate matrices
	for i in range(0, len(t)):      # for time steps
		for j in range(0, len(r)):  # for length of domain
			u[i][j] = pow(r[j],2)*S/(4*T*t[i])

			# Initialize series
			seq_plus = 0
			seq_minus = 0
			# Start loop to iterate W(u) series expansion
			for n in range(2,5,2):
				# Build series up to 11th power
				seq_plus = seq_plus + u[i][j]**(n+1)/((n+1)*math.factorial(n+1))
				seq_minus = seq_minus - u[i][j]**n/(n*math.factorial(n))
				n = n + 1

			# Build well function    
			
			Wu[i][j] = -0.577216-np.log(u[i][j])+u[i][j] + seq_plus+seq_minus
			# Seq[i][j] = seq_plus+seq_minus
			# Calculate drawdowns
			s[i][j] = Q*Wu[i][j]/(4*math.pi*T) + So
			j = j + 1                      
		i = i + 1

	# Remove errors due to model limitations (ie, drawdowns recurring at large
	# distances and (-) drawdowns)
	for i in range(0, len(t)):      # for time steps
		# Write the distance from the well at the top of the file
		if i == 0:
			for entry in r:
				fout.write('\t'+str(entry))
				# foutU.write('\t'+str(entry))
				# foutWu.write('\t'+str(entry))
				# foutS.write('\t'+str(entry))
			fout.write('\n')
			# foutU.write('\n')
			# foutWu.write('\n')
			# foutS.write('\n')
		
		# Write the time at beginning of line	
		fout.write(str(ttime[i])+'\t'+str(s[i][0])+'\t')
		# foutU.write(str(ttime[i])+'\t'+str(u[i][0])+'\t')
		# foutWu.write(str(ttime[i])+'\t'+str(Wu[i][0])+'\t')
		# foutS.write(str(ttime[i])+'\t'+str(Wu[i][0])+'\t')
		for j in range(1, len(r)):  # for length of domain
			# Remove data if it exceeds boundary
			if s[i][j] > s[i][j-1]:
				# s[i][j] = None
				s[i][j] = 0
				# print('s[i][j]= ',i,j,s[i][j])
			# Point it out if it doesn't make sense
			elif s[i][j] < 0:
				# print('s[i][j]= ',i,j,s[i][j])
				# print('s[i][j-1]= ',i,j-1,s[i][j-1])
				s[i][j] = 0
			fout.write(str(s[i][j])+'\t')
			# foutU.write(str(u[i][j])+'\t')
			# foutWu.write(str(Wu[i][j])+'\t')
			# foutS.write(str(Wu[i][j])+'\t')
		fout.write('\n')
		# foutU.write('\n')
		# foutWu.write('\n')
		# foutS.write('\n')
	fout.close()
	# foutU.close()
	# foutWu.close()
	# foutS.close()
	# Convert to a numpy array so time can be plotted
	s = np.array(s)
	return s
	
# ------------------------------------------------------
# INPUTS
# ------------------------------------------------------

# Define path
path = os.path.abspath(os.path.dirname(__file__))
# # Shorten path to one folder up
# path = path[:find_last(path,os.sep)]

# Read in well info
WellsIn = open(path+os.sep+'Wells.txt','r')
# Create Wells dictionary
Wells = {}
for line in WellsIn:
	x1 = find_nth(line, '\t', 1)
	x2 = find_nth(line, '\t', 2)
	Wells[line[:x1]] = [line[x1+1:x2], line[x2+1:-1]]
	
# Find distance to FCI observation wells
WellQ = Wells['R-01']
WellObs = Wells['I-01']
Robs = math.sqrt((float(WellQ[0])-float(WellObs[0]))**2 + (float(WellQ[1])-float(WellObs[1]))**2)

# Parameters to adjust
T = 622             # Transmissivity            [ft^2/day] 1000 for difference from average year
# T = K*b             # Transmissivity            [ft^2/day] 1000 for difference from average year
# S = Ss*b             # Storage coefficient       [-]
S = 0.0033            # Storage coefficient       [-]
So = 0             # Observed drawdown         [ft]

# Define all the discharges for steps
Q1 = 5					# Pumping rate              [gpm]
Q2 = 10					# Pumping rate              [gpm]
Q3 = 20					# Pumping rate              [gpm]
Q4 = 40					# Pumping rate              [gpm]
Q5 = 40				# Pumping rate              [gpm]

# Define array for radius  
rin = 5000       # Radius                    [ft]
r = int(math.ceil(10*math.log(rin,10)))
# r= 36
r = range((10*r)+1)        # All of this takes the input and
r = np.array(r)     # makes a log scale array up to the domain value
r = r/100.0
r = 10**r

# Define time array
ttime = 36 				# hours of testing
ttime = range(1,10*(ttime+1))
ttime = np.array(ttime)/10.0
t = ttime/float(24) # hours converted to days for calculations

# NTime = len(t)-1

# ------------------------------------------------------
## COMPUTATIONS
# ------------------------------------------------------

# Calculate drawdowns for each individual pumping step
s1 = Theis(Q1, S, T, r, t)
s2 = Theis(Q2-Q1, S, T, r, t)
s3 = Theis(Q3-Q2, S, T, r, t)
s4 = Theis(Q4-Q3, S, T, r, t)
s5 = Theis(Q5, S, T, r, t)

# Combine individual steps for total effect
# Create a space for the results
s = [[0 for _ in range(len(r))] for _ in range(len(t))]

# Open a loop to move through the matrix
for i in range(0, len(t)):      # for time steps
	for j in range(0, len(r)):  # for length of domain
		if i < 5:
			s[i][j] = s1[i][j]
		elif i >= 5 and i < 10:
			s[i][j] = s1[i][j] + s2[i-5][j]
		elif i >= 10 and i < 15:
			s[i][j] = s1[i][j] + s2[i-5][j] + s3[i-10][j]
		elif i >= 15 and i < 95:
			s[i][j] = s1[i][j] + s2[i-5][j] + s3[i-10][j] + s4[i-15][j]
		elif i >= 95:
			s[i][j] = s1[i][j] + s2[i-5][j] + s3[i-10][j] + s4[i-15][j] - s5[i-95][j]

# Remove errors due to model limitations (ie, drawdowns recurring at large
# distances and (-) drawdowns)
for i in range(0, len(t)):      # for time steps
	for j in range(1, len(r)):  # for length of domain
		# Remove data if it exceeds boundary
		if s[i][j] > s[i][j-1]:
			s[i][j] = None
		# Point it out if it doesn't make sense
		elif s[i][j] < 0:
			s[i][j] = 0
# Convert to a numpy array so time can be plotted
s = np.array(s)
				
# ------------------------------------------------------
## OUTPUTS
# ------------------------------------------------------

# For time in hours
times = [
			str((ttime[4]))+' hours', 
			str((ttime[9]))+' hours',
			str((ttime[14]))+' hours', 
			str((ttime[94]))+' hours',
			str((ttime[95]))+' hours', 
			str((ttime[109]))+' hours', 
			str((ttime[194]))+' hours'
			]

# Plot output
fig = plt.figure()

# Drawdown vs. distance plot
ax1 = fig.add_subplot(211)
ax1.plot(r,s[4],'.k', linewidth = .75)
ax1.plot(r,s[9],'-.k', linewidth = .75)
ax1.plot(r,s[14],'o-k', linewidth = .75, markersize=3, markeredgecolor='c', markeredgewidth='0.5', fillstyle='none')
ax1.plot(r,s[94],'-k', linewidth = .75)
ax1.plot(r,s[95],'.b', linewidth = .75)
ax1.plot(r,s[109],'-.b', linewidth = .75)
ax1.plot(r,s[194],'-.r', linewidth = .75)


# ax1.plot(r,s[4],'.k', linewidth = .75, markersize=3, markeredgecolor='b', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[9],'-.k', linewidth = .75, markersize=3, markeredgecolor='k', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[14],'-.-k', linewidth = .75, markersize=3, markeredgecolor='c', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[94],'-k', linewidth = .75, markersize=3, markeredgecolor='g', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[95],'.b', linewidth = .75, markersize=3, markeredgecolor='violet', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[109],'-.b', linewidth = .75, markersize=3, markeredgecolor='k', markeredgewidth='0.5', fillstyle='none')
# ax1.plot(r,s[194],'-.r', linewidth = .75, markersize=3, markeredgecolor='y', markeredgewidth='0.5', fillstyle='none')

# Now lets do stuff to the plot
plt.gca().invert_yaxis()    # Invert the y axis
plt.xlabel('Distance from pumping well [ft]', fontsize=12)
plt.ylabel('Drawdown [ft]', fontsize=12)
plt.legend(times, loc=0) # 4 is LR, 0 is best
# ax1.set_xscale('log')

# # Axis labels
# majorLocator   = MultipleLocator(r[len(r)-1]/4)
# majorFormatter = FormatStrFormatter('%d')
# minorLocator   = MultipleLocator(r[len(r)-1]/40)

# ax1.xaxis.set_major_locator(majorLocator)
# ax1.xaxis.set_major_formatter(majorFormatter)

# Drawdown vs. time plot
ax2 = fig.add_subplot(212)
# ax2.plot(ttime,s1[:,19],'o-r', linewidth = .75, markersize=3, markeredgecolor='r', markeredgewidth='0.5', fillstyle='none')
# ax2.plot(ttime,s5[:,19],'o-c', linewidth = .75, markersize=3, markeredgecolor='c', markeredgewidth='0.5', fillstyle='none')
ax2.plot(ttime,s[:,190],'o-k', linewidth = .75, markersize=3, markeredgecolor='b', markeredgewidth='0.5', fillstyle='none')

# Now lets do stuff to the plot

plt.gca().invert_yaxis()    # Invert the y axis
title = str('Drawdown: S = ')+str(S)+str(', T = ')+str(T)
fig.suptitle(title, fontsize=14) 
plt.xlabel('Time [hours]', fontsize=12)
plt.ylabel('Drawdown [ft]', fontsize=12)
Legend2 = ['Drawdown at '+str(round(r[19],1))+' ft',]
plt.legend(Legend2, loc=0) # 4 is LR, 0 is best

# # Axis labels
# majorLocator   = plt.MultipleLocator(5)
# majorFormatter = FormatStrFormatter('%d')
# #minorLocator   = MultipleLocator(rin/40)

# ax2.xaxis.set_major_locator(majorLocator)
# ax2.xaxis.set_major_formatter(majorFormatter)

#for the minor ticks, use no labels; default NullFormatter
#ax.xaxis.set_minor_locator(minorLocator)

'''
For analyzing effects on nearby wells

# Well screen info
top1 = 43
top2 = 58
top3 = 78

# Fill area to represent well screens
p.fill([0,r[len(r)-1],r[len(r)-1],0], [top1+10,top1+10,top1,top1],
       'b', alpha = 0.2, edgecolor='k')
p.fill([0,r[len(r)-1],r[len(r)-1],0], [top2+10,top2+10,top2,top2],
       'b', alpha = 0.2, edgecolor='k')
p.fill([0,r[len(r)-1],r[len(r)-1],0], [top3+10,top3+10,top3,top3],
       'b', alpha = 0.2, edgecolor='k')
# Project non interference drawdown
p.fill([0,r[len(r)-1],r[len(r)-1],0], [52,52,51.75,51.75],
       'k', alpha = 1)
'''
'''
# A second plot to debug

ax = fig.add_subplot(312)
ax.plot(r,u[0],'o-')
ax.plot(r,u[1], 'o-')
ax.plot(r,u[2], 'o-')
ax.plot(r,u[3], 'o-')
plt.legend(t, fontsize=12, loc=0)

ax = fig.add_subplot(313)
ax.plot(r,Wu[0],'o-')
ax.plot(r,Wu[1], 'o-')
ax.plot(r,Wu[2], 'o-')
ax.plot(r,Wu[3], 'o-')
plt.legend(t, fontsize=12, loc=0)
'''

plt.show()




