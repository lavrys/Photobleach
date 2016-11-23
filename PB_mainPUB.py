# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:56:41 2015
@author: Konstantinos Tsekouras
Copyright (C) 2015 Konstantinos Tsekouras 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import numpy as np

from KalafutPUB import KalafutC 
from LeffFinderPUB import LbarFind, PriorSlicer
from SeekerPUB import Slicer, mSICer
from OutputerPUB import Outputer

def exportData2(a1, a2): 
    poin  = len(a1)    
    dexf11 = open('DataY','w')
    dexf22 = open('Step2levels','w')
    for i in range(0,poin,1):
        gi1 = str(a1[i])
        dexf11.write(gi1)
        dexf11.write('\n')    
        gi2 = str(a2[i])
        dexf22.write(gi2)
        dexf22.write('\n')            
    dexf11.close()
    dexf22.close()        

#Preliminaries and choices
stat_flag = 3                                                                 #MANUALLY PICK How to get statistics depending on what you know
dropORpad = 0                                                                 #MANUALLY PICK If the size of your data set is not integer * your window size, set his to 0 to drop the extra data points at the end or to 1 to pad them                  
interruptflag = 0                                                             #MANUALLY PICK Set to 0 if you want the code to go through without interruption, or to 1 if you want to stop and check the result of your window choice in Step 2 
fluorupdateflag = 3                                                           #MANUALLY PICK Determines whether you should update your fluorophore numbers from Step 2 (2) or Step 3 (3)
priorsliceflag = 0                                                            #MANUALLY PICK Set to 0 if you want to have an estimate of λeffective from data or 1 if you have some idea of lbar and the last step location               
KVflag = 0                                                                    #MANUALLY PICK Set something other than 0 to give a slice of your data set around the last(in time) step to Kalafut 
zignal = np.loadtxt("TestDataSet", unpack=True)                            #Import your data set   
statzignal = np.array(zignal)
points = len(zignal)                                                          #Number of data points
windz = 50                                                                   #MANUALLY SET data set window size  
if points%windz==0:
    many_d = points//windz                                                    #Number of windows in data set 
else:
    many_d = points//windz + 1             
    if dropORpad==0:  
        zignal = np.array(zignal[:points-points%windz])    
    elif dropORpad==1:
        pad = np.array(zignal[points-windz+points%windz:])
        zignal = np.concatenate((zignal,pad))
print '# of windows:', many_d  
tignal = np.array(zignal) 
#Get or estimate statistics for background and single fluorophore
xignal = np.copy(tignal)
if KVflag==0:
    statsignal = statzignal[::-1]                                             #Invert the data set time-wise to streamline statistics_estimator    
else:
    byhand1 = 9800
    byhand2 = 250
    if byhand1+byhand2>=points:
        statsignal = np.copy(statzignal[byhand1-byhand2:])
        statsignal = statsignal[::-1]
    else:
        statsignal = np.copy(statzignal[byhand1-byhand2:byhand1+byhand2])
        statsignal = statsignal[::-1]          
if stat_flag==0:                                                              #Set stat_flag to 0 if you already know the mean and variance of background and single fluorophore
    mB, vB, mF, vF = 20.0, 0.001, 2.0, 0.2                                                                
    Try0 = KalafutC(statsignal)
    tzerom = Try0.tzero
elif stat_flag==1:                                                            #Set stat_flag to 1 if you already know the mean and variance of the background but not the single fluorophore
    mB, vB = 0.0, 10.0                                                        
    Try1 = KalafutC(statsignal)                                               #Call Kalafut C from Kalafooter to run the KV algorithm on the data set to get estimate of statistics you do not know
    Try1a = Try1.stats
    tzerom = Try1.tzero
    mF, vF = Try1a[2],Try1a[3]
elif stat_flag==2:                                                            #Set stat_flag to 2 if you already know the mean and variance of the single fluorophore but not the background
    mF, vF = 10.0, 10.0                                                       
    Try2 = KalafutC(statsignal)
    Try2a = Try2.stats
    tzerom = Try2.tzero    
    mB, vB = Try2a[0],Try2a[1]  
elif stat_flag==3:                                                            #Set stat_flag to 3 if you do not know the mean and variance of the background or the single fluorophore
    Try3 = KalafutC(statsignal)
    Try3a = Try3.stats  
    tzerom = Try3.tzero    
    mB, vB, mF, vF = Try3a[0], Try3a[1], Try3a[2], Try3a[3]
print 'KV completed',mB,vB,mF,vF,tzerom
signal = np.copy(xignal)                                                      #Make anew copy of the data set just for security   
signal = signal[::-1]                                                         #Invert the signal so that instead of fluorophores bleaching, you have fluorophores been "lit" 
elax = np.min(signal)
if elax <= 0.0: 
    lift = np.abs(elax) + 1.0                                                 #Remove negative values in the data set, adjust μB accordingly
    signal += lift
    mB += lift   
else:
    lift = 0.0                                                                
arraystats = np.zeros(4)                                                      #Array of statistics found
arraystats[0] = mB
arraystats[1] = vB
arraystats[2] = mF
arraystats[3] = vF    
fluor = 0
step_locations = np.zeros((many_d,windz))                                     #Precise step locations
stepsC = np.zeros(1)
slicelook = Slicer(signal,windz,0,arraystats)                                 #Call the Slicer to perform Step 2
if interruptflag == 1:
    exportData2(signal,slicelook.levelz)  
    raise ValueError
else:
    pass
fluorINI = 0                                                                  #Sets initial number of fluorophores in data set for calculating λeffective. Right now always set to 0, will add features later
if priorsliceflag==0:                                                         #If you want to find lbar, tzerom from the data                                                                    
    calba = LbarFind(signal, arraystats, tzerom, fluorINI)                                   
    lbar = calba.Lbar                                                
    zamm = PriorSlicer(signal, lbar, tzerom, fluorINI)   
    lefaray = zamm.leffarray                                                  
elif priorsliceflag==1:                                                       #Provide lbar, location of last step if you have an estimate somehow
    lbar = 0.0005           
    zamm = PriorSlicer(signal, lbar, tzerom, fluorINI)    
    lefaray = zamm.leffarray                                                  
for i in range(0,many_d,1):
    mSIClook = mSICer(signal, i, fluor, windz, arraystats, points, lefaray)   #Use the criterion                                                                    
    steps_found = mSIClook.SIClocs                                            #mSIC finding of steps in slice as numbers in array of window size
    stepsCi = np.zeros_like(steps_found)
    levelz = fluor
    for j in range(0,len(steps_found),1):                                     #Find active fluorophores at any moment  
        stepsCi[j] = levelz
        levelz += steps_found[j]
    if fluorupdateflag == 3:
        fluor = mSIClook.fluorOUT
    elif fluorupdateflag == 2 and i>0:
        fluor = slicelook.fosfor[i-1]
    stepsC = np.concatenate((stepsC,stepsCi))                
    for j in range (0,windz,1):
        step_locations[i][j] = steps_found[j]   
    print 'Step 3', i                                                         #prints the window just processed to give a measure of how long the code still has to go                          
stepsF = step_locations.flatten()                                             #Re-reverse arrays for output 
#stepsF = stepsF[::-1] 
stepsC = stepsC[1:]
#stepsC = stepsC[::-1]   
Out1 = Outputer(signal, stepsF, stepsC, arraystats)
Outputer.exportData(Out1, Out1.dataset, Out1.fstepsmeans, Out1.skinds, arraystats, 'OOO')
print ('\a')
input("\n\nPress Enter to exit.")













