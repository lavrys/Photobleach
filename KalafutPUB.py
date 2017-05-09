# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 08:47:29 2016
@author: Sina Jazani
@author: Konstantinos Tsekouras
Copyright (C) 2015 Sina Jazani, Konstantinos Tsekouras 
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
import math

class KalafutC(object):
    def __init__(self, array1):
        signal0 = array1 
        self.tzero = self.KalafutD(signal0)[1]
        self.stats = self.findstats(signal0, self.tzero)
    def uniq(self, input):
        output = []
        for x in input:
            if x not in output:
                output.append(x)
        return output 
    def KalafutD(self, signal):     
        threshold = 0.9                                                       #MANYALY SET Define a threshold to find the new steps
        jf=[0,len(signal)-1]
        jrhtest=10                                                            #Define a counter to be as a costraint in the while loop
        jignore=[]           
        while (jrhtest>0):    
            minSIC=pow(10,6)                                                  #Define a minimum SIC that is too huge and all the time compare and replace SIC with the lower one
            jpos=0
            for j2 in range(jf[0]+1,jf[len(jf)-1]-1):                         #For loop to find the new step position                         
                jtest=[None]*(len(jf))
                if jf.count(j2) != 1 and jignore.count(j2) != 1:                
                    jtest=jf[:]
                    jtest.extend([j2])
                    jtest.sort()
                    jtest=self.uniq(jtest)          
                    meannn=[None]*(len(jtest)-1)
                    tot=jtest[len(jtest)-1]
                    sigs=0
                    for hate in range(0,len(jtest)-1):
                        meannn[hate]=sum(signal[jtest[hate]:jtest[hate+1]])/(jtest[hate+1]-jtest[hate])
                        sigs+=(1./tot)*math.fsum(pow(signal[jtest[hate]:jtest[hate+1]]-meannn[hate],2))
                    if sigs==0:
                        sigs=1;
                    SIC=(len(jtest)*(math.log(tot)))+(tot*(math.log(sigs)))         
                    if SIC<=minSIC:      
                       minSIC=SIC
                       jpos=j2   
            meanl=0
            meanr=0
            if jf.count(jpos)!=1:                                             #Find the Mean of right and left side of the new step position           
                for s in range(0,len(jf)-1):
                    if jf[s]<=jpos and jf[s+1]>jpos and jpos!=0 and abs(jf[s]-jpos)>1 and abs(jf[s+1]-jpos)>1:
                        meanl=math.fsum(signal[jf[s]:jpos])/(jpos-jf[s])
                        meanr=math.fsum(signal[jpos:jf[s+1]])/(jf[s+1]-jpos)  
                        break           
            if abs(meanl-meanr)>threshold:                                    #Condition to accept the new step position
                jf.extend([jpos])
                jrhtest+=2
            else:
                if jignore.count(jpos) != 1:                                  #Put the rejected step position in the rejection array
                    jignore.extend([jpos])
                    jignore=self.uniq(jignore)
                    jrhtest-=1           
            jrhtest-=1       
            jf.sort()                                                         #Sort and "uniq" the step positions
            jf=self.uniq(jf)
            if len(jf)+len(jignore)>=len(signal):                             #A condition to break the while loop when the sum of the sizes of accepted and rejected point arrays is higher or equal that the signal size 
                break
            finalsignal=[[None]*(len(jf)-1),[None]*(len(jf)-1)]               #Create a final matrix to plot the results
            finalsignal[0][:]=range(len(signal)-1)     
            for mea in range(len(jf)-1):
                meann=(math.fsum(signal[jf[mea]:jf[mea+1]])/(jf[mea+1]-jf[mea]))
                finalsignal[1][jf[mea]:jf[mea+1]] = [meann]*(jf[mea+1]-jf[mea])
        return jf
    def findstats(self, data2, loc):                                          #loc is the location of the last step as found by Kalafut
        MIN = 1e-10 
        stat_window = 25                                                      #MANUALLY PICK a window size to calculate the mean and variance of the single fluorophore                                   
        stat_val = np.zeros(4)
        point = loc
        if point+stat_window<len(data2):
            limi = point + stat_window
        else:
            limi = 3
            print 'Insufficient data to acquire fluorophore statistics; results unreliable'   
        signal_sampleB = np.copy(data2[:point])                               #Use copy just for security
        signal_sampleF = np.copy(data2[point:limi])
        stat_val[0] = np.mean(signal_sampleB,dtype=np.float64)
        stat_val[1] = np.var(signal_sampleB,dtype=np.float64)
        stat_val[2] = np.mean(signal_sampleF,dtype=np.float64) - stat_val[0]
        stat_val[3] = np.var(signal_sampleF,dtype=np.float64) - stat_val[1]
        for i in range(0,4,1):
            if stat_val[i]<=MIN:
                stat_val[i] = MIN
                print i,'found value less than 1e-10, set to 1e-10: curate your dataset' 
        return stat_val      
      


    
    
    