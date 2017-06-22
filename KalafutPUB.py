# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 08:47:29 2016
@author: Sina Jazani
@author: Konstantinos Tsekouras
@author: Aditya Prakash
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
        #return list(set(input))
        output = []
        for x in input:
            if x not in output:
                output.append(x)
        return output 
    def KalafutD(self, signal):
        threshold = 0.9   
        tot = len(signal) - 1                                                    
        jf=[0,tot]
        jrhtest=1                                                             
        jignore=[]           
        while (jrhtest>0):    
            minSIC=pow(10,6)                                                  
            jpos=0
            SIC = []
            for j2 in range(1,len(signal)-2):
                if j2 not in jf and j2 not in jignore:               
                    jtest = sorted(list(set(jf + [j2])))
                    meannn=[None]*(len(jtest)-1)
                    meannn = [sum(signal[jtest[hate]:jtest[hate+1]])/(jtest[hate+1]-jtest[hate]) for hate in range(0,len(jtest)-1)]
                    sigs = max(0, np.sum([math.fsum(pow(signal[jtest[hate]:jtest[hate+1]]-meannn[hate],2))/tot for hate in range(0,len(jtest)-1)]))
                    if sigs==0:
                        sigs=1
                    SIC.append((len(jtest)*(math.log(tot)))+(tot*(math.log(sigs))))
            jpos = np.argmin(SIC) + 1
            minSIC = min(minSIC, min(SIC))   
            meanl=0
            meanr=0
            if jpos not in jf:
                for s in range(0,len(jf)-1):
                    if jf[s]<=jpos and jf[s+1]>jpos and jpos!=0 and abs(jf[s]-jpos)>1 and abs(jf[s+1]-jpos)>1:
                        meanl=math.fsum(signal[jf[s]:jpos])/(jpos-jf[s])
                        meanr=math.fsum(signal[jpos:jf[s+1]])/(jf[s+1]-jpos)  
                        break           
            if abs(meanl-meanr)>threshold:                                    
                jf.extend([jpos])
                jrhtest+=2
            else:
                if jpos not in jignore:
                    jignore.extend([jpos])
                    jignore = list(set(jignore))
                    jrhtest+=2           
            jrhtest-=1       
            jf = sorted(list(set(jf)))
            if len(jf)+len(jignore)>=len(signal):                             
                break
            finalsignal = []
            finalsignal.append(list(range(len(signal)-1)))
            finalsignal.append([math.fsum(signal[jf[mea]:jf[mea+1]])/(jf[mea+1]-jf[mea]) for mea in range(len(jf)-1) for _ in range(jf[mea+1]-jf[mea])])
        return jf

    def findstats(self, data2, loc):                                         
        MIN = 1e-10 
        stat_window = 25                                                      #MANUALLY PICK a window size to calculate the mean and variance of the single fluorophore                                   
        stat_val = np.zeros(4)
        point = loc
        if point+stat_window<len(data2):
            limi = point + stat_window
        else:
            limi = 3
            #print 'Insufficient data to acquire fluorophore statistics; results unreliable'
            print('Insufficient data to acquire fluorophore statistics; results unreliable')
        signal_sampleB = np.copy(data2[:point])                               #Use copy just for security
        signal_sampleF = np.copy(data2[point:limi])
        stat_val[0] = np.mean(signal_sampleB,dtype=np.float64)
        stat_val[1] = np.var(signal_sampleB,dtype=np.float64)
        stat_val[2] = np.mean(signal_sampleF,dtype=np.float64) - stat_val[0]
        stat_val[3] = np.var(signal_sampleF,dtype=np.float64) - stat_val[1]
        for i in range(0,4,1):
            if stat_val[i]<=MIN:
                stat_val[i] = MIN
                #print i,'found value less than 1e-10, set to 1e-10: curate your dataset'
                print('{0} found value less than 1e-10, set to 1e-10: curate your dataset'.format(i))
        return stat_val      
      


    
    
    
