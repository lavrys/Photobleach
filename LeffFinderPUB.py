# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:20:20 2015
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
import math


#This might just as well have been a function. 
#It finds an estimate of λbar which we need to calculate an array of local λeffectives in the next class

class LbarFind(object):                                                       #Gives a crude estimate of Lbar
    def __init__(self, dataP, arraystats, tzeroM, fluoMIN):                   #Give the data set, the statistics, the location of the bleach-to-background step and how many fluorophores are still lit at data set start (remember the data set is inversed)
        increm = 25                                                           #MANUALLY SET how many points to consider from the end of the inverted data set where fluorophore numbers are max
        burzum = dataP[len(dataP)-increm:]
        miuB = arraystats[0]
        miuF = arraystats[2]
        nee = np.ceil((np.mean(burzum, dtype=np.float64) - miuB)/miuF)        #Rough estimate for how many fluorophores you have in the "increm"-size end of the data set
        if fluoMIN==0:
            fluoMIN += 1                                                      #Make this into 1 just to avoid math issues
        if nee<=fluoMIN and nee>0:                                            #If your max fluorophore estimate is less than your starting one [SHOULD NOT HAPPEN!]
            self.Lbar = 10.0/(nee*(len(dataP)-tzeroM))                         
        elif nee<=fluoMIN and nee<=0:                                         #If your max fluorophore estimate is negative [SHOULD DEFINITELY NOT HAPPEN!]
            self.Lbar = 10.0/(fluoMIN*(len(dataP)-tzeroM))                    
        else:                                                                 #Normal case
            lin = np.arange(float(fluoMIN),float(nee))
            lin = lin**(-1.0)
            self.Lbar = 10.0/(np.sum(lin)*(len(dataP)-tzeroM))

#Here we calculate an array of local crude estimates of λeffective, the Poisson effective photobleach rate that is time-dependent 
class PriorSlicer(object):
    def __init__(self, dataP, lbarP, tzeroM, fluoMIN):                        #Give the data set, the lbar, the location of the bleach-to-background step and how many fluorophores are still lit at data set start (remember the data set is inversed)
        index_length = len(dataP) - tzeroM
        indexes = np.arange(index_length)
        j = fluoMIN + 1
        j_sum = 0
        while j_sum<index_length:
            index = math.ceil(1.0/(j*lbarP)) + j_sum
            if j==fluoMIN+1 and index>index_length:
                print "wrong lbar"
                raise ValueError
            j_sum = index
            j += 1
            for i in range(0,len(indexes),1):
                if indexes[i]==j_sum:
                    indexes[i]=0
        ptuple = np.where(indexes==0)
        self.pslice = np.array(ptuple[0])
        self.pslice += tzeroM
        k = 0
        self.leffarray = np.zeros_like(dataP)
        for i in range(0,len(dataP),1):
            if i in self.pslice:
                k += 1
            self.leffarray[i] = k*lbarP
                
