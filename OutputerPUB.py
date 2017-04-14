# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:01:16 2015
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
   
class Outputer(object):
    def __init__(self, array1, array2, array4, statz):                        #To create this class you need to specify arrays for dataset, found steps, an array that can be anything you want (true steps, means) and statistics.  
        self.dataset = np.copy(array1)                                        #Upon creation the class is initialized with the dwell means according to found steps
        self.foundsteps = np.copy(array2)
        self.fstepsmeans = self.findMeans(array1, array2)
        self.skinds = np.copy(array4)
    def exportData(self, a1, a2, a4, stz, KK):                                #Outputs files of the data, the found a step dwells means, plus the expected fluorescence given the found active fluorophore numbers and statistics 
        dexf1 = open('%s_data.txt' % KK,'w')
        dexf2 = open('%s_found_steps.txt' % KK,'w')
        dexf4 = open('%s_ active_fluors.txt' % KK,'w')
        points1 = len(a1)
        points4 = len(a4)
        for i in range(0,points1,1):
            gi1 = str(a1[i])
            dexf1.write(gi1)
            dexf1.write('\n')    
            gi2 = str(a2[i])
            dexf2.write(gi2)
            dexf2.write('\n')
        for i in range(0,points4,1):
            gi4 = str(a4[i])                                                  #str(a4[i]*stz[2]+stz[0]) 
            dexf4.write(gi4)
            dexf4.write('\n')
        maxfli = np.max(a4)
        nard = np.diff(a4)
        maxste = len(nard[nard!=0])
        print 'total number of fluorescence level changes: ' ,maxste
        print 'maximum fluorophores active: ' ,maxfli
        dexf1.close()
        dexf2.close()
        dexf4.close()
    def findMeans(self, dataz, stepz):                                        #Finds dwells means from the data given an array of equal length noting steps indices as ones and zeros
        pointz = len(dataz)
        galway = np.zeros(pointz)
        mark = 0
        for i in range(1,pointz,1):
            if stepz[i]!=0 and i!=pointz-1:
                for j in range(mark,i,1):
                    galway[j] = np.mean(dataz[mark:i],dtype=np.float64)
                mark = i
            if i==pointz-1 and stepz[i]!=0:  
                for j in range(mark,i,1):
                    galway[j] = np.mean(dataz[mark:i],dtype=np.float64)   
                galway[i] = dataz[i]
            if i==pointz-1 and stepz[i]==0:
                for j in range(mark,i+1,1):
                    galway[j] = np.mean(dataz[mark:i+1],dtype=np.float64)
        return galway
    def findMeans2(self, datax, stepx):                                       #Finds dwells means from the data given an array of step time indices    
        galwa = np.zeros_like(datax)
        marc = 0
        for i in range(0,len(stepx),1):
            for j in range(marc,int(stepx[i]+1),1):
                galwa[j] = np.mean(datax[marc:int(stepx[i]+1)],dtype=np.float64)
            marc = int(stepx[i]+1)
        for i in range(marc,len(datax),1):
            galwa[i] = np.mean(datax[marc:],dtype=np.float64)
        return galwa        

        
        
        
        
        
