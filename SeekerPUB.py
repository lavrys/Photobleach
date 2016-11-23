# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:24:06 2015
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
import itertools as it
import pandas as pd
from scipy import misc as mi
  
class mSICer(object):
    def __init__(self, dataa, interval, fluorR, winn, arraystats, NN, lefaray): 
        self.phi0 = 1.0
        self.phi1 = 1.0
        self.bigN = NN
        gamma = 0.5                                                           #MANUALLY SET a value for the Gamma0 cutoff.                                                   
        prior = 1                                                             #Set to 0 if you do not want to use a prior (only for testing) and 1 otherwise
        self.fluorOUT = 0
        self.leffarray = lefaray
        exte = float(len(dataa))
        laff = np.sum(self.leffarray)                                         #Just summing the array of λeffectives, as a crude way to estimate a general λeffective. This can be substantially improved, but is accurate enough for now.
        self.mB = arraystats[0]
        self.vB = arraystats[1]
        self.mF = arraystats[2]
        self.vF = arraystats[3]
        fluorL = fluorR    
        self.look_slice = np.copy(dataa[interval*winn:(interval+1)*winn])     #the interval in which we are going to look for steps, use copy for security                                               
        best_locs = np.zeros(winn)                                            #Temp store for the best step locations 
        best_locs2 = np.zeros(winn)
        best_mSIC2 = winn*np.log(fluorL*self.vF+self.vB)                      #Temp store for the best mSIC, start with value assuming there are NO events     
        for i in range(0,len(self.look_slice),1):
            best_mSIC2 += (self.look_slice[i] - fluorL*self.mF - self.mB)**2.0/(fluorL*self.vF+self.vB) #+ 2.0*np.sum(self.leffarray)
        all_cases = self.multichoose(winn, 2)                                 #find all possible cases with up to 2 change-points     
        look_multi = np.array(all_cases)                           
        for qq in range(0,len(look_multi),1):                      
            for ww in range(0,winn,1):
                if look_multi[qq][ww]==2:                       
                    look_multi[qq][ww] = 1                                      
        for j in range(0,len(look_multi),1):  
            roawx = np.copy(look_multi[j])                                    #Go over each possible permutations of step locations and list all possible overlap levels 
            cham = len(roawx[roawx!=0])
            dalist = self.dalistor(roawx, cham, winn)                  
            for kk in range(0, len(dalist),1):
                roaw = dalist[kk]
                if len(roaw[roaw!=0])==0:
                    pass
                else:
                    mSICs_case = 0                                            #mSIC for the particular step case 
                    fluorLin = fluorL
                    if self.negativeCheck(roaw,fluorLin)==0:                  #Check whether you ever fall below zero fluorophores      
                        intv = 0                                              #Index for slicing data 
                        for k in range(0,winn,1):
                            if roaw[k]!=0 and k<winn-1:
                                partial1 = np.copy(self.look_slice[intv:k+1]) #Use copy for security, calculate criterion value
                                n_p = len(partial1)                                       
                                intv = k+1
                                mSICs_case += n_p*np.log(fluorLin*self.vF+self.vB)
                                for q in range(0,len(partial1),1):
                                    mSICs_case += (partial1[q]  - fluorLin*self.mF - self.mB)**2.0/(fluorLin*self.vF+self.vB)
                                fluorLin += roaw[k] 
                        partial2 = np.copy(self.look_slice[intv:winn])
                        if len(partial2)!=0:
                            n_p = len(partial2)
                            mSICs_case += n_p*np.log(fluorLin*self.vF+self.vB)
                            for q in range(0,len(partial2),1):
                                mSICs_case +=  (partial2[q]  - fluorLin*self.mF - self.mB)**2.0/(fluorLin*self.vF+self.vB) 
                            if prior==0:
                                pass
                            elif prior==1:                                                                                     
                                dd1 = np.abs(roaw)
                                noo1 = int(np.max(dd1))
                                di1 = 0
                                for dd in range(1,noo1+1,1):
                                    di1 += np.log(mi.factorial(len(dd1[dd1==dd])))
                                kap1 = len(roaw[roaw!=0])
                                miua = np.sum(np.abs(roaw))                          
                                mSICs_case += 2.0*(di1+np.log(mi.factorial(miua-1))+kap1*np.log(exte)-kap1*np.log(laff)-np.log(mi.factorial(miua-kap1))-np.log(mi.factorial(kap1)))  
                                mSICs_case += 2.0*gamma*(float(miua-kap1+1)/float(kap1))+np.log(miua-kap1+2)+np.log(miua-kap1+1)   #the kap1*nplog(exte)=KlogN above is from the normalization of the lambda integral *(1/N)
                                mSICs_case -= 2.0*np.log(miua-kap1+2-(miua-kap1+1)*np.exp(-gamma/kap1))
                        else:
                            pass
                        if mSICs_case<best_mSIC2:
                            best_mSIC2 = mSICs_case
                            best_locs2 = roaw  
                            self.fluorOUT = fluorLin
                    else:
                        pass
        limit = 0                                                             #minimum number of extra steps to consider beyond what was found above  
        old_SIC = best_mSIC2
        old_locs = best_locs2
        new_locs = old_locs
        old_fluor = self.fluorOUT
        while limit<9:                                                        #MANUALLY SET maximum number of extra steps to consider beyond what was found above 
            for i in range(0,len(best_locs2),1):                              #Check various possible m's at each non-step location 
                ziplocs = np.zeros((6,len(best_locs2)))
                if best_locs2[i]==0:
                    for j in range(0,6,1):
                        ziplocs[j] = new_locs
                        ziplocs[j][i] = j+1
                        if ziplocs[j][i]>3:
                            ziplocs[j][i] += -7
                    for j in range(0,6,1):                                      
                        tamir = ziplocs[j]
                        mSICs_caseX = 0                         
                        fluorLin = fluorL
                        if self.negativeCheck(tamir,fluorLin)==0:             #Check whether you ever fall below zero fluorophores                        
                            intv = 0                                               
                            for k in range(0,winn,1):
                                if tamir[k]!=0 and k<winn-1:
                                    partial1 = np.copy(self.look_slice[intv:k+1])   #Use copy for security
                                    n_p = len(partial1)                                       
                                    intv = k+1
                                    mSICs_caseX += n_p*np.log(fluorLin*self.vF+self.vB)
                                    for q in range(0,len(partial1),1):
                                        mSICs_caseX += (partial1[q]  - fluorLin*self.mF - self.mB)**2.0/(fluorLin*self.vF+self.vB)
                                    fluorLin += tamir[k] 
                            partial2 = np.copy(self.look_slice[intv:winn])
                            if len(partial2)!=0:
                                n_p = len(partial2)
                                mSICs_caseX += n_p*np.log(fluorLin*self.vF+self.vB)
                                for q in range(0,len(partial2),1):
                                    mSICs_caseX +=  (partial2[q]  - fluorLin*self.mF - self.mB)**2.0/(fluorLin*self.vF+self.vB) 
                                if prior==0:
                                    pass
                                elif prior==1:                                                                                       
                                    dd2 = np.abs(tamir)
                                    noo2 = int(np.max(dd2))
                                    di2 = 0
                                    for de in range(1,noo2+1,1):
                                        di2 += np.log(mi.factorial(len(dd2[dd2==de])))
                                    kap2 = len(tamir[tamir!=0])
                                    miub = np.sum(np.abs(tamir))                                  
                                    mSICs_caseX += 2.0*(di2+np.log(mi.factorial(miub-1))+kap2*np.log(exte)-kap2*np.log(laff)-np.log(mi.factorial(miub-kap2))-np.log(mi.factorial(kap2)))
                                    mSICs_caseX += 2.0*gamma*(float(miub-kap2+1)/float(kap2))+np.log(miub-kap2+2)+np.log(miub-kap2+1)
                                    mSICs_caseX -= 2.0*np.log(miub-kap2+2-(miub-kap2+1)*np.exp(-gamma/kap2))                         
                            else:
                                pass
                            if mSICs_caseX<best_mSIC2:
                                best_mSIC2 = mSICs_caseX
                                best_locs2 = tamir   
                                self.fluorOUT = fluorLin                                                              
                        else:
                            pass 
                else:
                    pass
            limit += 1
            if best_mSIC2>=old_SIC:                                           #The extra step caused an increase in the criterion; stop further search (assumes smooth path to minimum of criterion value as a function of the number of steps, no ridges) 
                limit += 9
                best_mSIC = old_SIC                                           #IMPORTANT: If you are willing to accept a non-smooth path to minimum, you can set this to be equal to best_mSIC2 and set the value of th eprevious line to something smaller or zero. The code will be much slower. esp. if the value is set to zero. 
                best_locs = old_locs
                self.fluorOUT = old_fluor
            elif best_mSIC2<old_SIC:
                best_mSIC = best_mSIC2
                best_locs = best_locs2
                new_locs = best_locs2
        self.SIClocs = best_locs                                              #Final locations of events in the data set 
        self.SICbest = best_mSIC                                              #Final minimum SIC value                       
        self.SICmeans = self.findMeans2(self.look_slice, self.SIClocs)        #The means according to where the best steps are     
    def negativeCheck(self, Garray, Gnumber):                                 #Finds whether subtracting a number from an array ever results in a negative value in the array
        lista = np.zeros(len(Garray))
        for i in range(0,len(Garray),1):
            Gnumber += Garray[i]              
            lista[i] = Gnumber
        xalue = len(lista[lista<0])
        return xalue
    def findMeans2(self, dataz, stepz):                                       #Finds dwells means from the data given an array of equal length noting steps indices as ones and zeros
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
    def multichoose(self,n,k):                                               
        if k < 0 or n < 0: return "Error"
        if not k: return [[0]*n]
        if not n: return []
        if n == 1: return [[k]]
        return [[0]+val for val in self.multichoose(n-1,k)] + \
        [[val[0]+1]+val[1:] for val in self.multichoose(n,k-1)]   
    def dalistor(self, array_in, ksize, nsize):                               #Creates array with rows being possible step arrangements
        zelist = np.zeros(nsize)
        sym1, sym2, sym3, sym4, sym5, sym6 = '1', '2','3','4','5','6'         #MANUALLY PICK these set the various step sizes. Do not modify unless you are very sure what you are doing. 
        stringa = ksize*sym1+ksize*sym2+ksize*sym3+ksize*sym4+ksize*sym5+ksize*sym6
        kx = np.array([''.join(x) for x in it.permutations(stringa,ksize)])
        kz = np.zeros(len(kx))
        for i1 in range(0,len(kz),1):
            kz[i1] = float(kx[i1])
        ka = kz.flatten()
        tak = np.ones(ksize)
        for j1 in range(0, len(ka),1):
            vm = ka[j1]
            vz = np.zeros(ksize)
            for i1 in range(0,ksize,1):
                vz[ksize-(i1+1)] = (vm%10**(i1+1))//10**i1
                vm -= vz[ksize-(i1+1)]*10**i1
            for i1 in range(0, len(vz),1):
                if vz[i1]>3:
                    vz[i1] -= 7
            tak = np.vstack((tak, vz))
        tik = pd.DataFrame(tak)
        tik = tik.drop_duplicates().values
        tok = np.array(tik)
        fe = np.shape(tok)
        tek = np.zeros_like(tok, dtype=int)
        for i1 in range(0,fe[0],1):
            for j1 in range(0,fe[1],1):
                tek[i1][j1] = tok[i1][j1]  
        plax0 = np.copy(array_in)
        plax1 = np.array(np.nonzero(plax0))
        plax1  = plax1.flatten()
        for j1 in range(0,len(tek),1):
            plax2 = tek[j1]
            plax3 = np.zeros(nsize)
            for k1 in range(0,ksize,1):
                plax3[plax1[k1]] = plax2[k1]*plax0[plax1[k1]]
            zelist = np.vstack((zelist,plax3)) 
        zelist = zelist[1:]
        return zelist 

class Slicer(object):                                                         #Calculates initial (Step 2) estimate for step numbers and locations. Needs data, window size, how many fluors going in, stats.
   def __init__(self, array1, window1, degree1, arraystats):
       self.fluorr = degree1
       self.winnd = window1
       self.dataarray = array1
       self.evros = 100                                                       #MANUALLY SET value signifying how wide your local search is. Pick large numbers, it is inexpensive
       self.mB = arraystats[0]
       self.vB = arraystats[1]
       self.mF = arraystats[2]
       self.vF = arraystats[3]       
       self.levelz = np.zeros_like(self.dataarray)
       self.numb = int(len(self.dataarray)/self.winnd)                        #This ought to always be an integer anyway
       self.fosfor = np.zeros(self.numb)
       for i in range(0,self.numb,1):                                         #Sets each window to the fluorophore level found by localsearch  
           slicce = np.copy(self.dataarray[i*self.winnd:(i+1)*self.winnd])
           temprezstore = self.localsearch(slicce, int(self.fluorr), self.evros)
           self.levelz[i*self.winnd:(i+1)*self.winnd] = temprezstore
           self.fluorr = int((temprezstore[-1]-self.mB)/self.mF)
           self.fosfor[i] = self.fluorr
           print 'Step 2',i
   def localsearch(self, datarray, fluorr, aura):
       bonjovi = 0                                                            #Store results from each step of the local search
       seitan = 0                                                             #When this becomes 0 the local search is done
       saitan = 0                                                             #here just to allow the local search to run at least once
       leve = np.zeros_like(datarray)
       if fluorr>=aura:
           fluorFL = fluorr - aura                                            #what the fluorophore number is within the findL function, starting from the possible minimum
       elif fluorr==0:                                                        #if we have no fluorophores yet just look at the two next levels (three total)
           aura = 1
           fluorFL = 0
       elif fluorr>0 and fluorr<aura:
           aura = fluorr                                                      #if we have more than zero but fewer than our range fluorophores just look from 0 to double their number
           fluorFL = 0 
       fasma = 2*aura+1+fluorFL    
       while seitan!=0 or saitan==0:
           if saitan!=0: 
               fluorFL += seitan                                              #Set the fluorophore number to the new value for the next iteration
               if fluorFL<0:
                   raise ValueError('Congratulations! You have negative fluorophores: something must be amiss..')
           seitan  = 0                                                        #Reset the fluorophore number
           seitan += self.findL(datarray, fluorFL, fasma)                     #find fluorophore new number above previous threshold
           bonjovi += seitan                                                  #cumulative results
           saitan += 1
       leve += bonjovi*self.mF + self.mB           
       return leve         
   def findL(self, dataray, fluorL, fasma1):                                  #Find likelihood of slice belonging to fluorophore "group"
       mekos = len(dataray)
       pollap = np.zeros((fasma1,mekos))                                      #array for calculating the probabilities that a point belongs to a distribution 
       compar = np.zeros(fasma1)                                              #array to compare the probabilities that the window belongs to different distributions
       for i in range(0,fasma1,1):
           stdx = np.sqrt(self.vB + (fluorL+i)*self.vF)
           miux = self.mB + (fluorL+i)*self.mF  
           if stdx<=0 or miux<0:
               raise ValueError('Negative means or variances, something is very wrong')
           for j in range(0,mekos,1):
               pollap[i][j] = np.log(stdx*np.sqrt(2.0*np.pi)) + (dataray[j]-miux)**2.0/(2.0*stdx**2.0)    #find probability point belongs to a distribution            
       compar = np.sum(pollap, axis=1)                                        #Find probability window belongs to distribution
       group = np.argmin(compar)                                              #Find likeliest distribution
       return group     



       

           
 
    
