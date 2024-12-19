#!/usr/bin/python3

# File to process data from Bondville sites (From November 7 2024)
# General program flow:
#   First, read in 30METDAT data file, and assign corresponding MET LABELS
#   
# Usage: ./bv13analysis.py    
# Version .25, adding in plotting
#import statsmodels.api as sm
import sys        #Needed for command line args
import csv        #For reading in text files
import subprocess #System calls
import re         #Reg expressions
import glob       #For glob command to list files
from timeit import timeit 
import pandas as pd    #Runs on pandas 1.1.5
import numpy as np
#from windrose import WindroseAxes
import io         #For running StringIO through read_csv
import matplotlib.pyplot as plt
import warnings   #For filtering warnings like All-NaN slice encountered
#warnings.filterwarnings("error")
#metfile = sys.argv[1]                     # first command line argument is the year (2022)
#flxfile = sys.argv[2]                      # second command line arg is the day of year (101)
#fyear = iyear - 2000
fyear = 22
sb = 5.67e-8
#
#yrday = "%4d %d" % (iyear,idoy)
metfile = "BV13_180-220.30METDAT"
#
iyear = 2013
idoy = 180
yrday = "%4d %d" % (iyear,idoy)              # put into a string then
a2 = pd.to_datetime(yrday, format = '%Y %j') # convert to pd datetime 
start_datetime1 = a2.strftime('%Y-%m-%d 00:30:00')           # create start times for 30 minute data
#
df_met = pd.read_csv(metfile,delimiter= r"\s+",header=None) # read in flx data
df_metlab = pd.read_csv("./BV13_DAT.LABELS", delimiter=r"\s+",header=None) # read in the labels
met_cols = df_metlab[1].to_list()   # put the columns we want in a list and
df_met.columns = met_cols                                          # assign the column names to the input file
print (met_cols)
dt30min_index = pd.date_range(start_datetime1, periods=1968, freq="30min")         # daily pdtime index 30 min
df_met['Datetime'] = dt30min_index
df_met.set_index(df_met['Datetime'],inplace = True)                            # then make it the index
print (df_met)
#
df_met['sig_w'] = df_met['w2']**0.5
df_met['sig_u'] = df_met['u2']**0.5
df_met['sig_v'] = df_met['v2']**0.5
df_met['Tbulk'] = df_met['Ta']+273.18
df_met['ustar'] = (-df_met['uw'])**0.5
df_met['Ri_bulk'] = (df_met['Tair10']-df_met['Ts'])*9.8*10.0/(df_met['Tbulk']*df_met['wsp10']**2) #
#
df_met['itw'] = df_met['sig_w']/df_met['wsp10']     # compute the turbulence intensities using mean wind at 10 m
df_met['itu'] = df_met['sig_u']/df_met['wsp10']
df_met['itv'] = df_met['sig_v']/df_met['wsp10']
df_met['itustar'] = df_met['ustar']/df_met['wsp10']
#
df_met['C1_ice'] = 1.0003+df_met['pres']/1000000.0
df_met['C1_wat'] = 1.0007+df_met['pres']/1000000.0
df_met['es_surf_ice'] = df_met['C1_ice']*6.1115*np.exp(22.452*df_met['Ta']/(272.55+df_met['Ta']))
df_met['es_surf'] = df_met['C1_wat']*6.1121*np.exp(15.502*df_met['Ta']/(240.97+df_met['Ta']))
df_met['e_air_ice'] = df_met['C1_ice']*0.01*df_met['RH']*6.1115*np.exp(22.452*df_met['Ta']/(272.55+df_met['Ta']))
df_met['e_air'] = df_met['C1_wat']*0.01*df_met['RH']*6.1121*np.exp(15.502*df_met['Ta']/(240.97+df_met['Ta']))
df_met.loc[df_met['Ts'] < 0, 'es_surf'] = df_met.es_surf_ice
df_met.loc[df_met['Ta'] < 0, 'e_air'] = df_met.e_air_ice
df_met['factor'] = 1 - 0.378*df_met['e_air']/df_met['pres']
df_met['density'] = df_met['pres']*100*df_met['factor']/(287.04*(df_met['Ta']+273.15))  #  kg m-3
df_met['cp'] = 1004.7*(0.522*df_met['e_air']/df_met['pres'] + 1)      #specific heat moist air J kg-1 K -1
df_met['qs'] = df_met['es_surf']*0.622/(df_met['pres']-0.378*df_met['es_surf'])
df_met['qa'] = df_met['e_air']*0.622/(df_met['pres']-0.378*df_met['e_air'])
df_met['Lat'] = (2500 - 2.361*df_met['Ts'])*1000.0
df_met.loc[df_met['Ts'] < 0, 'Lat'] = 2834000.0
#df_tot['sig_w_mod'] = df_tot['wsp_3m']*0.07
#df_tot['sig_w_mod'] = df_tot['wsp_3m']*(0.15 - df_tot['Ri_bulk']*0.25)
#df_tot['u*_mod'] = df_tot['wsp_3m']*(0.1398 - df_tot['Ri_bulk']*0.25)
#df_tot['u*_mod'] = df_tot['wsp_3m']*0.063
df_met['alt_LE'] = (df_met['qs']-df_met['qa'])*df_met['density']*df_met['sig_w']*df_met['Lat']*0.06
#df_tot['u*_sqe'] = (df_tot['u*_mod']-df_tot['u*_3m'])**2
#
df_met['alt_H'] = -(df_met['Tair10']-df_met['Ts'])*df_met['sig_w']*df_met['density']*df_met['cp']*.06
#df_tot['u*_wsp5'] = df_tot['u*_3m']/(df_tot['wsp_3m'])
#df_tot['z0'] = 2.5/np.exp(0.4/df_tot['u*_wsp5'])
#df_tot['sig_w_wsp5'] = df_tot['sig_w_3']/df_tot['wsp_3m']
#df_tot['Rnet']=df_tot['Rg_in_avg']-df_tot['Rg_out_avg']+df_tot['Lw_in_avg']-df_tot['Lw_out_avg']
#df_tot['ghf']=(df_tot['ghflx_a']+df_tot['ghflx_b']+df_tot['ghflx_c'])/3.0 +(df_tot['stor_a']+
#    df_tot['stor_b']+df_tot['stor_c'])/3.0
#df_tot['GHF_ratio'] = df_tot['ghf']/df_tot['Rnet']
#print (df_tot['GHF_ratio'])
#df_tot_ratio = df_tot[df_tot['Rnet'] > 300]
#df_tot['alt_LE'] = df_tot['Rnet']-df_tot['H_3m']-df_tot['ghf']
newtot = df_met.loc['2013-06-29':'2013-08-08']
fx1 = newtot[(newtot['wsp10'] > 2) & (newtot['ustar'] > 0)]
#fx0_emm = fx0.resample('10D').mean()

#print (fx0_emm['emm'])
#print (newtot)
#fx1 = newtot[(newtot['qc_H_3'] < 1) & (newtot['wsp_10m'] > 2.5) & (newtot['qc_H_10'] < 1) & (newtot['u*_wsp5'] > 0)] # &
#fx2 = newtot[(newtot['qc_H_3'] < 1) & (newtot['wsp_10m'] > 2.5) & (newtot['u*_wsp5'] > 0) &
#             (newtot['Ri_bulk'] > -0.1) & (newtot['Ri_bulk'] < 0.1) & (newtot['stheta_5m'] < 20)]
#
fig = plt.figure(facecolor=(1,0,0,0.1),figsize=(24,12))
ax = fig.add_subplot(2,2,4)
plt.axis([-0.5,0.5,0,0.5])
plt.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
plt.xlabel('Ri',fontsize=20,fontweight = "bold")
plt.ylabel('sig_u/U',fontsize=16,fontweight = "bold")
plt.scatter(fx1['Ri_bulk'],fx1['itu'],color = 'red', edgecolor = 'black')
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.suptitle('Turbulence Statistics vs Stability  (Jul - Aug 2013)', fontsize=24, fontweight = "bold")
#Y = fx1['alt_H']
#X = fx1['H_3m']
#model = sm.OLS(Y, X).fit()
#predictions = model.predict(X) # make the predictions by the model
#print(predictions) 
# Print out the statistics
#print(model.summary())
#
plt.subplot(222)
plt.axis([-0.5,0.5,0,0.2])
plt.ylabel(r'$\sigma_w/\overline{U}_{10m}$',fontsize=20,fontweight = "bold")
plt.xlabel('Richardson Number',fontsize=16,fontweight = "bold")
plt.axvline(x = 0,color = 'b',label = 'axvline - full height')
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.scatter(fx1['Ri_bulk'],fx1['itw'],color='cyan', edgecolor = 'black')
#
plt.subplot(223)
plt.axis([-0.5,0.5,0,0.5])
plt.ylabel(r'$\sigma_v/\overline{U}_{10m}$',fontsize=20,fontweight = "bold")
plt.xlabel('Ri_bulk',fontsize=20)
#plt.axvline(x = 0,color = 'b',label = 'axvline - full height')
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.scatter(fx1['Ri_bulk'],fx1['itv'],color='blue', edgecolor = 'black')
#
plt.subplot(221)
plt.axis([-0.5,0.5,0,0.2])
plt.ylabel(r'$u_*/\overline{U}_{10m}$',fontsize = 16,fontweight="bold")
plt.xlabel('Richardson Number',fontsize=20,fontweight = "bold")
plt.axvline(x = 0,color = 'b',label = 'axvline - full height')
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.scatter(fx1['Ri_bulk'],fx1['itustar'],color='green', edgecolor = 'black')
#
#######################################################################################

#######################################################################################
fig2 = plt.figure(facecolor=(1,0,0,0.1),figsize=(12,12))
plt.suptitle('Simple Model Energy Fluxes- (Jan - Mar 2022)', fontsize=24, fontweight = "bold")
ax1 = fig2.add_subplot(2,1,1)
plt.axis([0,600,0,600])
plt.plot(ax1.get_xlim(), ax1.get_ylim(), ls="--", c=".3")
plt.ylabel(r'$mod-LE(W m^{-2})$',fontsize=20,fontweight = "bold")
plt.xlabel(r'$LE_(W m^{-2})$',fontsize=20,fontweight = "bold")
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.scatter(fx1['LE'],fx1['alt_LE'],color = 'red', edgecolor = 'black')
#
ax2 = plt.subplot(212)
plt.axis([-300,500,-300,500])
plt.plot(ax2.get_xlim(), ax2.get_ylim(), ls="--", c=".3")
plt.ylabel(r'$mod-H (W m^{-2})$',fontsize=20,fontweight = "bold")
plt.xlabel('H (W/m2)',fontsize=20,fontweight = "bold")
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.scatter(fx1['H'],fx1['alt_H'],color = 'red', edgecolor = 'black')
#

#fig, axes = plt.subplots(1, 2)  
#axes[0, 0].scatter(fx1['Ri_bulk'],fx1['u*_wsp10']) 
#axes[0, 1].scatter(fx1['Ri_bulk'],fx1['sig_w_wsp10']) 
#axes[1, 0].plot(x, y3, 'b--o') 
#axes[1, 1].plot(x, y4, 'r--o')
#plt.scatter(fx1['Ri_bulk'], fx1['u*_wsp10'])
plt.show()
