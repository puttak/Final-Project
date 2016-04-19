# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 18:27:18 2016

@author: samsung1
"""
import scipy as sp
import numpy as np
from scipy import integrate
import pylab as pl
import matplotlib.pyplot as plt

print 'Based on the capacity, enter the volume and influent flowrate values. \n '

V = raw_input('Volume of Tank:     (Default = 12000(m3))' '\n')
if V == '': V = 12000
else: V= float(V)

Q = raw_input('Inlet Sludge flow rate:     (Default = 300(m3/day))' '\n')
if Q == '': Q = 300.00
else: Q= float(Q)

"""The different Parameters defining our wastewater system based on Monod's Model"""  

print 'The parameters based on the Monods Constant Model can be entered upon the users microorganism system. \n '

Kdr = raw_input('Kdr:          (Default = 2 d^(-1))' '\n')
if Kdr == '': Kdr = 2
else: Kdr = float(Kdr)

Y = raw_input('Y:          (Default = 0.3 dimensionless)' '\n')
if Y == '': Y = 0.3
else: Y = float(Y)

Kh = raw_input('Kh:          (Default = 0.15 d^(-1))' '\n')
if Kh == '': Kh = 0.15
else: Kh = float(Y)

kas = raw_input('kas:          (Default = 8 g/g/day)' '\n')
if kas == '': kas= 8 
else: kas = float(kas)

Kas = raw_input('Kas:          (Default = 0.045 g/dm3)' '\n')
if Kas == '': Y = 0.045
else: Kas = float(Kas)

Yac = raw_input('Yac:          (Default = 0.2 g/g)' '\n')
if Yac == '': Yac = 0.2
else: Yac = float(Yac)

Da = raw_input('Da:          (Default = 0.1 d^(-1))' '\n')
if Da == '': Da = 0.1
else: Da = float(Da) 

fd= raw_input('fd:          (Default = 0.73 dimensionless)' '\n')
if fd == '': fd = 0.73
else: fd = float(fd) 

kvm = raw_input('kvm:          (Default = 6.2 g/g/day)' '\n')
if kvm == '': kvm = 6.2
else: kvm = float(kvm)

Kvm = raw_input('Kvm:          (Default = 0.045 g/dm3)' '\n')
if Kvm == '': Kvm = 0.045
else: Kvm = float(Kvm)

Dm= raw_input('Dm:          (Default = 0.015 d^(-1))' '\n')
if Dm == '': Dm = 0.015
else: Dm = float(Dm)

Ymc = raw_input('Ymc:          (Default = 0.057 g/g)' '\n')
if Ymc == '': Ymc = 0.057
else: Ymc = float(Ymc)

"""Kil = raw_input('Kil:          (Default = 1 g/dm3)' '\n')
if Kil == '': Kil = 1
else: Kil = float(Kil)

Ki2 = raw_input('Ki2:          (Default = 0.3 g/dm3)' '\n')
if Ki2 == '': Ki2 = 0.3
else: Ki2 = float(Ki2)

Ka= raw_input('Ka:          (Default = 10^(-4.5) d^(-1))' '\n')
if Ka == '': Ka = 0.3
else: Ka = float(Ka) """

print 'Enter the input variables. \n '

Xas0 = raw_input('Conc of AS biomass in influent:          (Default = 10(g/dm3))' '\n')
if Xas0 == '': Xas0 = 10.0
else: Xas0 = float(Xas0)

ps0 = raw_input('Conc of PS requiring hydrolysis influent:          (Default = 1(g/dm3))' '\n')
if ps0 == '': ps0 = 1
else: ps0 = float(ps0) 

ss0 = raw_input('Conc of soluble substrate in influent:          (Default = 5(g/dm3))' '\n')
if ss0 == '': ss0 = 5.0
else: ss0 = float(ss0)

xa0 = raw_input('Conc of acidogens influent: (Default = 0.05(g/dm3))' '\n')
if xa0 == '': xa0 = 0.05
else: xa0 = float(xa0)

xfa0 = raw_input('Conc of volatile fatty acid in influent:          (Default = 2(g/dm3))' '\n')
if xfa0 == '': xfa0 = 2
else: xfa0 = float(xfa0)

xm0 = raw_input('Conc of methanogens in influent:     (Default = 0.04(g/dm3))' '\n')
if xm0 == '': xm0 = 0.04
else: xm0= float(xm0)

m0 = raw_input('Conc of methane in influent:     (Default = 0.0(g/dm3))' '\n')
if m0 == '': m0 = 0.0
else: m0= float(m0)


""" Effluent Conditions Predefined """  

Xas1=11.31 #Conc of AS biomass in effluent in g/dm3
ps1=0 #Conc of PS biomass in effluent in g/dm3
ss1=0 #Conc of souluble substrate in effluent in g/dm3
xa1=0.1 #Conc of acidogens in effluent in g/dm3
xfa1=0 #Conc of volatile fatty acids in effluent in g/dm3
xm1=0.12 #Conc of methanogens in effluent in g/dm3
m1=0 ##Conc of methane  in effluent in g/dm3
theta=V/Q #day inverse


uo = np.array([Xas1,ps1,ss1,xa1,xfa1,xm1,m1])
iters = np.arange(0.0,15,0.1)

def sludge(u,t):
    dxasdt =(Xas0-u[0])/theta- Kdr*u[0]
    dpsdt= (ps0-u[1])/theta +(1-Y)*Kdr*u[0]-Kh*u[1]
    dssdt=(ss0-u[2])/theta+Kh*u[1]+ Y*Kdr*u[0]-(kas*u[2]*u[3]/(0.045+u[2]))
    dxadt=(xa0-u[3])/theta + Yac*(kas*u[2]*u[3])/(0.045+u[2])-Da*u[3]
    dxfadt=(xfa0-u[4])/theta + (kas*u[2]*u[3])/(0.045+u[2])- Yac*(kas*u[2]*u[3])/(0.045+(u[2]))-fd*Da*u[3]-kvm*u[4]*u[5]/(0.045+(u[4]))
    dxmdt=(xm0-u[5])/theta + (Ymc*kvm*u[3]*u[4]/(0.045+u[3]))-(Dm*u[5])
    dmdt=(m0-u[6])/theta + (kvm*u[4]*u[5]/(0.045+u[4]))-(Ymc*kvm*u[4]*u[5]/(Kvm+u[4]))
	
	
    return np.array([dxasdt,dpsdt,dssdt,dxadt,dxfadt,dxmdt,dmdt], dtype = float)


# The equations are integrated in vector form using odeint

results = integrate.odeint(sludge,uo,iters)

plt.title('Anaerobic Digester digesting Sewage Sludge(Monods Constant)')
plt.xlabel('Time (days)') # set x-axis label
plt.ylabel('Conc (g/L)') # set y-axis label
plt.plot(iters,results[:,0],'Blue',iters,results[:,1],'Green',iters,results[:,2],'Red',iters,results[:,3],'Turquoise',iters,results[:,4],'Magenta',iters,results[:,5],'Olive',iters,results[:,6],'Black')
plt.legend(['Activated Sludge Biomass', 'Particulate Substrate', 'Soluble Substrate','Acidogenic Biomass','Volatile fatty acid','Methanogens','Methane'], loc='')
fig_size = plt.rcParams["figure.figsize"]
fig_size[0]=10
fig_size[1]=10
plt.rcParams["figure.figsize"] = fig_size
plt.show()

