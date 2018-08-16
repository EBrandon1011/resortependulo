# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 08:37:32 2018

@author: brandon
"""
g=9.8   #Gravedad
l=1.0   #Longitud
kp=1.0    #Constante de amortiguamiento del péndulo
m=1.0

import numpy as np
import matplotlib.pyplot as plt

def dang(w):     #Dif. ángulo --> do/dt
    return w
    
def dw(ang,w):       #Dif. omega --> dw/dt
    return -g*np.sin(ang)/l-2*kp*w/(m*l**2)
    #return -g*ang/l-2*kp*w/(m*l**2)
    
#MÉTODO DE EULER
def Euler(ang0, w0, T, NI):
    #Condiciones iniciales
    ang=[]
    w=[]
    ang.append(ang0)
    w.append(w0)   
    h=T/NI     #Tamaño de paso
    i=1
    while i<NI:
        pw=w[i-1]+dw(ang[i-1],w[i-1])*h
        po=ang[i-1]+pw*h
        w.append(pw)
        ang.append(po)
        i+=1
    return ang, w
        
def cart(ang, NI):
    x=[]
    y=[]
    for i in range(NI):
        x.append(l*np.sin(ang[i]))
        y.append(-l*np.cos(ang[i]))
    return x, y
T=2     #Tiempo
n=500   #Iteraciones en el tiempo T
ang,w = Euler(np.pi/8, 0, T, n)
x,y=cart(ang,n)
#print(x,y)

#----------------------------------------------#

#TRAYECTORIA DEL PÉNDULO
plt.plot(x,y)

#----------------------------------------------#

#ANIMACIÓN DEL PÉNDULO

#plt.plot(x,y,'b--') 
#plt.grid(True)
#for i in range(100):
#    plt.xlim(-1,1)
#    plt.ylim(-1.1,0.1)
#    plt.plot(x,y,'b--')
#    plt.plot(0,0,'gs', markersize=20)
#    plt.plot(x[i],y[i],'bo', markersize=10)
#    plt.plot([0,x[i]],[0,y[i]],'y-')
#    plt.savefig('ft{}.png'.format(i))
#    plt.show()
#    plt.cla()
#plt.savefig('Pendulo.pdf')

#----------------------------------------------#

#ESPACIO DE FASES
#it=[]
#plt.xlabel(r'$\theta$',fontsize=15)
#plt.ylabel(r'$\dot{\theta}$', fontsize=15)
##plt.axvline(0, colors='g')
#plt.xlim(-1,1.8)
#for i in range(n+1):
#    it.append(i)
#plt.plot(ang,w,'#A52A2A', linewidth=2)
#plt.savefig('FasePenduloSin.pdf')
#plt.show()

#ÁNGULO CON RESPECTO AL TIEMPO
