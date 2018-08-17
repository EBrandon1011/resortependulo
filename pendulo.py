# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 08:37:32 2018

@author: brandon
"""
g=9.81  #Aceleración gravitacional.
k=8.3   #Constante de resorte [N/m].
br=0.724794  #Constante de amortiguamiento del resorte.  br=$\beta$
bp=0  #Constante de amorgiguamiento del péndulo.  bp=$\gamma$
M=0.200   #Masa del objeto - Resorte.
m=0.096   #Masa del objeto - Péndulo.
l=0.509   #Longitud del péndulo.

from numpy import cos, sin, pi
import matplotlib.pyplot as plt

#Definiendo las ecuaciones diferenciales dq_i/dt=f_i(q):

def dx(v):              #Dif. x --> dx/dt
    return v
    
def dv(x, v, ang, w):   #Dif. v --> dv/dt
    num=2*bp*cos(ang)*w-2*br*l*v+m*g*l*sin(ang)*cos(ang)-k*l*x+m*l**2*sin(ang)*w**2
    den=l*(M+m*(sin(ang))**2)
    return num/den
    
def dang(w):            #Dif. o --> do/dt
    return w
    
def dw(x, v, ang, w):   #Dif. w --> dw/dt
    num=-2*bp*m*w-2*bp*M*w+k*x*l*m*cos(ang) +2*br*l*m*v*cos(ang)-g*l*m**2*sin(ang) - g *l* m *M*sin(ang) - l**2* m**2* w**2*cos(ang)*sin(ang)
    den=m*l**2*(M+m*(sin(ang))**2) 
    return (1)*num/den
    
#APLICACIÓN DEL MÉTODO DE EULER
def Euler(x0, v0, ang0, w0, T, NI):
    #Condiciones iniciales
    x=[]
    v=[]
    ang=[]
    w=[]
    time=[]
    x.append(x0)
    v.append(v0)   
    ang.append(ang0)
    w.append(w0)
    time.append(0)
    h=T/NI     #Tamaño de paso
    n=1
    while n<NI:
        pw=w[n-1]+h*dw(x[n-1],v[n-1],ang[n-1],w[n-1])   #Obteniendo  w[n]
        pv=v[n-1]+h*dv(x[n-1],v[n-1],ang[n-1],w[n-1])   #Obteniendo  v[n]
        w.append(pw)
        v.append(pv)
        ang.append(ang[n-1]+h*dang(w[n-1]))
        x.append(x[n-1]+h*dx(v[n-1]))
        time.append(n*h)
        n+=1
    return x, v, ang, w, time

#Condiciones iniciales
x0=0.0
v0=0.0
ang0=10*pi/180
w0=0.0
T=10
N=T*1000
x,v,ang,w,t=Euler(x0, v0, ang0, w0, T, N)

#----------------------------------------------#

# TRAYECTORIAS
xp=[]   #Posición del péndulo
yp=[]
yr=[]
for i in range(N):
    ypi=(-1)*l*cos(ang[i])
    xpi=x[i]+l*sin(ang[i])
    yr.append(0)
    xp.append(xpi)
    yp.append(ypi)

plt.xlabel(r'$x$', fontsize=15)   
plt.ylabel(r'$y$', fontsize=15)
plt.xlim(-1.5,1.5)
plt.ylim(-0.7,0.7)  
plt.plot(0,0,'ro')  
plt.plot(x,yr,'#A52A2A',label=r'Resorte $M$')
plt.plot(xp,yp, '#006400', label=r'Péndulo $m$')
plt.legend()
plt.savefig('trayectoria_Amort_170g_20s.pdf')


#----------------------------------------------#

## ANIMACIÓN DEL RESORTE Y PÉNDULO

#i=0
#while i<N:
#    plt.xlim(-1,1)
#    plt.ylim(-1.1,0.1)
#    plt.plot(x[i],0,'#2F4F4F',marker='s', markersize=15)  #Carrito
#    plt.plot(xp[i],yp[i],'#191970', marker='o', markersize=10) #Péndulo
#    plt.plot([x[i],0],[xp[i],yp[i]],'#778899')
#    plt.savefig('pr{}.png'.format(i),dpi=150)
#    plt.show()#
#    plt.cla()
#    i+=10



#----------------------------------------------#

## X, ANGULO RESPECTO AL TIEMPO
#plt.ylim(-0.1,0.1)
#plt.xlabel(r'$t$', fontsize=15)
#plt.ylabel(r'$x$', fontsize=15)
#plt.plot(t,x, '#8B0000')
#plt.savefig('x_5g.pdf')

#----------------------------------------------#

## DIAGRAMA DE FASE
##tiempo=20 / iteraciones=100,000
##tiempo=100 / iteraciones=100,000
#plt.xlim(-0.25,0.25)
#plt.xlabel(r'$\theta$', fontsize=15)     
#plt.ylabel(r'$\omega$', fontsize=15) 
#plt.plot(ang,w, '#00008B') 
#plt.savefig('DFaseTW_10grados_Amort.pdf')
