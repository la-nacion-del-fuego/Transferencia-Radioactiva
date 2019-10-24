#!/usr/bin/env python3
import math
from astropy.modeling.blackbody import blackbody_lambda
from tqdm import tqdm

c=3e10 #cm/s
kB = 1.38e-16 #[ergk-1]

#temperature model [K]
def T(x):
    return 1e6

#density model [cm-3]
def n(x):
    return 1e7

#source function[erg/cm2 sec cm ster]
def S(x,wl):
    return blackbody_lambda(wl,T(x))

#opacity [cm-1]
def k(x, lw):
    nu = c/lw #frecuency
    #Ref Dulk (1985) eq. 21
    return 0.2*pow(n(x),2)*pow(T(x),-3/2)*pow(nu,-2)

#optical depth (adimensional)
def tau(dx,x, wl):
    return (dx/2.0)*(k(x-dx, wl) + k(x, wl))

def rayleigh(I, wl)
    return I*pow(wl, 4)/(2.0*c*kB)

N=6.96e3
I0 = 0.0 #[erg/cm2 sec cm ster]
nu = 1e8 #Hz
dx = 100e5 #[cm]
wl= c/nu #Armstrongs


layers = range(1,int(N+1))

I = I0
for i in tqdn(layers):
    x = float(i)*dx
    I = I*math.exp(-tau(dx,x, wl)) + S(x, wl)*(1-math.exp(-tau(dx,x, wl)))
    pass
print ("%e"%rayleigh(I.value,wl))
