from scipy.integrate import odeint
from numpy import linspace
from math import pi, sin
from matplotlib.pyplot import *

# class compartment:
#     """Single Compartment of model, consitent of update function with respective t, membrane voltage etc."""

#     def updateStep(self):
#         self.rec = odeint(self.ddt,self.d0,linspace(0,self.tstop,self.tstop/self.dt)) # Time is in milliseconds
        
#     def __init__(self,tstop, dt, dCompdt):
#         self.tstop = tstop
#         self.dt = dt
#         self.nsteps = tstop/dt
#         self.ddt = dCompdt

# class neuron:
#     """Composed of multiple compartments that may depend on each other"""

#     def __init__(self,comparts):

AM=22
BM=13
AH=0.5
BH=4
AN=2.2
BN=0.76
AQ=1.5
BQ=0.025

def ddt(V,t,p):
    # V[0] = dendrite voltage
    # V[1] = Soma Voltage
    # V[2] = m/dt
    # V[3] = h/dt
    # V[4] = n/dt
    # V[5] = q/dt

    # TODO: Implement smarter pulse paradihm implementation

    # if V[0]>0:
    #     a_m=AM
    #     b_m=0
    #     a_h=0
    #     b_h=BH
    #     a_n=AN
    #     b_n=0
    #     a_q=AQ
    #     b_q=0
    # else:
    a_m=0
    b_m=BM
    a_h=AH
    b_h=0
    a_n=0
    b_n=BN
    a_q=0
    b_q=BQ

    Iion = p.gNa * V[2]**3 * V[3]*(V[1]-p.ENa) +\
           p.gKf * V[4]**4 * (V[1]-p.EK) +\
           p.gKs * V[5]**2 * (V[1]-p.EK)
        
    dVdt = [
        (-p.Isyn_d-p.gld*(V[0]-p.El)-p.gc*(V[0]-V[1])+p.Iinj_d(t))/p.Cd,
        (-p.Isyn_s-p.gls*(V[1]-p.El)-p.gc*(V[1]-V[0])-Iion+p.Iinj_s(t))/p.Cs,
        a_m*(1-V[2])-b_h*V[2],
        a_h*(1-V[3])-b_h*V[3],
        a_n*(1-V[4])-b_n*V[4],
        a_q*(1-V[5])-b_q*V[5]
    ]
    print(dVdt[0])
    return dVdt


class Object(object):
    pass
p = Object()
    
p.ENa=120
p.gNa=30

p.EK=-10
p.gKf=4
p.gKs=16

p.El=0
Rmd=12350
Rms=1075
ld=5800
ls=79
rd=79/2
rs=50/2
p.gld= (2*pi*rd*ld)/(Rmd)
p.gls= (2*pi*rs*ls)/(Rms)
Ri=50
p.gc= 2 / ( ((Ri*ld)/(pi*rd**2)) + ((Ri*ls)/(pi*rs**2)) )

def inj(t):
    return (t>1000 and t<1500)*-800#*(sin(t/(2*pi*10)+1))

p.Iinj_d = inj
p.Iinj_s = lambda t: 0

p.Isyn_d=0
p.Isyn_s=0

dt=0.1
tstop=4000

Cm=70
p.Cd=2*pi*rd*ld*Cm
p.Cs=2*pi*rs*ls*Cm

sol = odeint(ddt,[0,0,0,0,0,0],linspace(0,tstop,int(tstop/dt)),args=(p,))

plot(sol[:,0])
plot(sol[:,1])
show()
