from scipy.integrate import odeint
from numpy import linspace
from math import pi, sin, exp
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
def inj(t):
    return ((t>1000 and t<1500) or (t>2000 and t<2500))*0.01#*(sin(t/(2*pi*10)+1))

AM=22
BM=13
AH=0.5
BH=4
AN=1.5
BN=0.1
AQ=1.5
BQ=0.025

def ddt(V,t,p):
    # V[0] = dendrite voltage
    # V[1] = Soma Voltage

    if not p.state and V[1]>0:
        p.state = True
        if (t-p.t0>0.6):
            # print(f"{t} : {p.t0}")
            p.state=False
        p.m0=p.m
        p.h0=p.h
        p.n0=p.n
        p.q0=p.q


    def gateVal(alpha,beta,v0):
        # public double getValueOn(double t) 
    	# value = v0 * Math.exp(-beta*(t - t0));
        # public double getValueOff(double t) 
    	# value = 1 + (v0 - 1) * Math.exp(-alpha*(t - t0));
        return (p.state)*(v0 * exp(-beta*(t - p.t0))) + \
            (not p.state)*(1 + (v0 - 1) * exp(-alpha*(t - p.t0)))

    p.m = gateVal(AM,BM,p.m0)
    p.h = gateVal(AH,BH,p.h0)
    p.n = gateVal(AN,BN,p.n0)
    p.q = gateVal(AQ,BQ,p.q0)
    
    Iion = p.gNa * p.m**3 * p.h*(V[1]-p.ENa) +\
           p.gKf * p.n**4 * (V[1]-p.EK) +\
           p.gKs * p.q**2 * (V[1]-p.EK)
        
    dVdt = [
        (-p.Isyn_d-p.gld*(V[0]-p.El)-p.gc*(V[0]-V[1])+p.Iinj_d(t))/p.Cd,
        (-p.Isyn_s-p.gls*(V[1]-p.El)-p.gc*(V[1]-V[0])-Iion+p.Iinj_s(t))/p.Cs
    ]
    # print(dVdt[0])
    return dVdt


class Object(object):
    pass
p = Object()

# Whether pulse is active
p.state=True
p.t0=0 # time of last pulse

# m,h,n,q values at pulse
p.m0=0
p.h0=0
p.n0=0
p.q0=0

# Sodium parameters
p.ENa=120 #equilibrium
p.gNa=30  #maximal conductance

# Potassium parameters
p.EK=-10  #equilibrium
p.gKf=4   #maximal fast conductance
p.gKs=16  #maximal slow conductance

# Leak parameters
p.El=0    # equilibrium leak

# Resistance of dendrite and soma
p.Rmd=12350 #Ohm cm^2
p.Rms=1075  #Ohm cm^2
# length of dendrite and soma
p.ld=5800   #um 5519-6789
p.ls=79     #um 77.5-82.5
## Convert units to cm
p.ld=p.ld/1000
p.ls=p.ls/1000
## End convert units

# Radius of dendrite and soma
p.rd=79/2   #um 38.75-41.75
p.rs=50/2   #um 20.75-31.75
## Convert units to cm
p.rd=p.ld/1000
p.rs=p.ls/1000
## End convert units

# maximal conductance of dendrite and soma
p.gld= (2*pi*p.rd*p.ld)/(p.Rmd)
p.gls= (2*pi*p.rs*p.ls)/(p.Rms)

p.Ri = 50 # Ohm cm : Cytoplasm resistance

# Rheobase of soma
p.rheo = 4.5 #nA: 3.5-6.5

# Calculate coupling parameter between soma and dendrite
p.gc = 2 / ( ((p.Ri*p.ld)/(pi*p.rd**2)) + ((p.Ri*p.ls)/(pi*p.rs**2)) )

# Calculate rn from conductances and coupling:
p.rn = 1/(p.gls + (p.gld*p.gc)/(p.gld+p.gc) )

# Calculate threshold of soma
p.threshold = p.rheo*p.rn
print(p.threshold)

p.Iinj_d = inj
p.Iinj_s = lambda t: 0

p.Isyn_d=0
p.Isyn_s=0

dt=0.1
tstop=4000

# Set capacitances
p.Cm=1 # Memrane specific capacitance
p.Cd=2*pi*p.rd*p.ld*p.Cm # dendrite capacitance
p.Cs=2*pi*p.rs*p.ls*p.Cm # soma capacitance

sol = odeint(ddt,[0,0],linspace(0,tstop,int(tstop/dt)),args=(p,))

plot(sol[:,0])
plot(sol[:,1])
show()
