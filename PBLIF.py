from math import pi, exp
import numpy as np

def inj(t):
    return ((t>1000 and t<1500) or (t>2000 and t<2500))*0.0001#*(sin(t/(2*pi*10)+1))

class PBLIF:
    """A Pulse Based Leaky Integrate and Fire model"""
    def __init__(self):
        # Initial State of variables
        # t, Vd, Vs
        self.t=[0]
        self.V=[[0,0]]
        # m,h,n,q
        self.m=[0]
        self.h=[0]
        self.n=[0]
        self.q=[0]
        # values for analytical solution to pulse based ODEs
        self.m0=0
        self.h0=0
        self.n0=0
        self.q0=0
        # Gate pulse parameters
        self.AM=22
        self.BM=13
        self.AH=0.5
        self.BH=4
        self.AN=1.5
        self.BN=0.1
        self.AQ=1.5
        self.BQ=0.025
        
        # Whether pulse is active
        self.state=False
        self.pulseState=True
        # Start of pulse
        self.t0=0
        
        # Sodium parameters
        self.ENa=120 # mV equilibrium
        self.gNa=30  # mS/cm^2maximal conductance
        
        # Potassium parameters
        self.EK=-10  # mV equilibrium
        self.gKf=4   # mS/cm^2 maximal fast conductance
        self.gKs=16  # mS/cm^2 maximal slow conductance
        
        # Leak parameters
        self.El=0    # mV equilibrium leak
        
        # Resistance of dendrite and soma
        self.Rmd=12350 #Ohm cm^2
        self.Rms=1075  #Ohm cm^2
        # length of dendrite and soma
        self.ld=5800   #um 5519-6789
        self.ls=79     #um 77.5-82.5
        ## Convert units to cm
        self.ld=self.ld/1000
        self.ls=self.ls/1000
        ## End convert units
        
        # Radius of dendrite and soma
        self.rd=79/2   #um 38.75-41.75
        self.rs=50/2   #um 20.75-31.75
        ## Convert units to cm
        self.rd=self.ld/1000
        self.rs=self.ls/1000
        ## End convert units
        
        # maximal conductance of dendrite and soma
        self.gld= (2*pi*self.rd*self.ld)/(self.Rmd)
        self.gls= (2*pi*self.rs*self.ls)/(self.Rms)
        
        self.Ri = 50 # Ohm cm : Cytoplasm resistance
        
        # Rheobase of soma
        self.rheo = 4.5/1000000 #nA/1000000: mA: 3.5-6.5
        
        # Calculate coupling parameter between soma and dendrite
        self.gc = 2 / ( ((self.Ri*self.ld)/(pi*self.rd**2)) + ((self.Ri*self.ls)/(pi*self.rs**2)) )
        
        # Calculate rn from conductances and coupling:
        self.rn = 1/(self.gls + (self.gld*self.gc)/(self.gld+self.gc) )
        
        # Calculate threshold of soma
        # self.threshold = self.rheo*self.rn
        self.threshold=1
        # print(self.threshold)
        
        self.Iinj_d = inj # Test input Current
        self.Iinj_s = lambda t: 0
        
        self.Isyn_d=0
        self.Isyn_s=0
        
        self.dt=0.1
        self.tstop=2000
        
        # Set capacitances
        self.Cm=1 # Memrane specific capacitance
        self.Cd=2*pi*self.rd*self.ld*self.Cm # dendrite capacitance
        self.Cs=2*pi*self.rs*self.ls*self.Cm # soma capacitance

    def ddt(self,slope,t,V):
        # V[0] = dendrite voltage
        # V[1] = Soma Voltage
        pulseState=True
        def gateVal(alpha,beta,v0,pulseState):
            # public double getValueOn(double t) 
    	    # value = v0 * Math.exp(-beta*(t - t0));
            # public double getValueOff(double t) 
    	    # value = 1 + (v0 - 1) * Math.exp(-alpha*(t - t0));
            ret = 0
            try:
                ret = (pulseState)*(v0 * exp(-beta*(t - self.t0))) + \
                    (not pulseState)*(1 + (v0 - 1) * exp(-alpha*(t - self.t0)))
            except OverflowError:
                print(f"{alpha}, {beta}, {v0}, {t}, {self.t0}")
                # if (t>900):
                # print(t-self.t0)
            return ret
        if slope==1:
            if ( self.state ):
                if(t - self.t0 > 0.6):
                    self.t0=t
                    self.m0=self.m
                    self.h0=self.h
                    self.n0=self.n
                    self.q0=self.q
                    self.state=False
                    # print("State change")
                    pulseState=True;
                else:
    	            pulseState=False;
            else:
    	        pulseState=True;
            
        if V[1]>self.threshold:
            self.t0=t
            self.m0=self.m
            self.h0=self.h
            self.n0=self.n
            self.q0=self.q
            self.state = True
        
        self.m = gateVal(self.AM,self.BM,self.m0,pulseState)
        self.h = gateVal(self.AH,self.BH,self.h0,not pulseState)
        self.n = gateVal(self.AN,self.BN,self.n0,pulseState)
        self.q = gateVal(self.AQ,self.BQ,self.q0,pulseState)
    
        Iion = self.gNa * self.m**3 * self.h*(V[1]-self.ENa) +\
               self.gKf * self.n**4 * (V[1]-self.EK) +\
               self.gKs * self.q**2 * (V[1]-self.EK)
        
        dVdt = np.array([
            (-self.Isyn_d-self.gld*(V[0]-self.El)-self.gc*(V[0]-V[1])+self.Iinj_d(t))/self.Cd,
            (-self.Isyn_s-self.gls*(V[1]-self.El)-self.gc*(V[1]-V[0])-Iion+self.Iinj_s(t))/self.Cs
        ])
    
        if (dVdt[0]>1e5):
            print(f"{dVdt[0]}, {Iion}")
            print(f"{self.m},{self.h},{self.n},{self.q}")
            print(f"{pulseState}")
            print("")
        return dVdt
        
    def rk4Step(self):
        k1 = self.dt * self.ddt(1,self.t[-1], self.V[-1])
        k2 = self.dt * self.ddt(2,self.t[-1] + 0.5 * self.dt, self.V[-1] + 0.5 * k1)
        k3 = self.dt * self.ddt(3,self.t[-1] + 0.5 * self.dt, self.V[-1] + 0.5 * k2)
        k4 = self.dt * self.ddt(4,self.t[-1] + self.dt, self.V[-1] + k3)
    
        V = self.V[-1] + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        
        t = self.t[-1] + self.dt

        self.V.append(V)
        self.t.append(t)

    
