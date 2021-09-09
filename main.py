import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p1 = PBLIF.PBLIF()
p1.record=True

p1.Iinj_d = lambda t: 100/1000000*(t>=250)*(t<=251)
p1.Iinj_s = lambda t: 10/1000000*(t>=150)*(t<=151)
p1.Iinj_a = lambda t: 50*(t>=50)*(t<=51)

p1.tstop = 2000

while (p1.t[-1]<=p1.tstop):
    p1.rk4Step()

print("\nDone!!")

subplot(211)
plot(p1.t,[item[1] for item in p1.V],color='#0000ff')
lim = ylim()
xlabel("Time (ms)")
ylabel("Membrane Voltage (mV) Soma")

subplot(212)
plot(p1.t,[item[0] for item in p1.V],color='#0000ff')
ylim(lim)
xlabel("Time (ms)")
ylabel("Membrane Voltage (mV) Dendrite")

suptitle("Membrane Voltage of Compartments")
tight_layout()
show()

subplot(221)
plot(p1.t,p1.m)
xlabel("Time (ms)")
ylabel("Gate Parameter m")
subplot(222)
plot(p1.t,p1.h)
xlabel("Time (ms)")
ylabel("Gate Parameter h")
subplot(223)
plot(p1.t,p1.n)
xlabel("Time (ms)")
ylabel("Gate Parameter n")
subplot(224)
plot(p1.t,p1.q)
xlabel("Time (ms)")
ylabel("Gate Parameter q")
suptitle("Time evolution of gate parameters")
tight_layout()
show()

subplot(131)
plot(p1.t,p1.INa)
xlabel("Time (ms)")
ylabel("Current Na (mA)")
subplot(132)
plot(p1.t,p1.IKf)
xlabel("Time (ms)")
ylabel("Current Kf (mA)")
subplot(133)
plot(p1.t,p1.IKs)
xlabel("Time (ms)")
ylabel("Current Ks (mA)")

suptitle("Somatic Ion Currents")
tight_layout()
show()
