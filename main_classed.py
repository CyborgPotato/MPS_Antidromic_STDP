import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p = PBLIF.PBLIF()
p.tstop=6000
while (p.t[-1]<=p.tstop):
    p.rk4Step()

print("\nDone!!")
# plot(p.V)
subplot(231)
plot(p.t,p.m,label='m')
subplot(232)
plot(p.t,p.h,label='h')
subplot(234)
plot(p.t,p.n,label='n')
subplot(235)
plot(p.t,p.q,label='q')
subplot(233)
plot(p.t,p.V)
show()
