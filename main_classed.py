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
subplot(331)
plot(p.t,p.m,label='m')
subplot(332)
plot(p.t,p.h,label='h')
subplot(334)
plot(p.t,p.n,label='n')
subplot(335)
plot(p.t,p.q,label='q')
subplot(333)
plot(p.t,[item[0] for item in p.V])
subplot(336)
plot(p.t,[item[1] for item in p.V])
subplot(337)
plot(p.t,p.INa)
subplot(338)
plot(p.t,p.IKf)
subplot(339)
plot(p.t,p.IKs)
# ylim(-10,100)
show()
