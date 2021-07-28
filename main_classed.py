import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p = PBLIF.PBLIF()
p.tstop=3000

while (p.t[-1]<=p.tstop):
    p.rk4Step()

print("\nDone!!")
# plot(p.V)
subplot(331)
title('m gate')
plot(p.t,p.m,label='m',color='#ff0000')
subplot(332)
title('h gate')
plot(p.t,p.h,label='h',color='#ff7800')
subplot(334)
title('n gate')
plot(p.t,p.n,label='n',color='#00ff28')
subplot(335)
title('q gate')
plot(p.t,p.q,label='q',color='#ff0028')
subplot(333)
title('Dendrite V')
plot(p.t,[item[0] for item in p.V],color='#0000ff')
subplot(336)
title('Soma V')
plot(p.t,[item[1] for item in p.V],color='#000000')
subplot(337)
title('Sodium Current')
plot(p.t,p.INa,color='#ff00ff')
subplot(338)
title('Fast Potassium')
plot(p.t,p.IKf,color='#780078')
subplot(339)
title('Slow Potassium')
plot(p.t,p.IKs,color='#783800')
# ylim(-10,100)
show()
