import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p1 = PBLIF.PBLIF()
p1.record=True
p1.recordCon=True
p2 = PBLIF.PBLIF()
p2.record=True

p1.tstop = 4000
p2.tstop = 4000

p2.Iinj_d = lambda t: 10/1000000
p1.Iinj_d = lambda t: 0#(t+500)/100000000
p1.gSyn=1/10000000
p1.connect(p2)

while (p1.t[-1]<=p1.tstop):
    p1.rk4Step()
    p2.rk4Step()

print("\nDone!!")

subplot(311)
plot(p1.t,[item[1] for item in p1.V],color='#0000ff')
subplot(312)
plot(p2.t,[item[1] for item in p2.V],color='#ff0000')
subplot(313)
plot(p3.t,[item[1] for item in p3.V],color='#ff0000')
show()
