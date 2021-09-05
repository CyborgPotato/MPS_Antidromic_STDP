import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p1 = PBLIF.PBLIF()
p1.record=True

p1.tstop = 2000

p1.Iinj_d = lambda t: 0#(t+500)/100000000

while (p1.t[-1]<=p1.tstop):
    p1.rk4Step()

print("\nDone!!")

plot(p1.t,[item[1] for item in p1.V],color='#0000ff')
show()
