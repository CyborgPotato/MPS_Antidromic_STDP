import PBLIF
import importlib
importlib.reload(PBLIF)
from matplotlib.pyplot import *

p1 = PBLIF.PBLIF()
p2 = PBLIF.PBLIF()

p1.Iinj_d = lambda t: 0
p1.connect(p2)

while (p1.t[-1]<=p1.tstop):
    p1.rk4Step()
    p2.rk4Step()

print("\nDone!!")

plot(p1.t,[item[1] for item in p1.V],color='#0000ff')
# plot(p2.t,[item[1] for item in p2.V],color='#ff0000')
show()
