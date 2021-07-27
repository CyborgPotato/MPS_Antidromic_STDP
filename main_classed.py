from PBLIF import PBLIF
from matplotlib.pyplot import *

p = PBLIF()

while (p.t[-1]<=p.tstop):
    p.rk4Step()

print("Done!!")
plot(p.V)
show()
