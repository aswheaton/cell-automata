import sys
import numpy as np
import time

from Cellular_Lattice import Cellular_Lattice

dynamic = str(sys.argv[1])
mode = str(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])

tic = time.clock()

try:

    if sys.argv[5] == "animate":
        simulation = Cellular_Lattice(size=(n,m), mode=mode)
        simulation.run(dynamic=dynamic, animate=True, max_iter=1000)

except IndexError:
    print("No animation argument!")
    pass

toc = time.clock()
print("Executed script in "+str(toc-tic)+" seconds.")
