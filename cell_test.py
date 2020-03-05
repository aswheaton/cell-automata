import sys
import numpy as np
import time
import matplotlib.pyplot as plt

from Cellular_Lattice import Cellular_Lattice
from SIRS_Lattice import SIRS_Lattice

dynamic = str(sys.argv[1])
mode = str(sys.argv[2])
n = int(sys.argv[3])
m = int(sys.argv[4])

tic = time.clock()

try:

    if sys.argv[5] == "animate":

        if dynamic == "conway":
            simulation = Cellular_Lattice(size=(n,m), mode=mode)
            simulation.run(dynamic=dynamic, animate=True, max_iter=1000)

        if dynamic == "SIRS":
            p1, p2, p3 = 0.8, 0.1, 0.01
            simulation = SIRS_Lattice(size=(n,m), mode=mode)
            simulation.run(dynamic=dynamic, animate=True, max_iter=1000,
                           p1=p1, p2=p2, p3=p3)

except IndexError:

    print("No animation argument, writing data to files!")

    if dynamic == "conway":

        equilibrium_steps = []
        for i in range(100):
            simulation = Cellular_Lattice(size=(n,m), mode=mode)
            equilibrium_steps.append(simulation.run(dynamic=dynamic, animate=False, max_iter=1000))
        plt.hist(equilibrium_steps)
        plt.savefig("plots/equilibrium_hist.png")

    elif dynamic == "SIRS":
        p1, p2, p3 = 0.8, 0.1, 0.01
        simulation = SIRS_Lattice(size=(n,m), mode=mode)
        simulation.run(dynamic=dynamic, animate=False, max_iter=1000,
                       p1=p1, p2=p2, p3=p3)

toc = time.clock()
print("Executed script in {} seconds.".format(toc-tic))
