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

if sys.argv[5] == "animate":

    if dynamic == "conway":
        simulation = Cellular_Lattice(size=(n,m), mode=mode)
        simulation.run(dynamic=dynamic, animate=True, max_iter=1000)

    if dynamic == "SIRS":
        p1, p2, p3 = 0.8, 0.1, 0.01
        simulation = SIRS_Lattice(size=(n,m), mode=mode,dynamic=dynamic,
                                  animate=True, max_iter=1000,p1=p1, p2=p2,
                                  p3=p3)
        simulation.run(max_iter=1000)

if sys.argv[5] == "histogram"

    equilibrium_steps = []
    for i in range(100):
        simulation = Cellular_Lattice(size=(n,m), mode=mode)
        equilibrium_steps.append(simulation.run(dynamic=dynamic, animate=False, max_iter=10000))
    plt.hist(equilibrium_steps)
    plt.xlabel("Sweeps to Equilibrium")
    plt.xlabel("Frequency")
    plt.savefig("plots/equilibrium_hist.png")

if sys.argv[5] == "phase":

    # Initialise probability domains.
    p1s = np.arange(0.0, 1.0, 0.025)
    p3s = np.arange(0.0, 1.0, 0.025)
    # Initialise phase space matrrices.
    phase_matrix = np.zeros((p1s.size, p3s.size))
    var_matrix = np.zeros((p1s.size, p3s.size))

    p2 = 0.5 # Fix the value of p2 for all simulations.

    for p1 in p1s:
        for p3 in p32:
            simulation = SIRS_Lattice(size=(n,m), mode=mode, dynamic=dynamic, animate=False, p1=p1, p2=p2, p3=p3)
            psis = []
            for sweep in range(1000):
                simulation.sweep()
                if sweep > 99:
                    psis.append(simulation.get_infected())
            # Dirty data-type change so math can be performed.
            psis = np.array(psis)
            # Get the phase-space location index.
            p1_index = np.where(p1s==p1)
            p3_index = np.where(p3s==p3)
            # Store the average infected fraction.
            phase_matrix[p1_index,p3_index] = np.mean(psis)/(n*m)
            var_matrix[p1_index,p3_index] = np.var(psis)/(n*m)
            # Write out relevant information.
            np.savetext(psis, "data/phase/psis_p1={}_p3={}.csv".format(p1,p3), delimiter=" ")

    # Write out the matrices.
    np.savetext(phase_matrix, "data/phase/phase_matrix.csv", delimiter=" ")
    np.savetext(var_matrix, "data/phase/var_matrix.csv", delimiter=" ")
    # Plot the phase diagram.
    plt.imshow(phase_matrix, origin='lower', cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Average Infected Fraction in p1-p3 Phase Space")
    plt.savefig("plots/phase_diagram.png")
    plt.clf()
    # Plot the variance in phase space.
    plt.imshow(var_matrix, origin='lower', cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Variance of Infected Sites in p1-p3 Phase Space")
    plt.savefig("plots/variance_diagram.png")
    plt.clf()

if sys.argv[5] = "slice":
    pass

toc = time.clock()
print("Executed script in {} seconds.".format(toc-tic))
