#! usr/bin/env/python3

import sys
import numpy as np
import time
import matplotlib.pyplot as plt

from Cellular_Lattice import Cellular_Lattice
from SIRS_Lattice import SIRS_Lattice

def bootstrap_error(full_data, statistic):
    statistics = []
    for sample_round in range(100):
        sample = []
        for i in range(len(full_data)):
            sample.append(full_data[np.random.randint(0,len(full_data))])
        statistics.append(statistic(sample))
    return(np.var(statistics)**0.5)

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
        try:
            p1, p2, p3 = sys.argv[6], sys.argv[7], sys.argv[8]
        except IndexError:
            p1, p2, p3 = 0.8, 0.1, 0.01
        simulation = SIRS_Lattice(size=(n,m), mode=mode,dynamic=dynamic,
                                  animate=True, max_iter=1000,p1=p1, p2=p2,
                                  p3=p3)
        simulation.run(max_iter=1000)

if sys.argv[5] == "speed":

    if dynamic == "conway":
        simulation = Cellular_Lattice(size=(n,m), mode=mode)
        simulation.run(dynamic=dynamic, animate=False, max_iter=1000)

if sys.argv[5] == "histogram":

    equilibrium_steps = []
    for i in range(100):
        print("\nRunning simulation {} of {}".format(i, 100))
        simulation = Cellular_Lattice(size=(n,m), mode=mode)
        equilibrium_steps.append(simulation.run(dynamic=dynamic, animate=False, max_iter=10000))
    # Write the steps to equilibrium for each simulation out to a file.
    np.savetxt("data/hist/sweeps_to_equ.csv", np.array(equilibrium_steps), delimiter=" ")
    # Plot a bin histogram of the steps to equilibrium and write it to a file.
    plt.hist(equilibrium_steps, bins=15)
    plt.xlabel("Sweeps to Equilibrium")
    plt.ylabel("Frequency")
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
        for p3 in p3s:
            print("Running simulation for p1={} and p3={}".format(p1, p3), end="\r"),
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
            np.savetxt("data/phase/psis_p1={}_p3={}.csv".format(p1,p3), psis, delimiter=" ")

    # Write out the matrices.
    np.savetxt("data/phase/phase_matrix.csv", phase_matrix, delimiter=" ")
    np.savetxt("data/phase/phase_variance_matrix.csv", var_matrix, delimiter=" ")
    # Plot the phase diagram.
    plt.imshow(phase_matrix.T, extent=[p1s[0],p1s[-1],p3s[0],p3s[-1]], cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Average Infected Fraction in p1-p3 Phase Space")
    plt.savefig("phase_diagram.png")
    plt.clf()
    # Plot the variance in phase space.
    plt.imshow(var_matrix.T, extent=[p1s[0],p1s[-1],p3s[0],p3s[-1]], cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Variance of Infected Sites in p1-p3 Phase Space")
    plt.savefig("variance_diagram.png")
    plt.clf()

if sys.argv[5] == "slice":

    # Initialise probability domains.
    p1s = np.arange(0.2, 0.5, 0.01)
    # Initialise phase space matrrices.
    variance = np.zeros(p1s.size)
    variance_err = np.zeros(p1s.size)

    p2, p3 = 0.5, 0.5 # Fix the value of p2, p3 for all simulations.

    for p1 in p1s:
        print("Running simulation for p1={}.".format(p1), end="\r"),
        simulation = SIRS_Lattice(size=(n,m), mode=mode, dynamic=dynamic, animate=False, p1=p1, p2=p2, p3=p3)
        psis = []
        for sweep in range(10000):
            simulation.sweep()
            if sweep > 99:
                psis.append(simulation.get_infected())
        # Dirty data-type change so math can be performed.
        psis = np.array(psis)
        # Get the phase-space location index.
        p1_index = np.where(p1s==p1)
        # Store the variance of the average infected fraction and errors.
        variance[p1_index] = np.var(psis)/(n*m)
        variance_err[p1_index] = bootstrap_error(psis, np.var)/(n*m)

        # Write out the psis for the particular value of p1.
        np.savetxt("data/slice/psis_p1={}.csv".format(p1), psis, delimiter=" ")
        # Plot slice in phase space with errors.
        plt.errorbar(p1s, variance, yerr=variance_err, elinewidth=1, capsize=3, barsabove=True)
        plt.xlabel("p1 (S -> I)")
        plt.ylabel("Variance Infected Fraction")
        plt.title("Variance Infected Fraction in p1 for p2,p3=0.5")
        plt.savefig("plots/slice_diagram.png")
        plt.clf()

    # Write out the variance of the infected fraction and error on the variance.
    variance_data = np.stack((p1s, variance, variance_err), axis=-1)
    np.savetxt("data/slice/variance.csv", variance_data, delimiter=" ")

if sys.argv[5] == "immunity":

    # Initialise immunity domain.
    immune_fracs = np.arange(0.0, 1.0, 0.01)
    # Fix probability domain.
    p1, p2, p3 = 0.5, 0.5, 0.5
    # Arrays for storing infected fraction and errors.
    inf_frac = np.zeros(immune_fracs.shape)
    inf_frac_err = np.zeros(immune_fracs.shape)

    for immune_frac in immune_fracs:
        # Get the infected fraction and an error for a value of immunity fraction.
        avg_inf_fracs = []
        for sim in range(5):
            print("Running simulation {} of 5 for frac={}.".format(sim, immune_frac), end="\r"),
            # Initialise the simulation object.
            simulation = SIRS_Lattice(size=(n,m), mode=mode, dynamic=dynamic, animate=False, p1=p1, p2=p2, p3=p3)
            # Replace random sites with immune members.
            for i in range(n):
                for j in range(m):
                    if np.random.rand() < immune_frac:
                        simulation.lattice[i,j] = 2
            # Run the simulation.
            psis = []
            for sweep in range(1000):
                simulation.sweep()
                if sweep > 99:
                    psis.append(simulation.get_infected())
            # Collect the average infected fraction for this run.
            avg_inf_fracs.append(np.mean(psis)/(n*m))
            # Write out the arrays.
            np.savetxt("data/immunity/psis_frac={}_sim{}.csv".format(immune_frac,sim), psis, delimiter=" ")
        # Get a mean and standard error for plotting.
        index = np.where(immune_fracs == immune_frac)
        inf_frac[index] = np.mean(avg_inf_fracs)
        inf_frac_err[index] = np.std(avg_inf_fracs) / len(avg_inf_fracs)**0.5

    # Write out the infected fracion vs. immune fraction data with error.
    immune_frac_data = np.stack((immune_fracs, inf_frac, inf_frac_err), axis=-1)
    np.savetxt("data/immunity/inf_frac_immune_frac.csv", immune_frac_data, delimiter=" ")
    # Plot slice in phase space with errors.
    plt.errorbar(immune_fracs, inf_frac, yerr=inf_frac_err, elinewidth=1, capsize=3, barsabove=True)
    plt.xlabel("Immunity Fraction")
    plt.ylabel("Average Infected Fraction")
    plt.title("Average Infected Fraction vs. Immunity Fraction")
    plt.savefig("plots/immunity_diagram.png")
    plt.clf()
    # Write out the infected fraction and errors.
    immune_frac_data = np.stack((immune_fracs, inf_frac, inf_frac_err), axis=-1)
    np.savetxt("data/immunity/inf_frac_immune_frac.csv", immune_frac_data, delimiter=" ")

toc = time.clock()
print("Executed script in {} seconds.".format(toc-tic))
