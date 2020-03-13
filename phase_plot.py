import numpy as np
import matplotlib.pyplot as plt

def main():
    n, m = 50, 50
    # Initialise probability domains.
    p1s = np.arange(0.0, 1.0, 0.025)
    p3s = np.arange(0.0, 1.0, 0.025)
    # Initialise phase space matrrices.
    phase_matrix = np.zeros((p1s.size, p3s.size))
    var_matrix = np.zeros((p1s.size, p3s.size))

    for p1 in p1s:
        for p3 in p3s:
            print("Loading data for p1={} and p3={}".format(p1, p3), end="\r"),
            try:
                psis = np.loadtxt("data/phase/psis_p1={}_p3={}.csv".format(p1,p3), delimiter=" ")
            except OSError:
                break
            # Get the phase-space location index.
            p1_index = np.where(p1s==p1)
            p3_index = np.where(p3s==p3)
            # Store the average infected fraction.
            phase_matrix[p1_index,p3_index] = np.mean(psis)/(n*m)
            var_matrix[p1_index,p3_index] = np.var(psis)/(n*m)
    print()

    # Plot the phase diagram.
    plt.imshow(phase_matrix.T, extent=[p1s[0],p1s[-1],p3s[0],p3s[-1]], cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Average Infected Fraction in p1-p3 Phase Space")
    plt.savefig("plots/phase_diagram.png")
    plt.show()
    plt.clf()
    # Plot the variance in phase space.
    plt.imshow(var_matrix.T, extent=[p1s[0],p1s[-1],p3s[0],p3s[-1]], cmap='viridis')
    plt.colorbar()
    plt.xlabel("p1 (S -> I)")
    plt.ylabel("p3 (R -> S)")
    plt.title("Variance of Infected Sites in p1-p3 Phase Space")
    plt.savefig("plots/variance_diagram.png")
    plt.show()
    plt.clf()
main()
