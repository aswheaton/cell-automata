import numpy as np

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
        print("Loading simulation {} of 5 for frac={}.".format(sim, immune_frac), end="\r"),
        psis = np.loadtxt("data/immunity/tot_inf_frac={}_sim{}.csv".format(immune_frac,sim), delimiter=" ")
        # Collect the average infected fraction for this run.
        avg_inf_fracs.append(np.mean(psis)/(n*m))

    # Get a mean and standard error for plotting.
    index = np.where(immune_fracs == immune_frac)
    inf_frac[index] = np.mean(avg_inf_fracs)
    inf_frac_err[index] = np.std(avg_inf_fracs) / len(avg_inf_fracs)**0.5

# Write out the infected fraction and errors.
immune_frac_data = np.stack((immune_fracs, inf_frac, inf_frac_err), axis=-1)
np.savetxt("data/immunity/immune_frac.csv", immune_frac_data, delimiter=" ")
