#! usr/bin/env/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Cellular_Lattice(object):

    def __init__(self, **kwargs):
        self.size = kwargs.get("size")
        self.mode = kwargs.get("mode")
        self.build()

    def build(self):
        """
        Creates a new class attribute of type ndarray and fills it with cell
        values, depending on the user specified mode.
        """
        if self.mode == "random":
            self.lattice = np.random.choice(a=[0,1], size=self.size)
        if self.mode == "glider":
            self.lattice = np.zeros(self.size, dtype=int)
            self.lattice[25:28,25:28] = np.array([[1,1,1],
                                                  [1,0,0],
                                                  [0,1,0]])
        if self.mode == "beehive":
            self.lattice = np.zeros(self.size, dtype=int)
            self.lattice[25:29,25:28] = np.array([[0,1,0],
                                                [1,0,1],
                                                [1,0,1],
                                                [0,1,0]])
        if self.mode == "blinker":
            self.lattice = np.zeros(self.size, dtype=int)
            self.lattice[25:28,25:28] = np.array([[0,1,0],
                                              [0,1,0],
                                              [0,1,0]])

    def bc(self, indices):
        """
        Determines if a pair of indices falls outside the boundary of the
        lattice and if so, applies a periodic (toroidal) boundary condition
        to return new indices.
        """
        return((indices[0]%self.size[0], indices[1]%self.size[1]))

    def get_neighbours(self, indices):
        """

        """
        n, m = indices
        neighbours = 0
        for i in [n-1, n, n+1]:
            for j in [m-1, m, m+1]:
                if self.lattice[self.bc((i,j))] == 1:
                    neighbours += 1
        if self.lattice[n,m] == 1:
            neighbours -= 1
        return(neighbours)

    def gen_next_lattice(self):

        new_lattice = np.zeros(self.size, dtype=int)

        if self.dynamic == "conway":
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    neighbours = self.get_neighbours((i,j))
                    # Condition for currently dead cells.
                    if self.lattice[i,j] == 0:
                        if neighbours == 3:
                            new_lattice[i,j] = 1
                        else:
                            new_lattice[i,j] = 0
                    # Condition for currently live cells.
                    elif self.lattice[i,j] == 1:
                        if neighbours == 2 or neighbours == 3:
                            new_lattice[i,j] = 1
                        else:
                            new_lattice[i,j] = 0

        return(new_lattice)

        if self.dynamic == "SIRS":
            print("No SIRS code implemented yet!")
            pass

    def weighted_mean_2D(self):
        x_sum, y_sum = np.sum(self.lattice, axis=0), np.sum(self.lattice, axis=1)
        x_avg = np.average(range(len(x_sum)), weights=x_sum)
        y_avg = np.average(range(len(y_sum)), weights=y_sum)
        return((x_avg, y_avg))

    def get_displacement(self, indices_a, indices_b):
        displacement = ((indices_b[0]-indices_a[0])**2 +
                        (indices_b[1]-indices_a[1])**2
                        )**0.5
        return(displacement)

    def remove_outliers(self, data):
        """
        Removes data falling outside two standard deviation of a dataset.
        """
        return data[abs(data - np.mean(data)) < 2 * np.std(data)]

    def step(self, *args):
        """
        Steps the simulation forward by attempting 1000 spin flips.
        Takes *args for call by animation.FuncAnimation instance.
        # TODO: Make make number of attempted spin flips configurable!
        # TODO: Determine the purpose of the trailing comma in return().
        """

        self.lattice = self.gen_next_lattice()

        if self.animate == True:
            self.image.set_array(self.lattice)
            return(self.image,)

    def run(self, **kwargs):
        """
        Sets up a figure, image, and FuncAnimation instance, then runs the
        simulation to the specified maximum number of iterations.
        # TODO: Utilise the Boolean animate attribute here, and implement
        datafile output.
        # TODO: Make the number of Metropolis trials more understandable.
        (Currently number of attempted flips is more than those specified by
        the user by a factor of 10^3 due to the way animate() method works.)
        """
        self.dynamic = kwargs.get("dynamic")
        self.max_iter = kwargs.get("max_iter")
        self.animate = kwargs.get("animate")

        if self.animate == True:
            self.figure = plt.figure()
            self.image = plt.imshow(self.lattice, animated=True)
            self.animation = animation.FuncAnimation(self.figure, self.step,
                                                    frames=self.max_iter,
                                                    repeat=False,
                                                    interval=50, blit=False
                                                    )
            plt.show()

        elif self.animate == False:

            self.com = np.zeros((self.max_iter, 2))
            self.disp = np.zeros(self.max_iter)
            self.live_cells = np.zeros(self.max_iter)

            self.com[0,:] = self.weighted_mean_2D()

            for step in range(self.max_iter):

                print("Step {} of {}".format(step, self.max_iter), end="\r"),
                self.step()

                self.com[step,:] = self.weighted_mean_2D()
                self.disp[step] = self.get_displacement(self.com[step-1,:], self.com[step,:])
                self.live_cells[step] = np.sum(self.lattice, dtype=int)

                if self.mode == "random":
                    if (self.live_cells[step] == self.live_cells[step-1]) and (self.live_cells[step-1] == self.live_cells[step-2]):
                        print("\nEquilibrium reached at step {}!".format(step))
                        break

            if self.mode == "glider":
                print()
                print("Max displacement: {}".format(np.amax(self.disp)))
                print("Mean displacement: {}".format(np.mean(self.remove_outliers(self.disp[:step]))))

            return(step)

    def exportAnimation(self, filename, dotsPerInch):
        """
        Exports the animation to a .gif file without compression. (Linux
        distributions with package "imagemagick" only. Files can be large!)
        # TODO: rename this for PEP8 compliance and add support for other
        image writing packages.
        """
        self.animation.save(filename, dpi=dotsPerInch, writer="imagemagick")
