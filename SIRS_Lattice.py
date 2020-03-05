#! usr/bin/env/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class SIRS_Lattice(object):

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
            self.lattice = np.random.choice(a=[-1,0,1], size=self.size)
        # Create empty lattice for storing next iteration.
        self.next_lattice = np.zeros(self.size, dtype=int)

    def bc(self, indices):
        """
            Determines if a pair of indices falls outside the boundary of the
            lattice and if so, applies a periodic (toroidal) boundary condition
            to return new indices.
        """
        return((indices[0]%self.size[0], indices[1]%self.size[1]))

    def get_neighbours(self, indices, **kwargs):
        """

        """
        n, m = indices
        id = kwargs.get("id")
        neighbours = 0
        for i in [n-1, n, n+1]:
            for j in [m-1, m, m+1]:
                if self.lattice[self.bc((i,j))] == id:
                    neighbours += 1
        if self.lattice[n,m] == id:
            neighbours -= 1
        return(neighbours)

    def gen_next_lattice(self):

        if self.dynamic == "conway":
            new_lattice = np.zeros(self.size, dtype=int)
            print("No Game of Life code implemented yet!")
            pass

        if self.dynamic == "SIRS":
            # Create new lattice for next sweep.
            new_lattice = self.lattice
            # Attempt to evolve N^2 sites on the lattice.
            for step in range(self.size[0] * self.size[1]):
                # Select a random site on the lattice.
                i = np.random.randint(0, self.size[0])
                j = np.random.randint(0, self.size[1])
                # Condition for susceptible state.
                if self.lattice[i,j] == -1:
                    if self.get_neighbours((i,j), id=0):
                        if np.random.rand() < self.p1:
                            new_lattice[i,j] = 0
                        else:
                            new_lattice[i,j] = -1
                # Condition for infected state.
                if self.lattice[i,j] == 0:
                    if np.random.rand() < self.p2:
                        new_lattice[i,j] = 1
                    else:
                        new_lattice[i,j] = 0
                # Condition for recovered state.
                if self.lattice[i,j] == 1:
                    if np.random.rand() < self.p3:
                        new_lattice[i,j] = -1
                    else:
                        new_lattice[i,j] = 1

        return(new_lattice)

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
        self.p1, self.p2, self.p3 = kwargs.get("p1"), kwargs.get("p2"), kwargs.get("p3")

        if kwargs.get("animate") == True:
            self.figure = plt.figure()
            self.image = plt.imshow(self.lattice, animated=True)
            self.animation = animation.FuncAnimation(self.figure, self.step,
                                                    frames=self.max_iter,
                                                    repeat=False,
                                                    interval=100, blit=False
                                                    )
            plt.show()

        elif kwargs.get("animate") == False:

            for step in range(self.max_iter):
                print("Step {} of {}".format(step, self.max_iter), end="\r"),
                self.step()

    def exportAnimation(self, filename, dotsPerInch):
        """
            Exports the animation to a .gif file without compression. (Linux
            distributions with package "imagemagick" only. Files can be large!)
            # TODO: rename this for PEP8 compliance and add support for other
            image writing packages.
        """
        self.animation.save(filename, dpi=dotsPerInch, writer="imagemagick")