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
            self.lattice[0:2,0:2] = np.array([[1,1,1],
                                              [1,0,0],
                                              [0,1,0]])
        # Create empty lattice for storing next iteration.
        self.next_lattice = np.zeros(self.size)
    def bc(self, indices):
        """
            Determines if a pair of indices falls outside the boundary of the
            lattice and if so, applies a periodic (toroidal) boundary condition
            to return new indices.
        """
        return((indices[0]%self.size[0], indices[1]%self.size[1]))

    def get_neighbours(self, indices):
        """
            Calculates and returns the change in energy from flipping a given
            of a given lattice site (n,m) due to spin interaction with its
            neighbours. Change in energy given by:

            E = - S_n,m * (S_n+1,m + Sn-1,m + S_n,m-1, + S_n,m+1)

            Four other lattice points enter this expression.
            # TODO: Verify this calculation!
        """
        n, m = indices
        neighbours = -self.lattice[n,m]
        for i in [n-1, n, n+1]:
            for j in [m-1, m, m+1]:
                neighbours += self.lattice[self.bc((i,j))]
        return(neighbours)

    def gen_next_lattice(self):

        if self.dynamic == "conway":
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    neighbours = self.get_neighbours((i,j))
                    if self.lattice[i,j] == 0:
                        if neighbours == 3:
                            self.next_lattice[i,j] = 1
                        else:
                            self.next_lattice[i,j] = 0
                    elif self.lattice[i,j] == 1:
                        if neighbours == 2 or neighbours == 3:
                            self.next_lattice[i,j] = 1
                        else:
                            self.next_lattice[i,j] = 0

        if self.dynamic == "sirs":
            print("No SIRS code implemented yet!")
            pass

    def sweep(self, *args):
        """
            Steps the simulation forward by attempting 1000 spin flips.
            Takes *args for call by animation.FuncAnimation instance.
            # TODO: Make make number of attempted spin flips configurable!
            # TODO: Determine the purpose of the trailing comma in return().
        """
        self.gen_next_lattice()
        self.lattice = self.next_lattice

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
        if kwargs.get("animate") == True:
            self.figure = plt.figure()
            self.image = plt.imshow(self.lattice, animated=True)
            # TODO: Make line wrapping PEP8 compliant, here.
            self.animation = animation.FuncAnimation(self.figure, self.sweep,
                                                    frames=self.max_iter,
                                                    repeat=False,
                                                    interval=1, blit=True
                                                    )
            plt.show()

        elif kwargs.get("animate") == False:
            f = open("dat/"+self.dynamic+"_"+str(self.temp)+".csv","w+")
            for sweep in range(self.max_iter):
                print("Sweep "+str(sweep)+" of "+str(self.max_iter)+" for T="+str(self.temp)+".", end="\r"),
                self.sweep()
                if sweep > 99 and sweep % 10 == 0:
                    f.write(str(self.total_energy())+", "+str(self.magnetization())+"\n")
            print("")
            f.close()


    def exportAnimation(self, filename, dotsPerInch):
        """
            Exports the animation to a .gif file without compression. (Linux
            distributions with package "imagemagick" only. Files can be large!)
            # TODO: rename this for PEP8 compliance and add support for other
            image writing packages.
        """
        self.animation.save(filename, dpi=dotsPerInch, writer="imagemagick")
