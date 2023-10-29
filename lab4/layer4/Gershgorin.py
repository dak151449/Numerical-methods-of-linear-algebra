import matplotlib.pyplot as plt
import numpy as np


class Gershgorin:

    def __init__(self, matrix):
        self.matrix = matrix

    def getCircles(self):
        circles = []
        for i in range(len(self.matrix)):
            a = self.matrix[i][i]
            r = sum(map(abs, self.matrix[i]))
            circles.append([a, r])
        return circles

    def ShowCrcles(self, circles, lambds):
        _, axes = plt.subplots()
        print(lambds)
        for i in lambds:
            plt.plot(i, 0, 'ro')
        l, r = 0, 0

        for i in circles:
            print(i)
            d = plt.Circle((i[0], 0), i[1], fill=False)
            radius = i[1]
            l = min(l, i[0] - 2*radius)
            r = max(r, i[0] + 2*radius)
            axes.add_patch(d)

        axes.set_aspect("equal", "box")
        plt.axis('scaled')
        x = np.linspace(l, r, 2000)
        axes.plot(x, 0*x)
        plt.title("Gershgorin circles")
        plt.show()
