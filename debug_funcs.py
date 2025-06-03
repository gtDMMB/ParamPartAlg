from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog

def nice_print_HP(HPs):
   for hp in HPs:
      da, dc, const = hp
      print(f"{da}x + {dc}y + {const} <= 0")

def nice_print_points(points, hull=True):
   if hull:
      hull = ConvexHull(points)
      points = points[hull.vertices]
   print_str = ",".join(str((float(p[0]),float(p[1]))) for p in points)
   print_str += "," + str((float(points[0][0]),float(points[0][1])))
   print(print_str)

def nice_graph(halfspaces, bp):
   fig = plt.figure()
   ax = fig.add_subplot(1, 1, 1, aspect='equal')
   xlim, ylim = (-10000, 10000), (-10000, 10000)
   ax.set_xlim((bp[0]-100,bp[0]+100))
   ax.set_ylim((bp[1]-100,bp[1]+100))
   x = np.linspace(-10000, 10000, 10000)
   print(halfspaces.shape)
   fmt = {"color": "b", "edgecolor": "b", "alpha": 1/halfspaces.shape[0]}
   for dx, dy, b in halfspaces:
      print(dx, dy, b)
      
      # fmt["hatch"] = sym
      if dy== 0:
         # if dx*xlim[0] :
         val_x = xlim[0] if (dx*xlim[0] + b) <= 9 else xlim[1]

         ax.axvline(-b/dx, color="black")
         xi = np.linspace(val_x, -b/dx)
         ax.fill_between(xi, ylim[0], ylim[1], **fmt)
      else:
         val_y = ylim[1] if (dx*xlim[0] + dy*100000 + b) <= 0 else xlim[0]
         ax.plot(x, (-b-dx*x)/dy, color="black")
         ax.fill_between(x, (-b-dx*x)/dy, val_y, **fmt)
   # x, y = zip(*hs.intersections)
   plt.show()
   # ax.plot(x, y, 'o', markersize=8)


if __name__ == "__main__":
   halfspaces = np.array([[-1, 0., 0.],
                        [0., -1., 0.],
                        [2., 1., -4.],
                        [-0.5, 1., -2.]])
   feasible_point = np.array([0.5, 0.5])
   nice_graph(halfspaces, (0,0))