import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

#custom colour map, default looks funny
cmap_bgr = colors.LinearSegmentedColormap.from_list("bgr", [
    (0.0, "darkred"),
    (0.2, "red"),
    (0.3, "orange"),
    (0.4, "yellow"),
    (0.5, "lightgreen"),
    (0.6, "cyan"),
    (0.7, "blue"),
    (1.0, "darkblue"),
])


def plotFile( filename ): 
    f = open(filename, 'r')
    values = [float(line) for line in f]
    f.close()

    x, y = np.mgrid[ slice(min(radii),max(radii)+dr,dr), 
              slice(min(angles), max(angles)+dtheta, dtheta) ]

    value_min = min(values)
    value_max = max(values)

    npvalues = np.array(values).reshape(len(radii),len(angles))

    plt.pcolor(x, y, npvalues, cmap=cmap_bgr, vmin=value_min, vmax=value_max)
    plt.title(filename)
    plt.xlabel("Radius [r_0]")
    plt.ylabel("Azimuth [theta]")
    plt.axis([x.min(), x.max(), y.min(), y.max()]);

f = open('r_project.data', 'r')
radii = [float(line) for line in f]
f.close()
dr = radii[1]-radii[0]

f = open('theta_project.data', 'r')
angles = [float(line) for line in f]
f.close()
dtheta = angles[1]-angles[0]

plt.subplot(2,2,1)
plotFile ('density_project.data')
plt.subplot(2,2,2)
plotFile ('force_mag.data')
plt.subplot(2,2,3)
plotFile ('force_r.data')
plt.subplot(2,2,4)
plotFile ('force_theta.data')

plt.show()
