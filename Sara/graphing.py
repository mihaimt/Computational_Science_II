import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

import os

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

def plot2(a,b):
	plt.subplot(1,2,1)
	plotFile(a)
	plt.subplot(1,2,2)
	plotFile(b)

def plotDiff2(a1, a2, b1, b2):
	plt.subplot(1,2,1)
	plotDiff(a1, a2)
	plt.subplot(1,2,2)
	plotDiff(b1, b2)

def plot4(a,b,c,d):
	plt.subplot(2,2,1)
	plotFile(a)
	plt.subplot(2,2,2)
	plotFile(b)
	plt.subplot(2,2,3)
	plotFile(c)
	plt.subplot(2,2,4)
	plotFile(d)

def plotFile(filename): 
    print "plotFile", filename
    f = open(filename, 'r')
    values = [float(line) for line in f]
    f.close()
    plotArray(values, os.path.basename(filename))

def plotDiff(file_a, file_b):
    print "plotDiff", file_a, file_b
    f = open(file_a, 'r')
    values_a = [float(line) for line in f]
    f.close()
    f = open(file_b, 'r')
    values_b = [float(line) for line in f]
    f.close()
    diff = [b_i - a_i for a_i, b_i in zip(values_a, values_b)]
    plotArray(diff, 'diff_' + os.path.basename(file_a))

def plotArray(values, name):
    x, y = np.mgrid[ slice(min(radii),max(radii)+dr,dr), 
              slice(min(angles), max(angles)+dtheta, dtheta) ]

    value_min = min(values)
    value_max = max(values)

    npvalues = np.array(values).reshape(len(radii),len(angles))

    #plt.pcolor(x, y, npvalues, cmap=cmap_bgr, vmin=value_min, vmax=value_max)
    plt.pcolor(x, y, npvalues, cmap=plt.get_cmap('jet'), vmin=value_min, vmax=value_max)
    plt.title(name + '\n[' + str(value_min) + ', ' + str(value_max) + ']')
    plt.xlabel("Radius [r_0]")
    plt.ylabel("Azimuth [theta]")
    plt.axis([x.min(), x.max(), y.min(), y.max()]);

def plotMass(level):
    f = open('mass_%.2d.txt' % level, 'r')
    masses = [float(line) for line in f]
    f.close()
    size = int(math.sqrt(len(masses)/2))
    x, y = np.mgrid[slice(0, size, 1), slice(0, size*2, 1)]
    value_min = min(masses)
    value_max = max(masses)
    print "plotMass", len(masses), size, value_min, value_max

    npvalues = np.array(masses).reshape(size, size*2)
    plt.pcolor(x, y, npvalues, cmap=cmap_bgr, vmin=value_min, vmax=value_max)
    plt.title('Mass Level %d' % level)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis([x.min(), x.max(), y.min(), y.max()]);

f = open('r_project.data', 'r')
radii = [float(line) for line in f]
f.close()
dr = radii[1]-radii[0]

f = open('theta_project.data', 'r')
angles = [float(line) for line in f]
f.close()
dtheta = angles[1]-angles[0]

# Show all local
#plot4('density_project.data', 'force_mag.data', 'force_r.data', 'force_theta.data')

# Show local
#plot2('force_r.data', 'force_theta.data')

path = '../dephil/data/'
suffix = '_pure_lvl0.data'
suffix2 = '_pure_lvl1.data'
# Show suffix
#plot2(path + 'radial' + suffix, path + 'angular' + suffix)

# Diff local vs local base
plotDiff2('force_r.data', 'data/force_r.data', 'force_theta.data', 'data/force_theta.data')

# Diff local vs dephil base
#plotDiff2('data/force_r.data', path + 'f_radial.data', 'data/force_theta.data', path + 'f_angular.data')

# Diff local vs. suffix
#plotDiff2('force_r.data', path + 'radial' + suffix, 'force_theta.data', path + 'angular' + suffix)

# Diff local vs dephil default
#plotDiff2('force_r.data', path + 'f_radial.data', 'force_theta.data', path + 'f_angular.data')

# Diff suffix vs suffix2
#plotDiff2(path + 'radial' + suffix, path + 'radial' + suffix2, path + 'angular' + suffix, path + 'angular' + suffix2)

plt.show()

#plt.savefig("forces.png")
#print "saved forces.png"

#for i in range(0,4):
#    plt.subplot(2,2,i+1)
#    plotMass(i)
#plt.show()
#plt.savefig("masses.png")
#print "saved masses.png"
