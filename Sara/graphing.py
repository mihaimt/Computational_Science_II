# -*- coding: utf-8 -*-

import argparse
import math
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import os
import sys
import time

piformat = mtick.FormatStrFormatter(u'%.1fπ')

# A bunch of command line flags to define which graphs to generate.
parser = argparse.ArgumentParser(description="Various plots.")
parser.add_argument("--prefix", help="Prefix for input files", default="run/default")
parser.add_argument("--polar", action="store_true", help="Make polar plots instead.")
parser.add_argument("-a", "--all", action="store_true", help="Plot all graphs")
parser.add_argument("-b", "--base", action="store_true", help="Plot the level0 graph")
parser.add_argument("-e", "--expt", action="store_true", help="Plot the higher level graph")
parser.add_argument("-d", "--diff", action="store_true", help="Plot the diff between level 0 and higher")
parser.add_argument("-d2", "--diff2D", action="store_true", help="Plot the 2d diff between level 0 and higher")
parser.add_argument("-d1", "--diff1D", action="store_true", help="Plot the 1d diff between level 0 and higher")
parser.add_argument("-l", "--levels", action="store_true", help="Debug plot the lookup distribution for one sample")
parser.add_argument("-m", "--masses", action="store_true", help="Debug plot all the masses.")
parser.add_argument("--summary", action="store_true", help="Plot the summary of density and forces.")
parser.add_argument("-s","--save", action="store_true", help="Save plots instead of displaying them")

FLAGS = parser.parse_args(sys.argv[1:])
print FLAGS

# 2 graphs from two files, next to each other (horizontally)
def plot2(a, name_a, b, name_b):
  fig = plt.figure()
  ax1 = fig.add_axes([0.05, 0.1, 0.35, 0.8], polar=FLAGS.polar)
  ax2 = fig.add_axes([0.55, 0.1, 0.35, 0.8], polar=FLAGS.polar)
  ax1_cbar = fig.add_axes([0.41, 0.1, 0.03, 0.8])
  ax2_cbar = fig.add_axes([0.91, 0.1, 0.03, 0.8])
  plotFile(a, name_a, fig, ax1, ax1_cbar)
  plotFile(b, name_b, fig, ax2, ax2_cbar)

# 2 graphs of the difference of 4 files
def plotDiff2(a1, a2, b1, b2):
  fig = plt.figure()
  ax1 = fig.add_axes([0.05, 0.1, 0.35, 0.8], polar=FLAGS.polar)
  ax2 = fig.add_axes([0.55, 0.1, 0.35, 0.8], polar=FLAGS.polar)
  ax1_cbar = fig.add_axes([0.41, 0.1, 0.03, 0.8])
  ax2_cbar = fig.add_axes([0.91, 0.1, 0.03, 0.8])
  plotDiff(a1, a2, 'Radial Force (%)', fig, ax1, ax1_cbar)
  plotDiff(b1, b2, 'Tangential Force (%)', fig, ax2, ax2_cbar)

# the 1D max value of the differences, 4 graphs.
def plotDiff2_1D(a1, a2, b1, b2):
  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
  plotDiff_1D(a1, a2, ax1, ax2, 'radial')
  plotDiff_1D(b1, b2, ax3, ax4, 'tangential')


# Read a file, plot it to ax.
def plotFile(filename, name, fig, ax, ax_cbar):
  print "plotFile", filename
  fp = open(os.path.join(FLAGS.prefix, filename), 'r')
  values = [float(line) for line in fp]
  fp.close()
  plotArray(values, name, fig, ax, ax_cbar, FLAGS.polar)

# Read two files, plot the percent difference to ax.
def plotDiff(file_a, file_b, name, fig, ax, ax_cbar):
  print "plotDiff", file_a, file_b
  fp = open(os.path.join(FLAGS.prefix, file_a), 'r')
  values_a = [float(line) for line in fp]
  fp.close()
  fp = open(os.path.join(FLAGS.prefix, file_b), 'r')
  values_b = [float(line) for line in fp]
  fp.close()
  spread = max(values_b) - min(values_b)
  # assumption, file_b is baseline
  diff = np.array([b_i - a_i for a_i, b_i in zip(values_a, values_b)])
  diff_percent = abs(diff) * 100.0 / spread
  mean_squared_error = (diff ** 2).mean()
  name = name + '\nMean Squared Error: %f' % mean_squared_error
  plotArray(diff_percent, name, fig, ax, ax_cbar, FLAGS.polar)

# Read two files, plot the max percent difference in both dimensions separately.
def plotDiff_1D(file_a, file_b, plt1, plt2, title_prefix):
  print "plotDiff_1D", file_a, file_b
  fp = open(os.path.join(FLAGS.prefix, file_a), 'r')
  values_a = [float(line) for line in fp]
  fp.close()
  fp = open(os.path.join(FLAGS.prefix, file_b), 'r')
  values_b = [float(line) for line in fp]
  fp.close()
  spread = max(values_b) - min(values_b)
  # assumption, file_b is baseline
  diff = [abs(b_i - a_i) * 100.0 / spread for a_i, b_i in zip(values_a, values_b)]

  npvalues = np.array(diff).reshape(len(radii),len(angles))
  # max values per row
  err_radial = map(max, npvalues)
  # max values per column
  err_tangential = map(max, npvalues.T)

  my_blue = '#39b3e6'

  plt1.set_title("Relative " + title_prefix + " error")
  plt1.set_xlim([-0.1, 2.1])
  plt1.set_xlabel('azimuth')
  plt1.set_ylabel('maximal difference')
  plt1.xaxis.set_major_formatter(piformat)
  plt1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f%%'))
  plt1.plot(angles_no_pi, err_tangential, color=my_blue, ls='-')
  plt1.grid()


  plt2.set_title("Relative " + title_prefix + " error")
  plt2.set_xlabel('radius')
  plt2.set_xlim([radii[0] - 0.03, radii[-1] + 0.03])
  plt2.set_ylabel('maximal difference')
  plt2.plot(radii, err_radial, color=my_blue, ls='-')
  plt2.grid()

# Read a mass file and plot it in a single chart.
def plotMass(level, fig, ax, ax_cbar, polar):
  fp = open(os.path.join(FLAGS.prefix, 'mass_%.2d.txt' % level), 'r')
  # Masses are negative (from -sigma term)
  masses = [-float(line) for line in fp]
  fp.close()

  plotArray(masses, 'Mass Level %d' % (level), fig, ax, ax_cbar, polar)

# Generalized plotting method. Autodetects higher levels, can switch to polar rendering if requested.
def plotArray(values, name, fig, ax, ax_cbar, polar):
  size = int(math.sqrt(len(values)/2))
  value_min = min(values)
  value_max = max(values)
  npvalues = np.array(values).reshape(size, size*2)

  level_mult = len(radii) / size
  dr_level = dr * level_mult
  dt_level = dtheta * level_mult

  if polar:
    t, r = np.mgrid[slice(min(angles),max(angles)+dt_level, dt_level),
                    slice(min(radii), max(radii)+dr_level, dr_level) ]
    cax = ax.pcolor(t, r, npvalues.T, cmap=plt.get_cmap('jet'), vmin=value_min, vmax=value_max)
    ax.set_xlim([t.min(), t.max()])
    ax.set_ylim([0, r.max()]);
    ax.get_yaxis().set_visible(False)
  else:
    dt_level = dt_level / np.pi
    r, t = np.mgrid[slice(min(radii),max(radii)+dr_level, dr_level),
                    slice(min(angles_no_pi), max(angles_no_pi)+dt_level, dt_level) ]
    cax = ax.pcolor(r, t, npvalues, cmap=plt.get_cmap('jet'), vmin=value_min, vmax=value_max)
    ax.set_xlim([r.min(), r.max()])
    ax.set_ylim([t.min(), t.max()]);
    ax.set_xlabel("r")
    ax.set_ylabel(u"θ")
    ax.yaxis.set_major_formatter(piformat)
  ax.set_title(name)

  if ax_cbar:
    fig.colorbar(cax, cax=ax_cbar)

# Always read the radii and angles.
fp = open('r_project.data', 'r')
radii = [float(line) for line in fp]
fp.close()
dr = radii[1]-radii[0]

fp = open('theta_project.data', 'r')
angles = [float(line) for line in fp]
fp.close()

# Divide the angles by pi so that the numbers are in interval [0, 2)
angles_no_pi = np.array(angles) / np.pi

dtheta = angles[1]-angles[0]

# Save to a file.
def save(filename, is2d=True):
  if is2d:
    suffix = '-polar' if FLAGS.polar else '-euklid'
    filename = filename.replace('.png', suffix + '.png')
  plt.savefig(os.path.join(FLAGS.prefix, filename))
  print "saved " + filename

# Save to a file if --save is set on the commandline, otherwise show the plot.
def saveOrShow(filename, is2d=True):
  if FLAGS.save:
    save(filename, is2d)
  else:
    plt.show()
  plt.close()

if (FLAGS.all or FLAGS.levels):
  fig = plt.figure()
  ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], polar=FLAGS.polar)
  plotFile('level.data', 'Lookup level', fig, ax, ax_cbar=None)
  saveOrShow("debug_level.png")


if (FLAGS.all or FLAGS.expt):
  plot2('force_r.data', 'Radial force', 'force_theta.data', 'Tangential Force')
  saveOrShow("forcesA.png")

if (FLAGS.all or FLAGS.base):
  plot2('level0/force_r.data', 'Radial force (L0)', 'level0/force_theta.data', 'Tangential Force (L0)')
  saveOrShow("forces0.png")

if FLAGS.all or FLAGS.diff or FLAGS.diff2D:
  plotDiff2('force_r.data', 'level0/force_r.data',
            'force_theta.data', 'level0/force_theta.data')
  saveOrShow("diff_2D.png")

if FLAGS.all or FLAGS.diff or FLAGS.diff1D:
  plotDiff2_1D('force_r.data', 'level0/force_r.data',
            'force_theta.data', 'level0/force_theta.data')
  saveOrShow("diff_1D.png", is2d=False)

if (FLAGS.all or FLAGS.masses):
  for i in range(0,100):
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.1, 0.75, 0.8], polar=FLAGS.polar)
    ax_cbar = fig.add_axes([0.83, 0.1, 0.03, 0.8])
    if not os.path.isfile(os.path.join(FLAGS.prefix, 'mass_%.2d.txt' % i)):
      break
    plotMass(i, fig, ax, ax_cbar, FLAGS.polar)
    save("mass_%.2d.png" % i)

if (FLAGS.all or FLAGS.summary):
  plot2('density_project.data', 'Density (input)', 'force_mag.data', 'Force magnitude')
  saveOrShow("summary.png")

