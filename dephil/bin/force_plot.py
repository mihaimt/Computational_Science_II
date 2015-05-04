#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:32:47 2015

@author: dephil
"""
import matplotlib.pyplot as plt
from numpy import zeros

max_level = 3

# define the paths to the files
r_path = "../data/r_project.data"
theta_path = "../data/theta_project.data"
density_path = "../data/density_project.data"
#force_r_path = "../data/f_radial.data"
#force_theta_path = "../data/f_angular.data"
force_lvl_r_paths = list()
force_lvl_theta_paths = list()
for i in range(max_level+1):
    r = "../data/f_radial_lvl"+str(i)+".data"
    t = "../data/f_angular_lvl"+str(i)+".data"
    force_lvl_r_paths.append(r)
    force_lvl_theta_paths.append(t)
for i in range(max_level+1):
    r = "../data/radial_osc_mass_lvl"+str(i)+".data"
    t = "../data/angular_osc_mass_lvl"+str(i)+".data"
    force_lvl_r_paths.append(r)
    force_lvl_theta_paths.append(t)
#for i in range(max_level+1):
#    r = "../data/radial_osc_force_lvl"+str(i)+".data"
#    t = "../data/angular_osc_force_lvl"+str(i)+".data"
#    force_lvl_r_paths.append(r)
#    force_lvl_theta_paths.append(t)

#diff_path = "../data/rad_m3_to_f3.data"


def read_path(path):
    """
    Reads data and stores it in a 1-D list/array as floats
    """
    # open data
    d_data = open(path, 'r')
    d_ = list()
    # first line read
    aline = d_data.readline()
    while aline:
        items = aline.split() # in case there are spaces
        d_.append(float(items[0])) # have python floats double precision ???
        aline = d_data.readline() # next line
    d_data.close()
    return d_

    
def to_matrix(ilist, N_r, N_theta):
    """
    Rearranges the olist 1-D to 2-D array
    """
    surface = zeros((N_theta,N_r))
    for i in range(N_r):
        for j in range(N_theta):
            surface[j][i] = ilist[j+N_theta*i]
    return surface


def plot_density(density, r, theta, name):
    """
    Plots the density on a r vs. theta plot
    """
    # rearrange density 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = to_matrix(density, N_r, N_theta)
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("density map")
    ax.set_xlabel(r'radius [$r_{0}$]')
    ax.set_ylabel('azimuth')
    # plot the 2-D array
    cax = ax.imshow(surface, aspect='auto', origin='lower',
              extent=[min(r), max(r), min(theta), max(theta)],
              cmap=plt.get_cmap('jet'),# interpolation='hamming', # change interpolation to None, if not wanted
              vmin=min(density), vmax=max(density))
    fig.colorbar(cax)
    fig.savefig(name)
    plt.close()


def plot_force_radial(force, r, theta, name):
    """
    Plots the force on a r vs. theta plot
    """
    # rearrange force 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = to_matrix(force, N_r, N_theta)
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("force map (radial component)")
    ax.set_xlabel('radius [r]')
    ax.set_ylabel('azimuth')
    # plot the 2-D array
    cax = ax.imshow(surface, aspect='auto', origin='lower',
              extent=[min(r)-.5*(r[1]-r[0]), max(r)-.5*(r[-1]-r[-2]), min(theta)-.5*(theta[1]-theta[0]), max(theta)-.5*(theta[-1]-theta[-2])],
              cmap=plt.get_cmap('jet'),# interpolation='hamming', # change interpolation to None, if not wanted
              vmin=min(force), vmax=max(force))
    fig.colorbar(cax)
    fig.savefig(name)
    plt.close()


def plot_force_angular(force, r, theta, name):
    """
    Plots the force on a r vs. theta plot
    """
    # rearrange force 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = to_matrix(force, N_r, N_theta)
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("force map (angular component)")
    ax.set_xlabel('radius [r]')
    ax.set_ylabel('azimuth')
    # plot the 2-D array
    cax = ax.imshow(surface, aspect='auto', origin='lower',
              extent=[min(r)-.5*(r[1]-r[0]), max(r)-.5*(r[-1]-r[-2]), min(theta)-.5*(theta[1]-theta[0]), max(theta)-.5*(theta[-1]-theta[-2])],
              cmap=plt.get_cmap('jet'),# interpolation='hamming', # change interpolation to None, if not wanted
              vmin=min(force), vmax=max(force))
    fig.colorbar(cax)
    fig.savefig(name)
    plt.close()


if __name__ == "__main__":

# read data
    r = read_path(r_path)
    theta = read_path(theta_path)
    density = read_path(density_path)
    
    #force_r = read_path(force_r_path)
    #force_theta = read_path(force_theta_path)

    force_lvl_r = list()
    force_lvl_theta = list()
    for s in force_lvl_r_paths:
        force_lvl_r.append(read_path(s))
    for s in force_lvl_theta_paths:
        force_lvl_theta.append(read_path(s))

#    diff = read_path(diff_path)

        
# plots
    #plot_density(density, r, theta, '../pictures/density.png')
    #plot_force_radial(force_r, r, theta, '../pictures/radial_force.png')
    #plot_force_angular(force_theta, r, theta, '../pictures/angular_force.png')

    for d in range(max_level+1):
        rstring = "../pictures/radial_force_lvl"+str(d)+".png"
        plot_force_radial(force_lvl_r[d], r, theta, rstring)
    for d in range(max_level+1):
        tstring = "../pictures/angular_force_lvl"+str(d)+".png"
        plot_force_angular(force_lvl_theta[d], r, theta, tstring)

    for d in range(max_level+1):
        rstring = "../pictures/osc_mass_radial_lvl"+str(d)+".png"
        plot_force_radial(force_lvl_r[1*(max_level+1)+d], r, theta, rstring)
    for d in range(max_level+1):
        tstring = "../pictures/osc_mass_angular_lvl"+str(d)+".png"
        plot_force_angular(force_lvl_theta[1*(max_level+1)+d], r, theta, tstring)

#    for d in range(max_level+1):
#        rstring = "./pictures/osc_force_radial_lvl"+str(d)+".png"
#        plot_force_radial(force_lvl_r[2*(max_level+1)+d], r, theta, rstring)
#    for d in range(max_level+1):
#        tstring = "./pictures/osc_force_angular_lvl"+str(d)+".png"
#        plot_force_angular(force_lvl_theta[2*(max_level+1)+d], r, theta, tstring)

#    plot_force_radial(diff, r, theta, '../pictures/diff_m3_to_f3.png')
