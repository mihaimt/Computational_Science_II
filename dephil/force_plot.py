#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:32:47 2015

@author: dephil
"""
import matplotlib.pyplot as plt
from numpy import zeros

# define the paths to the files
r_path = "./data/r_project.data"
theta_path = "./data/theta_project.data"
density_path = "./data/density_project.data"
force_r_path = "./data/f_radial.data"
force_theta_path = "./data/f_angular.data"

def read_path(path):
    """
    Reads r's and stores it in a 1-D list/array as floats
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

    
def density_matrix(density, N_r, N_theta):
    """
    Rearranges the density 1-D to 2-D array
    """
    surface = zeros((N_theta,N_r))
    for i in range(N_r):
        for j in range(N_theta):
            surface[j][i] = density[j+N_theta*i]
    return surface


def force_matrix(force, N_r, N_theta):
    """
    Rearranges the force 1-D to 2-D array
    """
    surface = zeros((N_theta,N_r))
    for i in range(N_r):
        for j in range(N_theta):
            surface[j][i] = force[j+N_theta*i]
    return surface


def plot_density(density, r, theta):
    """
    Plots the density on a r vs. theta plot
    """
    # rearrange density 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = density_matrix(density, N_r, N_theta)
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
    fig.savefig('density.png')


def plot_force_radial(force, r, theta):
    """
    Plots the force on a r vs. theta plot
    """
    # rearrange force 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = force_matrix(force, N_r, N_theta)
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
    fig.savefig('radial_force.png')


def plot_force_angular(force, r, theta):
    """
    Plots the force on a r vs. theta plot
    """
    # rearrange force 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    surface = force_matrix(force, N_r, N_theta)
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
    fig.savefig('angular_force.png')


if __name__ == "__main__":

    # read data
    r = read_path(r_path)
    theta = read_path(theta_path)
    density = read_path(density_path)
    force_r = read_path(force_r_path)
    force_theta = read_path(force_theta_path)
    
    # plots
    plot_density(density, r, theta)
    plot_force_radial(force_r, r, theta)
    plot_force_angular(force_theta, r, theta)
