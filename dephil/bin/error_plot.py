#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 10:23:15 2015

@author: dephil
"""

import matplotlib.pyplot as plt
from numpy import zeros, max
import os.path

max_level = 4
counter_files = 0
counter_plots = 0


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


def plot_maxerror(diff_force, r, theta, name):
    """
    Plots the density on a r vs. theta plot
    """
    # rearrange density 1-D array to 2-D array
    N_r = len(r)
    N_theta = len(theta)
    diffs = to_matrix(diff_force, N_r, N_theta)
    maxerrors = list()
    for j in range(N_theta):
        maxerrors.append(max(diffs[j][:]))
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("maximal error")
    ax.set_xlabel('azimuth')
    ax.set_ylabel('maximal difference')
    my_blue = '#39b3e6'
    ax.plot(theta, maxerrors, color=my_blue, ls='-')
    ax.grid()
    plt.savefig(name)
    plt.close()



if __name__ == "__main__":

# define paths to files in ./diff/
    r_path = "../data/r_project.data"
    theta_path = "../data/theta_project.data"

    diff_radial_paths = list()
    diff_angular_paths = list()

    for i in range(1, max_level+1):
        r = "../data/diff/radial_pure_diff"+str(i)+".data"
        t = "../data/diff/angular_pure_diff"+str(i)+".data"
        diff_radial_paths.append(r)
        diff_angular_paths.append(t)

    for i in range(1, max_level+1):
        r = "../data/diff/radial_refined_diff"+str(i)+".data"
        t = "../data/diff/angular_refined_diff"+str(i)+".data"
        diff_radial_paths.append(r)
        diff_angular_paths.append(t)


# read and plot data
    r = read_path(r_path)
    theta = read_path(theta_path)

    diff_radial = list()
    diff_angular = list()

    r_files_exist = [os.path.isfile(i) for i in diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]]
    t_files_exist = [os.path.isfile(i) for i in diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]]
    if (all(r_files_exist)) and (all(t_files_exist)):
        print diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]
        print diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]
        for s in diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]:
            diff_radial.append(read_path(s))
        for d in range(1, max_level+1):
            rstring = "../pictures/radial_pure_error"+str(d)+".png"
            plot_maxerror(diff_radial[counter_plots*(max_level)+d-1], r, theta, rstring)
        for s in diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]:
            diff_angular.append(read_path(s))
        for d in range(1, max_level+1):
            tstring = "../pictures/angular_pure_error"+str(d)+".png"
            plot_maxerror(diff_angular[counter_plots*(max_level)+d-1], r, theta, tstring)
        counter_plots = counter_plots + 1
        print "pure errors plotted"
    counter_files = counter_files + 1

    r_files_exist = [os.path.isfile(i) for i in diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]]
    t_files_exist = [os.path.isfile(i) for i in diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]]
    if (all(r_files_exist)) and (all(t_files_exist)):
        print diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]
        print diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]
        for s in diff_radial_paths[counter_files*max_level:((counter_files*max_level)+max_level)]:
            diff_radial.append(read_path(s))
        for d in range(1, max_level+1):
            rstring = "../pictures/radial_refined_error"+str(d)+".png"
            plot_maxerror(diff_radial[counter_plots*(max_level)+d-1], r, theta, rstring)
        for s in diff_angular_paths[counter_files*max_level:((counter_files*max_level)+max_level)]:
            diff_angular.append(read_path(s))
        for d in range(1, max_level+1):
            tstring = "../pictures/angular_refined_error"+str(d)+".png"
            plot_maxerror(diff_angular[counter_plots*(max_level)+d-1], r, theta, tstring)
        counter_plots = counter_plots + 1
        print "refined errors plotted"
    counter_files = counter_files + 1

