# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 21:16:03 2015

@author: omer
"""

from pylab import *
import os
# set the interactive mode of pylab to ON
ion()
# opens a new figure to plot into
fig_hndl = figure()
# make an empty list into which we'll append
# the filenames of the PNGs that compose each
# frame.
files=[]   
# filename for the name of the resulting movie
filename = 'animation'
number_of_frames = 10
for t in range(number_of_frames):
    # draw the frame
    imshow(rand(1000,1000))
    # form a filename
    fname = 'Animation/Level1/DFr_theta/DFrtheta%03d.png'%t
    # save the frame
    savefig(fname)
    # append the filename to the list
    files.append(fname)
# call mencoder
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o " + filename + ".mpg")
# cleanup
for fname in files: os.remove(fname):
    
ffmpeg()