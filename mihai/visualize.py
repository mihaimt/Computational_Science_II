from matplotlib import pyplot
import numpy
from sys import argv

def reading_1D(filename):
	fileb = open(str(filename), 'r')
	raw = fileb.readlines()
	a = numpy.zeros(len(raw))
	i = 0
	while i < len(raw):
		a[i] = float(raw[i][0:-1])
		i = i + 1

	fileb.close()
	return a

def colorbar(data, outpath, title):

#	figure(num=None, figsize=(8, 2), dpi=80, facecolor='w', edgecolor='k')

	ax1 = pyplot.subplot(111)

#	pyplot.figure(num=None, figsize=(8, 2), dpi=80, facecolor='w', edgecolor='k')


	import matplotlib as mpl
	cmap = mpl.cm.jet
	norm = mpl.colors.Normalize(vmin=min(data), vmax = max(data))

	cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,orientation='vertical')

	pyplot.savefig(outpath + "colorbar_"+title+".png")

	pyplot.cla()
	pyplot.clf()









def polar_plot(x, y, data, title, outpath):
	from matplotlib import pyplot

	ax = pyplot.subplot( projection="polar")
	ax.pcolormesh(y,x ,data_2D)
#	ax.set_rmax(max(x))
	ax.grid(True)

	ax.set_title(title, va='bottom')
	pyplot.savefig(outpath + title + "_polar.png")

def linear_plot(x, y, data, title, outpath):
	from matplotlib import pyplot

	import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	import numpy as np

	plt.figure()
	ax = plt.gca()
	

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.


        ax = pyplot.subplot()
	im =        ax.pcolormesh(y,x ,data_2D)
#        ax.set_rmax(max(x))



	ax.set_title(title, va='bottom', fontsize = 20)
        pyplot.ylabel("Radius  [$r_0$]", fontsize = 20)
        pyplot.xlabel(r"$\theta$", fontsize = 20)

        ax.grid(True)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)

	plt.colorbar(im, cax=cax)


        
        pyplot.savefig(outpath + title + "_linear.png")



# python visualize /home/ics/mihai/git/Computational_Science_II/Data/r_project.data /home/ics/mihai/git/Computational_Science_II/Data/theta_project.data /home/ics/mihai/git/Computational_Science_II_Closed/force_r.data Force_r


#script, x, y, map, outpath, title = argv

script, param = argv
param_f = open(str(param), 'r')
nr = int(param_f.readline()[0:-1])
print nr
nt = int(param_f.readline()[0:-1])
print nt
r_path = str(param_f.readline()[0:-1])
print r_path
t_path = str(param_f.readline()[0:-1])
print t_path
m_path = str(param_f.readline()[0:-1])
print m_path
o_path = str(param_f.readline()[0:-1])
print o_path
title = str(param_f.readline()[0:-1])
print title


param_f.close()


x =   reading_1D(r_path)
y =   reading_1D(t_path)
data = reading_1D(m_path)
data_2D = data.reshape(nr, nt)

r = x
theta = y


colorbar(data, str(o_path), str(title))
polar_plot(x,y, data_2D, str(title), o_path)
linear_plot(x,y, data_2D, str(title), o_path)






