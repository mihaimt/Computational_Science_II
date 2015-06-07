import numpy
from matplotlib import pyplot
import time
#nr = 128+4
#nt = 256

nr = 128 #128
nt = 256# 256

def reading_data(nr,nt):

        r_file = open("/home/ics/mihai/git/Computational_Science_II/Data/r_project.data", 'r')
        t_file = open("/home/ics/mihai/git/Computational_Science_II/Data/theta_project.data", 'r')
        d_file = open("/home/ics/mihai/git/Computational_Science_II_Open/acc_r_lvl0_new.data", 'r')
#	d_file = open("/home/ics/mihai/git/Computational_Science_II_Open/acc_r_project.data", 'r')
#	d_file = open("/home/ics/mihai/git/Computational_Science_II_Open/acc_r_lvl0_new.data", 'r')

#      r_file = open("/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/radius_l.data", 'r')
#       t_file = open("/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/theta_l.data", 'r')
#       d_file = open("/home/ics/mihai/git/Computational_Science_II_Open/New_attempt/mass_l.data", 'r')

        r = []

        for line in r_file.readlines():
                r = r + [float(line[0:-1])]

        r = numpy.array(r)
        r_file.close()
        dr = r[1]-r[0]
        t = []

        for line in t_file.readlines():
                t = t + [float(line[0:-1])]

        t = numpy.array(t)
        t_file.close()
        dt = t[1]-t[0]

        d = numpy.zeros((nr,nt))
        m = numpy.zeros((nr,nt))

        i = 0
        for line in d_file.readlines():
		
                d[i/nt, i%nt] = float(line[0:-1])
#                m[i/nt, i%nt] = float(line[0:-1])*r[i/nt]*dr*dt
                i = i + 1

        d_file.close()


        return r, t, d


def linear_plot(x,y, data):
        pyplot.imshow(data, interpolation='none', aspect = 'auto')#extent=[min(x),max(x),min(y), max(y)], aspect = 'auto')
	pyplot.colorbar(label = r"radial acceleration [AU]")
	pyplot.title("Level 0")
#        pyplot.colorbar(label = "density [DU]")
	pyplot.clim(-1.25,1.25)
	pyplot.xlabel(r"$\theta$ array")
	pyplot.ylabel(r"r array")
#        pyplot.show()
	pyplot.savefig("/home/ics/mihai/Desktop/Presentations/Computational_Science_II_June/level0_radial.png")


r, t, d = reading_data(nr, nt)
#print d[1::2,1::2].shape
#dbig = d[1::2,1::2]+d[0::2,0::2]+d[0::2,1::2]+d[1::2,0::2]

#pyplot.imshow(dbig*10**5, interpolation = 'none', aspect = 'auto')
#pyplot.colorbar()
#pyplot.show()

#for element in d:
#	print element


linear_plot(r, t, d)

