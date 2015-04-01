import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

## the data
N = 9
time = [409.2, 235, 249, 242, 22.8, 22.1, 19.3, 18.3, 16.8]


## necessary variables
ind = np.arange(N)                # the x locations for the groups
width = 0.35                      # the width of the bars

## the bars
rects1 = ax.bar(ind, time, width,color='black')



# axes and labels
ax.set_xlim(-width,len(ind)+width)
#ax.set_ylim(0,450)
ax.set_ylabel('cpu time [s]', fontsize = 20)
xTickMarks = ["raw", "precalc 1", "reduct", "precalc 2", "optimiz 1", "symmetry", "precalc 3", "optimiz 2", "memory" ]
ax.set_xticks(ind+width)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=45, fontsize=10)
plt.yscale('log')
i = 0
for t in time:
         
        ax.text(ind[i]+0.2, 1.05*t, '%d'%int(t)+' s',
                ha='center', va='bottom')
	i = i+1
## add a legend
plt.savefig("/home/ics/mihai/git/Computational_Science_II_Open/time_plot.png")

