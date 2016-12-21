#The purpose of this macro is to see how the evolution of a system varies with the number of planetesimals.
#All runs have the same initial conditions except for the number of planetesimals (same powerlaw, total disk mass, etc.). How is the evolution of the system affected by the number of planetesimals?

import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

#---------------------------------------------------------
#makes every simulation output at same time.
#If equal to 1, it takes the largest possible values common to all runs
#If equal to 2, the user supplies a time as the second argument to the program.
same_output_time = 1

#choice for what to plot
choice = sys.argv[3]

#Npl_cutoff - ignore all files with initial Np > N_cutoff
Npl_cutoff = int(sys.argv[2])

#Number of points to average the values over (i.e. the values oscilate)
N_trailing_avg = 30
#---------------------------------------------------------

def get_vars(filename,choice,same_output_time,t_min):
    time, dE, N, N_mini, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(f, delimiter=',', unpack=True)
    index = 0
    if same_output_time >= 1:
        index = max(np.where(time<t_min)[0])
    else:
        index = 1

    if choice == 'P':
        P = np.sqrt((a2**3)/(a1**3))
        var = np.mean(P[index-N_trailing_avg:index])
        var2 = 0
        name = 'Period Ratio (P2/P1)'
    if choice == 'N':
        Nratio = N/N[0]
        var = np.mean(Nratio[index-N_trailing_avg:index])
        var2 = 0
        name = 'N_i / N_tot'
    if choice == 'e':
        var = np.mean(e1[index-N_trailing_avg:index])
        var2 = np.mean(e2[index-N_trailing_avg:index])
        name = 'e'
    if choice == 'a':
        var = np.mean(a1[index-N_trailing_avg:index])
        var2 = np.mean(a2[index-N_trailing_avg:index])
        name = 'a'
    return time[index], var, var2, N[0], name

dir = sys.argv[1]
files = glob.glob(dir+'*sd*.txt')
N = len(files)
i=0
while i < N:    #just want the main .txt files
    f = files[i]
    string = f.split("_")
    try:
        Npl = float(string[-2].split("Np")[1])
    except:
        Npl = 10e10
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-2]=="eiasnapshot" or Npl>Npl_cutoff:
        files.remove(files[i])
        N -= 1
    else:
        i += 1

t_min = 1e15
if same_output_time == 1:   #find tmin
    for f in files:
        with open(f) as fh:
            for line in fh:
                pass
            t_min = min(t_min, float(line.split(",")[0]))
elif same_output_time == 2:
    t_min = float(sys.argv[2])

var = []
var2 = []
Np = []
t = []
for f in files:
    time, variable, variable2, N, name = get_vars(f,choice,same_output_time,t_min)
    var.append(variable)
    var2.append(variable2)
    t.append(time)
    Np.append(N)


im=plt.scatter(Np, var, c=t, cmap=cm.rainbow,lw=0,label='inner')
if sum(var2) > 0:
    plt.scatter(Np, var2, marker='+', color='black', label='outer')
plt.colorbar(im, label='elapsed time (yr)')
plt.legend(loc='lower left')
plt.xscale('log')
plt.ylabel(name)
plt.xlabel('Np')
plt.savefig(dir+'NpvEvolution-'+choice+'_Np%d'%Npl_cutoff+'.png')
#plt.show()

