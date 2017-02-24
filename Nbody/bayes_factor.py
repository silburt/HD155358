import numpy as np
import glob
import sys
from progress.bar import Bar
import matplotlib.pyplot as plt
import rebound

#initial args/settings
dir = sys.argv[1]
files = glob.glob(dir+'*_flatlnprob.npy')
n_params = 6  #x_s, x_t, y_s, y_t, phi, jitter2
glob_prob = np.empty(0)

#main loop
bar = Bar('Processing', max=len(files))
for f in files:
    lnL = np.load(f)[:,500:].reshape((-1))
    glob_prob = np.append(glob_prob, np.sum(np.exp(lnL)))
    
    bar.next()
bar.finish()

plt.plot(glob_prob)
plt.yscale('log')
plt.show()

print glob_prob/glob_prob[0]
print glob_prob/glob_prob[1]
print glob_prob/glob_prob[2]
