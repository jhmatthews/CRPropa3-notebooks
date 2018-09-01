
# coding: utf-8

# # Galactic backtracking
# The following setup shows how to use CRPropa for backtracking simulations.  
# In the JF12 model the Galaxy is a sphere of 20 kpc radius.
# For the magnetic field we are going to consider the regular component of the JF2012 model. T
# The large-scale (striated) and small-scale (turbulent) random 
# components can optionally be activated with the outcommented sections and a random seed can be set for reproducability.

# In[1]:

from crpropa import *
import matplotlib.pyplot as plt
import numpy as np
#import astropy

# magnetic field setup
# B = JF12Field()
# seed = 691342
# B.randomStriated(seed)
# B.randomTurbulent(seed)

# simulation: a sequence of simulation modules
sim = ModuleList()

# add propagator for rectalinear propagation
sim.add( SimplePropagation() )

# add interaction modules
sim.add( PhotoPionProduction(CMB) )
sim.add( PhotoDisintegration(CMB) )
sim.add( ElectronPairProduction(CMB) )
sim.add( NuclearDecay() )
sim.add( MinimumEnergy( 1 * EeV) )

obs = Observer()
obs.add( ObserverPoint() )  # observer at x = 0
sim.add(obs)
print sim

# event output# event 
output2 = TextOutput('events.txt', Output.Event1D)
obs.onDetection(output2)

# cosmic ray source
source = Source()
source.add( SourcePosition(100 * Mpc) )
source.add( SourceParticleType(nucleusId(1, 1)) )
source.add( SourcePowerLawSpectrum(1 * EeV, 200 * EeV, -1) )
print source



sim.setShowProgress(True)  # switch on the progress bar
sim.run(source, 10000)
output2.close()


data = np.genfromtxt('events.txt', names=True)
print 'Number of events', len(data)

logE0 = np.log10(data['E0']) + 18
logE  = np.log10(data['E']) + 18

plt.figure(figsize=(10, 7))
h1 = plt.hist(logE0, bins=25, range=(18, 20.5), histtype='stepfilled', alpha=0.5, label='At source')
h2 = plt.hist(logE,  bins=25, range=(18, 20.5), histtype='stepfilled', alpha=0.5, label='Observed')
plt.xlabel('log(E/eV)')
plt.ylabel('N(E)')
plt.legend(loc = 'upper left', fontsize=20)
plt.show()
