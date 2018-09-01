
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
B = JF12Field()
seed = 691342
# B.randomStriated(seed)
# B.randomTurbulent(seed)

# simulation setup
sim = ModuleList()
sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
obs = Observer()
obs.add(ObserverLargeSphere(Vector3d(0), 20 * kpc))
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)
print sim


## Backtracking including uncertainties
#The impact of the cosmic ray uncertainties backtracked directions can be investigated with a MC approach. 
# In the following, the cosmic ray energy and direction are varied within the statistical uncertainties before 
# backtracking.

R = Random()  # CRPropa random number generator
pid = - nucleusId(1,1)
meanEnergy = 60 * EeV
sigmaEnergy = 0.1 * meanEnergy  # 10% energy uncertainty
position = Vector3d(-8.5, 0, 0) * kpc

lons, lats = [], []
x,y = [], []

lons0 = np.arange(-130,195,65)*np.pi/180.0 + 1e-6
lats0 = np.arange(0,180,30)*np.pi/180.0  + 1e-6

# lons0 = np.array([(360.0-275.0)*np.pi/180.0])
# lats0 = np.array([(165.0)*np.pi/180.0])

for i in range(len(lons0)):
	for j in range(len(lats0)):

		direction = Vector3d()
		direction.setRThetaPhi(1, lats0[j], lons0[i])
		energy = 10.0 * EeV

		c = Candidate(ParticleState(pid, energy, position, direction))
		sim.run(c)

		d1 = c.current.getDirection()
		lons.append(d1.getPhi())
		lats.append(d1.getTheta())

		x.append(lons0[i])
		y.append(lats0[j])


special_lon = -(360.0-275.0)*np.pi/180.0
special_lat = (165.0)*np.pi/180.0
direction = Vector3d()
direction.setRThetaPhi(1, special_lat, special_lon)
c = Candidate(ParticleState(pid, energy, position, direction))
sim.run(c)

d1 = c.current.getDirection()
lons.append(d1.getPhi())
lats.append(d1.getTheta())
x.append(special_lon)
y.append(special_lat)

# ## (Optional) Plotting
# Finally we are plotting a skymap of the observed direction along with the distribution of
# directions at the galactic border.

# In[4]:

# Angle definitions:
# CRPropa uses
#   longitude (phi) [-pi, pi] with 0 pointing in x-direction
#   colatitude (theta) [0, pi] with 0 pointing in z-direction
# matplotlib expects
#   longitude [-pi, pi] with 0 = 0 degrees
#   latitude [pi/2, -pi/2] with pi/2 = 90 degrees (north)

# lat0 = np.pi/2 - lat0
lats = np.pi/2 - np.array(lats)
y = np.pi/2 - np.array(y)

plt.figure(figsize=(12,7))
plt.subplot(111, projection = 'mollweide')

import matplotlib
cmap = matplotlib.cm.get_cmap('Spectral')
indices = np.arange(len(lons))
norm = matplotlib.colors.Normalize(vmin=0, vmax=len(lons))
colors = cmap(norm(indices))

# for i in range(len(lons)):
# 	plt.plot([x.flatten()[i],lons[i]], [y.flatten()[i],lats[i]], c="k")

plt.gca().set_xticklabels(np.arange(180,-180,-30))
plt.scatter(-np.array(x), y, marker='o', c=colors, s=100)
plt.scatter(-np.array(lons), lats, marker='s', c=colors, alpha=1, s=100)

plt.scatter((360.0-240.16)*np.pi/180.0, -56.69*np.pi/180.0, marker='*', c=colors, s=100)
#plt.grid(True)
plt.savefig("farrar.png")
