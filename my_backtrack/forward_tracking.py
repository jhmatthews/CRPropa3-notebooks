
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
from astropy.coordinates import SkyCoord
from astropy import units as u

def CoLat2Lat(latitude):
	return (np.pi/2 - latitude)

def Lat2CoLat(latitude):
	return (np.pi/2 - latitude)

# magnetic field setup
B = JF12Field()
seed = 691342
B.randomStriated(seed)
B.randomTurbulent(seed)

# simulation setup
sim = ModuleList()
sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
obs = Observer()
obs.add(ObserverSmallSphere(Vector3d(-8.5, 0, 0) * kpc))
sim.add(obs)
print sim


## Backtracking including uncertainties
R = Random()  # CRPropa random number generator
pid = nucleusId(1,1)

# TA hotspot coordinates
lats0 = np.array([0,0,0,41.5,41.5,41.5,41.5,41.5,41.5,41.5,41.5]) * np.pi / 180.0
lons0 = np.array([-90,0,90,-135,-45,45,135,-135,-45,45,135]) * np.pi / 180.0

lats0 = CoLat2Lat(lats0)

for i in range(len(lats0)):
    position = Vector3d(20, 0, 0) * kpc
    meanDir = Vector3d()
    meanDir.setRThetaPhi(1, lats0[i], lons0[i])
    sigmaDir = 0.002  # directional uncertainty

    print (i)


    for j in range(10):

        print j
        energy = 10.0 * EeV
        direction = R.randVectorAroundMean(meanDir, sigmaDir)
        c = Candidate(ParticleState(pid, energy, position, direction))
        sim.run(c)
        d1 = c.current.getDirection()
        lons.append(d1.getPhi())
        lats.append(d1.getTheta())

        x.append(direction.getPhi())
        y.append(direction.getTheta())
        energies.append(energy/EeV)


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
lats = Lat2CoLat(np.pi/2 - np.array(lats))
y = Lat2CoLat(np.pi/2 - np.array(y))


observed = SkyCoord(l=x*u.radian, b=y*u.radian, frame='galactic')
source = SkyCoord(l=lons*u.radian, b=lats*u.radian, frame='galactic')



plt.figure(figsize=(12,7))
plt.subplot(111, projection = 'mollweide')

import matplotlib
cmap = matplotlib.cm.get_cmap('Spectral')
indices = np.arange(len(lons))
norm = matplotlib.colors.Normalize(vmin=0, vmax=len(lons))
colors = cmap(norm(indices))

# for i in range(len(lons)):
# 	plt.plot([x.flatten()[i],lons[i]], [y.flatten()[i],lats[i]], c="k")

plt.gca().set_xticklabels(np.arange(150,-180,-30))
# plt.scatter(-np.array(x), y, marker='o', c=colors, s=100)
# plt.scatter(-np.array(lons), lats, marker='s', c=colors, alpha=1, s=100)


plt.scatter(-observed.supergalactic.sgl.radian, observed.supergalactic.sgb.radian, alpha=0.5, marker='s', c=energies, s=80, cmap="Spectral")
plt.scatter(-source.supergalactic.sgl.radian, source.supergalactic.sgb.radian, alpha=0.8, marker='.', c=energies, cmap="Spectral")
plt.colorbar(shrink=0.5)
plt.savefig("ta_hotspot_backtrack.png", dpi=300)
