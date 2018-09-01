
# coding: utf-8

# # Galactic backtracking
# The following setup shows how to use CRPropa for backtracking simulations.  
# In the JF12 model the Galaxy is a sphere of 20 kpc radius.
# For the magnetic field we are going to consider the regular component of the JF2012 model. T
# The large-scale (striated) and small-scale (turbulent) random 
# components can optionally be activated with the outcommented sections and a random seed can be set for reproducability.

# In[1]:

def CoLat2Lat(latitude):
	return (np.pi/2 - latitude)

def Lat2CoLat(latitude):
	return (np.pi/2 - latitude)
def Convert2Plot(lon):
	lon[(lon > np.pi)] = lon[(lon > np.pi)] - (2.0*np.pi)
	lon[(lon <= np.pi)] = -lon[(lon <= np.pi)]
	return lon

# def flip_data()

from crpropa import *
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

# magnetic field setup
B = JF12Field()
seed = 691342
B.randomStriated(seed)
B.randomTurbulent(seed)

# simulation setup
sim = ModuleList()
sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
obs = Observer()
obs.add(ObserverLargeSphere(Vector3d(0), 20 * kpc))
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)
print sim


## Backtracking including uncertainties
R = Random()  # CRPropa random number generator
pid = - nucleusId(1,1)
meanEnergy = 1 * EeV
sigmaEnergy = 0.1 * meanEnergy  # 10% energy uncertainty
position = Vector3d(-8.5, 0, 0) * kpc

# TA hotspot coordinates
# 177.14661479433852 49.59816303592428
lat0 = (49.598)*np.pi/180.0
lon0 = (177.146)*np.pi/180.0
meanDir = Vector3d()
meanDir.setRThetaPhi(1, lat0, lon0)
sigmaDir = 0.02  # directional uncertainty
lons, lats = [], []
x,y = [], []
energies = []

ta_energy, ta_ra, ta_dec = np.loadtxt("ta_hotspot.dat", unpack=True, usecols=(1,2,3))
ta_hotspot = SkyCoord(ra=ta_ra*u.degree, dec=ta_dec*u.degree, frame='icrs')

energies = []
energies0 = []

for i in range(len(ta_energy)):

	print (ta_hotspot[i].galactic.l.degree)

	meanDir = Vector3d()
	if ta_hotspot[i].galactic.l.radian > np.pi: 
		l = ta_hotspot[i].galactic.l.radian - (2.0 * np.pi)
	else:
		l = ta_hotspot[i].galactic.l.radian


	meanDir.setRThetaPhi(1, np.pi/2 - ta_hotspot[i].galactic.b.radian, l)

	x.append(meanDir.getPhi())
	y.append(meanDir.getTheta())
	energies0.append(ta_energy[i])

	print (i, ta_energy[i])
		

	for j in range(10):
		energy = R.randNorm(10.0 * EeV, sigmaEnergy)
		#energy = ta_energy[i] * EeV
		#print (energy)
		#direction = Vector3d()
		#direction.setRThetaPhi(1, np.pi/2 - ta_hotspot[i].galactic.b.radian, float(ta_hotspot[i].galactic.l.radian))
		direction = R.randVectorAroundMean(meanDir, sigmaDir)


		c = Candidate(ParticleState(pid, energy, position, direction))
		sim.run(c)

		d1 = c.current.getDirection()
		lons.append(d1.getPhi())
		lats.append(d1.getTheta())
		energies.append(energy/EeV)

		# if j == 0:
		# 	x.append(direction.getPhi())
		# 	y.append(direction.getTheta())
		# 	energies0.append(energy/EeV)



# Plotting
lats = np.pi/2 - np.array(lats)
y = np.pi/2 - np.array(y)


observed = SkyCoord(l=x*u.radian, b=y*u.radian, frame='galactic')
source = SkyCoord(l=lons*u.radian, b=lats*u.radian, frame='galactic')



plt.figure(figsize=(12,7))
plt.subplot(111, projection = 'mollweide')
plt.title("Backtracking TA Events, Energies ~10 EeV")
import matplotlib
cmap = matplotlib.cm.get_cmap('Spectral')
indices = np.arange(len(lons))
norm = matplotlib.colors.Normalize(vmin=0, vmax=len(lons))
colors = cmap(norm(indices))

# for i in range(len(lons)):
# 	plt.plot([x.flatten()[i],lons[i]], [y.flatten()[i],lats[i]], c="k")

plt.gca().set_xticklabels(np.arange(150,-180,-30))


plt.scatter(Convert2Plot(observed.supergalactic.sgl.radian), observed.supergalactic.sgb.radian, alpha=0.5, marker='s', c=energies0, s=80, cmap="Spectral")
plt.scatter(Convert2Plot(source.supergalactic.sgl.radian), source.supergalactic.sgb.radian, alpha=0.5, marker='.', c=energies, cmap="Spectral")
plt.colorbar(shrink=0.5)
#plt.scatter((360.0-240.16)*np.pi/180.0, -56.69*np.pi/180.0, marker='*', c=colors, s=100)
#plt.grid(True)
plt.savefig("ta_event_backtrack_10EeV.png", dpi=300)
