
# coding: utf-8

# # Galactic Lensing of Simulated Cosmic Rays
# 
# Deflection in galactic magnetic fields can be accounted for using galactic lenses. To use the lenses efficiently, the ouput of a CRPropa simulation needs first to be filled in probability maps. These maps are then transformed by the lenses according to the rigidity of the particles. From the maps then finally new particles can be generated. 
# 
# ## Input Results from Extragalactic Simulation

# In[1]:

import crpropa
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import astropy_healpix as healpy
import matplotlib.pyplot as plt
from projection import *
from aab_subroutines import *
import powerlaw


def add_composition_run(sim, source, Zs, mass, names, fracs, gamma, rcut, num=1e4):
    composition = SourceComposition(1 * EeV, (rcut/1e18) * EeV, gamma)
    
    for i in range(len(Zs)):
        composition.add(mass[i], Zs[i],  fracs[i]/mass[i])  # H
    source.add( composition )
    
    # run simulation
    sim.setShowProgress(True)
    sim.run(source, int(num), True)

def CoLat2Lat(latitude):
    return (np.pi/2 - latitude)

def Lat2CoLat(latitude):
    return ( (np.pi/2.0) - latitude)

def Convert2Plot(lon):
    lon[(lon > np.pi)] = lon[(lon > np.pi)] - (2.0*np.pi)
    lon[(lon <= np.pi)] = -lon[(lon <= np.pi)]
    return lon

def initial_directions(M, source_l, source_b, weights, d, NN, energies=None):
    lons0 = []
    lats0 = []

    for isource in range(len(source_l)):
        l = source_l[isource]*np.pi/180.0 - (2.0*np.pi)
        b = Lat2CoLat(source_b[isource]*np.pi/180.0)
        #M.addParticles(Id, E, l, b, E0**-1)
        print (isource, weights[isource], l, b)
        meanDir = crpropa.Vector3d(-1, 0, 0)
        meanDir.setRThetaPhi(-1, b, l)
        kappa = 50.0


        ncrs = int((weights[isource] / total) * NN)
        #if ncrs>0: ncrs = 500
        #Emin = 10.0
        energies = powerlaw.Power_Law(xmin=10.0, xmax = 100.0, parameters=[2.2]).generate_random(ncrs)

        #print (energies)
        for i in xrange(ncrs):
            particleId = crpropa.nucleusId(1,1)
            energy = (energies[i]) * crpropa.EeV

            energy = 10.0 * crpropa.EeV
            galCenter = crpropa.Vector3d(-1,0,0)

            theta_k = 0.8 * (100.0*crpropa.EeV/energy) * np.sqrt(d[isource] / 10.0)
            kappa = (81.0 / theta_k)**2
            #print (kappa)


            momentumVector = crpropa.Random.instance().randFisherVector(meanDir, kappa)
            M.addParticle(particleId, energy, momentumVector) 
            
            lons0.append(momentumVector.getPhi())
            lats0.append(momentumVector.getTheta())

    lons0 = np.array(lons0)-np.pi
    lats0 = -CoLat2Lat(np.array(lats0))

    return (lons0, lats0)


scenario = sys.argv[1]
mylens = sys.argv[2]
fornax = sys.argv[3]

use_lens = False
if mylens == "JF12":
    use_lens = True


# read data from crpropa output into container. 
# The data is weighted with the source energy ~E**-1
M = crpropa.ParticleMapsContainer()

cont = dict()
ll, bb, d, flux, atta, attb, attc, cont["a"], cont["b"], cont["c"] = np.genfromtxt("agn.dat", usecols = np.arange(1,11), unpack=True)
lons0, lats0 = [], []


weights = cont[scenario]
updated = cont[scenario]
# updated[0] *= 7.4/3.9
# updated[1] *= 9.55/5.12

#print updated, weights
# add Fornax 
if fornax == "f":
    updated[-1] = cont[scenario][1]
    updated[1] = cont[scenario][1]/100.0 
    print updated, weights
    weights2 = updated

total = np.sum(weights)
NN = 100000

theta_k = 0.8 *  np.sqrt(20.0 / 10.0) * 10.0
kappa = (81.0 / theta_k)**2
print ("KAPPA:", kappa)
print d

lons0, lats0 = initial_directions(M, ll, bb, weights, d, NN)

# The lens can be downloaded here: https://crpropa.desy.de/ --> Additional downloads
lens_path = "../../JF12full_Gamale/lens.cfg"
lens = crpropa.MagneticLens(lens_path)
lens.normalizeLens()
if use_lens: M.applyLens(lens)

#stack all maps
# import healpy
# crMap = np.zeros(49152)
# for pid in M.getParticleIds():
#     energies = M.getEnergies(int(pid))
#     for i, energy in enumerate(energies):
#         crMap += M.getMap(int(pid), energy * crpropa.eV )
# #plot maps using healpy
# healpy.mollview(map=crMap, title='Lensed')
# plt.savefig('lensed_map.png')

pids, energies, lons, lats = M.getRandomParticles(50000)
lons = Convert2Plot(lons)

# # bin things up 
nb = 100
bins_b = np.arccos(np.linspace(1,-1 - (2.0/nb), num=nb))-np.pi/2
print (bins_b)
bins_b[-1] = 1.57079632679
bins_b[0] =  -np.pi/2
binsa1 = np.arccos(np.linspace(1,-1, num=nb/2))
binsa2 = 2.0*np.pi - binsa1
bins_a = np.linspace(-np.pi,np.pi, num=200)
bins = (bins_a, bins_b)


# # create histogram and interpolate
hist, x, y = np.histogram2d(lons, lats, bins = bins)
print (hist.shape, x[:-1].shape, y[:-1].shape)
from scipy.interpolate import interp2d 
interp_func = interp2d(x[:-1],y[:-1], hist.T, kind='cubic')

from scipy.stats import kde
# nbins=200
# x,y = lons, lats
# k = kde.gaussian_kde([lons, lats])
# xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
#plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
#



CMAP = "viridis"
# create a scatter plot of the particles
plt.figure(figsize=(10,4))
#plt.subplot(111, projection='mollweide')
plt.subplot(111, projection='custom_auger')
plt.pcolormesh( x[:-1],y[:-1],interp_func(x[:-1],y[:-1]), cmap=CMAP)

#plt.pcolormesh(xi, yi, np.zeros_like(xi), vmin=0, vmax=1)
#plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

print (np.max(x), np.max(y))
#plt.scatter(lons, lats, c=np.log10(energies), lw=0, alpha=0.5)
#plt.scatter(Convert2Plot(lons0), lats0, c="r", lw=0, alpha=0.5)
#plt.colorbar()
plt.gca().set_xticklabels(np.arange(150,-180,-30))
plt.gca().set_yticklabels(np.arange(-60,90,30))
#plt.scatter(Convert2Plot(ll*np.pi/180.0), bb*np.pi/180.0, c="w", lw=0, alpha=1, marker="*", s=100)
#plt.grid()
plt.savefig("crmap_{}_{}_{}.png".format(scenario, lens, fornax))


# plt.figure(figsize=(10,4))
# plt.subplot(111, projection='custom_auger')
# plt.pcolormesh(x[:-1],y[:-1],interp_func(x[:-1],y[:-1]), cmap=CMAP)
# #plt.scatter(lons, lats, c=np.log10(energies), lw=0, alpha=1)
# plt.colorbar()
# plt.gca().set_xticklabels(np.arange(150,-180,-30))
# #plt.gca().set_yticklabels(np.arange(-60,90,30))
# plt.scatter(Convert2Plot(ll*np.pi/180.0), bb*np.pi/180.0, c="k", lw=0, alpha=1, marker="*", s=100)
# plt.grid()
# plt.savefig('scatterd_particles_custom.png')
# plt.figure(figsize=(20,8))
# plt.subplot(111, projection='mollweide')
# hist, x, y = np.histogram2d(-lons0, -np.pi/2+lats0, bins = bins)
# print (hist.shape, x[:-1].shape, y[:-1].shape)
# from scipy.interpolate import interp2d 
# interp_func = interp2d(x[:-1],y[:-1], hist.T, kind='cubic')
# plt.pcolormesh(x[:-1],y[:-1],interp_func(x[:-1],y[:-1]), cmap="Spectral")
# #plt.show()
