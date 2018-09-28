import matplotlib
matplotlib.use("TkAgg")
from crpropa import *
import numpy as np 
from aab_subroutines import * 
import matplotlib.pyplot as plt

def plot_auger_data():
    logE, n, dn1, dn2 = np.loadtxt("auger.dat", unpack=True)

    # convert from m^{-2} / s to km^{-2} yr^{-1}
    conversion = 1e6 * 3.154e+7 * 1
    Elinear = 10.0**logE
    y = Elinear**2 * n * conversion
    yerr1 = Elinear**2*dn1 * conversion
    yerr2 = Elinear**2*dn2 * conversion
    pl.errorbar(logE, y, yerr=(yerr1, yerr2), fmt="o", label="PAO data", c="k")
    pl.semilogy()
    
    return (Elinear, y)

def add_composition_run(sim, source, Zs, mass, names, fracs, gamma, rcut, num=1e4):
    composition = SourceComposition(1 * EeV, (rcut/1e18) * EeV, gamma)
    
    for i in range(len(Zs)):
        composition.add(mass[i], Zs[i],  fracs[i]/mass[i])  # H
    source.add( composition )
    
    # run simulation
    sim.setShowProgress(True)
    sim.run(source, int(num), True)

Zs = [1,2,7,14,26]
mass = [1,4,14,28,56]
names = ["H", "He", "N", "Si", "Fe"]
fracs = [0.064, 0.467, 0.375, 0.094, 0.0]
rcut = 10.0**18.66
gamma = 0.93

fname = "modela.txt"
sim, source, output = setup_sim(fname, maxdistance=1000)
add_composition_run(sim, source, Zs, mass, names, fracs, 1, 1e20, num=1e5)
output.close()

pl.figure(figsize=(10,7))
Elinear, y = plot_auger_data()
make_plot_observed2(fname, bins="auger", norm=None, rcut=rcut, gamma=gamma, y=y)
pl.xlim(18,20.5)
pl.ylim(5e35,1e38)
pl.show()