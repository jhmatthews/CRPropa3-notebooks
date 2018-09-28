import matplotlib.pyplot as pl
import numpy as np 
from crpropa import * 
from scipy.optimize import minimize


def smooth(x,window_len=5,window='hanning'):

	'''smooth data x by a factor with window of length window_len'''

	if x.ndim != 1:
		raise ValueError("smooth only accepts 1 dimension arrays.")

	if x.size < window_len:
		raise ValueError("Input vector needs to be bigger than window size.")

	if window_len<3:
		return x

	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]


	if window == 'flat': #moving average
		w = np.ones(window_len,'d')
	else:  
		w = eval('np.'+window+'(window_len)')

	y=np.convolve(w/w.sum(),s,mode='same')

	return y[window_len:-window_len+1]

def make_plot_input(Zs, names, fracs, gamma, rcut, weights):
	E = np.logspace(18,21,num=1000)
	mysum = 0.0
	norm = E ** 3
	norm = 1.0
	for i in range(len(Zs)):
		flux = dn_dE(E, gamma=gamma, R_cut=rcut, Z=Zs[i], f_A=fracs[i], J0=1) * weights[i]
		J = norm * flux 
		pl.plot(E, J / J[0], label=names[i])
		mysum += flux

	J = norm * mysum
	pl.plot(E, J / J[0], c="k", linewidth=3) 
	pl.loglog()
	pl.legend(loc=4)
	pl.ylim(1e-5,1)
	pl.xlim(1e18,10.0**(20.5))

# simulation setup
def setup_sim(fname, EBL_model=IRB_Gilmore12, maxdistance=1000):
	sim = ModuleList()
	sim.add( SimplePropagation(1*kpc, 10*Mpc) )
	sim.add( Redshift() )	

	# add attenuation processes 	
	sim.add( PhotoPionProduction(CMB) )
	sim.add( PhotoPionProduction(EBL_model) )
	sim.add( PhotoDisintegration(CMB) )
	sim.add( PhotoDisintegration(EBL_model) )
	sim.add( ElectronPairProduction(CMB) )
	sim.add( ElectronPairProduction(EBL_model) )

	# add nuclear decay and minimum energy floor
	sim.add( NuclearDecay() )
	sim.add( MinimumEnergy( 1 * EeV) )

	# observer and output
	obs = Observer()
	obs.add( ObserverPoint() )
	output = TextOutput(fname, Output.Event1D)
	obs.onDetection( output )
	sim.add( obs )
	# source
	if maxdistance > 0:
		source = Source()
		source.add( SourceUniform1D(1 * Mpc, maxdistance * Mpc) )
		source.add( SourceRedshift1D() )
	else:
		source = None

	return sim, source, output


def make_plot_observed(fname, bins=None, norm = None):

	d = pl.genfromtxt(fname, names=True, invalid_raise=False)

	# observed quantities
	Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
	lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))

	if bins == None:
		lEbins = pl.arange(18, 20.51, 0.01)  # logarithmic bins
		lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	else:
		lEbins = pl.arange(17.5, 20.2, 0.1)  # logarithmic bins
	
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths

	# identify mass groups
	EE = (10.0 ** lEcens) ** 3
	J1 = EE * pl.histogram(lE, bins=lEbins)[0] / dE
	J = np.zeros_like(J1)
	Jn = np.zeros([len(J1, 58)])
	for iz in range(1,57):
		idx = (Z == iz)
		J  += pl.histogram(lE[idx], bins=lEbins)[0] / dE


	# idx1 = A == 1
	# idx2 = (A > 1) * (A < 5) 
	# idx3 = (A > 4) * (A < 23) 
	# idx4 = (A > 22) * (A < 39) 

	# # calculate spectrum: J(E) = dN/dE, we want E^3 J(E)
	# EE = (10.0 ** lEcens) ** 3

	#lum = np.sum(J * dE)
	# J  = EE * pl.histogram(lE, bins=lEbins)[0] / dE
	# J1 = EE * pl.histogram(lE[idx1], bins=lEbins)[0] / dE
	# J2 = EE * pl.histogram(lE[idx2], bins=lEbins)[0] / dE
	# J3 = EE * pl.histogram(lE[idx3], bins=lEbins)[0] / dE
	# J4 = EE * pl.histogram(lE[idx4], bins=lEbins)[0] / dE

	# normalize the histograms 
	# if norm == None:
	# 	norm = 1.0 / J[0]
	# else:
	# 	norm = norm / J[10]
	norm = 1.0 / J[0]

	J1 *= norm
	# J2 *= norm
	# J3 *= norm
	# J4 *= norm
	J *= norm

	pl.plot(lEcens, smooth(J),  color='k', linewidth=3, label="Total")
	pl.plot(lEcens, smooth(J1), color='#e84118', label='A = 1')
	pl.plot(lEcens, smooth(J2), color="#7f8fa6", label='A = 2-4')
	pl.plot(lEcens, smooth(J3), color="#44bd32", label='A = 5-22')
	pl.plot(lEcens, smooth(J4), color="#00a8ff", label='A = 23-38')

	pl.legend(fontsize=20, frameon=True, loc=3)
	pl.semilogy()
	#pl.ylim(1e-5)
	pl.grid()
	pl.ylabel('$E^3 J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	#pl.savefig('sim1D_spectrum.png')



def norm_func(norm, y1, y2):
	if (norm < 0):
		return (-norm*1e300)
	x = np.sum( ((norm * y1) - y2)**2 / y2)
	print (x, norm)
	return np.sum( ((norm * y1) - y2)**2 / y2)


def make_plot_observed2(fname, bins=None, norm = None, gamma=2, rcut=1e20, y=None):

	d = np.genfromtxt(fname, names=True, invalid_raise=False)

	# observed quantities
	Z = np.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
	A = np.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
	lE = np.log10(d['E']) + 18  # energy in log10(E/eV))

	Einit = d['E'] * 1e18
	Efinal = d['E0'] * 1e18
	weights = np.zeros_like(Efinal)
	for i in range(len(weights)):
		weights[i] = renorm(Efinal[i], Einit[i], Z[i], rcut, gamma)

	if bins == None:
		lEbins = np.arange(18, 20.51, 0.01)  # logarithmic bins
		lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	else:
		lEbins = np.arange(17.5, 20.25, 0.1)  # logarithmic bins
	
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths

	# identify mass groups
	EE = (10.0 ** lEcens) ** 3
	J1 = EE * np.histogram(lE, bins=lEbins, weights=weights)[0] / dE
	J = np.zeros_like(J1)
	J_n = [np.zeros_like(J1) for i in range(0,27)]
	for iz in range(1,27):
		idx = (Z == iz)
		J_n[iz] = np.histogram(lE[idx], bins=lEbins, weights=weights[idx])[0] / dE
		J  += J_n[iz]

	J *= EE

	select = (lEcens > 19.5)
	norm = np.max(y[select]) / np.max(J[select])

	J *= norm
	J1 *= norm
	for i in range(1,27):
		J_n[i] *= norm * EE
	# J2 *= norm
	# J3 *= norm

	pl.plot(lEcens, smooth(J),  color='k', linewidth=3, label="Total")
	pl.plot(lEcens, smooth(J_n[1]), label='Z = 1')
	pl.plot(lEcens, smooth(J_n[2]), label='Z = 2')
	pl.plot(lEcens, smooth(J_n[7]), label='Z = 7')
	pl.plot(lEcens, smooth(J_n[14]), label='Z = 14')
	pl.plot(lEcens, smooth(J_n[26]), label='~ = 26')
	#pl.plot(lEcens, smooth(J3), color="#44bd32", label='A = 5-22')
	#pl.plot(lEcens, smooth(J4), color="#00a8ff", label='A = 23-38')

	pl.legend(fontsize=20, frameon=True, loc=3)
	pl.semilogy()
	#pl.ylim(1e-5)
	pl.grid()
	pl.ylabel('$E^3 J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	#pl.savefig('sim1D_spectrum.png')

def dn_dE(E, gamma=2, R_cut=1e19, Z=1, f_A=1, J0=1):
	x = f_A * J0 * (E/1e18)**-gamma
	x *= fcut(E, Z, R_cut)
	return x

def fcut(E, Z, R_cut):
	if isinstance(E, float):
		if E < Z * R_cut:
			return 1
		else:
			return np.exp(1.0 - E/Z/R_cut)
	else:
		x = np.ones_like(E)
		select = (E > Z * R_cut)
		x[select] *= np.exp(1.0 - E[select]/Z/R_cut)
	return x

def renorm(Efinal, Einit, Z, R_cut, gamma):
	weight = Einit**(1-gamma) * fcut(Efinal, Z, R_cut)
	return weight




