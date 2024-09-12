import matplotlib.pyplot as plt

from ocelot.gui import *
from ocelot import *
from thz_match_Functions import *
from pylab import *
import numpy as np
import cmath
from scipy.special import jv, yv
import scipy.integrate as integrate
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
from ocelot.cpbd.match import *
from scipy.special import genlaguerre
from scipy.special import factorial
import seaborn as sns

cLight = 299792458
er = 5.63 # 
a = 1e-3

mirror_radiusSmall = 50e-3
# mirror_radiusSmall = 60e-3
mirror_radiusBig = 75e-3
window_radius = 25e-3
# mirror_radius = 37e-3


energy = 16 # GeV
gamma = energy/(0.511e-3)
# M_squared = emittance / (1/(2*k))

Resonance_Frequency = 300e9
omega = 2 * np.pi * Resonance_Frequency
lambda_thz = 2*np.pi*cLight/omega
k = 2*np.pi/lambda_thz
extension = "300GHz-short"

sigma_0 = 0.00072221
emittance = 254.171e-6
distance=7e-2
distance=6e-2

beam = Beam()
beam.E = 16.0

beam.emit_x = emittance
beam.emit_y = emittance

beam.beta_x = sigma_0**2 /emittance
beam.beta_y = sigma_0**2 /emittance
beam.alpha_x = 0.0
beam.alpha_y = 0.0
beam.disp_x = 0.0
beam.disp_y = 0.0
beam.disp_dx = 0.0
beam.disp_dy = 0.0



# d_lengths = [distance, 0.9, 2, 3.2, 2.6, 3.38, 2.15, 3.61, 2.5, 2, 1.5]
# d_lengths = [distance, 0.9, 2, 2.931, 2.8546203, 3.0960211, 1.9183118, 3.0856586, 3.0067, 2.3429, 1.2]
d_lengths = [distance, 0.9, 1.240, 2.993, 1.0]

Dists = [Drift(l=length, eid=f'd{ii+1:.0f}') for ii, length in enumerate(d_lengths)]
QuadLength = 0.001
Quads = [Quadrupole(l=QuadLength, k1=0.1, eid=f'q{ii+1:.0f}') for ii in range(len(d_lengths)-1)]
# Quads = [Quadrupole(l=QuadLength, k1=0.1, eid=f'q{ii+1:.0f}') for ii in range(len(d_lengths)-1)]

m1 = Monitor(eid="start")
mSmall = Monitor(eid="Small")
mBig = Monitor(eid="Big")
mWindow = Monitor(eid="Window")
mAny = Monitor(eid="Any")

Any_Radius = 2*mirror_radiusBig

# Monitors = [m1, mSmall, mBig, mWindow]
MonitorDictionary = {m1:mirror_radiusSmall, mSmall:mirror_radiusSmall, mBig:mirror_radiusBig, mWindow:window_radius, mAny:Any_Radius}
dbzTotal = (m1,Dists[0],mAny,Quads[0],Dists[1],mSmall,Quads[1],Dists[2],mBig,Quads[2],Dists[3],mBig,Quads[3],
            Dists[4],mWindow)
varsTotal = Quads
print(varsTotal)

MSquareMax = 5.5
# MSquareMax = 10
SigmaLevel = 2
latTotal, tw0, quads = OptimizeLattice(Dists,Quads,MonitorDictionary, dbzTotal, varsTotal, beam, MSquareMax, SigmaLevel, k, InitialLattice='Parallel')

MonitorDictionary.update({mAny:mirror_radiusSmall})
print(Dists[0].l)

FixedQuads = False
FixedQuadsList = []
quads_tot = quads

lat = MagneticLattice(latTotal.sequence)
tws = twiss(lat, Twiss(beam), nPoints = 1000)

# print(max_betaBig)
PosList, fList, sizeList, betaList = PlotLattice(lat,beam, MonitorDictionary,k,MSquareMax,SigmaLevel,extension,quads_tot,FixedQuadsList=[],mWindow=mWindow)
print(PosList)
print(sizeList)


# 3 THz part
Dists[0].l = 0.7
Dists[0].l = 0.3
varsTotal = [Quads[0]]
sigma_0 = 0.000726211512749662
emittance = 25.526727393361995e-6

Resonance_Frequency = 3e12
omega = 2 * np.pi * Resonance_Frequency
lambda_thz = 2*np.pi*cLight/omega
k = 2*np.pi/lambda_thz
# M_squared = emittance / (1/(2*k))
extension = "3THz-short"

beam = Beam()
beam.E = 16.0

beam.emit_x = emittance
beam.emit_y = emittance

beam.beta_x = sigma_0**2 /emittance
beam.beta_y = sigma_0**2 /emittance
beam.alpha_x = 0.0
beam.alpha_y = 0.0
beam.disp_x = 0.0
beam.disp_y = 0.0
beam.disp_dx = 0.0
beam.disp_dy = 0.0

MonitorDictionary.update({mAny:Any_Radius})


latTotal, tw0, quads = OptimizeLattice(Dists,Quads,MonitorDictionary,dbzTotal, varsTotal, beam, MSquareMax, SigmaLevel, k, InitialLattice='FOFO',Initialize=False)

MonitorDictionary.update({mAny:mirror_radiusSmall})

lat = MagneticLattice(latTotal.sequence)
tws = twiss(lat, Twiss(beam), nPoints = 1000)

# print(max_betaBig)
PosList, fList, sizeList, betaList = PlotLattice(lat,beam, MonitorDictionary,k,MSquareMax,SigmaLevel,extension,quads_tot,FixedQuadsList=[],mWindow=mWindow)
print(PosList)


# 30 THz part
sigma_0 = 723.37e-6
emittance = 2.544e-6

Resonance_Frequency = 30e12
omega = 2 * np.pi * Resonance_Frequency
lambda_thz = 2*np.pi*cLight/omega
k = 2*np.pi/lambda_thz
# M_squared = emittance / (1/(2*k))
extension = "30THz-short"

beam = Beam()
beam.E = 16.0

beam.emit_x = emittance
beam.emit_y = emittance

beam.beta_x = sigma_0**2 /emittance
beam.beta_y = sigma_0**2 /emittance
beam.alpha_x = 0.0
beam.alpha_y = 0.0
beam.disp_x = 0.0
beam.disp_y = 0.0
beam.disp_dx = 0.0
beam.disp_dy = 0.0

# print(max_betaBig)
PosList, fList, sizeList, betaList = PlotLattice(lat,beam, MonitorDictionary,k,MSquareMax,SigmaLevel,extension,quads_tot,FixedQuadsList=[],mWindow=mWindow)
print(PosList)