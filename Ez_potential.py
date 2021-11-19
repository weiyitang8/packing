#!/usr/bin/env python

#Implementation from 
#Ez, a Depth-dependent Potential for Assessing the Energies of Insertion of Amino Acid Side-chains into Membranes: Derivation and Applications to Determining the Orientation of Transmembrane and Interfacial Helices
#
#http://www.sciencedirect.com/science/article/pii/S0022283606012095
#

import numpy as np
from prody import *


#parameters
#equation 2 from paper
Ez_potential = {}
Ez_potential['ALA'] = [0, -0.29, 10.22, 4.67]
Ez_potential['ASP'] = [0, 1.19, 14.25, 8.98]
Ez_potential['GLU'] = [0, 1.30, 14.66, 4.16]
Ez_potential['PHE'] = [0, -0.65, 19.67, 7.12]   # mm mod, since this value looks to high and I only get TRP 
Ez_potential['GLY'] = [0, -0.01, 13.86, 6.00]
Ez_potential['HIS'] = [0, 0.75, 12.26, 2.77]
Ez_potential['ILE'] = [0, -0.56, 14.34, 10.69]
Ez_potential['LYS'] = [0, 1.66, 11.11, 2.09]
Ez_potential['LEU'] = [0, -0.64, 17.34, 8.61]
Ez_potential['MET'] = [0, -0.28, 18.04, 7.13]
Ez_potential['ASN'] = [0, 0.89, 12.78, 6.28]
Ez_potential['PRO'] = [0, 0.83, 18.09, 3.53]
Ez_potential['GLN'] = [0, 1.21, 10.46, 2.59]
Ez_potential['ARG'] = [0, 1.55, 9.34, 4.68]
Ez_potential['SER'] = [0, 0.10, 13.86, 6.00]
Ez_potential['THR'] = [0, 0.01, 13.86, 6.00]
Ez_potential['VAL'] = [0, -0.47, 11.35, 4.97]
#equation 3
Ez_potential['TRP'] = [1, -0.65, 11.65, 7.20] # mm mod to -.65, since this value looks to high and I only get TRP 
Ez_potential['TYR'] = [1, -0.42, 13.04, 6.20]

def compute_Ez_energy(positions):
    """
    Parameters:
        positions: list of residue name (three letter code), position in membrane
    Returns:
        Ez: Ez energy
    """
    Ez = 0
    for r,z in positions:
        Ez += compute_Ez(r,z)
    return Ez

def compute_Ez_energy_multistate(positions):
    Ez = 0
    Ez1 = 0
    Ez2 = 0

    for r,z in positions:
        Ez += compute_Ez(r,z)
        Ez1 += compute_Ez(r,z + 5)
        Ez2 += compute_Ez(r,z - 5)
    return Ez,Ez1,Ez2

def compute_Ez(r, z):
    if r not in Ez_potential:
        return 0
    equation,E0,zm,n = Ez_potential[r]

    if equation == 0:
        dEz =  E0/(1+(z/zm)**n)
    elif equation == 1:
        dEz = E0*np.exp(-(z-zm)**2/2/n**2)

    return dEz
