"""
This file contains fucntions for calculating the bond stretching and bond bending parameters.
The bond torsion parameters are not yet implemented.
These methods are used in `write_forcefield.py`.
"""

import numpy as np
import unyt as u


def equil_bond_distance(atom1, atom2):
    """
    From the paper:
    ---------------
    Rij = Ri + Rj - delta
    delta = 0.01 Ang = 0.001 nm
    """
    return atom1.get("bond_r").to("nm").value + atom2.get("bond_r").to("nm").value - 0.001


def bond_energy(atom1, atom2):
    """
    From the paper:
    ---------------
    K_ij = 700kcal/mol for all single bonds
    K_ij(n) = n(K_ij) for bonds of order n
    """
    # Determine n from atom 1 and atom 2?
    n = 1
    return (n * 700 * u.kcal).to("kJ").value
    

def equil_bond_angle(central_atom):
    """
    From the paper:
    ---------------
    Equilibrium angle is determined only by middle atom
    E_ijk = 0.5 * C_ijk[cos(theta_ijk) - cos0(theta_j)]^2
    """
    return central_atom.get("theta").to("rad").value
    

def angle_energy(central_atom):
    """
    From the paper:
    ---------------
    - E_ijk = 0.5 * C_ijk[cos(theta_ijk) - cos0(theta_j)]^2
    - C_ijk = K_ijk / sin(theta_j)^2
    - K_ijk = 100kcal/mol for all angles
    - The equilibrium angle is determined by the middle atom (j/atom2)
    
    If the angle has linear geometry: (theta_j = 180 degrees) then the functional form changes
    In this work, we will use E_ijk for linear angles as well, but replace C with K (can't divide by sin(180)^2.
    """
    theta_j = central_atom.get("theta").to("rad")
    K = (100 * u.kcal).to("kJ")
    if np.round(np.sin(theta_j), 4) != 0:
        C = K / (np.sin(theta_j)**2)
        return C.value / 2
    else:
        return K.value / 2


def equil_torsion_angle(atom1, atom2, atom3, atom4):
    """
    From the paper:
    ---------------
    E_ijkl = 0.5 * Vjk{1 - cos[n_jk(phi_ijkl - phi_jk)]}
    
    n: Periodicity
    phi_ijkl: Dihedral angle
    phi_jk: Equilibrium dihedral angle
    V_jk: Force constant

    V_jk, n and phi_jk only depend on atom 2 and aotm 3
    """
    raise NotImplementedError
