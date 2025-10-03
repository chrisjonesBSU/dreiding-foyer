"""
This file contains fucntions for calculating the bond stretching and bond bending parameters.
The bond torsion parameters are not yet implemented.
These methods are used in `write_forcefield.py`.
"""

import numpy as np
import unyt as u
from atom_types import dreiding_atom_types

halogens = ["F", "Cl", "Br", "I"]


def equil_bond_distance(atom1, atom2):
    """
    From the paper:
    ---------------
    Rij = Ri + Rj - delta
    delta = 0.01 Ang = 0.001 nm
    """
    return atom1.get("bond_r").to("nm").value + atom2.get("bond_r").to("nm").value - 0.001


borderDict = {"-":1, "=":2, "#":3,":":1.5}
def bond_energy(bond_order):
    """
    From the paper:
    ---------------
    K_ij = 700kcal/mol for all single bonds
    K_ij(n) = n(K_ij) for bonds of order n
    n = bond_order
    """
    if isinstance(bond_order, str):
        bond_order = borderDict[bond_order]
    return (bond_order * 700 * u.kcal).to("kJ").value

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
    K = (100 * u.kcal).to("kJ") # kcal/mol/rad**2
    if np.round(np.sin(theta_j), 4) != 0:
        C = K / (np.sin(theta_j)**2)
        return C.value / 2
    else:
        return K.value / 2
    
def harmonic_cos_angle_energy():
    """
    - E_ijk = 0.5 * C_ijk[cos(theta_ijk) - cos0(theta_j)]^2
    - C_ijk = K_ijk / sin(theta_j)^2
    - K_ijk = 100kcal/mol for all angles

    Problem with using harmonic_cos form is when theta_j=0. This can 
    be handled in LAMMPS https://docs.lammps.org/angle_cosine.html
    but not HOOMD-Blue or GROMACS.
    """
    pass

def equil_torsion_angle(atomi, bondij, atomj, bondjk, atomk, bondkl, atoml):
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
    # TODO: validate fix if statements
    if atomj[-2:] == "_3" and atomk[-2:] == "_3" and bondjk == "-": # case a
        return {"n":3*u.dimensionless, "V_jk":2*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-a"}
    elif not atomi[-2:] in ["_2", "_R"]  and atomj[:-2] in ["_2", "_R"] and atomk[-2:] == "_3" and bondjk == "-": # case j exception to case b
        return {"n":3*u.dimensionless, "V_jk":2*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-j"}
    elif atomj[-2:] in ["_2", "_R"] and atomk[-2:] == "_3" and bondjk == "-": # case b
        return {"n":6*u.dimensionless, "V_jk":1*u.kcal/u.mol, "phi_jk":0*u.degree, "name":"case-b"}
    elif atomj[-2:] == "_2" and atomk[-2:] == "_2" and bondjk == "=": # case c
        return {"n":2*u.dimensionless, "V_jk":45*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-c"}
    elif atomj[-2:] == "_R" and atomk[-2:] == "_R" and bondjk==":": # case d
        return {"n":2*u.dimensionless, "V_jk":25*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-d"}
    elif atomj[-2:] == "_R" and atomk[-2:] == "_R" and bondjk == "-": # casef checks before casee
        return {"n":2*u.dimensionless, "V_jk":10*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-f"}
    elif (atomj[-2:] == "_2" or atomj[-2:] == "_R") and (atomk[-2:] == "_2" or atomk[-2:] == "_R") and bondjk == "-": # case e
        return {"n":2*u.dimensionless, "V_jk":5*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-e"}
    elif (atomj[-2:] == "_1" or dreiding_atom_types[atomj].get("is_metal") or dreiding_atom_types[atomj].get("element") in halogens) or (
        atomk[-2:] == "_1" or dreiding_atom_types[atomk].get("is_metal") or dreiding_atom_types[atomk].get("element") in halogens
        ): # case g
        return None
    elif atomj[-2:] == "_3" and atomk[-2:] == "_3" and atomj.split("_")[0] in ["O", "S", "Se"] and atomk.split("_")[0] in ["O", "S", "Se"]: # case h
        return {"n":2*u.dimensionless, "V_jk":2*u.kcal/u.mol, "phi_jk":90*u.degree, "name":"case-h"}
    elif atomj[-2:] == "_3" and atomj.split("_")[0] in ["O", "S", "Se"] and atomk[-2:] in ["_R", "_2"]: # case i
        return {"n":2*u.dimensionless, "V_jk":5*u.kcal/u.mol, "phi_jk":180*u.degree, "name":"case-i"}
    print(f"NotImplementedError{(atomi, atomj, atomk, atoml, bondij, bondjk, bondkl)}")
    return None
