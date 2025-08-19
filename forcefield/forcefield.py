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


def create_forcefield_xml(output_file, atom_types_dict, bonds, angles):
    ForceField = ET.Element(
        "ForceField",
        name="mBuild_Dreiding",
        version="0.1.0",
        combining_rule="geometric"
    )

    # Add AtomTypes
    AtomTypes = ET.SubElement(ForceField, "AtomTypes")
    for atom_type, vals in atom_types_dict.items():
        ET.SubElement(
            AtomTypes,
            "Type",
            **{
                "name": atom_type,
                "class": atom_type,
                "element": vals["element"],
                "mass": str(vals["mass"]),
                "def": vals["_def"],
                "desc": vals["desc"],
                "doi": vals["doi"]
            }
        )
    
    # Add HarmonicBondForce
    HarmonicBondForce = ET.SubElement(ForceField, "HarmonicBondForce")
    for bond in bonds:
        atom1_dict = dreiding_atom_types[bond[0]]
        atom2_dict = dreiding_atom_types[bond[1]]
        eq_length = equil_bond_distance(atom1=atom1_dict, atom2=atom2_dict)
        bond_k = bond_energy(atom1=atom1_dict, atom2=atom2_dict)
        
        ET.SubElement(
            HarmonicBondForce,
            "Bond",
            **{
                "class1": bond[0],
                "class2": bond[1],
                "length": str(float(eq_length)),
                "k": str(float(bond_k)),
            }
        )

    # Add HarmonicAngleForce
    HarmonicAngleForce = ET.SubElement(ForceField, "HarmonicAngleForce")
    for angle in angles:
        atom2_dict = dreiding_atom_types[angle[1]]
        theta = equil_bond_angle(central_atom=atom2_dict)
        angle_k = angle_energy(central_atom=atom2_dict)

        if theta == 0:
            continue
        
        ET.SubElement(
            HarmonicAngleForce,
            "Angle",
            **{
                "class1": angle[0],
                "class2": angle[1],
                "class3": angle[2],
                "angle": str(float(theta)),
                "k": str(float(angle_k)),
            }
        )

    # Add NonbondedForce
    NonbondedForce = ET.SubElement(
        ForceField, "NonbondedForce",
        coulomb14scale="0.0",
        lj14scale="1.0"
    )
    for atom_type, vals in atom_types_dict.items():
        ET.SubElement(
            NonbondedForce,
            "Atom",
            **{
                "type": atom_type,
                "charge": "0.0",
                "sigma": str(float(vals["vdw_r"].to("nm"))),
                "epsilon": str(float(vals["D"].to("kJ"))),
            }
        )
    tree = ET.ElementTree(ForceField)
    ET.indent(tree, space="  ", level=0)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
