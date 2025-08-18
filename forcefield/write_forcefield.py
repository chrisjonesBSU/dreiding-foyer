import unyt as u
import ele
import numpy as np

import xml.etree.ElementTree as ET

from atom_types import dreiding_atom_types
from bonding_rules import bonding_rules
from forcefield import (
        equil_bond_distance,
        bond_energy,
        equil_bond_angle,
        angle_energy
)


def write_xml(fname="dreiding.xml", exclude_metals=True, likely_hood_limit=1):
    all_bonds = _get_all_bonds(exclude_metals=exclude_metals)
    all_angles = _get_all_angles(likely_hood_limit=likely_hood_limit, exclude_metals=exclude_metals)

    _create_forcefield_xml(
        output_file=fname,
        atom_types_dict=dreiding_atom_types,
        bonds=all_bonds,
        angles=all_angles
    )


def _create_forcefield_xml(output_file, atom_types_dict, bonds, angles):
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
                "overrides": vals["overrides"],
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


def _parse_bond_types(atom_type, exclude_metals):
    this_atom_dict = dreiding_atom_types[atom_type]
    element = dreiding_atom_types[atom_type]["element"]
    allowed_bonds = bonding_rules[element]["allowed_bonds"]
    
    bonds_list = set()
    for bonding_element in allowed_bonds:
        for other_atom_type, vals in dreiding_atom_types.items():
            if exclude_metals and vals["is_metal"]:
                continue
            if vals["element"] == bonding_element:
                sorted_bonds = sorted([atom_type, other_atom_type])
                bonds_list.add((sorted_bonds[0], sorted_bonds[1]))
    return bonds_list


def _get_all_bonds(exclude_metals=True):
    master_bonds_list = set()
    for atom_type in dreiding_atom_types:
        if exclude_metals and dreiding_atom_types[atom_type]["is_metal"]:
            continue
        atom_type_bonds = _parse_bond_types(atom_type, exclude_metals=exclude_metals)
        for bond_type in atom_type_bonds:
            master_bonds_list.add(bond_type)
    return master_bonds_list


def _parse_angle_types(central_atom, exclude_metals):
    """Returns all possible angles where central_atom is the middle atom"""
    central_dict = dreiding_atom_types[central_atom]
    central_element = central_dict["element"]
    allowed_bonds = bonding_rules[central_element]["allowed_bonds"]
    # Get all atom types that can bond to the central atom
    bonded_atom_types = []
    for element in allowed_bonds:
        for atom_type, vals in dreiding_atom_types.items():
            if vals["element"] == element:
                if exclude_metals and vals["is_metal"]:
                    pass
                else:
                    bonded_atom_types.append(atom_type)
    # All unique pairs of bonded atom types for angles
    angles = set()
    for i, atom1 in enumerate(bonded_atom_types):
        for atom2 in bonded_atom_types[i:]:  # include i: to allow same atom type twice
            # Sort the outer atoms so (A,B,C) == (B,A,C)
            sorted_outer = sorted([atom1, atom2])
            angles.add((sorted_outer[0], central_atom, sorted_outer[1]))
    return angles


def _get_all_angles(likely_hood_limit=1, exclude_metals=True):
    master_angles_list = set()
    for atom_type, vals in dreiding_atom_types.items():
        if exclude_metals and vals["is_metal"]:
            continue
        if bonding_rules[vals["element"]]["middle_bond"] >= likely_hood_limit:
            atom_type_angles = _parse_angle_types(atom_type, exclude_metals=exclude_metals)
            for angle_type in atom_type_angles:
                master_angles_list.add(angle_type)
    return master_angles_list
