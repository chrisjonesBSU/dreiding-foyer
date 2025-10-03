import unyt as u
import ele
import numpy as np

import xml.etree.ElementTree as ET

from bonding_rules import bonding_rules
from forcefield import (
        equil_bond_distance,
        bond_energy,
        equil_bond_angle,
        angle_energy,
        equil_torsion_angle
)
def _get_all_atomtypes(atom_typesDict):
    """Add to dictionary sigma parameter from Ro
    
    Notes
    -----
    See Here for details on Dreiding Sigma conversion
    https://docs.lammps.org/pair_lj.html#coefficients
    """
    for key in atom_typesDict:
        # D0 = epsilon
        atom_typesDict[key]["sigma"] = atom_typesDict[key]["vdw_r"]*2**(-1/6)
    return atom_typesDict

def write_xml(fpath, exclude_metals=True, likely_hood_limit=1):
    from atom_types import dreiding_atom_types
    all_atom_types = _get_all_atomtypes(dreiding_atom_types)
    all_bonds = _get_all_bonds(all_atom_types, exclude_metals=exclude_metals)
    all_angles = _get_all_angles(all_atom_types, likely_hood_limit=likely_hood_limit, exclude_metals=exclude_metals)
    all_dihedrals = _get_all_dihedrals(all_atom_types, likely_hood_limit=likely_hood_limit, exclude_metals=exclude_metals)

    _create_forcefield_foyer_xml(
        output_file=fpath+"/foyer-dreiding.xml",
        atom_typesDict=all_atom_types,
        bonds=all_bonds,
        angles=all_angles
    )
    _create_forcefield_gmso_xml(
        output_file=fpath+"/gmso-dreiding.xml",
        atom_typesDict=all_atom_types,
        bonds=all_bonds,
        angles=all_angles,
        dihedrals=all_dihedrals,
    )


def _create_forcefield_foyer_xml(output_file, atom_typesDict, bonds, angles):
    ForceField = ET.Element(
        "ForceField",
        name="mBuild_Dreiding",
        version="0.1.0",
        combining_rule="geometric"
    )

    # Add AtomTypes
    AtomTypes = ET.SubElement(ForceField, "AtomTypes")
    for atom_type, vals in atom_typesDict.items():
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
        atom1_dict = atom_typesDict[bond[0]]
        atom2_dict = atom_typesDict[bond[2]]
        eq_length = equil_bond_distance(atom1=atom1_dict, atom2=atom2_dict)
        bond_k = bond_energy(bond[1])
        
        ET.SubElement(
            HarmonicBondForce,
            "Bond",
            **{
                "class1": bond[0],
                "class2": bond[2],
                "length": str(float(eq_length)),
                "k": str(float(bond_k)),
            }
        )

    # Add HarmonicAngleForce
    HarmonicAngleForce = ET.SubElement(ForceField, "HarmonicAngleForce")
    for angle in angles:
        atom2_dict = atom_typesDict[angle[1]]
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
    for atom_type, vals in atom_typesDict.items():
        ET.SubElement(
            NonbondedForce,
            "Atom",
            **{
                "type": atom_type,
                "charge": "0.0",
                "sigma": str(float(vals["sigma"].to("nm"))),
                "epsilon": str(float(vals["D"].to("kJ"))),
            }
        )
    tree = ET.ElementTree(ForceField)
    ET.indent(tree, space="  ", level=0)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)

def _create_forcefield_gmso_xml(output_file, atom_typesDict, bonds, angles, dihedrals):
    ForceField = ET.Element(
        "ForceField",
        name="GMSO_Dreiding",
        version="0.1.0",
    )
    # Add metadata
    Metadata = ET.SubElement(ForceField, "FFMetaData", combining_rule="geometric", electrostatics14Scale="0.5", nonBonded14Scale="0.5")
    ET.SubElement(Metadata, "Units", energy="kJ", distance="nm", mass="amu", charge="elementary_charge")

    # Add AtomTypes
    AtomTypes = ET.SubElement(ForceField, "AtomTypes", expression="4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)")
    ET.SubElement(AtomTypes, "ParametersUnitDef", **{"parameter": "epsilon", "unit":"kJ/mol"})
    ET.SubElement(AtomTypes, "ParametersUnitDef", **{"parameter": "sigma", "unit":"nm"})
    # <ParametersUnitDef parameter="epsilon" unit="kJ/mol"/>
    # <ParametersUnitDef parameter="sigma" unit="nm"/>
    for atom_type, vals in  atom_typesDict.items():
        atype = ET.SubElement(
            AtomTypes,
            "AtomType",
            **{
                "name": atom_type,
                "atomclass": atom_type,
                "element": vals["element"],
                "mass": str(vals["mass"]),
                "definition": vals["_def"],
                "description": vals["desc"],
                "overrides": vals["overrides"],
                "doi": vals["doi"],
            }
        )
        params = ET.SubElement(atype, "Parameters")
        ET.SubElement(params, "Parameter", **{"name":"sigma", "value":str(round(float(vals["sigma"].to("nm")), 7))})
        ET.SubElement(params, "Parameter", **{"name":"epsilon", "value":str(round(float(vals["D"].to("kJ")), 9))})
    
    # Add HarmonicBondForce
    HarmonicBondForce = ET.SubElement(ForceField, "BondTypes", expression="0.5*k*(r - r_eq)**2")
    ET.SubElement(HarmonicBondForce, "ParametersUnitDef", **{"parameter": "k", "unit":"kJ/mol/nm**2"})
    ET.SubElement(HarmonicBondForce, "ParametersUnitDef", **{"parameter": "r_eq", "unit":"nm"})
    bkeys = set()
    for bond in bonds:
        atom1_dict = atom_typesDict[bond[0]]
        atom2_dict = atom_typesDict[bond[2]]
        bond_order = bond[1]
        eq_length = equil_bond_distance(atom1=atom1_dict, atom2=atom2_dict)
        bond_k = bond_energy(bond_order)
        
        btype = ET.SubElement(
            HarmonicBondForce,
            "BondType",
            name="HarmonicBondType",
            classes = "".join(bond),
        )

        params = ET.SubElement(btype, "Parameters") #TODO: Check units of k and r_eq
        ET.SubElement(params, "Parameter", name="k", value=str(float(bond_k)))
        ET.SubElement(params, "Parameter", name="r_eq", value=str(float(eq_length)))

    # Add HarmonicAngleForce
    # TODO: Check angle units
    HarmonicAngleForce = ET.SubElement(ForceField, "AngleTypes", expression="0.5*k*(theta - theta_eq)**2")
    ET.SubElement(HarmonicAngleForce, "ParametersUnitDef", **{"parameter": "k", "unit":"kJ/mol/radian**2"})
    ET.SubElement(HarmonicAngleForce, "ParametersUnitDef", **{"parameter": "theta_eq", "unit":"radian"})
    for angle in angles:
        atom2_dict = atom_typesDict[angle[1]]
        theta = equil_bond_angle(central_atom=atom2_dict)
        angle_k = angle_energy(central_atom=atom2_dict)

        if theta == 0:
            continue
        angletype = ET.SubElement(
            HarmonicAngleForce,
            "AngleType",
            name="HarmonicAngleType",
            classes= "*~"+angle[1]+"~*", 
        )

        params = ET.SubElement(angletype, "Parameters") #TODO: Check units of k and theta_eq
        ET.SubElement(params, "Parameter", name="k", value=str(float(angle_k)))
        ET.SubElement(params, "Parameter", name="theta_eq", value=str(float(theta)))

    # Add PeriodicTorsionForce
    PeriodicTorsionForce = ET.SubElement(ForceField, "DihedralTypes", expression="k * (1 + cos(n * phi - phi_eq))")
    ET.SubElement(PeriodicTorsionForce, "ParametersUnitDef", **{"parameter": "k", "unit":"kJ/mol"})
    ET.SubElement(PeriodicTorsionForce, "ParametersUnitDef", **{"parameter": "n", "unit":"dimensionless"})
    ET.SubElement(PeriodicTorsionForce, "ParametersUnitDef", **{"parameter": "phi_eq", "unit":"degree"})
    for dihedral in dihedrals:
        paramsDict = equil_torsion_angle(*dihedral)
        if not paramsDict: # assume torsion energy is 0
            continue
        classes = "".join(dihedral)
        
        btype = ET.SubElement(
            PeriodicTorsionForce,
            "DihedralType",
            name="PeriodicTorsion-"+paramsDict["name"],
            classes=classes
        )

        params = ET.SubElement(btype, "Parameters") #TODO: Check units of k and r_eq
        ET.SubElement(params, "Parameter", name="k", value=str(float(paramsDict["V_jk"])))
        ET.SubElement(params, "Parameter", name="n", value=str(float(paramsDict["n"])))
        ET.SubElement(params, "Parameter", name="phi_eq", value=str(float(paramsDict["phi_jk"])))
    
    tree = ET.ElementTree(ForceField)
    ET.indent(tree, space="  ", level=0)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)

def _parse_bond_types(atom_typesDict, atom_type, exclude_metals):
    element = atom_typesDict[atom_type]["element"]
    allowed_bonds = bonding_rules[element]["allowed_bonds"]
    
    bonds_list = set()
    for bonding_element in allowed_bonds:
        for other_atom_type, vals in atom_typesDict.items():
            if exclude_metals and vals["is_metal"]:
                continue
            if vals["element"] == bonding_element:
                sorted_bonds = sorted([atom_type, other_atom_type])
                for border in ("-", "=", "#", ":"):
                    bonds_list.add((sorted_bonds[0], border, sorted_bonds[1]))
    return bonds_list


def _get_all_bonds(atom_typesDict, exclude_metals=True):
    master_bonds_list = set()
    for atom_type in atom_typesDict:
        if exclude_metals and atom_typesDict[atom_type]["is_metal"]:
            continue
        atom_type_bonds = _parse_bond_types(atom_typesDict, atom_type, exclude_metals=exclude_metals)
        for bond_type in atom_type_bonds:
            master_bonds_list.add(bond_type)
    return master_bonds_list


def _parse_angle_types(central_atom, exclude_metals):
    """Returns all possible angles where central_atom is the middle atom"""
    angles = set((("*", central_atom, "*", "~", "~"),)) # for now only use wildcards
    return angles


def _get_all_angles(atom_typesDict, likely_hood_limit=1, exclude_metals=True):
    master_angles_list = set()
    for atom_type, vals in atom_typesDict.items():
        if exclude_metals and vals["is_metal"]:
            continue
        if bonding_rules[vals["element"]]["middle_bond"] >= likely_hood_limit:
            atom_type_angles = _parse_angle_types(atom_type, exclude_metals=exclude_metals)
            for angle_type in atom_type_angles:
                master_angles_list.add(angle_type)
    return master_angles_list

def _get_all_dihedrals(atom_typesDict, likely_hood_limit=1, exclude_metals=True):
    master_dihedrals_list = set()
    for atom_type in atom_typesDict:
        if exclude_metals and atom_typesDict[atom_type]["is_metal"]:
            continue
        atom_type_bonds = _parse_bond_types(atom_typesDict, atom_type, exclude_metals=exclude_metals) # all sets of middle bond pairs
        for bond_type in atom_type_bonds:
            atomj, bondjk, atomk = bond_type
            master_dihedrals_list.add(("*", "~", atomj, bondjk, atomk, "~", "*"))
    return master_dihedrals_list


if __name__ == "__main__":
    from pathlib import Path
    import mbuild as mb
    import gmso
    from gmso.parameterization import apply


    path = Path(__file__).parent / "../xmls"
    
    write_xml(str(path))

    # test gmso loading forcefield
    ff = gmso.ForceField(str(path/"gmso-dreiding.xml"))

    # DIHEDRALS TESTING
    from collections import Counter
    test_molecules = ["CC", "CC(=O)O", "C1CCCCC1CN", "C=C", "C1(=CC=CC=C1)C2=CC=CC=C2", "C1=CC=CC=C1CC1=CC=CC=C1"]
    test_dihedrals = [{"case-a":9}, {"case-b":1, "case-c":1}]
    for mol, dih in zip(test_molecules, test_dihedrals):
        top = mb.load(mol, smiles=True).to_gmso()
        typed_top = apply(top, ff, identify_connections=True, ignore_params=["dihedral", "improper"]) # gmso
        dihedralsCounter = Counter([dihedral.dihedral_type.name for dihedral in typed_top.dihedrals])
        for key in dih:
            assert dihedralsCounter["PeriodicTorsion-"+key] == dih[key]

    # HOOMD TESTING
