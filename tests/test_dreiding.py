import pytest
from pathlib import Path
import mbuild as mb
import unyt as u
import gmso
from gmso.parameterization import apply

# Temporary code
path = Path(__file__).parent / "../xmls"
import sys
sys.path.append("../forcefield")
from write_forcefield import write_xml 
write_xml(str(path), exclude_metals=False)

def print_dihedral_type(dihedral):
    print([atom.atom_type.name for atom in dihedral.connection_members])


class TestDreiding:
    def test_load_xml(self):
        # test gmso loading forcefield
        ff = gmso.ForceField(str(path/"gmso-dreiding.xml"))
        assert ff

    @pytest.fixture
    def ff(self):
        return gmso.ForceField(str(path/"gmso-dreiding.xml"))

    def test_atom_type_style(self, ff):
        atype = ff.atom_types["H_"]
        assert atype.name == "H_"
        assert atype.mass == 1.008
        assert str(atype.expression) == "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        assert atype.parameters["sigma"] == 0.2846421 * u.nm
        assert atype.parameters["epsilon"] == 0.0635968 * u.kJ/u.mol

    def test_bond_type_style(self, ff):
        btype = ff.bond_types["C_3-H_"]
        btype2 = ff.get_potential("bond_type", "H_-C_3")[0]
        assert btype is btype2
        assert btype.name == "HarmonicBondType"
        assert str(btype.expression) == "0.5*k*(r - r_eq)**2"
        assert btype.parameters["k"] == 2928.8 * u.Unit("kJ/mol/nm**2")
        assert btype.parameters["r_eq"] == 0.109 * u.nm

    def test_angle_type_style(self, ff):
        angtype = ff.angle_types["*~C_3~*"]
        angtype2 = ff.get_potential("angle_type", "*~C_3~*")[0]
        assert angtype is angtype2
        assert angtype.name == "HarmonicAngleType"
        assert str(angtype.expression) == "0.5*k*(theta - theta_eq)**2"
        assert angtype.parameters["k"] == 235.34935916425985 * u.Unit("kJ/mol/radian**2")
        assert angtype.parameters["theta_eq"] == 1.9106293854507126 * u.radian

    def test_dihedral_type_style(self, ff):
        from collections import Counter
        test_molecules = ["CC", "CC(=O)O", "C1=CC=CC=C1CN", "C=C", "C1(=CC=CC=C1)C2=CC=CC=C2", "C1=CC=CC=C1CC1=CC=CC=C1"]
        test_dihedrals = [
            {"case-a":9}, {"case-b":4, "case-j":4}, {"case-d":24, "case-a":6, "case-j":6}, {"case-c":4}, 
            {"case-d":48, "case-f":4}
        ]
        for mol, dih in zip(test_molecules, test_dihedrals):
            top = mb.load(mol, smiles=True).to_gmso()
            typed_top = apply(top, ff, identify_connections=True, ignore_params=["dihedral", "improper"]) # gmso
            dihedralsCounter = Counter([dihedral.dihedral_type.name for dihedral in typed_top.dihedrals])
            for key in dih:
                assert dihedralsCounter["PeriodicTorsion-"+key] == dih[key]

    def test_improper_type_style(self):
        pass

    def test_typing(self, ff):
        """['N[C@@H](CCCNC(=N)N)C(=O)O', 'C1C=CC=N1', 'C1=CC=NCC1', 'C1=CC(=O)N=CC1']"""
        with open("smiles_list.txt") as f:
            smiles_list = [line.strip() for line in f if line.strip()]
        passed_atom_typing = []
        failed_smiles = []
        failed_topologies = []
        errors = []
        for smiles in set(smiles_list):
            try:
                comp = mb.load(smiles, smiles=True)
            except:
                print(f"{smiles=} does not work via RDKit")
                continue
            gmso_top = comp.to_gmso()
            try:
                apply(gmso_top, ff, identify_connections=True, ignore_params=["dihedral", "improper"]) # gmso
                passed_atom_typing.append(smiles)
            except Exception as e:
                failed_smiles.append(smiles)
                failed_topologies.append(gmso_top)
                errors.append(e)
                print(smiles)
                print(e)
                print("------------------------")
        assert len(errors) == 0, failed_smiles

    def test_run_hoomd(self, ff):
        pass

    def test_hoomd_structure(self):
        pass

    # HOOMD TESTING