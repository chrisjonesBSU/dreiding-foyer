import foyer
import mbuild as mb

from rdkit import Chem

# Check each SMILES string
with open("smiles_list.txt") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

valid_smiles = []
failed = []
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        try:
            mb.load(smiles, smiles=True)
            valid_smiles.append(smiles)
        except:
            failed.append(smiles)
    else:
        failed.append(smiles)


from ../forcefield/write_forcefield.py import write_xml 

write_xml(fname="dreiding-test.xml", exclude_metals=False, likely_hood_limit=1)

passed_atom_typing = []
failed_atom_typing = []
dreiding_ff = foyer.Forcefield("dreiding-test.xml")
for smiles in set(valid_smiles):
    comp = mb.load(smiles, smiles=True)
    try:
        typed_comp = dreiding_ff.apply(comp, assert_dihedral_params=False)
        passed_atom_typing.append(smiles)
    except Exception as e:
        print(smiles)
        print(e)
        failed_atom_typing.append(smiles)
        print()


with open('failed_smiles.txt', 'w') as f:
    for smiles in set(failed_atom_typing):
        f.write(smiles + '\n')
