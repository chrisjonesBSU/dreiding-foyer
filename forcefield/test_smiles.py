import foyer
import mbuild as mb
import forcefield_utilities as ff_utils
from gmso.parameterization import apply
from gmso import GMSOLogger

from rdkit import Chem

import warnings
warnings.filterwarnings("ignore")

logger = GMSOLogger()
logger.print_level("error")
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


from write_forcefield import write_xml 

write_xml(fname="dreiding-test.xml", exclude_metals=False, likely_hood_limit=1)


ffxml_loader = ff_utils.FoyerFFs()
dreiding_foyer = ffxml_loader.load("dreiding-test.xml")
dreiding_gmso = dreiding_foyer.to_gmso_ff()
dreiding_gmso.to_xml("dreiding-test-gmso.xml", overwrite=True)

passed_atom_typing = []
failed_atom_typing = []
for smiles in set(valid_smiles):
    comp = mb.load(smiles, smiles=True)
    gmso_top = comp.to_gmso()
    try:
        typed_comp = apply(gmso_top, dreiding_gmso, identify_connections=True, ignore_params=["dihedral", "improper"])
        passed_atom_typing.append(smiles)
    except Exception as e:
        failed_atom_typing.append(smiles)
        print(smiles)
        print(e)
        print("------------------------")


if len(failed_atom_typing) > 0:
    print("Logging failed SMILES strings to failed_smiles.txt")
    with open('failed_smiles.txt', 'w') as f:
        for smiles in set(failed_atom_typing):
            f.write(smiles + '\n')
else:
    print("All tested SMILES strings passed atom typing.")
