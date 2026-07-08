from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def docking_prep(clean_receptor_path, center_x, center_y, center_z):
    mol = Chem.SDMolSupplier('././data/processed/best_rutin.sdf', removeHs=False)
    mk_prep = MoleculePreparation()
    molsetup_list = mk_prep(mol)

    molsetup = molsetup_list[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)

