from scripts.ligand_prep import ligand_prep
from scripts.receptor_prep import receptor_prep

if __name__ == '__main__':

    # Ligand Preparation - Pubchem CID 5280805 - Rutin
    smiles = 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O'
    # try:
    #     ligand_prep(smiles)
    #     print('\n Successfully prepared ligand!')
    # except ValueError:
    #     print(f'Error: {smiles} invalid.')


    # Receptor Preparation

    with open("4AL0.pdb", "r") as f:
        het_names = set()
        for line in f:
            if line.startswith("HETSYN"):
                parts = line.split()
                if len(parts) >= 2:
                    het_names.add(parts[1])
            elif line.startswith("HETNAM"):
                parts = line.split()
                if len(parts) >= 2:
                    het_names.add(parts[1])
        print("Found chemical entities in file:", het_names)

    protein_name = '4AL0'
    file_path = './4AL0.pdb'
    native_ligand = 'GSH'
    try:
        clean_path, (cx, cy, cz) = receptor_prep(protein_name, file_path, target_ligand_name=native_ligand)
        print('\nSuccessfully prepared receptor!')
        print(f"Receptor File: {clean_path}")
        print(f"AutoDock Vina Grid Parameters:")
        print(f"  --center_x {cx:.3f}")
        print(f"  --center_y {cy:.3f}")
        print(f"  --center_z {cz:.3f}")
        print(f"  --size_x 22.0 --size_y 22.0 --size_z 22.0")
    except FileNotFoundError:
        print(f"\nError: '{file_path}' not found! Please download it from rcsb.org/structure/4AL0")
