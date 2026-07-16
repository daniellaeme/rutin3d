'''
This master script coordinates Stage 1 to Stage 4.
It enforces directory structures, processes geometries, executes Vina,
and generates the visual PDB complex automatically.
'''


import os
from scripts.ligand_prep import prep_ligand
from scripts.receptor_prep import prep_receptor
from scripts.docking import prep_ligand_pdbqt, prep_receptor_pdbqt, run_vina_docking
from scripts.generate_complex import extract_best_pose_and_merge


if __name__ == '__main__':
    # Initialize workspace structured subdirectories
    dir_ligand_prep = os.path.join('.', 'data', 'processed', 'ligand_prep')
    dir_receptor_prep = os.path.join('.', 'data', 'processed', 'receptor_prep')
    dir_docking_prep = os.path.join('.', 'data', 'processed', 'docking')
    dir_results = os.path.join('.', 'results')

    # Enforce directory creation and that it matches active workspace exactly
    for directory in [dir_ligand_prep, dir_receptor_prep, dir_docking_prep, dir_results]:
        os.makedirs(directory, exist_ok=True)

    # 1
    # Ligand Preparation - Pubchem CID 5280805 - Rutin
    smiles = 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O'
    sdf_output = os.path.join(dir_ligand_prep, "best_rutin.sdf")

    try:
        # Generate conformers and save directly to ligand_prep directory
        prep_ligand(smiles, output_path=sdf_output)
        print('\nSuccessfully prepared ligand!')
    except Exception as e:
        print(f'Error: {e}')

    # 2
    # Receptor Preparation
    raw_pdb_path = os.path.join('.', 'data', 'raw', '4AL0.pdb')

    # Check for chemical entities
    try:
        with open(raw_pdb_path, 'r') as f:
            het_names = set()
            for line in f:
                if line.startswith('HETSYN') or line.startswith('HETNAM'):
                    parts = line.split()
                    if len(parts) >= 2:
                        het_names.add(parts[1])
            print('Found chemical entities in file:', het_names)
    except FileNotFoundError:
        print(f'\nError: "{raw_pdb_path}" not found! Please download it from rcsb.org/structure/4AL0')
        exit(1)

    protein_name = '4AL0'
    native_ligand = 'GSH'
    clean_receptor_dest = os.path.join(dir_receptor_prep, 'receptor_clean.pdb')

    try:
        # Process and write clean receptor straight to receptor_prep subdirectory
        clean_path, (cx, cy, cz) = prep_receptor(
            protein_name=protein_name,
            file_path=raw_pdb_path,
            target_ligand_name=native_ligand,
            output_path=clean_receptor_dest
        )

        print('\nSuccessfully prepared receptor!')
        print(f'Receptor File: {clean_path}')
        print(f'AutoDock Vina Grid Parameters:')
        print(f'  --center_x {cx:.3f}')
        print(f'  --center_y {cy:.3f}')
        print(f'  --center_z {cz:.3f}')
        print(f'  --size_x 22.0 --size_y 22.0 --size_z 22.0')

        # Define precise inputs and destinations within respective processed folders
        lig_pdbqt_dest = os.path.join(dir_docking_prep, 'rutin_prepared.pdbqt')
        rec_pdbqt_dest = os.path.join(dir_docking_prep, 'receptor_prepared.pdbqt')
        docked_poses_dest = os.path.join(dir_docking_prep, 'rutin_docked_poses.pdbqt')
        docking_log_dest = os.path.join(dir_results, 'vina_docking.log')
        merged_complex_dest = os.path.join(dir_results, 'rutin_mpges1_complex.pdb')

        # 3
        # Convert files to PDBQT format
        lig_pdbqt = prep_ligand_pdbqt(sdf_output, lig_pdbqt_dest)
        rec_pdbqt = prep_receptor_pdbqt(clean_path, rec_pdbqt_dest)

        # 4
        # Run Docking Simulation using calculated coordinates
        # Both outputs (docked poses and vina logs) are kept inside docking_prep
        best_score = run_vina_docking(
            receptor_pdbqt=rec_pdbqt,
            ligand_pdbqt=lig_pdbqt,
            center_coords=[cx, cy, cz],
            output_poses_path=docked_poses_dest,
            log_path=docking_log_dest
        )

        print(f'\nSimulation successful. Top binding score: {best_score:.2f} kcal/mol')
        print(f'Results stored under: {dir_docking_prep}')

        extract_best_pose_and_merge(
            receptor_path=clean_path,
            docked_pdbqt_path=docked_poses_dest,
            output_complex_pdb=merged_complex_dest
        )
        print(f'Operational Complex saved directly to: {merged_complex_dest}')
        print('Standalone PML script ready. Load via PyMOL command line: @style_complex.pml')

    except Exception as e:
        print(f'\nPipeline failed during execution: {str(e)}')