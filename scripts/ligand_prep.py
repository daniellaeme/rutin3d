import os
from rdkit import Chem
from rdkit.Chem import AllChem


def prep_ligand(smiles, output_path=os.path.join('.', 'data', 'processed', 'ligand_prep', 'best_rutin.sdf')):
    """
    This function takes in the smiles string, adds hydrogen molecules, performs ETKDG v3 Embedding (N=50),
    MMFF Optimization, Calculates the Energies, Sorts according to the threshold and performs Heavy Atom Pruning.
    Saves the final structure directly to the specified ligand_prep destination.
    :param smiles:
    :return: unique_confs
    """

    # Load SMILES and add Hydrogen
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('SMILES is not valid!')
    mol_h = Chem.AddHs(mol)

    # Configure ETKDGv3 parameters
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useBasicKnowledge = True     # Enforces basic chemical rules

    # Embed Conformers
    num_confs = 50
    cids = AllChem.EmbedMultipleConfs(mol_h, num_confs, params)

    # Optimise properties and Calculate Energies (Get MMFF94 Force Field)
    optm = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=1000, mmffVariant='MMFF94')     # Returns a list of (converged_status, energy) tuples
    conf_energies = []

    # Map IDs to their calculated energies
    for i, cid in enumerate(cids):
        status, energy = optm[i]                                                 # Get results for this specific conformer

        if status == 0:                                                          # Status == 0 means successfully converged
            conf_energies.append((cid, energy))
        else:
            print(f"Conformer {cid} did not converge within the iteration limit.")

    # Energy Threshold Filtering (Threshold = 5.0 kcal/mol)
    sorted_confs = sorted(conf_energies, key=lambda x: x[1])            # List sorted in ascending order based on energy value.
    if not sorted_confs:
        raise ValueError('No valid conformations found')

    min_energy = sorted_confs[0][1]                                                     # Global minimum conformation
    threshold = 5.0
    allowed_confs = [c for c in sorted_confs if c[1] - min_energy <= threshold]     # Conformers with a relative energy to the global minimum conformation <= 5.0 kcal/mol

    print(f'{len(allowed_confs)} conformations <= {threshold} found')

    # RMSD Pruning
    unique_confs = []
    if allowed_confs:
        unique_confs.append(allowed_confs[0][0])                            # Add lowest energy

        for i in range(1, len(allowed_confs)):
            cid_test = allowed_confs[i][0]
            is_unique = True

            # Calculate RMSD
            for cid_ref in unique_confs:
                rmsd = AllChem.GetBestRMS(mol_h, mol_h, cid_test, cid_ref)  # Accounts for chemical symmetries, understands compounds existing with the same spatial shape and prevents incorrect calculations.
                if rmsd < 1.0:                                              # If the calculated RMSD to all accepted reference conformers is greater than the threshold (1.0 Angstrom), the conformer is structurally unique.
                    is_unique = False
                    break

            if is_unique:
                unique_confs.append(cid_test)

    print(f'Found {len(unique_confs)} unique conformations')

    # Save the best conformer
    best_conf_id = unique_confs[0]  # Retrieve the conformer with the absolute lowest energy after RMSD pruning.
    print(f'Best conformer ID: {best_conf_id}')

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if os.path.exists(output_path):
        print(f'f"File "{output_path}" already exists! Skipping save to prevent overwrite."')
    else:
        writer = Chem.SDWriter(output_path)
        writer.write(mol_h, best_conf_id)
        writer.close()
        print(f'Best conformer saved to {output_path}!')

    return unique_confs
