import os
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import openbabel
import subprocess
# from vina import Vina

def prep_ligand_pdbqt(input_sdf_path, output_pdbqt_path):
    """
    Reads a 3D SDF ligand, processes its torsional tree with Meeko,
    calculates partial charges, and exports a valid PDBQT file.
    """
    print("\nStarting Ligand PDBQT Preparation using Meeko")

    if not os.path.exists(input_sdf_path):
        raise FileNotFoundError(f"Error: Input SDF file not found at {input_sdf_path}")

    supplier = Chem.SDMolSupplier(input_sdf_path, removeHs=False)
    if len(supplier) == 0 or supplier[0] is None:
        raise ValueError("Error: Could not parse a valid molecule from the SDF file.")

    mol = supplier[0]
    print(f"[1/3] Successfully loaded 3D ligand from: {input_sdf_path}")

    # Initialize the Meeko preparer
    preparator = MoleculePreparation()
    prepared_setups = preparator.prepare(mol)
    mol_setup = prepared_setups[0]
    print("[2/3] Torsional tree mapped. Non-polar hydrogens merged.")

    # Write coordinates to PDBQT format
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setup)
    if not is_ok:
        raise RuntimeError(f"Meeko writing failed: {error_msg}")


    # Safely save to file without overwriting your original SDF
    with open(output_pdbqt_path, 'w') as f:
        f.write(pdbqt_string)

    print(f"[3/3] Prepared ligand saved to: {output_pdbqt_path}")

    return output_pdbqt_path


def prep_receptor_pdbqt(input_pdb_path, output_pdbqt_path):
    """
    Uses OpenBabel to add polar hydrogens and convert a PDB receptor to PDBQT format.
    """
    print("\nStarting Receptor PDBQT Preparation using OpenBabel")

    if not os.path.exists(input_pdb_path):
        raise FileNotFoundError(f"Error: Input PDB file not found at {input_pdb_path}")

    # Initialize the OpenBabel conversion engine
    ob_conv = openbabel.OBConversion()
    ob_conv.SetInAndOutFormats('pdb', 'pdbqt')
    mol = openbabel.OBMol()

    # Read the clean PDB file
    if not ob_conv.ReadFile(mol, input_pdb_path):
        raise IOError(f"Could not read input file: {input_pdb_path}")

    # Add polar hydrogens ONLY (True, False, 7.4)
    # This aligns the receptor with Vina's united-atom scoring parameters
    mol.AddHydrogens(True, False, 7.4)
    if not ob_conv.WriteFile(mol, output_pdbqt_path):
        raise IOError(f"Could not write output file: {output_pdbqt_path}")
    print("[1/2] Polar hydrogens added to receptor side-chains.")

    # Write coordinates to PDBQT format
    if not ob_conv.WriteFile(mol, output_pdbqt_path):
        raise IOError(f"Could not write prepared receptor file: {output_pdbqt_path}")

    # READ the file OpenBabel created
    with open(output_pdbqt_path, "r") as f:
        lines = f.readlines()

    # STRIP out ligand specific tags that break rigid receptors
    bad_tags = ["ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"]
    clean_lines = [line for line in lines if not any(tag in line for tag in bad_tags)]

    # WRITE the clean file back out
    with open(output_pdbqt_path, "w") as f:
        f.writelines(clean_lines)
    print(f"[2/2] Prepared receptor stripped of invalid tags and saved to: {output_pdbqt_path}")

    return output_pdbqt_path


def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center_coords, box_size=22.0):
    """
    Runs the standalone Vina executable using Python's subprocess module
    and extracts the binding affinity energy table.
    """
    print("\nStarting AutoDock Vina Standalone Executable")

    output_poses_path = 'rutin_docked_poses.pdbqt'
    log_path = 'vina_docking.log'

    # Build the exact command line arguments for your Windows .\vina.exe
    cmd = [
        r".\vina.exe",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center_coords[0]),
        "--center_y", str(center_coords[1]),
        "--center_z", str(center_coords[2]),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", "32",
        "--num_modes", "9",
        "--out", output_poses_path
    ]

    print("[1/4] Command line arguments structured.")
    print(f"[2/4] Bounding box centered at: {center_coords} with size {box_size}")
    print("[*] Launching external Vina process. Exhaustiveness set to 32...")

    # Run the compiled Vina binary and redirect outputs to a log file
    with open(log_path, 'w', encoding='utf-8') as log_file:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=0)
        # result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, text=True) #  capture_output=True,

        while True:
            char = process.stdout.read(1)
            if not char:
                break  # Process finished streaming

            print(char, end="", flush=True)  # Print immediately without waiting for \n
            log_file.write(char)  # Write to backup log file

        process.wait()
    # Re-read the log to show the output in terminal so it doesn't look frozen
    # with open(log_path, "r") as log_file:
    #     print(log_file.read())

    if process.returncode != 0:
        raise RuntimeError(f"Vina failed with error:\n{process.stderr}")

    print("[3/4] Simulation finished successfully.")
    print("[4/4] Output poses written to disk.")

    # Parse the log file to extract the affinity table and the best score
    print("\nDOCKING RESULTS")
    best_score = 0.0

    # Read from the log file
    with open(log_path, "r", encoding="utf-8") as log_file:
        for line in log_file:
            if "   1 " in line:  # Finds the line starting with mode 1
                try:
                    best_score = float(line.split()[1])  # Grabs the second item (-7.434)
                except (IndexError, ValueError):
                    pass
                break

    # for line in result.stdout.splitlines():
    #     if "   1 " in line:  # Finds the line starting with mode 1
    #         best_score = float(line.split()[1])  # Grabs the second item (-7.442)
    #         break
    # found_table = False
    #
    # with open(log_path, "r") as log_file:
    #     for line in log_file:
    #         # Check for the Vina result header
    #         if "mode |   affinity" in line or "-----+" in line:
    #             found_table = True
    #             print(line.strip())
    #             continue
    #
    #         # Print and capture the score lines
    #         if found_table and line.strip() and line.strip()[0].isdigit():
    #             print(line.strip())
    #             parts = line.split()
    #             # The first parsed score line represents the top binding pose (index 1)
    #             if best_score == 0.0:
    #                 best_score = float(parts[1])
    #         elif found_table and not line.strip():
    #             # Stop parsing once the empty space below the table is hit
    #             found_table = False

    return best_score


# def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center_coords, box_size=22.0):
#     """
#     Instantiates AutoDock Vina, pre-computes interaction grids,
#     and runs a 32-thread Monte Carlo conformational docking search.
#     """
#     print("\nStarting AutoDock Vina Docking Engine")
#
#     # Instantiate the engine using the standard Vina empirical scoring function
#     v = Vina(sf_name='vina')
#
#     # Load input PDBQT files
#     v.set_receptor(receptor_pdbqt)
#     v.set_ligand_from_file(ligand_pdbqt)
#     print("[1/4] Receptor and Ligand files loaded into docking memory.")
#
#     # Define grid box dimensions
#     grid_size = [box_size, box_size, box_size]
#     print(f"[2/4] Bounding box centered at: {center_coords} with size {grid_size}")
#
#     # Calculate grid-based interaction potentials
#     v.compute_vina_maps(center=center_coords, size=grid_size)
#     print("[3/4] Interaction affinity grid maps successfully pre-calculated.")
#
#     # Launch global Monte Carlo search
#     print("[*] Launching Monte Carlo search. Exhaustiveness set to 32...")
#     v.dock(exhaustiveness=32, n_poses=9)
#
#     # Write output poses to disk
#     output_poses_path = 'rutin_docked_poses.pdbqt'
#     v.write_poses(output_poses_path, n_poses=9, overwrite=True)
#
#     # Print quantitative scoring table
#     print("\nDOCKING RESULTS")
#     energies = v.energies()
#     for idx, energy in enumerate(energies):
#         print(f" Pose {idx + 1:2d} | Free Energy of Binding (ΔG): {energy[0]:6.2f} kcal/mol")
#
#     return energies[0][0]

