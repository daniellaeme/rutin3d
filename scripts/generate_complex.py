import os

def extract_best_pose_and_merge(receptor_path, docked_pdbqt_path, output_complex_pdb):
    """
    Extracts Model 1 from a Vina PDBQT file, cleans its structural syntax,
    and merges it with the receptor PDB to build a single operational complex.
    Writes output complex to the target directory.
    """
    best_pose_lines = []
    with open(docked_pdbqt_path, 'r') as f:
        in_model_1 = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_model_1 = True
                continue
            if line.startswith("ENDMDL") and in_model_1:
                break
            if in_model_1:
                # Keep only atomic coordinates (ATOM/HETATM)
                if line.startswith(("ATOM", "HETATM")):
                    # Truncate the Vina specific columns (charges/atom types)
                    # to keep standard PDB formatting compliance
                    best_pose_lines.append(line[:66] + "\n")

    if not best_pose_lines:
        raise RuntimeError("Could not isolate Model 1 from the docked PDBQT file.")

    print(f"[1/2] Extracted {len(best_pose_lines)} ligand atom lines from Pose 1.")

    # Read the cleaned receptor structure
    with open(receptor_path, 'r') as f:
        receptor_lines = f.readlines()

    # Combine them into a singular structural complex file
    with open(output_complex_pdb, 'w') as f:
        # Write protein atoms first
        for line in receptor_lines:
            if line.startswith(("ATOM", "HETATM", "TER")):
                f.write(line)

        f.write("TER\n")  # Tell structural parsers the protein chain ended

        # Write the ligand atoms as a distinct hetero-group
        for line in best_pose_lines:
            f.write(line)

        f.write("END\n")

    print(f"[2/2] Success! Publication-ready complex saved to: {output_complex_pdb}")

    return






