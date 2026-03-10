#!/usr/bin/env python

import os
import glob
import argparse
import subprocess
import string
import json
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBExceptions import PDBIOException
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
import freesasa
import tempfile
from pathlib import Path
from pathlib import Path
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

# path to this file: .../af3-pipeline/src/process_result.py
this_dir = os.path.dirname(os.path.abspath(__file__))

# go up to af3-pipeline/
project_root = os.path.abspath(os.path.join(this_dir, os.pardir))

# vendor/ folder is here
vendor_dir = os.path.join(project_root, "vendor")

default_pdockq_script  = os.path.join(vendor_dir, "pdockq.py")
default_prodigy_script = os.path.join(vendor_dir, "prodigy", "src", "prodigy_prot", "Modified_predict_IC.py")

def collect_confidence(input_folder):
    """
    Scan *folder* (recursively) for files whose names contain
    'summary_confidences' and build a DataFrame with columns:
    [iptm, ptm, ranking_score, min_pae].

    The DataFrame index is the cleaned filename, e.g.
    'fold_5b4x_all_summary_confidences_0.json' → 'fold_5b4x_all_0'.
    """
    rows = []
    pattern = os.path.join(input_folder, "**", "*summary_confidences*.json")

    for path in glob.glob(pattern, recursive=True):
        with open(path) as f:
            data = json.load(f)

        iptm  = float(data.get("iptm","nan"))
        ptm   = float(data.get("ptm","nan"))
        rank  = float(data.get("ranking_score","nan"))

        flat  = [v for row in data.get("chain_pair_pae_min", []) for v in row]
        min_pae = float(min(flat)) if flat else float("nan")

        filename = os.path.splitext(os.path.basename(path))[0]
        key = filename.replace("_summary_confidences", "_model")

        rows.append(
            dict(
                FileName=key,
                Json_path=path,
                iptm=iptm,
                ptm=ptm,
                ranking_score=rank,
                min_pae=min_pae,
            )
        )

    df = pd.DataFrame(
        rows,
        columns=["FileName", "Json_path", "iptm", "ptm", "ranking_score", "min_pae"]
    )
    return df

def convert_cif_to_pdb(input_folder, output_folder):
    """
    Convert all .cif files in the input_folder to .pdb files
    and save them in the output_folder. Also generate an Excel
    summary (summary.xlsx) with file paths.
    """
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Gather all .cif files from the input folder
    cif_files = glob.glob(os.path.join(input_folder, "*.cif"))
    
    # Data list to store file info
    data_list = []
    
    # Initialize the Biopython MMCIF Parser and PDBIO
    parser = MMCIFParser(QUIET=True)
    pdb_writer = PDBIO()

    # Loop through each .cif file
    for cif_file in cif_files:
        # Extract the file name without extension
        base_name = os.path.splitext(os.path.basename(cif_file))[0]
        
        # Parse the MMCIF file into a structure object
        structure = parser.get_structure(base_name, cif_file)
        
        pdb_name = f"{base_name}.pdb"
        # Define the output PDB file path
        pdb_file = os.path.join(output_folder, pdb_name)
        
        # Write the structure to PDB format
        pdb_writer.set_structure(structure)
        pdb_writer.save(pdb_file)
        
        # Store information in list
        data_list.append({
            "FileName": base_name,
            "CIF_Location": cif_file,
            "PDB_Location": pdb_file,
            "PDB_File": pdb_name
        })
    
    # Create a DataFrame from the collected data
    df = pd.DataFrame(data_list, columns=["FileName", "CIF_Location", "PDB_Location", "PDB_File"])
    
    return df

def convert_cif_to_pdb_big(input_folder, output_folder):
    """
    Convert all .cif files in the input_folder to .pdb files
    and save them in the output_folder. If a chain’s ID would
    exceed single-character PDB format (e.g. 'AA'), just omit
    that chain completely.
    """

    os.makedirs(output_folder, exist_ok=True)

    # Gather all .cif files from input_folder
    cif_files = glob.glob(os.path.join(input_folder, "*.cif"))

    data_list = []
    parser = MMCIFParser(QUIET=True)
    pdb_writer = PDBIO()

    # Single-character chain IDs that are valid in standard PDB
    valid_chain_ids = [chr(x) for x in range(ord('A'), ord('Z') + 1)]

    for cif_file in cif_files:
        base_name = os.path.splitext(os.path.basename(cif_file))[0]

        # Parse the structure
        try:
            structure = parser.get_structure(base_name, cif_file)
        except Exception as e:
            print(f"Error parsing {cif_file}: {e}")
            continue

        # Build a new empty Structure to hold "cleaned" models/chains
        cleaned_structure = Structure(structure.id)

        for model in structure:
            # Make a new empty Model for this model ID
            new_model = Model(model.id)

            chain_idx = 0
            for chain in model.get_chains():
                # Only keep a chain if we have a valid one-character ID left
                if chain_idx < len(valid_chain_ids):
                    new_chain_id = valid_chain_ids[chain_idx]
                    # Make a copy of the chain
                    new_chain = Chain(new_chain_id)
                    new_chain.level = "C"  # Just to be explicit

                    # Copy over residues/atoms from old chain
                    for residue in chain:
                        new_chain.add(residue.copy())

                    # Add new_chain to the new_model
                    new_model.add(new_chain)

                    chain_idx += 1
                else:
                    # Skip extra chains if we run out of valid IDs
                    print(f"Skipping chain {chain.id} in {cif_file} because it exceeds PDB chain ID limit.")

            # Add the new_model to the cleaned structure
            cleaned_structure.add(new_model)

        # Write out the cleaned structure
        pdb_name = f"{base_name}.pdb"
        pdb_file = os.path.join(output_folder, pdb_name)
        try:
            pdb_writer.set_structure(cleaned_structure)
            pdb_writer.save(pdb_file)
        except PDBIOException as e:
            print(f"Error writing to {pdb_file}: {e}")
            continue

        # Track info for summary
        data_list.append({
            "FileName": base_name,
            "CIF_Location": cif_file,
            "PDB_Location": pdb_file,
            "PDB_File": pdb_name
        })

    # Create and return a DataFrame summary
    df = pd.DataFrame(data_list,
                      columns=["FileName", "CIF_Location", "PDB_Location", "PDB_File"])
    return df

def run_prodigy(output_folder, temperature=None, script_path=None):
    """
    For each PDB file in output_folder, run 'predict_IC.py' with the specified temperature.
    """

    if script_path is None:
        raise ValueError("Path to PRODIGY script must be provided via --prodigy_script")

    pdb_files = glob.glob(os.path.join(output_folder, "*.pdb"))
    
    results = []
    column_names = [
        "No. of intermolecular contacts",
        "No. of charged-charged contacts",
        "No. of charged-polar contacts",
        "No. of charged-apolar contacts",
        "No. of polar-polar contacts",
        "No. of apolar-polar contacts",
        "No. of apolar-apolar contacts",
        "Percentage of apolar NIS residues",
        "Percentage of charged NIS residues",
        "Predicted binding affinity (kcal.mol-1)",
        "Temperature",
        "Predicted dissociation constant (M)"
    ]
    
    for pdb_file in pdb_files:
        # IMPORTANT: Pass --temperature as a named argument
        command = [
            "python",
            script_path,
            "--temperature",
            str(temperature),
            pdb_file
        ]
        
        try:
            proc_output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode("utf-8")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {pdb_file}:\n{e.output.decode('utf-8')}")
            continue
        
        for line in proc_output.splitlines():
            line = line.strip()
            if line.startswith("output:"):
                data_str = line.split("output:")[1].strip()
                values = data_str.split(",")
                # e.g. "78,18,2,21,0,13,24,36.60,26.35,-12.2,25.0,1.1e-09"
                
                # We expect 12 values
                if len(values) == 12:
                    try:
                        values = [float(x) for x in values]
                    except ValueError:
                        pass
                    results.append([os.path.basename(pdb_file)] + values)
                else:
                    print(f"Warning: 'output:' line for {pdb_file} does not have 12 elements.")
    
    df_columns = ["PDB_File"] + column_names
    df = pd.DataFrame(results, columns=df_columns)
    return df

def run_pdockq(output_folder, script_path=None):
    """
    For each PDB file in output_folder, run 'pdockq.py' and capture
    lines of the form 'pDockQ= 0.608 ,PPV= 0.9400192'.
    
    Returns a DataFrame with columns:
    ['PDB_File', 'pDockQ', 'PPV'].
    """

    if script_path is None:
        raise ValueError("Path to pDockQ script must be provided via --pdockq_script")
    
    # Find all .pdb files in the output folder
    pdb_files = glob.glob(os.path.join(output_folder, "*.pdb"))
    
    # Prepare a list to store results
    results = []
    
    # Loop over each pdb file
    for pdb_file in pdb_files:
        # Command to run pdockq with the --pdbfile argument
        command = [
            "python",
            script_path,
            "--pdbfile",
            pdb_file
        ]
        
        try:
            # Run pdockq and capture the output
            proc_output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode("utf-8")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {pdb_file}:\n{e.output.decode('utf-8')}")
            continue
        
        pDockQ_val = None
        PPV_val = None
        
        for line in proc_output.splitlines():
            line = line.strip()
            if line.startswith("pDockQ="):
            
                parts = line.replace(" ", "").split(",")  # => ["pDockQ=0.608", "PPV=0.9400192"]
                if len(parts) == 2:
                    # Parse pDockQ
                    if parts[0].startswith("pDockQ="):
                        pDockQ_str = parts[0].split("pDockQ=")[1]
                        pDockQ_val = float(pDockQ_str)
                    # Parse PPV
                    if parts[1].startswith("PPV="):
                        PPV_str = parts[1].split("PPV=")[1]
                        PPV_val = float(PPV_str)
                
                # Once we parse them, we can break from the loop
                break
 
        # If we got values, append them to results
        if pDockQ_val is not None and PPV_val is not None:
            results.append([os.path.basename(pdb_file), pDockQ_val, PPV_val])
        else:
            print(f"Warning: Could not parse pDockQ/PPV in output for {pdb_file}.")

    # Create the DataFrame
    df = pd.DataFrame(results, columns=["PDB_File", "pDockQ", "PPV"])
    return df

def compute_total_sasa(pdb_path):
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    return result.totalArea()

def write_chain_subset_pdb(input_pdb, output_pdb, chains):
    chains = set(chains)
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21].strip()
                if chain_id in chains:
                    fout.write(line)
        fout.write("END\n")

def compute_sasa_metrics(pdb_path, chain_a="A", chain_b="B"):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        pdb_a = tmpdir / f"{chain_a}.pdb"
        pdb_b = tmpdir / f"{chain_b}.pdb"
        pdb_ab = tmpdir / f"{chain_a}_{chain_b}.pdb"

        write_chain_subset_pdb(pdb_path, pdb_a, [chain_a])
        write_chain_subset_pdb(pdb_path, pdb_b, [chain_b])
        write_chain_subset_pdb(pdb_path, pdb_ab, [chain_a, chain_b])

        sasa_a = compute_total_sasa(pdb_a)      # A alone
        sasa_b = compute_total_sasa(pdb_b)      # B alone
        sasa_ab = compute_total_sasa(pdb_ab)    # A+B complex

        expose_a, expose_b = compute_chain_sasa_in_complex(
            pdb_ab, chain_a=chain_a, chain_b=chain_b
        )

    delta_a = sasa_a - expose_a
    delta_b = sasa_b - expose_b

    bsa_total = delta_a + delta_b
    bsa = bsa_total / 2.0

    return {
        "SASA_A": sasa_a,
        "SASA_B": sasa_b,
        "SASA_AB": sasa_ab,
        "Expose_A": expose_a,
        "Expose_B": expose_b,
        "Delta_A": delta_a,
        "Delta_B": delta_b,
        "BSA_total": bsa_total,
        "BSA": bsa,
    }

def compute_chain_sasa_in_complex(pdb_path, chain_a="A", chain_b="B"):
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)

    expose_a = 0.0
    expose_b = 0.0

    for atom_index in range(structure.nAtoms()):
        chain = structure.chainLabel(atom_index).strip()
        area = result.atomArea(atom_index)

        if chain == chain_a:
            expose_a += area
        elif chain == chain_b:
            expose_b += area

    return expose_a, expose_b

def run_freesasa(output_folder, chain_a="A", chain_b="B"):
    pdb_files = glob.glob(os.path.join(output_folder, "*.pdb"))
    results = []

    for pdb_file in pdb_files:
        try:
            m = compute_sasa_metrics(pdb_file, chain_a=chain_a, chain_b=chain_b)
            results.append([
                os.path.basename(pdb_file),
                m["SASA_A"],
                m["SASA_B"],
                m["SASA_AB"],
                m["Expose_A"],
                m["Expose_B"],
                m["Delta_A"],
                m["Delta_B"],
                m["BSA_total"],
                m["BSA"],
            ])
        except Exception as e:
            print(f"Error running FreeSASA on {pdb_file}: {e}")

    return pd.DataFrame(
        results,
        columns=[
            "PDB_File",
            "SASA_A", "SASA_B", "SASA_AB",
            "Expose_A", "Expose_B",
            "Delta_A", "Delta_B",
            "BSA_total", "BSA"
        ]
    )

def residue_to_one_letter(residue):
    resname = residue.get_resname().strip().upper()
    return protein_letters_3to1.get(resname, "X")

def residue_sort_key(residue):
    _, resseq, icode = residue.id
    return (resseq, icode.strip())

def residue_label(residue):
    _, resseq, icode = residue.id
    aa = residue_to_one_letter(residue)
    icode = icode.strip()
    return f"{resseq}{icode}{aa}" if icode else f"{resseq}{aa}"

def get_chain_residues(chain):
    return [res for res in chain if is_aa(res, standard=False)]

def get_interface_labels(pdb_path, chain_a="A", chain_b="B", radius=5.0):
    """
    Return two sorted lists of interface residue labels:
    labels_A, labels_B
    Example label: 123K or 101AY (if insertion code exists)
    """
    pdb_path = Path(pdb_path)
    structure = PDBParser(QUIET=True).get_structure("complex", str(pdb_path))

    model = next(structure.get_models())
    chainA = model[chain_a]
    chainB = model[chain_b]

    residues_A = get_chain_residues(chainA)
    residues_B = get_chain_residues(chainB)

    atoms_A = [atom for res in residues_A for atom in res.get_atoms()]
    atoms_B = [atom for res in residues_B for atom in res.get_atoms()]

    ns_A = NeighborSearch(atoms_A)
    ns_B = NeighborSearch(atoms_B)

    interface_A = set()
    interface_B = set()

    for res in residues_A:
        for atom in res.get_atoms():
            hits = ns_B.search(atom.coord, radius, level="R")
            hits = [r for r in hits if r.get_parent().id == chain_b and is_aa(r, standard=False)]
            if hits:
                interface_A.add(res)
                break

    for res in residues_B:
        for atom in res.get_atoms():
            hits = ns_A.search(atom.coord, radius, level="R")
            hits = [r for r in hits if r.get_parent().id == chain_a and is_aa(r, standard=False)]
            if hits:
                interface_B.add(res)
                break

    labels_A = [residue_label(res) for res in sorted(interface_A, key=residue_sort_key)]
    labels_B = [residue_label(res) for res in sorted(interface_B, key=residue_sort_key)]

    return labels_A, labels_B

def run_interface_residue(output_folder, chain_a="A", chain_b="B", radius=5.0):
    pdb_files = glob.glob(os.path.join(output_folder, "*.pdb"))
    results = []

    for pdb_file in pdb_files:
        try:
            labels_A, labels_B = get_interface_labels(
                pdb_file,
                chain_a=chain_a,
                chain_b=chain_b,
                radius=radius
            )

            results.append({
                "PDB_File": os.path.basename(pdb_file),
                "Interface_A": ",".join(labels_A),
                "Interface_B": ",".join(labels_B),
                "N_Interface_A": len(labels_A),
                "N_Interface_B": len(labels_B),
            })

        except Exception as e:
            print(f"Error finding interface residues for {pdb_file}: {e}")
            results.append({
                "PDB_File": os.path.basename(pdb_file),
                "Interface_A": None,
                "Interface_B": None,
                "N_Interface_A": None,
                "N_Interface_B": None,
            })

    return pd.DataFrame(
        results,
        columns=["PDB_File", "Interface_A", "Interface_B", "N_Interface_A", "N_Interface_B"]
    )

def main():
    parser = argparse.ArgumentParser(description="Convert .cif files to .pdb format using Biopython, then run Prodigy.")
    parser.add_argument("--input_folder", "-i", required=True, help="Path to the input folder containing .cif files.")
    parser.add_argument("--output_folder", "-o", required=True, help="Path to the output folder for .pdb files and summary.")
    parser.add_argument("--temperature", "-t", required=False, default= 25.0, help="Temperature (in °C) to use in PRODIGY calculation.")
    parser.add_argument("--big", action="store_true", default=False,  help="add this when .cif chain ID going to exceed 'Z'")
    parser.add_argument("--prodigy_script", required=False, default=default_prodigy_script,help="Path to Prodigy Modified_predict_IC.py script.")
    parser.add_argument("--pdockq_script", required=False, default=default_pdockq_script,help="Path to pDockQ script.")
    parser.add_argument("--chain_a", default="A", help="First chain for interface residue detection")
    parser.add_argument("--chain_b", default="B", help="Second chain for interface residue detection")
    parser.add_argument("--interface_radius", type=float, default=5.0, help="Distance cutoff in Å for interface residues")
    args = parser.parse_args()
    
    # Step 1: Collect AF3 Confidence
    confidence_df = collect_confidence(args.input_folder)

    # Step 2: Convert CIF -> PDB
    if not args.big:
        pdb_sum_df = convert_cif_to_pdb(args.input_folder, args.output_folder)
    else:
        pdb_sum_df = convert_cif_to_pdb_big(args.input_folder, args.output_folder)
    
        
    # Step 3: Run pDockQ on each PDB file in output_folder
    pdockq_results_df = run_pdockq(args.output_folder, args.pdockq_script)

    # Step 4: Run freeSasa on each PDB file in output_folder using Chain A, Chain B, Complex
    freesasa_results_df = run_freesasa(args.output_folder, args.chain_a, args.chain_b)

    # Step 5: Run Prodigy on each PDB file in output_folder, using the specified or default temperature
    prodigy_results_df = run_prodigy(args.output_folder, args.temperature, args.prodigy_script)

    # Step 6: Calculate interface residues
    interface_df = run_interface_residue(
        args.output_folder,
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        radius=args.interface_radius
    )
    
    # Step 7: Combine three df and output as one
    summary_df = (
        pdb_sum_df
        .merge(confidence_df, on="FileName", how="outer")
        .merge(pdockq_results_df, on="PDB_File", how="outer")
        .merge(prodigy_results_df, on="PDB_File", how="outer")
        .merge(freesasa_results_df, on="PDB_File", how="outer")
        .merge(interface_df, on="PDB_File", how="outer")
    )

    summary_file = os.path.join(args.output_folder, "summary.xlsx")
    summary_df.to_excel(summary_file, index=False)
    print(f"Conversion complete! Summary saved at: {summary_file}")

if __name__ == "__main__":
    main()

