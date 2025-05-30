#!/usr/bin/env python

import os
import glob
import argparse
import subprocess
import string
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBExceptions import PDBIOException
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

# path to this file: .../af3-pipeline/src/process_result.py
this_dir = os.path.dirname(os.path.abspath(__file__))

# go up to af3-pipeline/
project_root = os.path.abspath(os.path.join(this_dir, os.pardir))

# vendor/ folder is here
vendor_dir = os.path.join(project_root, "vendor")

default_pdockq_script  = os.path.join(vendor_dir, "pdockq.py")
default_prodigy_script = os.path.join(vendor_dir, "prodigy_prot", "src", "prodigy_prot", "Modified_predict_IC.py")

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
    # If you also want digits, you could do something like:
    # valid_chain_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")

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
                        # Biopython 1.80+ has .copy(), or you can replicate manually
                        new_chain.add(residue.copy())
                    # Alternatively:
                    # for residue in chain:
                    #     new_chain.add(residue.copy())

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

    # Path to your 'predict_IC.py' script
    #prodigy_script = r"C:\Users\mingz\Desktop\Alphafold\AF3\prodigy\src\prodigy_prot\predict_IC.py"
    
    # If no temperature was provided, default to 25°C
    if temperature is None:
        temperature = 25
    else:
        # Convert the incoming temperature (string) to a float if possible
        temperature = float(temperature)

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

    # Path to your 'pdockq.py' script
    #pdockq_script = r"C:\Users\mingz\Desktop\Alphafold\AF3\pdockq.py"
    
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
        
        # Look for the line containing 'pDockQ='
        # Example line: "pDockQ= 0.608 ,PPV= 0.9400192"
        pDockQ_val = None
        PPV_val = None
        
        for line in proc_output.splitlines():
            line = line.strip()
            if line.startswith("pDockQ="):
                # Attempt to parse "pDockQ= 0.608 ,PPV= 0.9400192"
                # Remove any commas or extra spaces as needed
                # e.g. split by spaces or by commas
                # Alternatively, you can do a more robust parse:
                
                # If the line always looks like "pDockQ= X ,PPV= Y"
                # then we can do a simple approach:
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

def main():
    parser = argparse.ArgumentParser(description="Convert .cif files to .pdb format using Biopython, then run Prodigy.")
    parser.add_argument("--input_folder", required=True, help="Path to the input folder containing .cif files.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder for .pdb files and summary.")
    parser.add_argument("--temperature", required=False, help="Temperature (in °C) to use in PRODIGY calculation.")
    parser.add_argument("--big", action="store_true", default=False,  help="add this when .cif chain ID going to exceed 'Z'")
    parser.add_argument("--prodigy_script", required=False, default=default_prodigy_script,help="Path to Prodigy Modified_predict_IC.py script.")
    parser.add_argument("--pdockq_script", required=False, default=default_pdockq_script,help="Path to pDockQ script.")
    args = parser.parse_args()
    
    # Step 1: Convert CIF -> PDB
    if not args.big:
        pdb_sum_df = convert_cif_to_pdb(args.input_folder, args.output_folder)
    else:
        pdb_sum_df = convert_cif_to_pdb_big(args.input_folder, args.output_folder)
    
    # Step 2: Run pDockQ on each PDB file in output_folder
    pdockq_results_df = run_pdockq(args.output_folder, args.pdockq_script)

    # Step 3: Run Prodigy on each PDB file in output_folder, using the specified or default temperature
    prodigy_results_df = run_prodigy(args.output_folder, args.temperature, args.prodigy_script)

    # Step 4: Combine three df and output as one
    summary_df = (pdb_sum_df  .merge(pdockq_results_df, on="PDB_File", how="outer")
                                    .merge(prodigy_results_df, on="PDB_File", how="outer"))
    summary_file = os.path.join(args.output_folder, "summary.xlsx")
    summary_df.to_excel(summary_file, index=False)
    print(f"Conversion complete! Summary saved at: {summary_file}")

if __name__ == "__main__":
    main()


#python process_result.py --input_folder /Users/mac/Desktop/AF3/folds_2025_05_22_07_31/5b4x_all --output_folder /Users/mac/Desktop/AF3/folds_2025_05_22_07_31/5b4x_all/pdb --temperature 37
