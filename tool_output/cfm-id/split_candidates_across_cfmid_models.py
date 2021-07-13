import os
import argparse
import pandas as pd
import numpy as np
import sqlite3

from tqdm import tqdm
from rdkit.Chem.AllChem import MolFromSmiles, MolToInchiKey


DEFAULT_PRETRAINED_MODEL_IDX = 0


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("candidate_database")
    arg_parser.add_argument("training_molecules_fn__positive")
    arg_parser.add_argument("training_molecules_fn__negative")
    arg_parser.add_argument("output_directory")
    args = arg_parser.parse_args()

    # Read in training molecules and compute InChIKey

    df_train = {
        "positive": pd.read_csv(
            args.training_molecules_fn__positive, sep=" ", usecols=[1, 2], names=["smiles", "cv_fold"], header=None
        ),
        "negative": pd.read_csv(
            args.training_molecules_fn__negative, sep=" ", usecols=[1, 2], names=["smiles", "cv_fold"], header=None
        )
    }
    for ion in df_train.keys():
        df_train[ion]["inchikey1"] = [MolToInchiKey(MolFromSmiles(smi)).split("-")[0] for smi in df_train[ion]["smiles"]]

    # There is a candidate set for each CFM-ID model
    candidates = {
        "positive": [set() for _ in range(10)],
        "negative": [set() for _ in range(10)]
    }

    # Track which model was used for which spectrum
    df_spec2model = []

    mb_db = sqlite3.connect(args.candidate_database)
    try:
        # Get all spectrum ids and the corresponding InChIKey(1)s
        rows = mb_db.execute(
            "SELECT accession, inchikey1, precursor_type FROM scored_spectra_meta"
            "   INNER JOIN molecules m on m.cid = scored_spectra_meta.molecule"
        ).fetchall()

        for acc, ikey1, ptype in tqdm(rows, desc="Process spectra"):
            # Determine ionization time
            if ptype.endswith("+"):
                ion = "positive"
            elif ptype.endswith("-"):
                ion = "negative"
            else:
                raise ValueError("Invalid precursor type: '%s'" % ptype)

            # Check for the spectrum, whether it is used for the CFM-ID training and if yes in which fold
            try:
                idx = df_train[ion]["inchikey1"].tolist().index(ikey1)
                cv_fold = df_train[ion].iloc[idx]["cv_fold"]
            except ValueError:
                cv_fold = DEFAULT_PRETRAINED_MODEL_IDX  # we use model for fold 0 as default

            # Get the candidates for the current spectrum
            for ikey1_cnd, smi_cnd in mb_db.execute(
                "SELECT inchikey1, smiles_can FROM candidates_spectra "
                "   INNER JOIN molecules m ON m.cid = candidates_spectra.candidate"
                "   WHERE spectrum IS ?", (acc, )
            ):
                # Add the molecule and its smiles representation to prediction list for the current model
                candidates[ion][cv_fold] |= {(ikey1_cnd, smi_cnd)}

            df_spec2model.append((acc, ikey1, cv_fold, ion))

    finally:
        mb_db.close()

    # Write out which model is used for which spectrum
    pd.DataFrame(df_spec2model, columns=["accession", "inchikey1", "model", "ionization"]) \
        .to_csv(os.path.join(args.output_directory, "spec2model.tsv"), sep="\t", index=False)

    # Write out the model specific candidate sets
    for ion in ["positive", "negative"]:
        for cv_fold in tqdm(range(10), desc="Write out candidate files (%s)" % ion):
            if len(candidates[ion][cv_fold]) > 0:
                with open(
                        os.path.join(args.output_directory, "candidates__cv=%d__ion=%s.txt" % (cv_fold, ion)), "w"
                ) as ofile:
                    for ikey1, smi in candidates[ion][cv_fold]:
                        ofile.write("%s %s\n" % (ikey1, smi))
