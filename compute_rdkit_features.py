####
#
# The MIT License (MIT)
#
# Copyright 2021 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####
import sqlite3
import argparse
import numpy as np
import more_itertools as mit

from joblib.parallel import Parallel, delayed
from sklearn.impute import SimpleImputer

from rdkit.Chem import Descriptors
from rdkit import __version__ as rdkit_version

from rosvm.feature_extraction.featurizer_cls import FeaturizerMixin
from rosvm import __version__ as rosvm_version

from utils import get_backup_db


def get_descriptors(data_batch):
    # Get the list of descriptors available in RDKit
    l_rdkit_desc = sorted(Descriptors.descList)

    # Matrix storing all descriptor values
    X = np.zeros((len(data_batch), len(l_rdkit_desc)))

    cids = []
    for idx, (cid, smi) in enumerate(data_batch):
        # Get RDKit mol-objects
        mol = FeaturizerMixin.sanitize_mol(smi)

        # Compute descriptors
        for jdx, (_, dfun) in enumerate(l_rdkit_desc):
            X[idx, jdx] = dfun(mol)

        cids.append(cid)

    # Find and replace np.inf - values
    X[np.isinf(X)] = np.nan

    # How many molecules did have problems with the descriptor computation?
    print("%d / %d molecules with nan-descriptors." % (np.sum(np.any(np.isnan(X), axis=1)).item(), len(data_batch)))

    # Impute missing values
    X = SimpleImputer(copy=False).fit_transform(X).astype(str)

    # Convert descriptor matrix to list of strings
    dvals = [",".join(x_i) for x_i in X]

    return cids, dvals


if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument(
        "--batch_size",
        type=int,
        default=2**14,
        help="Size of the batches in which the descriptors are computed and inserted to the DB."
    )
    arg_parser.add_argument(
        "--n_jobs",
        type=int,
        default=4
    )
    args = arg_parser.parse_args()

    # Open connection to database
    conn = get_backup_db(args.massbank_db_fn, overwrite=True, postfix="with_descriptors")
    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        # Create meta and data tables
        with conn:
            desc_name = "bouwmeester"

            conn.execute(
                "CREATE TABLE IF NOT EXISTS descriptors_meta("
                "   name        VARCHAR NOT NULL PRIMARY KEY,"
                "   length      INTEGER NOT NULL,"
                "   mode        VARCHAR NOT NULL,"
                "   timestamp   VARCHAR NOT NULL,"
                "   library     VARCHAR NOT NULL" 
                ")"
            )

            conn.execute(
                "INSERT OR REPLACE INTO descriptors_meta "
                "   VALUES ('%s', %d, 'real', DATETIME('now', 'localtime'), 'rosvm: %s, RDKit: %s')"
                % (desc_name, len(Descriptors.descList), rosvm_version, rdkit_version)
            )

            conn.execute(
                "CREATE TABLE IF NOT EXISTS descriptors_data__%s("
                "   molecule    INTEGER NOT NULL PRIMARY KEY,"
                "   desc_vals   VARCHAR NOT NULL,"
                "   FOREIGN KEY (molecule) REFERENCES molecules(cid)"
                ")" % desc_name
            )

        # Load all molecular descriptors
        print("Load molecules ...")
        data = conn.execute("SELECT cid, smiles_iso FROM molecules ORDER BY random()").fetchall()

        # Compute descriptors in batches
        print("Compute descriptors ...")
        batches = list(mit.chunked(range(len(data)), args.batch_size))
        res = Parallel(n_jobs=args.n_jobs, verbose=11)(
            delayed(get_descriptors)([data[idx] for idx in batch])
            for batch in batches[:4]
        )

        # Insert descriptors
        print("Insert descriptors ...")
        with conn:
            for cids, dvals in res:
                conn.executemany(
                    "INSERT OR REPLACE INTO descriptors_data__%s VALUES (?, ?)" % desc_name,
                    zip(cids, dvals)
                )

    finally:
        conn.close()
