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
import logging
import argparse
import more_itertools as mit

from joblib.parallel import delayed, Parallel

from rdkit.Chem.EState.Fingerprinter import FingerprintMol as EStateFingerprinter
from rdkit.Chem import MolFromSmiles, SanitizeFlags, SanitizeMol
from rdkit import __version__ as rdkit_version


# ================
# Setup the Logger
LOGGER = logging.getLogger("Compute EState Fingerprints")
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

CH = logging.StreamHandler()
CH.setLevel(logging.INFO)

FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)

LOGGER.addHandler(CH)
# ================


def get_fp_meta_data(fingerprint_definition: str):
    name = fingerprint_definition
    fp_type = "estate"
    fp_mode = "real" if fingerprint_definition == "estate_idc" else "count"
    param = ": ".join(["molecule_presentation", args.molecule_representation])
    library = ": ".join(["RDKit", rdkit_version])
    length = 79
    is_folded = 0
    hash_keys = None

    return name, fp_type, fp_mode, param, library, length, is_folded, hash_keys


def get_fingerprints(batch):
    cids = []
    cnts_bits = []
    cnts_vals = []
    idc_bits = []
    idc_vals = []

    # Computer fingerprints and convert into index-strings
    for cid, smi in batch:
        # Get RDKit mol-object: that might fails
        mol = MolFromSmiles(smi)

        if mol is None:
            LOGGER.warning("[cid = %d] Cannot get mol-object for '%s'. Retry using 'sanitize=False'." % (cid, smi))

            # Approach taking from th RDKit cookbook:
            # http://rdkit.org/docs/Cookbook.html#explicit-valence-error-partial-sanitization
            mol = MolFromSmiles(smi, sanitize=False)
            mol.UpdatePropertyCache(strict=False)
            san_ret = SanitizeMol(
                mol,
                SanitizeFlags.SANITIZE_FINDRADICALS |
                SanitizeFlags.SANITIZE_KEKULIZE |
                SanitizeFlags.SANITIZE_SETAROMATICITY |
                SanitizeFlags.SANITIZE_SETCONJUGATION |
                SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True
            )

            if san_ret != SanitizeFlags.SANITIZE_NONE:
                LOGGER.error("[cid = {}] Cannot get mol-object for '{}'. Manual sanitization failed: {}."
                             .format(cid, smi, san_ret))
                continue

            if mol is None:
                LOGGER.error("[cid = %d] Cannot get mol-object for '%s'." % (cid, smi))
                continue

        try:
            res = EStateFingerprinter(mol)

            cids.append(cid)

            cnts_bits.append(",".join(["%d" % i for i, c in enumerate(res[0]) if c > 0]))
            cnts_vals.append(",".join(["%d" % c for c in res[0] if c > 0]))

            idc_bits.append(",".join(["%d" % i for i, r in enumerate(res[1]) if r != 0]))
            idc_vals.append(",".join(["%f" % r for r in res[1] if r != 0]))
        except RuntimeError as err:
            LOGGER.error("[cid = %d] Cannot compute estate-fingerprint for '%s'." % (cid, smi))
            LOGGER.error(err)

    return cids, cnts_bits, cnts_vals, idc_bits, idc_vals


if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("--n_jobs", type=int, default=4,
                            help="Number of parallel jobs used to compute the fingerprints.")
    arg_parser.add_argument("--batch_size", type=int, default=1024,
                            help="Size of the batches in which the fingerprints are computed and inserted to the DB.")
    arg_parser.add_argument("--molecule_representation", type=str, default="smiles_can",
                            choices=["smiles_can", "smiles_iso"])
    args = arg_parser.parse_args()

    # Open connection to database
    conn = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        # Insert fingerprint meta-data
        with conn:
            for name, fp_type, fp_mode, param, library, length, is_folded, hash_keys in [
                get_fp_meta_data(fp_def) for fp_def in ["estate_idc", "estate_cnt"]
            ]:
                conn.execute("INSERT OR REPLACE INTO fingerprints_meta "
                             "  VALUES (?, ?, ?, ?, DATETIME('now', 'localtime'), ?, ?, ?, ?, ?)",
                             (name, fp_type, fp_mode, param, library, length, is_folded, hash_keys, None))

        # Get all molecular candidate structures from the DB
        rows = conn.execute("SELECT cid, %s FROM molecules" % args.molecule_representation).fetchall()

        # Insert separate table for all the fingerprints
        with conn:
            for name in ["estate_idc", "estate_cnt"]:
                conn.execute("CREATE TABLE IF NOT EXISTS fingerprints_data__%s("
                             "  molecule    INTEGER NOT NULL PRIMARY KEY,"
                             "  bits        VARCHAR NOT NULL,"
                             "  vals        VARCHAR,"
                             "  FOREIGN KEY (molecule) REFERENCES molecules(cid))" % name)

        # Compute fingerprints
        res = Parallel(n_jobs=args.n_jobs, backend="multiprocessing")(
            delayed(get_fingerprints)(batch) for batch in mit.chunked(rows, args.batch_size)
        )

        # Insert fingerprints
        with conn:
            for cids, cnts_bits, cnts_vals, idc_bits, idc_vals in res:
                conn.executemany(
                    "INSERT INTO fingerprints_data__estate_idc VALUES (?, ?,?)",
                    zip(cids, idc_bits, idc_vals)
                )
                conn.executemany(
                    "INSERT INTO fingerprints_data__estate_cnt VALUES (?, ?, ?)",
                    zip(cids, cnts_bits, cnts_vals)
                )

        # Create index on the molecules
        with conn:
            for name in ["estate_idc", "estate_cnt"]:
                conn.execute(
                    "CREATE INDEX IF NOT EXISTS fpd__molecule__%s ON fingerprints_data__%s(molecule)" % (name, name)
                )

    finally:
        conn.close()

