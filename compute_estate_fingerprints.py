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
import more_itertools as mit

from joblib.parallel import delayed, Parallel

from rdkit.Chem.EState.Fingerprinter import FingerprintMol as EStateFingerprinter
from rdkit.Chem import MolFromSmiles
from rdkit import __version__ as rdkit_version


def get_fp_meta_data(fingerprint_definition: str):
    if fingerprint_definition == "estate_idc":
        name = fingerprint_definition
        fp_type = "estate"
        fp_mode = "real"
        param = "molecule_presentation=%s" % args.molecule_representation
        library = ":".join(["RDKit", rdkit_version])
        length = 79
        is_folded = 0
        hash_keys = None
    elif fingerprint_definition == "estate_cnt":
        name = fingerprint_definition
        fp_type = "estate"
        fp_mode = "count"
        param = "molecule_presentation=%s" % args.molecule_representation
        library = ":".join(["RDKit", rdkit_version])
        length = 79
        is_folded = 0
        hash_keys = None
    else:
        raise ValueError("Invalid fingerprint definition: '%s'" % fingerprint_definition)

    return name, fp_type, fp_mode, param, library, length, is_folded, hash_keys


def get_fingerprints(batch):
    cids = []
    fp_cnt = []
    fp_idc = []

    # Computer fingerprints and convert into index-strings
    for cid, smi in batch:
        cids.append(cid)

        res = EStateFingerprinter(MolFromSmiles(smi))

        fp_cnt.append(",".join(["%d:%d" % (i, c) for i, c in enumerate(res[0]) if c > 0]))
        fp_idc.append(",".join(["%d:%f" % (i, r) for i, r in enumerate(res[1]) if r != 0]))

    return cids, fp_cnt, fp_idc


if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("--n_jobs", type=int, default=4,
                            help="Number of parallel jobs used to compute the fingerprints.")
    arg_parser.add_argument("--batch_size", type=int, default=1024,
                            help="Size of the batches in which the fingerprints are computed and inserted to the DB.")
    arg_parser.add_argument("--molecule-representation", type=str, default="smiles_can",
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
                             "  VALUES (?, ?, ?, ?, DATETIME('now', 'localtime'), ?, ?, ?, ?)",
                             (name, fp_type, fp_mode, param, library, length, is_folded, hash_keys))

        # Get all molecular candidate structures from the DB
        rows = conn.execute("SELECT cid, smiles_can FROM molecules").fetchall()

        # Compute fingerprints
        res = Parallel(n_jobs=args.n_jobs, backend="multiprocessing")(
            delayed(get_fingerprints)(batch) for batch in mit.chunked(rows, args.batch_size)
        )

        # Insert fingerprints
        with conn:
            for cids, fps_cnt, fps_idc in res:
                conn.executemany(
                    "INSERT INTO fingerprints_data VALUES (?, ?, ?)", zip(cids, len(cids) * ["estate_cnt"], fps_cnt)
                )
                conn.executemany(
                    "INSERT INTO fingerprints_data VALUES (?, ?, ?)", zip(cids, len(cids) * ["estate_idc"], fps_idc)
                )

    finally:
        conn.close()

