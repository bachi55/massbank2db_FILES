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

from rdkit import __version__ as rdkit_version

from rosvm.feature_extraction.featurizer_cls import CircularFPFeaturizer
from rosvm import __version__ as rosvm_version


# ================
# Setup the Logger
LOGGER = logging.getLogger("Compute Circular Fingerprints")
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

CH = logging.StreamHandler()
CH.setLevel(logging.INFO)

FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)

LOGGER.addHandler(CH)
# ================


def get_fp_meta_data(fprinter: CircularFPFeaturizer):
    fp_type = fprinter.fp_type
    fp_mode = fprinter.fp_mode
    name = "__".join([fp_type, fp_mode, args.key_training_set])  # e.g. ECFP__count__gt_structures
    param = ", ".join(["{}: {}".format(k, v) for k, v in fprinter.get_params(deep=False).items()])
    param += ", molecule_representation: %s" % args.molecule_representation
    library = ", ".join(["{}: {}".format(p, v) for p, v in [("rosvm", rosvm_version), ("RDKit", rdkit_version)]])
    length = fprinter.get_length()
    is_folded = int(fprinter.fp_mode == "binary_folded")

    if hasattr(fprinter, "freq_hash_set_"):
        hash_keys = ",".join(["%d" % k for k in fprinter.freq_hash_set_.keys()])
    else:
        hash_keys = None

    return name, fp_type, fp_mode, param, library, length, is_folded, hash_keys


def get_fingerprints(batch, fprinter):
    cids, smis = map(list, zip(*batch))
    return cids, fprinter.transform(smis)


if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("key_training_set", type=str, choices=["gt_structures", "random_candidates"])
    arg_parser.add_argument("fp_type", type=str, choices=["ECFP", "FCFP"])
    arg_parser.add_argument("--n_jobs", type=int, default=4,
                            help="Number of parallel jobs used to compute the fingerprints.")
    arg_parser.add_argument("--batch_size", type=int, default=2**14,
                            help="Size of the batches in which the fingerprints are computed and inserted to the DB.")
    arg_parser.add_argument("--molecule_representation", type=str, default="smiles_iso",
                            choices=["smiles_can", "smiles_iso"])
    args = arg_parser.parse_args()

    # Open connection to database
    conn = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        # Load the SMILES used to train the frequent circular fingerprint hashes
        if args.key_training_set == "random_candidates":
            rows = conn.execute(
                "SELECT %s FROM molecules ORDER BY random() LIMIT 250000" % args.molecule_representation
            ).fetchall()
            min_subs_freq = 250
        else:  # gt_structures
            rows = conn.execute(
                "SELECT m.%s FROM scored_spectra_meta "
                "   INNER JOIN molecules m on scored_spectra_meta.cid = m.cid" % args.molecule_representation
            ).fetchall()
            min_subs_freq = 25

        LOGGER.info("Number of SMILES used for training (%s): %d." % (args.key_training_set, len(rows)))

        # Train the circular fingerprinter
        fprinter = CircularFPFeaturizer(
            fp_type=args.fp_type, only_freq_subs=True, min_subs_freq=min_subs_freq, n_jobs=args.n_jobs, radius=3,
            use_chirality=True, output_format="sparse_string"
        ).fit([row[0] for row in rows])
        LOGGER.info("Size of frequent hash set: %d" % len(fprinter))

        # Insert fingerprint meta-data
        with conn:
            name, fp_type, fp_mode, param, library, length, is_folded, hash_keys = get_fp_meta_data(fprinter)
            conn.execute(
                "INSERT OR REPLACE INTO fingerprints_meta "
                "  VALUES (?, ?, ?, ?, DATETIME('now', 'localtime'), ?, ?, ?, ?)",
                (name, fp_type, fp_mode, param, library, length, is_folded, hash_keys)
            )

        # Get all molecular candidate structures from the DB
        rows = conn.execute("SELECT cid, %s FROM molecules" % args.molecule_representation).fetchall()

        # Compute and insert fingerprints
        for batch in mit.chunked(rows, args.batch_size):
            cids, fps = get_fingerprints(batch, fprinter)

            with conn:
                conn.executemany(
                    "INSERT OR REPLACE INTO fingerprints_data VALUES (?, ?, ?)", zip(cids, len(cids) * [name], fps)
                )
                conn.executemany(
                    "INSERT OR REPLACE INTO fingerprints_data VALUES (?, ?, ?)", zip(cids, len(cids) * [name], fps)
                )

    finally:
        conn.close()

