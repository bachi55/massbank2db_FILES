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
import tarfile
import os
import pandas as pd
import re
import glob
import logging
import traceback

from massbank2db.db import MassbankDB

# ================
# Setup the Logger
LOGGER = logging.getLogger("Import CSI:FingerID Scores")
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

CH = logging.StreamHandler()
CH.setLevel(logging.INFO)

FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)

LOGGER.addHandler(CH)
# ================


def create_tables(conn: sqlite3.Connection):
    conn.execute("CREATE TABLE IF NOT EXISTS scored_spectra_meta ("
                 "  accession           VARCHAR PRIMARY KEY,"
                 "  dataset             VARCHAR NOT NULL," 
                 "  original_accessions VARCHAR NOT NULL,"
                 "  precursor_mz        REAL NOT NULL,"
                 "  precursor_type      VARCHAR NOT NULL,"
                 "  cid                 INTEGER NOT NULL,"
                 "  retention_time      REAL NOT NULL,"
                 "  FOREIGN KEY (cid) REFERENCES molecules(cid),"
                 "  FOREIGN KEY (dataset) REFERENCES datasets(name))")

    conn.execute("CREATE TABLE IF NOT EXISTS merged_accessions ("
                 "  massbank_accession  VARCHAR PRIMARY KEY,"
                 "  merged_accession    VARCHAR NOT NULL,"
                 "  FOREIGN KEY (massbank_accession) REFERENCES spectra_meta (accession),"
                 "  FOREIGN KEY (merged_accession) REFERENCES scored_spectra_meta(accession))")

    conn.execute("CREATE TABLE IF NOT EXISTS scoring_methods ("
                 "  name        VARCHAR PRIMARY KEY,"
                 "  method      VARCHAR NOT NULL,"
                 "  description VARCHAR)")

    conn.execute("CREATE TABLE IF NOT EXISTS spectra_candidate_scores ("
                 "  spectrum        VARCHAR NOT NULL,"
                 "  candidate       INTEGER NOT NULL,"
                 "  scoring_method  VARCHAR NOT NULL,"
                 "  dataset         VARCHAR NOT NULL,"
                 "  score           REAL NOT NULL,"
                 "  FOREIGN KEY (spectrum) REFERENCES scored_spectra_meta(accession),"
                 "  FOREIGN KEY (candidate) REFERENCES molecules(cid),"
                 "  FOREIGN KEY (dataset) REFERENCES datasets(name),"
                 "  FOREIGN KEY (scoring_method) REFERENCES scoring_methods(name),"
                 "  PRIMARY KEY (spectrum, candidate, scoring_method))")

    conn.execute("CREATE TABLE IF NOT EXISTS fingerprints_meta ("
                 "  name        VARCHAR PRIMARY KEY,"
                 "  type        VARCHAR NOT NULL,"
                 "  mode        VARCHAR NOT NULL,"
                 "  param       VARCHAR,"
                 "  timestamp   VARCHAR NOT NULL,"
                 "  library     VARCHAR NOT NULL,"
                 "  length      INTEGER NOT NULL,"
                 "  is_folded   BOOLEAN NOT NULL,"
                 "  hash_keys   VARCHAR,"
                 "  CONSTRAINT type_mode_combination UNIQUE (type, mode, param))")

    conn.execute("CREATE TABLE IF NOT EXISTS fingerprints_data ("
                 "  molecule    INTEGER NOT NULL,"
                 "  name        VARCHAR NOT NULL,"
                 "  fingerprint VARCHAR NOT NULL,"
                 "  FOREIGN KEY (molecule) REFERENCES molecules(cid),"
                 "  FOREIGN KEY (name) REFERENCES fingerprints_met(name),"
                 "  PRIMARY KEY (molecule, name))")


def create_indices(conn: sqlite3.Connection):
    conn.execute("CREATE INDEX IF NOT EXISTS scc__spectrum ON spectra_candidate_scores(spectrum)")
    conn.execute("CREATE INDEX IF NOT EXISTS scc__molecule ON spectra_candidate_scores(candidate)")
    conn.execute("CREATE INDEX IF NOT EXISTS scc__scoring_method ON spectra_candidate_scores(scoring_method)")
    conn.execute("CREATE INDEX IF NOT EXISTS scc__dataset ON spectra_candidate_scores(dataset)")

    conn.execute("CREATE INDEX IF NOT EXISTS ma__merged_accession ON merged_accessions(merged_accession)")

    conn.execute("CREATE INDEX IF NOT EXISTS fpd__molecule ON fingerprints_data(molecule)")
    conn.execute("CREATE INDEX IF NOT EXISTS fpd__name ON fingerprints_data(name)")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("idir")
    arg_parser.add_argument("score_tgz_fn")
    arg_parser.add_argument(
        "--pubchem_db_fn",
        default="/home/bach/Documents/doctoral/projects/local_pubchem_db/db_files/pubchem_01-02-2021.sqlite",
        type=str,
        help="Filepath of the PubChem database.")
    args = arg_parser.parse_args()

    conn_mb = sqlite3.connect(args.massbank_db_fn)
    conn_pc = sqlite3.connect(args.pubchem_db_fn)

    prefix_pattern = re.compile(r"[A-Z]{2,3}")

    try:
        with conn_mb:
            conn_mb.execute("PRAGMA foreign_keys = ON")

        with conn_mb:
            create_tables(conn_mb)

            conn_mb.execute("INSERT OR IGNORE INTO fingerprints_meta "
                            "   VALUES (?, ?, ?, ?, DATETIME('now', 'localtime'), ?, ?, ?, ?)",
                            ("sirius_fps", "UNION", "binary", None, "CDK: None", 3047, False, None))

            conn_mb.execute("INSERT OR IGNORE INTO scoring_methods VALUES (?, ?, ?)",
                            ("sirius__sd__correct_mf", "SIRIUS: 4, CSI:FingerID: 1.4.5",
                             "Dr. Kai Dührkop ran CSI:FingerID to predict the candidate scores in a structure disjoint "
                             "(sd) fashion. That means, the ground truth molecular structured associated with the "
                             "Massbank spectra have not been used for training. The correct molecular formula was used "
                             "to construct the candidate sets."))

        spectra_idx = 0
        with tarfile.open(os.path.join(args.idir, args.score_tgz_fn), "r:gz") as score_archive:
            for entry in score_archive:
                if entry.isreg():
                    # Found a spectrum ...
                    spectra_idx += 1
                    spec_id = os.path.basename(entry.name).split(".")[0]  # e.g. AC01111385
                    LOGGER.info("Process spectrum %05d: id = %s" % (spectra_idx, spec_id))

                    # ==================================
                    # Find meta-information for spectrum
                    meta_info = None

                    ds_pref = re.findall(prefix_pattern, spec_id)[0]  # e.g. AC
                    for d in glob.iglob(os.path.join(args.idir, "%s*" % ds_pref)):
                        _df = pd.read_csv(os.path.join(d, "spectra_summary.tsv"), sep="\t")
                        try:
                            row_idx = _df["accession"].to_list().index(spec_id)
                            meta_info = _df.iloc[row_idx, :]
                            meta_info["dataset"] = ds_pref
                            break
                        except ValueError:
                            # Spectrum ID was not found in the dataset
                            continue

                    if meta_info is None:
                        # Fail here
                        raise RuntimeError("[%s] Meta-data for spectrum not found." % spec_id)
                    # ==================================

                    # =================================================================================
                    # Load candidates and with CSI:FingerID scores and combine with PubChem information
                    cands = pd.read_csv(score_archive.extractfile(entry), sep="\t")  # type: pd.DataFrame
                    cands = pd.merge(
                        left=cands,
                        right=pd.read_sql(
                            "SELECT * FROM compounds WHERE InChIKey_1 IN %s" % MassbankDB._in_sql(cands["inchikey"]),
                            conn_pc
                        ),
                        left_on="inchikey", right_on="InChIKey_1", how="inner", suffixes=("__SIRIUS", "__PUBCHEM")
                    )

                    if meta_info["pubchem_id"] not in cands["cid__PUBCHEM"]:
                        raise RuntimeError("[%s] Correct molecular structure (cid = %d, inchikey = %s) is not in the "
                                           "candidate set."
                                           % (spec_id, meta_info["pubchem_id"], meta_info["inchikey"]))

                    _unq_mf = cands["molecular_formula__PUBCHEM"].unique()
                    if meta_info["molecular_formula"] not in _unq_mf:
                        raise RuntimeError("[{}]Correct molecular formula is not on the candidate set: {} not in {}."
                                           .format(spec_id, meta_info["molecular_formula"], _unq_mf))

                    if len(_unq_mf > 1):
                        LOGGER.warning("[{}] There is more than one molecular formula in the candidate set: {}."
                                       .format(spec_id, _unq_mf))
                    # =================================================================================

                    # ===========================
                    # Insert new data into the DB
                    with conn_mb:
                        # Meta-information about the scored spectra
                        conn_mb.execute("INSERT OR IGNORE INTO scored_spectra_meta VALUES (?, ?, ?, ?, ?, ?, ?)",
                                        (
                                            meta_info["accession"], meta_info["dataset"],
                                            meta_info["original_accessions"], meta_info["precursor_mz"],
                                            meta_info["precursor_type"], meta_info["pubchem_id"],
                                            meta_info["retention_time"]
                                        ))

                        # Molecule structures associated with the candidates
                        conn_mb.executemany("INSERT OR IGNORE INTO molecules VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                            [
                                                (row["cid__PUBCHEM"], row["InChI__PUBCHEM"], row["InChIKey__PUBCHEM"],
                                                 row["InChIKey_1__PUBCHEM"], row["InChIKey_2__PUBCHEM"],
                                                 row["SMILES_ISO__PUBCHEM"], row["SMILES_CAN__PUBCHEM"],
                                                 row["exact_mass__PUBCHEM"], row["monoisotopic_mass__PUBCHEM"],
                                                 row["molecular_formula__PUBCHEM"], row["xlogp3"])
                                                for _, row in cands.iterrows()
                                            ])

                        # CSI:FingerID candidate scores
                        conn_mb.executemany(
                            "INSERT INTO spectra_candidate_scores VALUES (?, ?, ?, ?, ?)",
                            [
                                (meta_info["accession"], row["cid__PUBCHEM"], meta_info["dataset"],
                                 "sirius__sd__correct_mf", row["score__SIRIUS"])
                                for _, row in cands.iterrows()
                            ]
                        )

                        # Insert information about the merged Massbank accessions and their new ids
                        conn_mb.executemany("INSERT INTO merged_accessions VALUES (?, ?)",
                                            [
                                                (acc, meta_info["accession"])
                                                for acc in meta_info["original_accessions"].split(",")
                                            ])

                        # CSI:FingerID fingerprints as index strings
                        conn_mb.executemany(
                            "INSERT OR IGNORE INTO fingerprints_data VALUES (?, ?, ?)",
                            [
                                (row["cid__PUBCHEM"],
                                 ",".join(map(str, [idx for idx, fp in enumerate(row["fingerprint"]) if fp == "1"])),
                                 "sirius_fps")
                                for _, row in cands.iterrows()
                            ]
                        )
                    # ===========================

        with conn_mb:
            create_indices(conn_mb)

    except RuntimeError as err:
        traceback.print_exc()
        LOGGER.error(err)
    finally:
        conn_mb.close()
        conn_pc.close()

