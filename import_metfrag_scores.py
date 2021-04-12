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
import gzip
import os
import pandas as pd
import glob
import logging
import traceback
import numpy as np

# ================
# Setup the Logger
LOGGER = logging.getLogger("Import MetFrag Scores")
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

CH = logging.StreamHandler()
CH.setLevel(logging.INFO)

FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)

LOGGER.addHandler(CH)
# ================


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("idir", help="Input directory containing the MetFrag scores.")
    args = arg_parser.parse_args()

    sqlite3.register_adapter(np.int64, int)

    conn_mb = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn_mb:
            conn_mb.execute("PRAGMA foreign_keys = ON")

        with conn_mb:
            ms2scorer_name = os.path.split(args.idir.rstrip(os.path.sep))[-1]  # e.g. tool_output/metfrag/ --> metfrag
            if ms2scorer_name == "metfrag":
                description = "The MetFrag software (in version 2.4.5) was used to calculate in-silico MS2 scores for " \
                              "predefined candidate sets. The 'FragmenterScore' was used as MS2 score, abs. mass dev. " \
                              "was set to 0.0001, rel. mass dev. was set to 5, and the maximimum tree depth was 2. If " \
                              "multiple collision energies where available, the corresponding normalized (maximum " \
                              "intensity equals 100) MS2 spectra have been merged using the 'mzClust_hclust' function " \
                              "from XCMS."
            elif ms2scorer_name == "metfrag__norm_after_merge":
                description = "The MetFrag software (in version 2.4.5) was used to calculate in-silico MS2 scores for " \
                              "predefined candidate sets. The 'FragmenterScore' was used as MS2 score, abs. mass dev. " \
                              "was set to 0.0001, rel. mass dev. was set to 5, and the maximimum tree depth was 2. If " \
                              "multiple collision energies where available, the corresponding MS2 spectra have been " \
                              "merged using the 'mzClust_hclust' function from XCMS. Spectra normalization was performed " \
                              "after the merging."
            else:
                raise ValueError("Invalid MetFrag input directory: '%s'.")

            conn_mb.execute(
                "INSERT OR REPLACE INTO scoring_methods VALUES (?, ?, ?)",
                (ms2scorer_name, "MetFrag: 2.4.5", description)
            )

        # Find scored candidate lists
        l_score_fns = sorted(glob.glob(os.path.join(args.idir, "*", "*.csv.gz")))

        for idx, score_fn in enumerate(l_score_fns):
            spec_id = os.path.basename(score_fn).split(os.path.extsep)[0]  # e.g. AC01111385
            dataset = score_fn.split(os.sep)[-2]  # e.g. AC_001
            LOGGER.info(
                "Process spectrum %05d / %05d: id = %s, dataset = %s" % (idx + 1, len(l_score_fns), spec_id, dataset)
            )

            # =====================================================================================
            # Load candidates and with MetFrag scores and combine with all stereo-isomers in the DB
            with gzip.open(score_fn) as score_file:
                cands = pd.read_csv(score_file)  # type: pd.DataFrame

            with gzip.open(score_fn.replace("csv.gz", "cands.gz")) as original_cands_file:
                orig_cands = pd.read_csv(original_cands_file, sep="|")  # type: pd.DataFrame

            if len(cands) < len(orig_cands):
                LOGGER.warning(
                    "[%s] Less candidates with MetFrag score than in the original candidate set: %d < %d" %
                    (spec_id, len(cands), len(orig_cands))
                )

            if len(cands) == 0:
                LOGGER.warning("[%s] Empty candidate set." % spec_id)
                continue
            else:
                LOGGER.info("[%s] Number of candidates (MetFrag): %d" % (spec_id, len(cands)))

            # Here we add the stereo-isomers for each 2D structure
            cands = pd.merge(
                left=cands,
                right=pd.read_sql(
                    "SELECT * FROM molecules "
                    "   WHERE cid in ("
                    "      SELECT candidate FROM candidates_spectra WHERE spectrum IS '%s')" % spec_id,
                    con=conn_mb
                ),
                left_on="InChIKey1", right_on="inchikey1", how="inner"
            )
            LOGGER.info("[%s] Number of candidates (with added stereo-configurations): %d" % (spec_id, len(cands)))

            _gt_cid = conn_mb.execute(
                "SELECT molecule FROM scored_spectra_meta WHERE accession IS ?", (spec_id, )
            ).fetchone()[0]

            if _gt_cid not in cands["cid"].to_list():
                LOGGER.error(
                    "[%s] Correct molecular structure (cid = %d) is not in the candidate set." % (spec_id, _gt_cid)
                )
                continue
            # =====================================================================================

            # ===========================
            # Insert new data into the DB
            with conn_mb:
                # MetFrag candidate scores
                conn_mb.executemany(
                    "INSERT INTO spectra_candidate_scores VALUES (?, ?, ?, ?, ?)",
                    [
                        (spec_id, row["cid"], ms2scorer_name, dataset, row["FragmenterScore"])
                        for _, row in cands.iterrows()
                    ]
                )
            # ===========================

    except RuntimeError as err:
        traceback.print_exc()
        LOGGER.error(err)
    finally:
        conn_mb.close()
