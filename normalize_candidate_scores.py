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

from tqdm import tqdm

if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("ms2scorer", help="MS2-scoring method for which the scores should be normalized.")
    args = arg_parser.parse_args()

    conn = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        # Read meta-information of the MS2-scorer
        meta_info = list(conn.execute("SELECT * FROM scoring_methods WHERE name IS ?", (args.ms2scorer, )).fetchone())
        print(meta_info)

        # Insert meta-information for the normalized scores
        ms2scorer_name = meta_info[0] = (meta_info[0] + "__norm")
        with conn:
            conn.execute("INSERT OR REPLACE INTO scoring_methods VALUES (?, ?, ?)", tuple(meta_info))

        # Iterate over all spectra in the DB
        for spectrum, dataset in tqdm(conn.execute(
                "SELECT accession, dataset FROM scored_spectra_meta"
                "   INNER JOIN datasets d on d.name = scored_spectra_meta.dataset"
        ), total=conn.execute("SELECT count(*) FROM scored_spectra_meta").fetchone()[0]):
            # Load candidate scores for the specific spectrum
            res = conn.execute(
                "SELECT candidate, score FROM spectra_candidate_scores"
                "   WHERE spectrum IS ? AND scoring_method IS ?",
                (spectrum, args.ms2scorer)
            ).fetchall()

            if len(res) == 0:
                continue

            candidates, scores = zip(*res)

            # Normalize scores
            scores = np.array(scores)
            if args.ms2scorer.startswith("sirius"):
                # SIRIUS provides log-scores
                scores = scores - np.max(scores)
            else:
                # MetFrag
                max_score = np.max(scores)
                if max_score == 0:
                    scores = np.ones_like(scores)
                else:
                    scores = scores / np.max(scores)

            # Insert normalized scores
            with conn:
                conn.executemany(
                    "INSERT OR REPLACE INTO spectra_candidate_scores VALUES ('%s', ?, '%s', '%s', ?)" % (
                        spectrum, ms2scorer_name, dataset
                    ),
                    zip(candidates, scores)
                )
    finally:
        conn.close()




