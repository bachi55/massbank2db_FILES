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
import os


if __name__ == "__main__":
    # Read CLI arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument(
        "--massbank_record_dir", default="/home/bach/Documents/doctoral/data/MassBank-data_upstream"
    )
    args = arg_parser.parse_args()

    conn = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        res = conn.execute(
            "SELECT dataset, contributor, accession, original_accessions FROM scored_spectra_meta "
            "   INNER JOIN datasets d ON d.name = scored_spectra_meta.dataset")

        n_total = {}
        n_erroneous = {}

        for idx, (ds, contr, acc, orig_accs) in enumerate(res):
            all_orig_accs_still_in_mb = True
            for orig_acc in orig_accs.split(","):
                if not os.path.exists(os.path.join(args.massbank_record_dir, contr, orig_acc + ".txt")):
                    all_orig_accs_still_in_mb = False
                    break

            if not all_orig_accs_still_in_mb:
                try:
                    n_erroneous[ds] += 1
                except KeyError:
                    n_erroneous[ds] = 1

                with conn:
                    conn.execute("DELETE FROM scored_spectra_meta WHERE accession IS ?", (acc, ))

            try:
                n_total[ds] += 1
            except KeyError:
                n_total[ds] = 1

        print(n_total)
        print(n_erroneous)

    finally:
        conn.close()




