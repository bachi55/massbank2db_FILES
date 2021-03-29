####
#
# The MIT License (MIT)
#
# Copyright 2020, 2021 Eric Bach <eric.bach@aalto.fi>
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

import os
import sqlite3
import tqdm
import argparse
import pandas as pd

from massbank2db.db import MassbankDB
from massbank2db.spectrum import MBSpectrum
from massbank2db.exceptions import UnsupportedPrecursorType


def write_ms_and_cand_files(ds, n_compounds):
    # Create output directory
    odir = os.path.join(args.odir, args.output_format.lower(), ds)
    os.makedirs(odir, exist_ok=True)

    # Handle different output formats
    if args.output_format.lower() == "sirius":
        merge_peak_lists = False
    elif args.output_format.lower() == "metfrag":
        merge_peak_lists = True
    else:
        raise ValueError("Unsupported output format: '%s'" % args.output_format)

    # Spectra summary
    df_summary = []

    with MassbankDB(args.massbank_db_fn) as mb_db:
        for idx, (mol, specs, cands) in tqdm.tqdm(enumerate(
                mb_db.iter_spectra(dataset=ds, grouped=True, return_candidates=args.return_candidates,
                                   pc_dbfn=args.pubchem_db_fn)),
                total=n_compounds):

            # Merge spectra meta-information (e.g. precursor-mz, retention-time) and peaks (if desired) into a single
            # spectrum.
            spec = MBSpectrum.merge_spectra(specs, merge_peak_lists=merge_peak_lists)

            if args.output_format.lower() == "sirius":
                # Output in SIRIUS (.ms) format
                output = spec.to_sirius_format(molecular_candidates=cands, add_gt_molecular_formula=True)
            elif args.output_format.lower() == "metfrag":
                # Output in MetFrag format
                try:
                    output = spec.to_metfrag_format(molecular_candidates=cands,
                                                    **{"MetFragScoreWeights": [1.0],
                                                       "MetFragScoreTypes": ["FragmenterScore"],
                                                       "LocalDatabasePath": "./",
                                                       "ResultsPath": "./",
                                                       "NumberThreads": 4,
                                                       "PeakListPath": "./"})
                except UnsupportedPrecursorType as err:
                    print(err)
                    output = {spec.get("accession") + ".failed": str(err)}
            else:
                raise ValueError("Unsupported output format: '%s'" % args.output_format)

            # Write output
            for k, v in output.items():
                if v is None:
                    continue

                with open(os.path.join(odir, k), "w") as ofile:
                    ofile.write(v)

            # Write out spectra summary
            df_summary.append([spec.get("accession"), spec.get("dataset"), ",".join(spec.get("original_accessions")),
                               spec.get("precursor_mz"), spec.get("precursor_type"),
                               spec.get("pubchem_id"),
                               spec.get("retention_time")])

    pd.DataFrame(df_summary, columns=["accession", "dataset", "original_accessions",
                                      "precursor_mz", "precursor_type",
                                      "pubchem_id",
                                      "retention_time"]) \
        .to_csv(os.path.join(odir, "spectra_summary.tsv"), sep="\t", index=False)


def parse_cli_arguments():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("output_format",
                            help="Set the output format for the desired In-silico tool, e.g. SIRIUS or MetFrag.",
                            choices=["sirius", "SIRIUS", "metfrag", "MetFrag", "METFRAG"])
    arg_parser.add_argument("odir", help="Path to the output directory.")

    arg_parser.add_argument(
        "--pubchem_db_fn",
        default="/home/bach/Documents/doctoral/projects/local_pubchem_db/db_files/pubchem_01-02-2021.sqlite",
        type=str,
        help="Filepath of the PubChem database.")
    arg_parser.add_argument(
        "--datasets",
        nargs="*",
        type=str,
        default=None,
        help="Names of the datasets to process. If not set all datasets are processed.")
    arg_parser.add_argument(
        "--return_candidates",
        type=str,
        choices=[None, "mf", "mz_gt", "mz_measured"],
        default=None,
        help="If set, the parameter controls the candidate set definition.")

    return arg_parser.parse_args()


if __name__ == "__main__":
    # Read command line interface arguments
    args = parse_cli_arguments()

    # Get list of all dataset in the DB
    con = sqlite3.connect(args.massbank_db_fn)
    ds_ncmp = con.execute("SELECT name, num_compounds FROM datasets").fetchall()
    con.close()

    # Process the datasets
    for idx, (ds, ncmp) in enumerate(ds_ncmp):
        if args.datasets is None or ds in args.datasets:
            print("\nProcess: %s (%d/%d) with %d compounds" % (ds, idx + 1, len(ds_ncmp), ncmp))
            write_ms_and_cand_files(ds, ncmp)

