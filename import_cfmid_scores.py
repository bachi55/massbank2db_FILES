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
import tarfile
import os
import pandas as pd
import logging
import traceback
import numpy as np

from typing import List, Dict

from matchms import Spectrum
from matchms.similarity import ModifiedCosine
from matchms.filtering import normalize_intensities

from massbank2db.utils import get_precursor_mz




# ================
# Setup the Logger
LOGGER = logging.getLogger("Import CFM-ID Scores")
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

CH = logging.StreamHandler()
CH.setLevel(logging.INFO)

FORMATTER = logging.Formatter('[%(levelname)s] %(name)s : %(message)s')
CH.setFormatter(FORMATTER)

LOGGER.addHandler(CH)
# ================


def parse_cfm_spectrum(spec_str: List[str], meta_data: Dict):
    """

    :param spec_str:
    :param meta_data:
    :return:
    """
    # Read in the spectra file in and parse into three (3) matchms.Spectrum objects
    mzs = []
    ints = []
    i = 0
    for row in spec_str:
        if row == "\n":
            continue
        elif row.startswith("energy"):
            mzs.append([])
            ints.append([])
        else:
            _mz, _int = map(float, row.strip().split(" "))
            mzs[-1].append(_mz)
            ints[-1].append(_int)

    assert len(mzs) == 3, "CFM-ID predicts the spectra for three (3) collision energies."

    specs = [Spectrum(np.array(_mzs), np.array(_ints), meta_data) for _mzs, _ints in zip(mzs, ints)]

    return specs


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("idir", help="Input directory containing the predicted MS2 spectra.")
    arg_parser.add_argument(
        "--build_unittest_db",
        action="store_true",
        help="Use this option to create a smaller database (subset of SIRIUS scores) that we can use for unittests."
    )
    args = arg_parser.parse_args()

    sqlite3.register_adapter(np.int64, int)

    conn_mb = sqlite3.connect(args.massbank_db_fn)

    try:
        with conn_mb:
            conn_mb.execute("PRAGMA foreign_keys = ON")

        with conn_mb:
            ms2scorer_name = os.path.split(args.idir.rstrip(os.path.sep))[-1]  # e.g. tool_output/metfrag/ --> metfrag
            description = \
                "The CFM-ID software (in version 2.0) was used to predict the MS2 spectra for all candidates " \
                "associated with the MS2 spectra in Massbank. We use the single energy (se) pre-trained models from " \
                "the online repository for the spectra prediction. The predictions are done in a structure disjoint " \
                "fashion. That means, the pre-trained model for each Massbank spectrum (and its candidates) is chosen " \
                "so that the ground truth structure was not part of the training. The similarity between the predicted " \
                "and measured MS2 spectra is computed using the modified cosine similarity (matchms Python package, " \
                "version 0.9.2). As in the original CFM-ID implementation, the similarity is computed for each energy " \
                "separately and the scores are summed up. If multiple collision energies where available, the " \
                "corresponding normalized (maximum intensity equals 100) MS2 spectra have been merged using the " \
                "'mzClust_hclust' function from XCMS prior to the similarity scoring."

            conn_mb.execute(
                "INSERT OR REPLACE INTO scoring_methods VALUES (?, ?, ?)",
                (ms2scorer_name, "CFM-ID: 2.0", description)
            )

        # Load the list spectra for which candidate spectra have been predicted
        df = pd.read_csv(os.path.join(args.idir, "spec2model.tsv"), sep="\t")

        for idx, (spec_id, _, model_id, ion_mode) in df.iterrows():
            dataset, precursor_mz = conn_mb.execute(
                "SELECT dataset, precursor_mz FROM scored_spectra_meta WHERE accession IS ?", (spec_id, )
            ).fetchone()
            LOGGER.info(
                "Process spectrum %05d / %05d: id = %s, dataset = %s" % (idx + 1, len(df), spec_id, dataset)
            )

            if args.build_unittest_db \
                    and not conn_mb.execute("SELECT accession FROM scored_spectra_meta WHERE accession IS ?", (spec_id, )) \
                    .fetchall():
                # In the unittest DB we only include MetFrag scores for spectra that are already in the DB
                continue

            # =======================================================================================
            # Load the measured spectrum
            with open(os.path.join(
                    args.idir, dataset, "%s.txt" % spec_id
            )) as spec_file:
                specs_measured = parse_cfm_spectrum(spec_file.readlines(), {"precursor_mz": precursor_mz})

            # Get all molecule IDs of the candidates belonging to the current MB spectrum
            cur = conn_mb.execute(
                "SELECT distinct(inchikey1), monoisotopic_mass FROM candidates_spectra "
                "   INNER JOIN molecules m on m.cid = candidates_spectra.candidate"
                "   WHERE spectrum IS ?", (spec_id, )
            )
            for ikey1, monoisotopic_mass in cur:
                try:
                    # Load and parse the predicted candidate spectra
                    with open(os.path.join(
                        args.idir, "predicted_spectra", "cv=%d__ion=%s" % (model_id, ion_mode), "%s.log" % ikey1
                    )) as spec_file:
                        specs_pred = parse_cfm_spectrum(
                            spec_file.readlines(),
                            {
                                "precursor_mz": get_precursor_mz(
                                    monoisotopic_mass, "[M+H]+" if ion_mode == "positive" else "[M-H]-"
                                )
                            }
                        )

                    # Compute the spectra similarity
                    sim = np.mean([
                        ModifiedCosine().pair(specs_measured[en], normalize_intensities(specs_pred[en])).item()[0]
                        for en in range(3)
                    ])


                except FileNotFoundError:
                    LOGGER.warning(
                        "Cannot find predicted spectra for mol=%s, spec=%s, model_id=%d, ion_mode=%s"
                        % (ikey1, spec_id, model_id, ion_mode)
                    )
            # =======================================================================================

            # # Open the archive containing the predicted candidate spectra for the current MB spectrum
            # with tarfile.open(
            #     os.path.join(args.idir, "predicted_spectra", "cv=%d__ion=%s.tar.gz" % (model_id, ion_mode)),
            #     mode="r:gz"
            # ) as pred_cand_spectra_archive:
            #     for idx, pred_cand_spec in enumerate(pred_cand_spectra_archive):
            #         print(idx + 1)
            #
            #         # The molecule identifier is given by the predicted spectrum's filename
            #         ikey1 = os.path.basename(pred_cand_spec.name).split(os.extsep)[0]
            #
            #         # Consider only the relevant candidate spectra
            #         if ikey1 not in candidate_set:
            #             continue
            #         else:
            #             print("found")
            #
            #         # Read the spectrum
            #         spec_str = pred_cand_spectra_archive.extractfile(os.extsep.join([ikey1, "log"])).readlines()
            #         print(spec_str)
            #
            #
            #         # for l in spec_str:
            #         #     if l.startswith("energy"):







            # # =====================================================================================
            # # Load candidates and with MetFrag scores and combine with all stereo-isomers in the DB
            # with gzip.open(score_fn) as score_file:
            #     cands = pd.read_csv(score_file)  # type: pd.DataFrame
            #
            # with gzip.open(score_fn.replace("csv.gz", "cands.gz")) as original_cands_file:
            #     orig_cands = pd.read_csv(original_cands_file, sep="|")  # type: pd.DataFrame
            #
            # if len(cands) < len(orig_cands):
            #     LOGGER.warning(
            #         "[%s] Less candidates with MetFrag score than in the original candidate set: %d < %d" %
            #         (spec_id, len(cands), len(orig_cands))
            #     )
            #
            # if len(cands) == 0:
            #     LOGGER.warning("[%s] Empty candidate set." % spec_id)
            #     continue
            # else:
            #     LOGGER.info("[%s] Number of candidates (MetFrag): %d" % (spec_id, len(cands)))
            #
            # # Here we add the stereo-isomers for each 2D structure
            # cands = pd.merge(
            #     left=cands,
            #     right=pd.read_sql(
            #         "SELECT * FROM molecules "
            #         "   WHERE cid in ("
            #         "      SELECT candidate FROM candidates_spectra WHERE spectrum IS '%s')" % spec_id,
            #         con=conn_mb
            #     ),
            #     left_on="InChIKey1", right_on="inchikey1", how="inner"
            # )
            # LOGGER.info("[%s] Number of candidates (with added stereo-configurations): %d" % (spec_id, len(cands)))
            #
            # _gt_cid = conn_mb.execute(
            #     "SELECT molecule FROM scored_spectra_meta WHERE accession IS ?", (spec_id, )
            # ).fetchone()[0]
            #
            # if _gt_cid not in cands["cid"].to_list():
            #     LOGGER.error(
            #         "[%s] Correct molecular structure (cid = %d) is not in the candidate set." % (spec_id, _gt_cid)
            #     )
            #     continue
            # # =====================================================================================
            #
            # # ===========================
            # # Insert new data into the DB
            # with conn_mb:
            #     # MetFrag candidate scores
            #     conn_mb.executemany(
            #         "INSERT INTO spectra_candidate_scores VALUES (?, ?, ?, ?, ?)",
            #         [
            #             (spec_id, row["cid"], ms2scorer_name, dataset, row["FragmenterScore"])
            #             for _, row in cands.iterrows()
            #         ]
            #     )
            # # ===========================

    except RuntimeError as err:
        traceback.print_exc()
        LOGGER.error(err)
    finally:
        conn_mb.close()
