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
import time

import pandas as pd
import logging
import traceback
import numpy as np

from typing import List, Dict, Union

from matchms import Spectrum
from matchms.similarity import ModifiedCosine
from matchms.filtering import normalize_intensities

from massbank2db.utils import get_precursor_mz
from massbank2db.spectrum import MBSpectrum

from utils import get_backup_db


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


def parse_cfm_spectrum(spec_str: List[str], meta_data: Dict, return_mbspectrum: bool = False) \
        -> Union[List[MBSpectrum], List[Spectrum]]:
    """
    Function that parses an CFM-ID spectrum representation (text-file) into a list of spectra

    :param spec_str: list of strings, lines of the CFM-ID spectrum file

    :param meta_data: dictionary, meta data stored in the spectrum object

    :param return_mbspectrum: boolean, indicating whether the spectrum object should be of type massbank2db.MBSpectrum
        rather than matchms.Spectrum.

    :return: list of spectrum objects, each energy (low, med and high) is stored in a separate spectrum object
    """
    # Read in the spectra file in and parse into three (3) MBSpectrum objects
    mzs = []
    ints = []
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

    # Create spectrum objects
    if return_mbspectrum:
        specs = []
        for _mzs, _ints in zip(mzs, ints):
            specs.append(MBSpectrum())
            specs[-1].set_mz(_mzs)
            specs[-1].set_int(_ints)

            for k, v in meta_data.items():
                specs[-1].set(k, v)
    else:
        specs = [Spectrum(np.array(_mzs), np.array(_ints), meta_data) for _mzs, _ints in zip(mzs, ints)]

    return specs


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("massbank_db_fn", help="Filepath of the Massbank database.")
    arg_parser.add_argument("idir", help="Input directory containing the measured MS2 spectra and meta-data.")

    arg_parser.add_argument(
        "--idir_pred_spec",
        help="Input directory containing the predicted MS2 spectra.",
        default=None,
        type=str
    )
    arg_parser.add_argument(
        "--build_unittest_db",
        action="store_true",
        help="Use this option to create a smaller database (subset of SIRIUS scores) that we can use for unittests."
    )
    args = arg_parser.parse_args()

    if args.idir_pred_spec is None:
        args.idir_pred_spec = args.idir

    sqlite3.register_adapter(np.int64, int)

    conn = get_backup_db(args.massbank_db_fn, exists="reuse", postfix="with_cfm_id")

    try:
        with conn:
            conn.execute("PRAGMA foreign_keys = ON")

        with conn:
            ms2scorer_name = os.path.split(args.idir.rstrip(os.path.sep))[-1]  # e.g. tool_output/cfm-id/ --> cfm-id
            ms2scorer_name = [ms2scorer_name + "__summed_up_sim__norm", ms2scorer_name + "__merge_pred_spec__norm"]
            description = [
                "The CFM-ID software (in version 2.0) was used to predict the MS2 spectra for all candidates \
                associated with the MS2 spectra in Massbank. We use the single energy (se) pre-trained models from \
                the online repository for the spectra prediction. The predictions are done in a structure disjoint \
                fashion. That means, the pre-trained model for each Massbank spectrum (and its candidates) is chosen \
                so that the ground truth structure was not part of the training. The similarity between the predicted \
                and measured MS2 spectra is computed using the modified cosine similarity (matchms Python package, \
                version 0.9.2). As in the original CFM-ID implementation, the similarity is computed for each energy \
                separately and the scores are summed up. If the measured MassBank entry was available at multiple \
                collision energies (i.e. we have multiple measured spectra for the same compound), MS2 spectra have \
                been merged using the 'mzClust_hclust' function from XCMS prior to the similarity scoring.",
                "The CFM-ID software (in version 2.0) was used to predict the MS2 spectra for all candidates \
                associated with the MS2 spectra in Massbank. We use the single energy (se) pre-trained models from \
                the online repository for the spectra prediction. The predictions are done in a structure disjoint \
                fashion. That means, the pre-trained model for each Massbank spectrum (and its candidates) is chosen \
                so that the ground truth structure was not part of the training. The similarity between the predicted \
                and measured MS2 spectra is computed using the modified cosine similarity (matchms Python package, \
                version 0.9.2). We merged the MS2 spectra for the three (3) collision energies CFM-ID is predicting. \
                If the measured MassBank entry was available at multiple collision energies (i.e. we have multiple \
                measured spectra for the same compound), MS2 spectra have been merged using the 'mzClust_hclust' \
                function from XCMS. The similarity between predicted and measured spectrum was computed on the merged \
                spectra."
            ]

            conn.executemany(
                "INSERT OR IGNORE INTO scoring_methods VALUES (?, ?, ?)",
                zip(ms2scorer_name, ["CFM-ID: 2.0"] * len(description), description)
            )

        # Load the list spectra for which candidate spectra have been predicted
        df = pd.read_csv(os.path.join(args.idir, "spec2model.tsv"), sep="\t")

        for idx, (spec_id, _, model_id, ion_mode) in df.iterrows():
            # if idx < 6550:
            #     continue

            dataset, accession, precursor_mz = conn.execute(
                "SELECT dataset, accession, precursor_mz FROM scored_spectra_meta WHERE accession IS ?", (spec_id, )
            ).fetchone()
            LOGGER.info(
                "Process spectrum %05d / %05d: id = %s, dataset = %s" % (idx + 1, len(df), spec_id, dataset)
            )

            if args.build_unittest_db \
                    and not conn.execute("SELECT accession FROM scored_spectra_meta WHERE accession IS ?", (spec_id,)) \
                    .fetchall():
                # In the unittest DB we only include MetFrag scores for spectra that are already in the DB
                continue

            # =======================================================================================
            # Load the measured spectrum and score it against the candidates
            with open(os.path.join(
                    args.idir, dataset, "%s.txt" % spec_id
            )) as spec_file:
                # We only use the first energy level. Remember, that we use the same (merged) spectrum from MassBank for
                # each collision energy anyway.
                spec_measured = parse_cfm_spectrum(spec_file.readlines(), {"precursor_mz": precursor_mz})[0]

            # Track the candidate scores
            ikey1s = []
            scores__summed_up_sim = []
            scores__merge_pred_spec = []

            # Get all molecule IDs of the candidates belonging to the current MB spectrum
            cur = conn.execute(
                "SELECT inchikey, monoisotopic_mass FROM candidates_spectra "
                "   INNER JOIN molecules m on m.cid = candidates_spectra.candidate"
                "   WHERE spectrum IS ?",
                (spec_id, )
            )

            _t_load = 0.0
            _t_merge = 0.0
            _t_sim_1 = 0.0
            _t_sim_2 = 0.0

            for ikey, monoisotopic_mass in cur:
                try:
                    # Load and parse the predicted candidate spectra
                    s = time.time()
                    with open(os.path.join(
                            args.idir_pred_spec,
                            "predicted_spectra",
                            "cv=%d__ion=%s" % (model_id, ion_mode),
                            "%s.log" % ikey
                    )) as spec_file:
                        specs_pred = parse_cfm_spectrum(
                            spec_file.readlines(),
                            {
                                "precursor_mz": get_precursor_mz(
                                    monoisotopic_mass, "[M+H]+" if ion_mode == "positive" else "[M-H]-"
                                ),
                                "accession": accession
                            },
                            return_mbspectrum=True
                        )
                    _t_load += (time.time() - s)

                    # Compute the spectra similarity (summed up similarity)
                    s = time.time()
                    scores__summed_up_sim.append(
                        np.sum([
                            ModifiedCosine().pair(
                                spec_measured,
                                normalize_intensities(
                                    Spectrum(
                                        np.array(spec.get_mz()),
                                        np.array(spec.get_int()),
                                        spec.get_meta_information()
                                    )
                                )
                            ).item()[0]
                            for spec in specs_pred
                        ])
                    )
                    _t_sim_1 += (time.time() - s)

                    # Merge the collision energies of the predicted spectra
                    s = time.time()
                    spec_pred = MBSpectrum.merge_spectra(
                        specs_pred, normalize_peaks_before_merge=False, rt_agg_fun=None
                    )
                    _t_merge += (time.time() - s)

                    # Convert the MBSpectrum into a matchms.Spectrum
                    spec_pred = Spectrum(
                        np.array(spec_pred.get_mz()), np.array(spec_pred.get_int()), spec_pred.get_meta_information()
                    )

                    ikey1s.append(ikey)

                    # Compute the spectra similarity (merged)
                    s = time.time()
                    scores__merge_pred_spec.append(ModifiedCosine().pair(spec_measured, spec_pred).item()[0])
                    _t_sim_2 += (time.time() - s)
                except FileNotFoundError:
                    LOGGER.warning(
                        "Cannot find predicted spectra for mol=%s, spec=%s, model_id=%d, ion_mode=%s"
                        % (ikey, spec_id, model_id, ion_mode)
                    )
                except AssertionError:
                    LOGGER.warning(
                        "Something went wrong while parsing the predicted spectrum: %s, %s, %s"
                        % (ikey, model_id, ion_mode)
                    )

            if len(scores__summed_up_sim) == 0:
                LOGGER.warning("[%s] Empty candidate set." % spec_id)
                continue

            # Normalize the scores per candidate set
            scores__summed_up_sim = np.array(scores__summed_up_sim)
            scores__merge_pred_spec = np.array(scores__merge_pred_spec)

            max_s = np.max(scores__summed_up_sim)
            if max_s > 0:
                scores__summed_up_sim /= max_s
            else:
                scores__summed_up_sim = np.ones_like(scores__summed_up_sim)

            max_s = np.max(scores__merge_pred_spec)
            if max_s > 0:
                scores__merge_pred_spec /= max_s
            else:
                scores__merge_pred_spec = np.ones_like(scores__merge_pred_spec)

            cands = pd.DataFrame({
                "inchikey1": ikey1s,
                "cfmid_score__summed_up_sim": scores__summed_up_sim.tolist(),
                "cfmid_score__merge_pred_spec": scores__merge_pred_spec.tolist()
            })

            LOGGER.debug("Loading the spectra file took: %.4fs" % (_t_load / len(cands)))
            LOGGER.debug("Computing the spectra similarity took (summed up): %.4fs" % (_t_sim_1 / len(cands)))
            LOGGER.debug("Merging the spectra took: %.4fs" % (_t_merge / len(cands)))
            LOGGER.debug("Computing the spectra similarity took (merged): %.4fs" % (_t_sim_2 / len(cands)))
            # =======================================================================================

            # =====================================================================================
            if len(cands) == 0:
                LOGGER.warning("[%s] Empty candidate set." % spec_id)
                continue
            else:
                LOGGER.info("[%s] Number of scored candidates (CFM-ID): %d" % (spec_id, len(cands)))

            # Here we add the stereo-isomers for each 2D structure
            cands = pd.merge(
                left=cands,
                right=pd.read_sql(
                    "SELECT * FROM molecules "
                    "   WHERE cid in ("
                    "      SELECT candidate FROM candidates_spectra WHERE spectrum IS '%s')" % spec_id,
                    con=conn
                ),
                on="inchikey1", how="inner"
            )
            LOGGER.info("[%s] Number of candidates (with added stereo-configurations): %d" % (spec_id, len(cands)))

            _gt_cid = conn.execute(
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
            with conn:
                # CFM-ID candidate scores (summed up similarities)
                conn.executemany(
                    "INSERT INTO spectra_candidate_scores VALUES (?, ?, ?, ?, ?)",
                    [
                        (spec_id, row["cid"], ms2scorer_name[0], dataset, row["cfmid_score__summed_up_sim"])
                        for _, row in cands.iterrows()
                    ]
                )

                # CFM-ID candidate scores (merged predicted spectra)
                conn.executemany(
                    "INSERT INTO spectra_candidate_scores VALUES (?, ?, ?, ?, ?)",
                    [
                        (spec_id, row["cid"], ms2scorer_name[1], dataset, row["cfmid_score__merge_pred_spec"])
                        for _, row in cands.iterrows()
                    ]
                )
            # ===========================
    except RuntimeError as err:
        traceback.print_exc()
        LOGGER.error(err)
    finally:
        conn.close()
