# MassBank DB

## "with_only_stereo_splits" (v8, 08.10.2021)

- added lc-ms data splits where all ground truth molecules in the
  test sets do have a stereo annotation

## "with_2D_FCFP" (v7, 24.09.2021)

- added FCFP fingerprints (count and binarized) calculated suing the 
  canonical SMILES strings

## "with_cfm_id" (v6, 24.09.2021)

- CFM-ID scores have been added (normalized)
- We compute the similarity between the predicted and measured spectra 
  using the modified cosine similarity:
  * compute the similarity for each predicted energy separately and sum up
    the consine similarities

## "with_classyfire" (v5, 24.08.2021)

- added the Classyfire molecule class classification for the groud
  truth structures of the MassBank entries (not candidates)

## "with_descriptors" (v4, ??.07.2021)

- added Bouwmeester et al. (2019) molecular descriptors from RDKit
- descriptors are added for 'smiles_iso' and 'smiles_can'

## "with_test_splits" (v3, 12.07.2021)

- random splits test splits have been generated for all datasets
  in the DB
- splits were generated with "rt > 3 * column_dead_time"
- splits added to the 'lcms_data_splits' table

## "only_normalized_scores" (v2, 26.05.2021)

- SIRIUS and MetFrag scores are normalized within each candidate set
  - SIRIUS: score - max_score (log-space)
  - MetFRag: score / max_score

## "with_metfrag" (v1, ??.??.????)

- still contains the un-normalized MS2 scores 
- 
