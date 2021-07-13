# MassBank DB

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
