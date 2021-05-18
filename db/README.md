# Database Generation

The tandem mass spectra (MS/MS) and their corresponding meta-data, such as ground truth annotation and retention 
time (RT), are extracted from the official [Massbank (MB)](https://github.com/MassBank/MassBank-data/) repository 
([release 2020.11](https://github.com/MassBank/MassBank-data/releases/tag/2020.11)). From that starting point, 
several steps where needed to generate the final evaluation database: 

1) Parse, standardize and cleanup the mass spectrometry data of the MB release 2020.11
2) Extract MS/MS data in the input formats of different in-silico spectra annotation tools, such as CSI:FingerID 
3) Organize the molecular candidate sets in the DB and add structure information from Pubchem 
4) Add molecule features used for Machine Learning, such as molecular fingerprints, to the DB. 
5) ....

In the following we give a detailed description if these steps.

## (1) Parsing and Cleaning up MB 2020.11

We use the '[massbank2db](https://github.com/bachi55/massbank2db)' Python package (in version 0.6.1) to build the 
initial database containing the MS/MS spectra plus meta-data. All data in the 




# Rebuilding the DB

## Import SIRIUS scores

SIRIUS scores serve as reference for candidate sets.

```bash
python import_sirius_scores.py \
  db/massbank__2020.11__v0.6.1.sqlite \
  tool_output/sirius \
  scores.tar.gz \
  --build_unittest_db
```

Outputs "DB_FILE.sqlite"

## Get Candidate Sets for MetFrag

Candidate sets for MetFrag are generated from the SIRIUS candidates sets.

```bash
python get_metfrag_candidates.py \
  db/DB_FILE.sqlite \
  tool_output/metfrag \
  --gzip
```

## Import MetFrag Scores

```bash
python import_metfrag_scores.py \
  db/DB_FILE.sqlite \
  tool_output/metfrag__norm_after_merge
```

## Compute Circular Fingerprints

```bash
python compute_circular_fingerprints.py \
  db/DB_FILE.sqlite \
  FCFP \
  --batch_size=100000 \
  --min_subs_freq=100
```

## Convert circular counting fingerprints to binary

```bash
python convert_to_binary_fingerprints.py \
  db/DB_FILE.sqlite \
  FCFP
```