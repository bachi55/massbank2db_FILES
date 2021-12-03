# Building the MassBank data

Here we describe the scripts used to generate / load the:

- Candidate sets 
- MS² matching scores (CFM-ID, MetFrag and SIRIUS)
- Molecular fingerprints and descriptors
- ...

into our SQLite database (DB). We start with the initial SQLite DB: ```db/massbank__2020.11__v0.6.1.sqlite```. This DB 
was generated as described in the Methods "Pre-processing pipeline for raw MassBank records" of the manuscript based on
the [MassBank release 2020.11](https://github.com/MassBank/MassBank-data/releases/tag/2020.11) using the ```massbank2db``` 
package (version 0.6.1).

## Generating the molecular candidate sets

The candidate sets where generated using the [SIRIUS software](https://bio.informatik.uni-jena.de/software/sirius/) by
Dr. Kai Dührkop (developer of SIRIUS). SIRIUS uses [PubChem](https://pubchem.ncbi.nlm.nih.gov/) as molecular structure 
DB and returns the candidate sets limited to molecules with the ground truth molecular formula of the particular MassBank
spectrum. It is important to note, that neither the GUI nor the CLI tool of SIRIUS was used for the candidate set and 
MS² score generation. Instead, the non-public internal SIRIUS library was used which allows the score prediction in a 
structure disjoint fashion. That means, for each MassBank spectrum a CSI:FingerID (prediction backend of SIRIUS) model
was used, that was _not_ trained using its ground truth structure. This setting was chosen to prevent overfitting. 

### Generating the "SIRIUS ready" ms-files

The following script call was used to generate the "SIRIUS ready" ms-files: 
```
python create_insilico_tool_outputs_for_database.py db/massbank__2020.11__v0.6.1.sqlite sirius tool_output
```
A directory ```tool_output/sirius``` will be created with sub-directory for each MassBank group (see Methods 
"Pre-processing pipeline for raw MassBank records") containing the ms-files (```*.ms```) for each group of original 
MassBank accessions (see Methods "Pre-computing the MS² matching scores"). For example "AU22543794" in "AU_001" relates
to the original MassBank accessions "[AU300907](https://massbank.eu/MassBank/RecordDisplay?id=AU300907)", "AU300908",
"AU300909", "AU300910" and "AU300911". The file ```tool_output/sirius/AU_001/AU22543794.ms``` can be directly loaded into
the SIRIUS software.

### Importing the SIRIUS candidates and MS² scores

By calling: 
```bash
python import_sirius_scores.py db/massbank__2020.11__v0.6.1.sqlite tool_output/sirius SIRIUS_PREDICTIONS.TAR.GZ \
  --pubchem_db_fn=/PATH/TO/LOCAL_PUBCHEM.SQLITE \
  --acc_to_be_removed_fn=grouped_accessions_to_be_removed.txt
```
a copy of our initial SQLite DB is generated (```db/massbank__with_sirius.sqlite```) and the following information is 
added to the database:

- all (spectrum, candidate)-pairs generated by SIRIUS
- all (spectrum, candidate, MS² scores) for SIRIUS
- enriched candidate sets (see Methods "Generating the molecular candidate sets")
- (optional, ```--include_sirius_fps```) the binary fingerprints for each candidate as used by SIRIUS

#### Candidate set enrichment

As SIRIUS does not return scores for stereoisomers we need to them manually to the candidate sets. For that, we perform 
an inner merge on first InChIKey part (e.g. **FMGSKLZLMKYGDP**-UHFFFAOYSA-N) of each candidate of the candidate set 
between the candidates provided by SIRIUS and a local copy of PubChem: 

| ROW | InChIKey | SIRIUS MS² score |
| --- | --- | --- |
| 1 | JYGXADMDTFJGBT | 0.3 |
| 2 | OIGNJSKKLXVSLS | 0.1 |
| 3 | OMFXVFTZEKFJBZ | 0.5 |

becomes

| IDX | InChIKey | SIRIUS MS² score |
| --- | --- | --- |
| 1.A | JYGXADMDTFJGBT-UHFFFAOYSA | 0.3 |
| 1.B | JYGXADMDTFJGBT-KZUKIWJVSA | 0.3 |
| 1.C | JYGXADMDTFJGBT-LHXCVJCSSA | 0.3 |
| 2.A | OIGNJSKKLXVSLS-UHFFFAOYSA | 0.1 |
| 2.B | OIGNJSKKLXVSLS-RDEQXLMJSA | 0.1 |
| 3.A | OMFXVFTZEKFJBZ-UHFFFAOYSA | 0.5 |
| 3.B | OMFXVFTZEKFJBZ-NCQLXYLKSA | 0.5 |

after the merge.

#### Removal of records associated with the [#152 pull-request in MassBank](https://github.com/MassBank/MassBank-data/pull/152)

As described in the Methods "Pre-processing pipeline for raw MassBank records" we remove a couple of MassBank records 
related to the "LU" datasets which where reported to have issues. For that, we compare which original "LU*" accessions
where removed from MassBank between release 2020.11 (our baseline) and release 2021.3. We list our internal accession IDs
in the file ```grouped_accessions_to_be_removed.txt```. Entries in this list will not be imported to ```massbank__with_sirius.sqlite```
database and hence are not part of our experiments. 

## MetFrag MS² scores

### Generating the scores using the MetFrag binary 



## ```compute_circular_fingerprints.py```

- 

## ```generate_train_test_splits.py```

Script to sample LC-MS2 data from the (MS2, RT)-tuples in the database. For each dataset (e.g. AU_001) we 
repeatedly sample (MS2, RT)-tuples, approximately 50 each time, which constitute the LC-MS2 data. This procedure 
ensures that ground truth annotation are available for the LC-MS2 data.

## ```compute_rdkit_features.py```

Script to compute all RDKit molecular descriptors for all candidate in the database. We use these features to build the 
retention time prediction model.