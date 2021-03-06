# Input Files for different In-silico Structure Identification Tools

Directory contains the tandem mass-spectra (MS2) from the [Massbank DB](../db/massbank__2020.11__v0.6.1.sqlite). MS2 
spectra from the same dataset, corresponding to the same ground truth structure and measured under the same 
conditions have been combined. Depending on the in-silico tool, e.g. SIRIUS or MetFrag, different combining 
strategies are used:
- **SIRIUS**: MS2 are separately added to a single ms-file. The software performs the spectra merging.
- **MetFrag**: MS2 spectra are merged, using hierarchical clustering, into a single spectrum. 

## Molecular Candidates

The candidate structures for each MS2 spectrum is defined using the ground truth molecular formula. That means, only 
those molecular structures, extract from PubChem, are used that have the same molecular formula as the ground truth 
structures. The SIRIUS software was used to create the candidate lists. 

SIRIUS outputs a single candidate for all associated stereo-isomers: "VUQLHQFKACOHNZ-UHFFFAOYSA-N", 
"VUQLHQFKACOHNZ-LURJTMIESA-N" and "VUQLHQFKACOHNZ-ZCFIWIBFSA-N" --> "VUQLHQFKACOHNZ". While adding the candidate 
scores to the database, we use a [local PubChem copy](https://github.com/bachi55/local_pubchem_db) to retrieve all 
stereo isomers for a candidate. To each of them the same fragmenter score is assigned, as MS2 cannot resolve 
stereo-chemistry.

## In-silico Tools

### SIRIUS

Dr. Kai Dührkop ran CSI:FingerID to predict the candidate scores in a structure disjoint (sd) fashion. That means, 
the ground truth molecular structured associated with the Massbank spectra have not been used for training.
The correct molecular formula was used to construct the candidate sets.

### MetFrag



#### 

The [MetFrag software](https://ipb-halle.github.io/MetFrag/projects/metfragcl/) (in version 2.4.5) was used to 
calculate in-silico MS2 scores for predefined candidate sets. The 'FragmenterScore' was used as MS2 score, abs. mass 
dev. was set to 0.0001, rel. mass dev. was set to 5, and the maximum tree depth was 2. If multiple collision energies 
where available, the corresponding normalized (maximum intensity equals 100) MS2 spectra have been merged using the '
[mzClust_hclust](https://github.com/sneumann/xcms/blob/master/src/mzClust_hclust.c)' function from [XCMS](https://github.com/sneumann/xcms).

### MetFrag