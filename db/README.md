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
