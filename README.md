## ```compute_circular_fingerprints.py```

- 

## ```generate_train_test_splits.py```

Script to sample LC-MS2 data from the (MS2, RT)-tuples in the database. For each dataset (e.g. AU_001) we 
repeatedly sample (MS2, RT)-tuples, approximately 50 each time, which constitute the LC-MS2 data. This procedure 
ensures that ground truth annotation are available for the LC-MS2 data.

## ```compute_rdkit_features.py```

Script to compute all RDKit molecular descriptors for all candidate in the database. We use these features to build the 
retention time prediction model.