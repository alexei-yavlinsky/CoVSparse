# CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data

Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp

University College London

[![DOI](https://zenodo.org/badge/1230811791.svg)](https://doi.org/10.5281/zenodo.20053707)  

## How to run this pipeline

Ensure that the `GISAID metadata_tsv.tar.xz` is located in the `data/` directory

Generate a filtered set of sequences for each protein:
`python preprocess.py`

Generate sparse matrices for each protein:
`python sparse_matrix.py`

Generate DMS features
`python dms.py`

Generate single AA variant features:
`python sav_features.py`

Train the XGBoost model:
`python predict.py`

Generate the SHAP plot and ranked feature importance:
`python feature_analysis.py`
