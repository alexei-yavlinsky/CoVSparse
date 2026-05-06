########################################################################################################################
#                                                                                                                      #
#                      CoVSparse - a suite of Python tools for working with GISAID SARS-CoV-2 data                     #
#                                                                                                                      #
#Dr. Alexei Yavlinsky, Dr. Marina Escalera Zamudio, Cedric Tan, Ali Demirci, Prof. Francois Balloux, Dr. Lucy van Dorp #
#                                                                                                                      #
#                                              University College London                                               #
#                                                                                                                      #
########################################################################################################################

# The code below is based on the code accompanying the following publication:
# Intrahost dynamics, together with genetic and phenotypic effects predict the success of viral mutations
# Cedric CS Tan, Marina Escalera-Zamudio, Alexei Yavlinsky, Lucy van Dorp, Francois Balloux
# bioRxiv 2024.10.18.619070; doi: https://doi.org/10.1101/2024.10.18.619070

# Companion to predict.py that performs SHAP-based feature-importance analysis on the SARS-CoV-2 mutation-frequency model.
# Loads the hyperparameters saved by predict.py, refits an XGBoost regressor on the full feature set for the chosen protein
# and region, then uses shap.TreeExplainer to compute per-feature SHAP values. Outputs a SHAP summary plot as a PDF and prints
# features ranked by mean absolute SHAP value. A region can be selected via a command-line index into the protein's region list
# defaulting to the whole protein.

import matplotlib
matplotlib.use('Agg')

from xgboost import XGBRegressor
import shap
import pandas as pd
import joblib
import numpy as np
import matplotlib.pyplot as plt

from predict import filter_region
from preprocess import protein_regions
import sys

def f_imp(protein,region):
	params = joblib.load(f'xgb_params_{protein}_{region}.pkl') #Load the protein-specific XGBoost hyperparameters saved in predict.py
	
	# Reinitialize the model with the hyperparameters
	model = XGBRegressor(**params)
	# Load the single AA variant features for the given protein
	X = pd.read_csv(f'data/features/{protein}.csv',low_memory=False)
	X = X.rename(columns={'Unnamed: 0':'mutation'})
	X = filter_region(X,protein,region) #Filter by protein region if necessary
	X = X.drop(columns=['mutation']) #Drop the index column
	# Features and target
	y = np.log10(X['count']+1) #Create the target variable (the same as in predict.py)
	X = X.drop(columns=['count']) #Create the feature matrix (the same as in predict.py)
	model.fit(X,y) #Fit the model
	# Initialize SHAP Explainer
	explainer = shap.TreeExplainer(model)
	
	# Compute SHAP values
	shap_values = explainer(X)
	
	# Plot feature importance
	shap.summary_plot(shap_values, X)
	plt.savefig(f"shap_summary_{protein}_{region}.pdf")
	

	# Rank features by importance
	feature_importance = np.abs(shap_values.values).mean(axis=0)
	ranked_features = pd.DataFrame({
	    'Feature': X.columns,
	    'Importance': feature_importance
	})
	ranked_features = ranked_features.sort_values(by='Importance', ascending=False)
	
	print(ranked_features)
	
#Optional function for displaying the histogram of log mutation counts
def count_data():
	X = pd.read_csv('data/features/Spike.csv',low_memory=False)
	X = X.drop(columns=['Unnamed: 0'])
	y = np.log10(X['count']+1)
	plt.hist(y, bins=100, edgecolor='black')  # bins define the number of bins
	plt.title('Histogram of Counts')
	plt.xlabel('Log(mutation_count)')
	plt.ylabel('Number of samples')
	plt.show()
	
def main():
    protein = "Spike"
    region_list = list(protein_regions[protein].keys())
    if len(sys.argv) < 2:
        region = "all"
    else:
        region = region_list[int(sys.argv[1])]
    print(region)
    f_imp(protein,region)


if __name__ == '__main__':
    main()
