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

# This script trains an XGBoost regressor to predict the observed frequency (log-count) of SARS-CoV-2 single AA variant mutations
# in a given protein from precomputed per-mutation features, optionally restricted to a structural protein region. Hyperparameters
# (max_depth, n_estimators, colsample_bytree) are tuned by 10-fold grid search, and generalisation performance is estimated
# by 10x10 nested cross-validation, reporting mean R^2 and MAE across outer folds. The chosen hyperparameters from a
# grid search on the full dataset are saved for downstream SHAP-based feature-importance analysis.

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import shap
from sklearn.model_selection import GridSearchCV, RepeatedKFold, cross_validate, KFold, cross_val_score, train_test_split
from sklearn.metrics import make_scorer, explained_variance_score, r2_score, mean_absolute_error
from xgboost import XGBRegressor
import joblib
from preprocess import protein_regions
import re


# Helper function to check if a mutation site is in a given region
def in_region(s,start,end):
	m = re.search(r'\d+', s)
	if m:
		site = int(m.group())
		return start <= site <= end
	return False

# Function to filter out mutations that do not belong to a particular region
def filter_region(df,protein,region):
	if region == 'all':
		return df
	else:
		protein_region = protein_regions[protein][region]
		start = protein_region[0]
		end = protein_region[1]
		return df[df['mutation'].apply(lambda x: in_region(x, start, end))]

# Load the single AA variant features for a given protein, filtering by region if necessary, and return the data frames for regression
def load_preprocess(protein,region):
	df = pd.read_csv(f'data/features/{protein}.csv',low_memory=False)
	df = df.rename(columns={'Unnamed: 0':'mutation'}) #rename index column to 'mutation'
	df = filter_region(df,protein,region)
	df.to_csv('filtered.csv',index=True) #save the resulting data frame as a sanity check
	df = df.drop(columns=['mutation']) #drop the index column
	# Features (X) and target (y)
	X = df.drop(columns=['count']) # remove mutation count from the features as that is what we are predicting
	y = np.log10(df['count']+1)#set the target to be the log of mutation count
	return X,y

def optimise_evaluate(X, y):
    np.random.seed(66) #set the random seed for reproducibility
    regressor = XGBRegressor(tree_method='hist', device='cpu') #our regression model is an XGBoost regressor that histogram-based and does not use the GPU

    # Hyperparemeter optimisation using grid search
    n_estimators = range(100, 1000, 200) 
    max_depth = range(1, 10, 1)
    colsample_bytree = np.linspace(0.1, 1, 10)
    param_grid = dict(max_depth=max_depth,
                      n_estimators=n_estimators,
                      colsample_bytree=colsample_bytree,
                      n_jobs = [1])

    inner_cv = KFold(n_splits=10, shuffle=True) #Inner cross-validation using 10 folds
    outer_cv = KFold(n_splits=10, shuffle=True) #Outer cross-validation using 10 folds

    # Inner cross-validation
    model = GridSearchCV(regressor,
                         param_grid,
                         scoring='neg_mean_squared_error',
                         n_jobs=64,
                         cv=inner_cv,
                         verbose=3)

    model.fit(X, y) #Fit the model on the full dataset while optimising the hyperparameters, we will use these hyperparameters to re-fit the model for SHAP visualisation later
    best_params = model.best_params_
    print(best_params)

    #Define custom evaluation metrics
    expl_var = make_scorer(explained_variance_score)
    rsquared = make_scorer(r2_score)
    mae = make_scorer(mean_absolute_error)

    scoring = {'r2': rsquared,
               'mae': mae}

    
    # Nested cross-validation gives an unbiased estimate of generalisation performance
    # because hyperparameter selection (inner CV) never sees the outer test fold.
    # For each outer train/test split:
    #   Inner CV runs on the training portion to select hyperparameters and fit the model
    #   The fitted model is scored on the held-out test portion (MAE and R^2)
    # We report the mean across outer folds as the headline estimate.
    outer_results = cross_validate(model, X=X, y=y, cv=outer_cv, scoring=scoring, n_jobs=4)
    outer_results = pd.DataFrame(outer_results)

    return outer_results, best_params

   
def main():
	protein = "Spike"
	#region = "SARS-CoV-2_Spike_S1_RBD"
	region = "all"
	X,y = load_preprocess(protein,region)
	raw_results, raw_params = optimise_evaluate(X, y)
	res = pd.DataFrame(raw_results).mean()[['test_r2', 'test_mae']] #Average error and R^2 estimates across all outer folds to produce headline estimates
	with open(f'xgb_results_{protein}_{region}.txt','w') as f:
		f.write(repr(res))
	print(res)
	joblib.dump(raw_params, f'xgb_params_{protein}_{region}.pkl') #Save the best hyperparameters that were chosen while fitting the model on the full dataset, so we can re-fit the model for SHAP visualisation later
	

if __name__ == '__main__':
    main()

