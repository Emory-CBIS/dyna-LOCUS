#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 08:55:37 2022

@author: scarlett
"""

# import modules
import os
import pyreadr
import numpy as np
from sklearn.model_selection import ParameterGrid
from sklearn.decomposition import DictionaryLearning

# import environment variable
iteration = int(os.environ["iter"])
print(type(iteration))

# parameter grid
param_grid = {'type': [1], 'n_subj': [20, 50], 'sd': np.arange(0, 3.5, 0.5), 'n_sims': range(1, 101)}
param_grid_list = list(ParameterGrid(param_grid))

# parameters
type_s = param_grid_list[iteration]['type']
n_subj = param_grid_list[iteration]['n_subj']
sd = param_grid_list[iteration]['sd']
n_sims = param_grid_list[iteration]['n_sims']

# set working directory
os.chdir("/projects/guo_lab/neurostat/project/dFC/simulation/sim_Y/type" + str(type_s))
# load data
if sd == 0 or sd == 1 or sd == 2 or sd == 3:
	Y = pyreadr.read_r("Y_n" + str(n_subj) + "_sd" + str(int(sd)) + "_" + str(n_sims) + ".RData")
else:
	Y = pyreadr.read_r("Y_n" + str(n_subj) + "_sd" + str(sd) + "_" + str(n_sims) + ".RData")
# Y = pyreadr.read_r('/Users/scarlett/Dropbox/dFC/simulation/data/sim_Y/Y_n10_sd4_1.RData')
# extract the data frame from result dictionary
Y = Y['Y'].to_numpy()
# remove the mean of each column of Y
mean_Y = Y.mean(axis=1)
Y_demean = Y - mean_Y[:, np.newaxis]

# init dictionary learning model 
dict_learner = DictionaryLearning(n_components = 8, tol = 1e-4, alpha = 1)
# fit the learner
dict_learner.fit(np.transpose(Y_demean))

# get latent source matrix
S = np.transpose(dict_learner.transform(np.transpose(Y_demean)))
# get loading matrix
A = np.dot(Y, np.dot(np.transpose(S), np.linalg.inv(np.dot(S, np.transpose(S)))))

# save as csv file
os.chdir("/projects/guo_lab/neurostat/project/dFC/simulation/decomp/dl/type" + str(type_s))
np.savetxt("S_n" + str(n_subj) + "_sd" + str(sd) + "_" + str(n_sims) + ".csv",  S, delimiter = ",")
np.savetxt("A_n" + str(n_subj) + "_sd" + str(sd) + "_" + str(n_sims) + ".csv",  A, delimiter = ",")

