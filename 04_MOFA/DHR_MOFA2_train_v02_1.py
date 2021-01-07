
######################################
## Train MOFA+ model with DHR data  ##
######################################

from mofapy2.run.entry_point import entry_point
import pandas as pd
import io
import requests # to download the online data

###############
## Load data ##
###############

# Two formats are allowed for the input data:

# Option 1: a nested list of matrices, where the first index refers to the view and the second index refers to the group.
#           samples are stored in the rows and features are stored in the columns.
# 			Missing values must be filled with NAs, including samples missing an entire view

# datadir = "/Users/ricard/data/mofaplus/test"
# views = ["0","1"]
# groups = ["0","1"]
# data = [None]*len(views)
# for m in range(len(views)):
#     data[m] = [None]*len(groups)
#     for g in range(len(groups)):
#         datafile = "%s/%s_%s.txt.gz" % (datadir, views[m], groups[g])
#         data[m][g] = pd.read_csv(datafile, header=None, sep=' ')

# Option 2: a data.frame with columns ["sample","feature","view","group","value"]
#           In this case there is no need to have missing values in the data.frame,
#           they will be automatically filled in when creating the corresponding matrices

# use  "Autoantibodies" "GCMS" "LCMS.neg" "LCMS.pos" "microb.16S" "PEA"  (default alphabetical order)

file = "/home/francesco.marabita/Projects/DHR/out/MOFA2_v02/DHR_MOFA2_v02_input.txt"
data = pd.read_csv(file, sep="\t")

lik = ["bernoulli","gaussian","gaussian","gaussian", "gaussian", "gaussian"]

###########################
## Initialise MOFA model ##
###########################


## (1) initialise the entry point ##
ent = entry_point()


## (2) Set data options ##
# - scale_groups: if groups have significantly different ranges, it is good practice to scale each group to unit variance
# - scale_views: if views have significantly different ranges, it is good practice to scale each view to unit variance
ent.set_data_options(
	scale_groups = False, 
	scale_views = False
)


## (3, option 1) Set data using the nested list of matrices format ##
# views_names = ["view1","view2"]
# groups_names = ["groupA","groupB"]

# samples_names nested list with length NGROUPS. Each entry g is a list with the sample names for the g-th group
# - if not provided, MOFA will fill it with default samples names
# samples_names = (...)

# features_names nested list with length NVIEWS. Each entry m is a list with the features names for the m-th view
# - if not provided, MOFA will fill it with default features names
# features_names = (...)

# ent.set_data_matrix(data, 
# 	views_names = views_names, 
# 	groups_names = groups_names, 
# 	samples_names = samples_names,   
# 	features_names = features_names
# )

# (3, option 2) Set data using a long data frame
ent.set_data_df(data, likelihoods = lik)


## (4) Set model options ##
# - factors: number of factors. Default is K=10
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli"). 
# 		Default is None, and they are infered automatically
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_factors: use ARD prior in the factors? (TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (TRUE if using multiple views)

# Simple (using default values)
# ent.set_model_options()

# Advanced (using personalised values)
ent.set_model_options(
	factors = 20,
	spikeslab_weights = True,
	ard_factors = True, 
	ard_weights = True
)


## (5) Set training options ##
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow". 
#		For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training.
# 		Default is None, inactive factors are not dropped during training
# - gpu_mode: use GPU mode? this needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed

# Simple (using default values)
# ent.set_train_options()

# Advanced (using personalised values)
ent.set_train_options(
	iter = 10000, 
	convergence_mode = "medium",
	dropR2= 0.02,
	gpu_mode = False, 
	verbose = False
)


## (6, optional) Set stochastic inference options##
# This is only recommended with very large sample size (>1e4)
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.25 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)

# Simple (using default values)
# ent.set_stochastic_options()

# Advanced (using personalised values)
# ent.set_stochastic_options(batch_size=0.5, learning_rate=0.75, forgetting_rate=0.5, start_stochastic=10)


####################################
## Build and train the MOFA model ##
####################################

# Build the model 
ent.build()

# Run the model
ent.run()

##################################################################
## (Optional) do dimensionality reduction from the MOFA factors ##
##################################################################

# ent.umap()
# ent.tsne()

####################
## Save the model ##
####################

outfile = "/home/francesco.marabita/Projects/DHR/out/MOFA2_v02/DHR_MOFA2_v02_train_1.hdf5"

# - save_data: logical indicating whether to save the training data in the hdf5 file.
# this is useful for some downstream analysis in R, but it can take a lot of disk space.
ent.save(outfile, save_data=True)
