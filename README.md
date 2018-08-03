# MultGraphModels
# Author: Elin Shaddox
# Contact: elin@rice.edu
#
# The provided Matlab files for Bayesian inference of multiple graphical models are 
# associated with the following publication:
#
# Shaddox, E., Stingo, F., Peterson, C.B., Jacobson, S., Cruickshank-Quinn, C., Kechris, K., # Bowler, R. and Vannucci, M. (2016). A Bayesian Approach for Learning Gene Networks 
# Underlying Disease Severity in COPD. Statistics in Biosciences, accepted.
#
# These scripts rely on functions from the Matlab code for G-wishart sampling provided by
# Hao #Wang at https://msu.edu/~haowang/ and are associated with the following publications
#
# Associated publications:
# H. Wang, Scaling It Up: Stochastic Search Structure Learning in Graphical Models Bayesian # Analysis 10 (2015): 351-377
#
# Wang, H. and Li, S. (2012). Efficient Gaussian graphical model determination
# under G-Wishart prior distributions. Electronic Journal of Statistics.
# 6: 168—198.
#
# Please cite all publications if you use this code. Thanks!
#
# OVERVIEW OF FILES ————————————————————————————————————————————
#
# Example_multiple_graphs_SSVS.m
# ——————————————————————————————————
# Basic example of running the MCMC sampler and generating results summaries on a simple 
# setting with 3 groups with identical dependence structure
#
# MCMC_multiple_graphs_SSVS_final.m
# ——————————————————————————————————
# Code for running the MCMC sampler
#
# calc_mrf_C.m
# ——————————————————————————————————
# Helper function for calculating the normalizing constant for the MRF prior
#
# generate_sim1_input..m
# ——————————————————————————————————
# Script to generate matrices similar to those used as input to the first Simulation
#
#
