#  cohsmix

An R package called CohsMix for Covariates in Hidden Structure Using
Mixture models provides variational EM algorithms for the three proposed CohsMix models. 
Theses algorithms allow to detect hidden structure and estimate the model parameters.

Generally we can distinguish three features:

* Simulation: Simulates data according to the desired model parameters: number of nodes, class proportions, 
connectivity and covariates values.

* Estimation: Estimates the parameters and the hidden structure. If the number of classes of the network is 
unknown, the ICL algorithm can also be solicited by the estimation function. If the latter is used, the parameter 
estimates would be provided for all classes. Therefore, the choice of number of classes is recommended 
but not required.

* Representation: Plots the network with the estimated structure. Given the number of group strategy, a 
representation of the ICL evolution is available


Further details concerning the models (and estimation stategies) can be found in:
- Clustering based on random graph model embedding vertex features (2010): http://www.sciencedirect.com/science/article/pii/S0167865510000413
- Model Based approaches for uncovering Web structures (2010): http://stat.genopole.cnrs.fr/_media/publications/zanghi.pdf

An application using the model and estimation strategy approach: http://constellations.labs.exalead.com/
