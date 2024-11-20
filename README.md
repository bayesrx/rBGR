# rBGR
rBGR is an R package that implements a flexible Bayesian model to construct heterogeneous graphs under non-normal continuous data. For more details, please see Yao et al. (2023+)  Robust Bayesian Graphical Regression Models for Assessing Tumor Heterogeneity in Proteomic Networks on [arXiv](https://arxiv.org/abs/2310.18474).

# Manual
The `rBGR` accommodates the non-normality by the random scales and builds graphs through graphical regressions. The coefficients of graphical regression incorporate the subject-specific information and encode the graph edges by the zero coefficients. Due to the formulation of graphical regression, `rBGR` obtains posterior samples of coefficients by Gibbs sampler and infers the random scales by the Metropolis-Hasting algorithm. We refer to more details in Yao et al. (2023+) Section 4. 

## Main MCMC Function
In the package, we wrap the MCMC algorithm in a function of `rBGR_mcmc_Int()`. Given proteomic expression data and the subject-specific information (immune component abundance in the Main Paper), `rBGR` regresses one variable on the rest of the variables. Given the covariates, users need to execute the same function on all variables to construct the whole graph. Fortunately, this algorithm can be run parallelly to reduce the computation time. 
 
The MCMC function `rBGR_mcmc_Int()` takes several arguments with different options for users to control the algorithm. To run the MCMC, users are required to specify the data by the following three arguments: (i) regressand by the argument $Y$, (ii) regressor by the argument $X$, and (iii) the covariate by the argument $U$. The function `rBGR_mcmc_Int()` also allows additional options to control the model, such as $N$ for the number of iterations, `burnin` for the number of iterations to be discarded, and `seed_` for the initial seeds. 

The function `rBGR_mcmc_Int()` produces the posterior samples of (i) two components of coefficients $\xi_{j,k,h}$ and $\eta_{j,k,h}$, (ii) threshold parameter $t_j$ and random scales $d_{ij}$. Given the covariates, users can obtain the posterior coefficients by $\alpha_{j,k,h}=\eta_{j,k,h}\xi_{j,k,h}$ and edges by $\beta_{j,k}( \mathbf{X_i})=\theta_{j,k}(\mathbf{X_i}) \lvert \theta_{j,k}(\mathbf{X_i}) \rvert > t_j  $, where $\theta_{j,k}(\mathbf{X_i})=\sum_{h=1}^q \alpha_{j,k,h}X_{ih}$.

 
## Summarizing Posterior Samples
Once we obtain posterior samples of coefficients and edges, users can symmetrize and summarize the results to obtain undirected graphs in two different levels: population ($\alpha_{j,k,h}$) and individual ( $\beta_{j,k}( \mathbf{X_i})$ ) levels. We offer codes in the package to demonstrate our symmetrizing and summarizing algorithm with the data used in Yao et al. (2023+).

For the population-level graph, users can run `pstSmpExt_pop.R` to extract the population-level information and visualize the results through the `PlotRes_popLevel.R`. Similarly, in the individual-level graph, we first obtain the symmetrized edges by executing `pstSmpExt_ind.R` and then visualize the results by `PlotRes_indQuantile.R`. Currently, we use the 5, 25, 50, 75, and 95-th percentiles of the covariates as five different individuals, as shown in Yao et al. (2023+).
