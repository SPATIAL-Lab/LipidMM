# LipidMM

## basic

This is a multi-source mixing model in a Bayesian framework that provides vegetation reconstruction using plant wax n-alkane chain length distribution and δ13C of multiple chains. The model consists of priors that include user-defined source groups and their associated parametric distributions of n-alkane concentration and δ13C. The mixing process involves newly defined mixing fractions such as fractional leaf mass contribution (FLMC) that can be used in vegetation reconstruction, and fractional source contribution to a specific n-alkane homologue (FSCn).

## Case studies

The code presented here involves two case studies with distinct sets of priors, to illustrate how the model can be applied to the interpretation of sedimentary n-alkane records. One involves published long-chain n-alkane records from lake surface sediments of Lake Qinghai, China [https://doi.org/10.1016/j.orggeochem.2015.03.017](Liu_et_al.,_2015). The other involves published long-chain n-alkane records from lake surface sediments along a vegetation gradient in Cameroon, western Africa [https://doi.org/10.1016/j.gca.2014.07.004](Garcin_et_al.,_2014). Both case studies use published data from the literature.

## Data

The "data" folder consists of pulished empirical records of plant wax n-alkanes selected for the case studies. These records are used to inform the prior distributions used in the Bayesian framework.

## Software requirement

The code is developed calling the standalone JAGS (Just Another Gibbs Sampler) program, which is required before the code can be run. The JAGS program can be downloaded via this link [https://sourceforge.net/projects/mcmc-jags/](https://sourceforge.net/projects/mcmc-jags/). Please make sure to download the version that is appropriate for your operating system.
To call JAGS from RStudio, the R packages "rjags" and "R2jags" are also required.
Other R packages specified in the file headers help to visualize data.
