# LipidMM: A multi-source mixing model for plant wax lipids using Bayesian framework

## Basic
This is a multi-source mixing model in a Bayesian framework that provides vegetation reconstruction using plant wax n-alkane chain length distribution and δ13C of multiple chains. The model consists of priors that include user-defined source groups and their associated parametric distributions of n-alkane concentration and δ13C. The mixing process involves newly defined mixing fractions such as fractional leaf mass contribution (FLMC) that can be used in vegetation reconstruction, and fractional source contribution to a specific n-alkane homologue (FSCn). These mixing fractions can be used in the interpretation of n-alkane records. 3 case studies are given to demonstrate the utility of the model framework in the analysis of sedimentary n-alkane mixtures. The first two case studies demonstrate model utility in disentangling the mixing scenario of multiple plant sources. The third case study demonstrates model utility in estimating δ2H of mean annual precipitation with the additional model component that involves n-alkane δ2H.

## Case studies
The code presented here involves three case studies with distinct sets of priors, to illustrate how the model can be applied to the interpretation of sedimentary n-alkane records. 
Case Study 1 involves published long-chain n-alkane records from lake surface sediments of Lake Qinghai, China [(W. Liu et al., 2015)](https://doi.org/10.1016/j.orggeochem.2015.03.017). CS1 evaluates chain length distribution and δ13C of 3 n-alkane chains.
Case Study 2 involves published long-chain n-alkane records from lake surface sediments along a vegetation gradient in Cameroon, western Africa [(Garcin et al., 2014)](https://doi.org/10.1016/j.gca.2014.07.004). CS2 evaluates chain length distribution and δ13C of 3 n-alkane chains. 
Case Study 3 involves published long-chain n-alkane records from a marine core off the Zambezi River Mouth [(Y. V. Wang et al., 2013)](https://doi.org/10.1016/j.gca.2012.10.016). CS3 evaluates chain length distribution, δ13C, and δ2H of multiple n-alkane chains.
All case studies use published data from the literature.

## Data
The "data" folder consists of compiled empirical records of modern plant wax n-alkanes for the case studies. These records are used to inform the prior distributions used in the Bayesian framework. It also contains published data on the chronology and sedimentary sourcing patterns for interpretation of results in Case Study 3.

## Software requirements
The code is developed in R calling the standalone JAGS (Just Another Gibbs Sampler) program, which is required before the code can be run. The JAGS program can be downloaded via [this link](https://sourceforge.net/projects/mcmc-jags/). Please make sure to download the version that is appropriate for your operating system.

To call JAGS from RStudio, the R packages "rjags" and "R2jags" are also required.

Other R packages specified in the file headers help to visualize data.

## Related Publication
in revision: D. Yang and G.J. Bowen, Vegetation reconstruction using plant wax n-alkane chain length distribution and δ13C of multiple chains: A multi-source mixing model using Bayesian framework, [Climate of the Past (preprint)](https://egusphere.copernicus.org/preprints/2022/egusphere-2022-23/).

## How to cite this software
D. Yang (2022), LipidMM: A multi-source mixing model for plant wax lipids using Bayesian framework, [https://github.com/SPATIAL-Lab/LipidMM/](https://github.com/SPATIAL-Lab/LipidMM/).