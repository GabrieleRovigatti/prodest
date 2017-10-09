---
title: 'Production Function Estimation in R: The prodest Package'
tags:
  - prodest
  - R software
  - production functions
  - productivity
authors:
 - name: Gabriele Rovigatti
   orcid: 0000-0003-4021-3070
   affiliation: 1
affiliations:
 - name: University of Chicago Booth School of Business
   index: 1
date: 10 July 2017
bibliography: paper.bib
---

# Summary

The Total Factor Productivity (TFP) - also called Multi-factor productivity - measures the change in output that cannot be accounted for 
by changes in the amounts of input.
The R package prodest provides functions for TFP estimation following the most widely-known methodologies using the 
control function approach. Focusing on Value Added production functions, it estimates the two--steps models presented by Olley--Pakes 
(1996) [@olley_etal96] and Levinshon--Petrin (2003) [@levinsohn_etal03], as well as their correction proposed by Ackerberg--Caves--Frazer 
(2015) [@ackerberg_etal15]. The system GMM framework proposed by Wooldridge (2009) [@wooldridge_09] is also implemented in two slightly 
different versions. 
Dealing with standard Cobb-Douglas technology in a panel framework, all methods assume that the productivity term evolves according to a
first-order Markov process and that a proxy variable exists - i.e., a function of state variables and productivity - invertible and
monotonically increasing in productivity. Exploiting these features and with different choices of the proxy variables, the methods yield 
consistent estimates of labor and capital inputs parameters, allowing for an immediate computation of TFP. 
The prodest package features also the Data Generating Process used by Ackerberg--Caves--Frazer (2015) [@ackerberg_etal15] and 
allows for the simulation of datasets according to several measurement errors and random shock variances. It can be used by 
practitioners for both running Monte Carlo simulations and benchmarking estimate results. 

# References
