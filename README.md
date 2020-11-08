# catTSclassifier
R code for “Classification of Categorical Time Series Using the Spectral Envelope and Optimal Scalings” by Li, Bruce, and Cai (2020)

## Dependencies: 
Code was developed using R version 4.0.2 ("Taking Off Again"), so code
may not function properly on older versions of R.  The packages listed
in the demo file must also be installed prior to use: astsa, abind.

Code has been tested and developed for data of the following dimensions:
<= 500 time points (T), <= 20 training time series, <= 100 testing time series.
Larger datasets may require more RAM for processing. 

## Description:
1) '**functions.r**': contains all necessary functions to run simulations.
2) '**demo_case1.r**': demonstrates the proposed method for Case 1 in the simulation studies.
3) '**demo_case2.r**': demonstrates the proposed method for Case 1 in the simulation studies.
4) '**demo_case3.r**': demonstrates the proposed method for Case 1 in the simulation studies.
5) '**seql executables (folder)**': contains necessary C++ implementation of the sequence learner  (Ifrim and Wiuf, 2011).
