# RIS phase optimization for Near-Field 5G Positioning: CRLB Minimization

This article addresses near-field localization using Reconfigurable Intelligent Surfaces (RIS) in 5G systems, where Line-of-Sight (LOS) between the base station and the user is obstructed. We propose a RIS phase optimization method based on the minimization of the Cramér-Rao Lower Bound (CRLB) for position estimation. The main contributions of the article are: (1) the derivation of the CRLB to optimize the RIS configuration; (2) the application of the mentioned framework in near-field considering reflective RIS; and (3) a set of simulation experiments showing the accuracy improvement of the proposed RIS phase optimization with respect to the state-of-the-art methods.

## Folder structure

```
└── /src
    ├── /results_generation     # Generate results for each plot
    └── ⁄optimization_functions # Needed functions
```

## Basic Run
The codes that reproduce the results are organize in four folders, one for each plot: 
* ***LocError_SNR20_Nris40*** -> **RMSE vs ITERATIONS** & **SNR vs ITERATIONS**

    Executing SNR20_Nris40_main.m you can obtain the result (long waiting time)

    For obtaining directly the plot and values open *SNR20_Nr40.mat*

 * ***LocError_SNR20diffNris*** -> **RMSE vs ITERATIONS**

    First one need to execute *SNR20_diffNris_10to100_main.m* and then *SNR20_diffNris_plot*


* ***LocErrordiffSNR_Nris40*** -> **RMSE vs SNR**

    First one need to execute *diffSNR_Nris40_main.m* and then *diffSNR_Nris40_plot*
 
* ***SNRvsCRLBdiffNris*** -> **CRLB vs SNR** 
    Execute *CRBvsSNRdiffNris*
