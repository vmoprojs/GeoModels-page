---
title: "GeoModels Packakge"
date: ""
output: md_document
---

The GeoModels package provides a set of procedures for simulation, estimation and prediction of spatio-temporal random fields.




<!--
<a href="https://www.buymeacoffee.com/samanyougarg"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" target="_blank"></a>

## Live Demo
## [Hanuman](https://samanyougarg.com/hanuman)
![Hanuman](/Screenshots/hanuman.jpg "Hanuman Preview")
-->

## Features

The main features of the package are:

-   The type of data that can be modeled are:
        -     spatial data
        -     space time data with spatial location sites possibly changing over time
        -     bivariate spatial data with (possibly different) spatial location sites

-   Data can be defined on euclidean space or on a sphere of arbitrary radius

-   The random fields can have the following marginal distributions:
    -   Gaussian
    -   Skew-Gaussian
    -   Student T
    -   Gamma
    -   Weibull
    -   LogGaussian
    -   Binomial
    -   Negative binomial
    -   Wrapped-Gaussian for directional data
-   Parametric models for both regression and dependence analysis through covariance models
-   Parametric (bivariate) spatial and spatiotemporal covariance models, including Matern, Generalized Wendland, Gneiting model, bivariate Matern

-   Estimation methods:
    -   Pairwise likelihood (optional parallel computation with OpenCL)
    -   Full likelihood (when feasible)
-   Optimal (linear) prediction


## Webpage

-   Visit our website for more information: [link](https://vmoprojs.github.io/GeoModels-page/)

-   Please report any bugs, suggestions and/or improvements: [link](https://github.com/vmoprojs/GeoModels-page/issues)