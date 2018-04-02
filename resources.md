---
layout: page
title: Resources and Download
permalink: /resources/
---
<figure class="ampstart-image-with-heading  m0 relative mb4">
<amp-img src="{{ site.baseurl }}assets/images/about.jpg" width="600" height="400" layout="responsive" alt="" class="mb3"></amp-img>
<figcaption class="absolute right-0 bottom-0 left-0">
</figcaption>
</figure>

The GeoModels provides a set of procedures for simulation,  estimation  and prediction of spatio-temporal random fields.


The main features of the package are:


* The type of data that can be modeled are:
    *   spatial data
    *   space time data with fixed location spatial sites
    *   space time data with dynamical location spatial sites
    *   bivariate spatial data with fixed location spatial sites
    *   bivariate spatial data with different location spatial sites
* Data can be defined on euclidean space or on a sphere of arbitrary radius
* The random fields can have the following marginal distributions:
    *   On the real line
        *   Gaussian
        *   Skew-Gaussian
        *   Student T
        *   Two Piece T
        *   Two Piece Gaussian
        *   Logistic
    *   On the positive real line
        *   Gamma
        *   Weibull
        *   LogLogistic (under development)
        *   LogGaussian
    *   On positive natural numbers
        *   binomial
        *   negative binomial
    *   On bounded support
        *   Uniform  (under development)
        *   Wrapped-Gaussian for directional data
* Parametric models for both regression and dependence analysis through covariance models
* Parametric (bivariate) spatial and   spatiotemporal covariance models, including Matern, Generalized Wendland, Gneiting model, bivariate Matern
* Estimation methods:  Pairwise likelihood,  full likelihood  (when feasible)
* Optimal (linear) prediction
