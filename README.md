---
title: "GeoModels Packakge"
date: ""
output: md_document
---

GeoModels: Procedures for Gaussian and Non Gaussian Geostatistical (Large) Data Analysis




<!--
<a href="https://www.buymeacoffee.com/samanyougarg"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" target="_blank"></a>

## Live Demo
## [Hanuman](https://samanyougarg.com/hanuman)
![Hanuman](/Screenshots/hanuman.jpg "Hanuman Preview")
-->

# GeoModels Package

The **GeoModels** package provides a set of procedures for simulation, estimation and prediction of spatio-temporal random fields.

## Main Features

### Types of Data Modeled

- **Spatial data**
- **Space-time data** with spatial location sites possibly changing over time
- **Bivariate spatial data** with (possibly different) spatial location sites

### Data Domain

- Euclidean space
- Sphere of arbitrary radius

### Marginal Distributions of Random Fields

#### Continuous distributions (supported on the whole real line)
- Gaussian  
- Skew-Gaussian  
- Student T  
- Logistic  
- Sinh-arcsinh  
- Two-piece  

#### Positive continuous distributions
- Gamma  
- Weibull  
- LogGaussian  
- LogLogistic  

#### Discrete distributions
- Binomial  
- Negative Binomial  
- Poisson  

#### Circular distributions
- Wrapped-Gaussian  

### Parametric Models

- For **regression** and **dependence analysis** through covariance models

### Covariance Models

- Parametric (bivariate) spatial and spatiotemporal models, including:
  - Matern  
  - Generalized Wendland  
  - Gneiting model  
  - Bivariate Matern  

### Estimation Methods

- **Composite likelihood** based on pairs (optional parallel computation with OpenCL)
- **Full likelihood** (when feasible)

### Prediction

- **Optimal (local) linear prediction** (kriging)



## Webpage

-   Visit our website for more information: [link](https://vmoprojs.github.io/GeoModels-page/)

-   Please report any bugs, suggestions and/or improvements: [link](https://github.com/vmoprojs/GeoModels-page/issues)

-   CRAN: https://doi.org/10.32614/CRAN.package.GeoModels