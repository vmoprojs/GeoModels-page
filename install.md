---
layout: page
title: Installation Instructions
permalink: /install/
---


We currently are loaded in Github only. This means that for `GeoModels` installation you will need to previously install [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package if you do not have it installed yet:

`
install.packages(devtools)
library(devtools)
`

`devtools` lets you install pacakges from github since they need to be installed from source code.


We have developed two GeoModels version, one *standard* version and one that uses the `OpenCL` framework for parallel computing. The standard version can be installed in any operating system: Windows, OSX and Linux,

`
install_github("vmoprojs/GeoModels")
library(GeoModels)
`

and you are good to go. 

The OpenCL `GeoModels` version is currently supported for **OSX only**. It is installed with this code:

`
install_github("vmoprojs/GeoModels-OCL")
library(GeoModels)
`

