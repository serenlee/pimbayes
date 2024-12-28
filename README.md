# R Package pimbayes

#### Comprehensive tools for Bayesian inference with the partially identified model


### How to use `pimbayes`

1. Install JAGS

[Source download](https://sourceforge.net/projects/mcmc-jags/files/) and [homepage](https://mcmc-jags.sourceforge.io/)

2. Install package
``` R
# install.packages("devtools")
library(devtools)
install_github("https://github.com/serenlee/pimbayes")
```

3. Setup julia and julia packages
   
A function `setup` includes the installation of necessary julia packages 
``` R
setup()
```

If a machine does not have `julia`, please use `installJulia = TRUE`.
``` R
setup(installJulia = TRUE)
```

It is also possible to specify the path of the installed julia binary.

``` R
setup(JULIA_HOME = "custom_julia_path")
```
