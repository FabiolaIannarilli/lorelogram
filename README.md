The __lorelogram__ package provides functions to explore and graphically describe correlation structures in binary data.
The package can be downloaded in R (RCoreTeam 2018) via:

    if (!require("devtools")) install.packages("devtools", repos = "http://cran.us.r-project.org", 
                                               dependencies = "Imports")
    devtools::install_github("FabiolaIannarilli/lorelogram")
    library(lorelogram)
    
As an alternative, users can download the package and simultaneously build the vignette using the following lines of code. However, 
this will take longer time to run.

    if (!require("devtools")) install.packages("devtools", 
                                           repos = "http://cran.us.r-project.org", 
                                           dependencies = "Imports")
    devtools::install_github("FabiolaIannarilli/lorelogram", build_vignettes=TRUE) 
    library(lorelogram)
