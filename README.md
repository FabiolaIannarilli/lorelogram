The __lorelogram__ package provides functions to explore and graphically describe correlation structures in binary data.
The package can be downloaded in R (RCoreTeam 2018) via:

    if (!require("devtools")) install.packages("devtools", repos = "http://cran.us.r-project.org", 
                                               dependencies = "Imports")
    devtools::install_github("FabiolaIannarilli/lorelogram")
    library(lorelogram)
 
The package's vignette is available at: <https://fabiolaiannarilli.github.io/lorelogram/articles/lorelogram.html> ('Get started' tab).
    
Alternatively, to install the package and simultaneously build the vignette, use the following lines of code. This will require several minutes.

    if (!require("devtools")) install.packages("devtools", 
                                           repos = "http://cran.us.r-project.org", 
                                           dependencies = "Imports")
    devtools::install_github("FabiolaIannarilli/lorelogram", build_vignettes=TRUE) 
    library(lorelogram)
