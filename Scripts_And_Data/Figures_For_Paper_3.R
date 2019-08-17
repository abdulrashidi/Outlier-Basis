# # -------------------------------------------------------------------------------
# FIGURES FOR PAPER
# COMMON - READ DATA
# EX1    - LOF COMPARISON
# EX2    - KNN COMPARISON
# EX3    - iFOREST COMPARISON
# -------------------------------------------------------------------------------

# CHANGE FOLDER TO DIRECTORY OF SCRIPTS AND DATA IN THE LINE BELOW
folder <- paste( getwd(), "/paper/", sep="")
source(paste(folder, "Functions_For_Paper.R", sep=""))
library("ggplot2")
library("tidyverse")

# -------------------------------------------------------------------------------
# COMMON - READ DATA
# -------------------------------------------------------------------------------

dat <- read.csv(paste(folder, "Data_For_Section_5.csv", sep=""))

# -------------------------------------------------------------------------------
# EX1    - LOF COMPARISON
# -------------------------------------------------------------------------------

out <- draw_nemenyi_plots(dat, method=1)
out$friedman
draw_density_plots(dat, method=1)


# -------------------------------------------------------------------------------
# EX2    - KNN COMPARISON
# -------------------------------------------------------------------------------
out <- draw_nemenyi_plots(dat, method=2)
out$friedman
draw_density_plots(dat,method=2)


# -------------------------------------------------------------------------------
# EX3    - iFOREST COMPARISON
# -------------------------------------------------------------------------------
out <- draw_nemenyi_plots(dat, method=3)
out$friedman
draw_density_plots(dat,method=3)
