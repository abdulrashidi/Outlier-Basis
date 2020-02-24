# # -------------------------------------------------------------------------------
# FIGURES FOR PAPER
# COMMON - READ DATA
# EX1    - LOF COMPARISON
# EX2    - KNN COMPARISON
# EX3    - iFOREST COMPARISON
# -------------------------------------------------------------------------------

# CHANGE FOLDER TO DIRECTORY OF SCRIPTS AND DATA IN THE LINE BELOW
folder <- paste( getwd(), "/Supp_Mat/", sep="")
source(paste(folder, "Functions_For_Paper.R", sep=""))
library("ggplot2")
library("tsutils")
library("tidyr")

# -------------------------------------------------------------------------------
# COMMON - READ DATA
# -------------------------------------------------------------------------------
# folder <- paste(getwd(), "/Cluster_Output/EX_5/N3_2/", sep="" )
dat <- read.csv(paste(folder, "Data_For_Section_5.csv", sep=""))


# -------------------------------------------------------------------------------
# EX1    - LOF COMPARISON
# -------------------------------------------------------------------------------
out <- draw_nemenyi_plots(dat, method=1)
out$nemenyi$fpval

draw_box_plots(dat, method=1)

# -------------------------------------------------------------------------------
# EX2    - KNN COMPARISON
# -------------------------------------------------------------------------------
out <- draw_nemenyi_plots(dat, method=2)
out$nemenyi$fpval

draw_box_plots(dat, method=2)

# -------------------------------------------------------------------------------
# EX3    - iFOREST COMPARISON
# -------------------------------------------------------------------------------
out <- draw_nemenyi_plots(dat, method=3)
out$nemenyi$fpval

draw_box_plots(dat, method=3)
