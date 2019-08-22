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
library("tsutils")
library("tidyr")

# -------------------------------------------------------------------------------
# COMMON - READ DATA
# -------------------------------------------------------------------------------

dat <- read.csv(paste(folder, "Data_For_Section_5.csv", sep=""))

# -------------------------------------------------------------------------------
# EX1    - LOF COMPARISON
# -------------------------------------------------------------------------------
#savepdf("LOF_Nemenyi_1")
out <- draw_nemenyi_plots(dat, method=1)
#dev.off()
out$friedman
#savepdf("LOF_BoxPlot")
draw_box_plots(dat, method=1)
#dev.off()



# -------------------------------------------------------------------------------
# EX2    - KNN COMPARISON
# -------------------------------------------------------------------------------
#savepdf("KNN_Nemenyi_1")
out <- draw_nemenyi_plots(dat, method=2)
#dev.off()
out$friedman
#savepdf("KNN_BoxPlot")
draw_box_plots(dat, method=2)
#dev.off()


# -------------------------------------------------------------------------------
# EX3    - iFOREST COMPARISON
# -------------------------------------------------------------------------------
#epdf("IFOREST_Nemenyi_1")
out <- draw_nemenyi_plots(dat, method=3)
#dev.off()
out$friedman
#savepdf("IFOREST_BoxPlot")
draw_box_plots(dat, method=3)
#dev.off()

