# Introduction to this script -----------
#the goal of this script is to identify differentially expressed transcripts (DETs)
#you should already know what pairwise comparisons are most important to you

# Load packages -----
library(sleuth) 

# set-up your model ----
#use relevel function to set a different intercept
groups <- relevel(groups, "untreated.B6")
design <- model.matrix(~groups)
#design <- model.matrix(~0+groups)
design

# Import Kallisto transcript counts into R using Sleuth ----
# construct sleuth object
mySleuth <- sleuth_prep(targets, design, target_mapping = Tx) 
head(mySleuth$data)

# fit a linear model to the data ----
so <- sleuth_fit(mySleuth, design)

# take a look at your models
models(so)
design_matrix(so)

so <- sleuth_wt(so, which_beta = "groupsNAMEHERE")
sleuth_live(so)

# You can also view the same graphs ----
plot_pca(so, pc_x=1L, pc_y=2L, 
         units="tpm", 
         color_by = groups, point_size = 5)

plot_sample_heatmap(so)

