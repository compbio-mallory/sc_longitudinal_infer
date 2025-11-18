#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("LACE")

library("LACE")

args = commandArgs(trailingOnly=TRUE)

load(args[1])
print(final_df)
inference1 = LACE(D=final_df) # running LACE with default values
print(inference1)
save(inference1,file=args[2])
longitudinal.tree = longitudinal.tree.plot(inference = inference1,
                                           clone_labels = NULL,
                                           filename=args[3],
                                           legend_position = "topright")
