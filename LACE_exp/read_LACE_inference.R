# This file gets the required values from the LACE inference. 
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
load(args[1])
C_inf = inference1['C']
clones_summary = inference1['clones_summary']

#print(C_inf)
#print(clones_summary)
#str(C_inf)
#str(clones_summary)

# Loop through the R lists and save them as a dictionary or json
#for (exps in C_inf){	
#	for (exp in exps){
#		str(exp)
#		for (i in 1:length(exp)){
#			print(exp[[i]])
#		}
#	}
#}

C_JSON=toJSON(C_inf,pretty=TRUE,auto_unbox=TRUE)
#C_JSON
write(C_JSON,args[2])

clonesJSON=toJSON(clones_summary,pretty=TRUE,auto_unbox=TRUE)
#clonesJSON
write(clonesJSON,args[3])
