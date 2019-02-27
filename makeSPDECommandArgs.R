# make the command arguments file for runAllSPDE.R
spdeCommandArgs = list(list(tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE), 
                       list(tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
                       list(tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE), 
                       list(tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
                       list(tausq=0, test=TRUE,  includeUrbanRural=FALSE), 
                       list(tausq=0, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
                       list(tausq=0, test=TRUE,  includeUrbanRural=TRUE), 
                       list(tausq=0, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
                       list(tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE), 
                       list(tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
                       list(tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE), 
                       list(tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
                       list(tausq=0, test=FALSE,  includeUrbanRural=FALSE), 
                       list(tausq=0, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
                       list(tausq=0, test=FALSE,  includeUrbanRural=TRUE), 
                       list(tausq=0, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE))
save(spdeCommandArgs, file="spdeCommandArgs.RData")
