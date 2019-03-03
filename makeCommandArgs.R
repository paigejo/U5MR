# make the command arguments file for runAllSPDE.R
makeSpdeCommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                            includeClustEffectVec=c(FALSE, TRUE), urbanEffectVec=c(FALSE, TRUE), 
                            testVec=c(FALSE, TRUE)) {
  spdeCommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        for(i4 in 1:length(includeClustEffectVec)) {
          includeClustEffect = includeClustEffectVec[i4]
          
          for(i5 in 1:length(urbanEffectVec)) {
            urbanEffect = urbanEffectVec[i5]
            
            for(i6 in 1:length(testVec)) {
              test = testVec[i6]
              
              spdeCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, 
                                          includeClustEffect=includeClustEffect, urbanEffect=urbanEffect, 
                                          test=test)
              i = i + 1
            }
          }
        }
      }
    }
  }
  
  save(spdeCommandArgs, file="spdeCommandArgs.RData")
}
# spdeCommandArgs = list(list(tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE), 
#                        list(tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE), 
#                        list(tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(tausq=0, test=TRUE,  includeUrbanRural=FALSE), 
#                        list(tausq=0, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(tausq=0, test=TRUE,  includeUrbanRural=TRUE), 
#                        list(tausq=0, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE), 
#                        list(tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE), 
#                        list(tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(tausq=0, test=FALSE,  includeUrbanRural=FALSE), 
#                        list(tausq=0, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(tausq=0, test=FALSE,  includeUrbanRural=TRUE), 
#                        list(tausq=0, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        
#                        list(margVar=0, tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE), 
#                        list(margVar=0, tausq=0.1^2, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE), 
#                        list(margVar=0, tausq=0.1^2, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0, test=TRUE,  includeUrbanRural=FALSE), 
#                        list(margVar=0, tausq=0, test=TRUE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0, test=TRUE,  includeUrbanRural=TRUE), 
#                        list(margVar=0, tausq=0, test=TRUE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE), 
#                        list(margVar=0, tausq=0.1^2, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE), 
#                        list(margVar=0, tausq=0.1^2, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0, test=FALSE,  includeUrbanRural=FALSE), 
#                        list(margVar=0, tausq=0, test=FALSE,  includeUrbanRural=FALSE, includeClustEffect=FALSE), 
#                        list(margVar=0, tausq=0, test=FALSE,  includeUrbanRural=TRUE), 
#                        list(margVar=0, tausq=0, test=FALSE,  includeUrbanRural=TRUE, includeClustEffect=FALSE))
# save(spdeCommandArgs, file="spdeCommandArgs.RData")

# tausq=0.1^2, big=TRUE, test
# make the command arguments file for runAllSPDE.R
makeDirectCommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                                 bigVec=c(FALSE, TRUE), testVec=c(FALSE, TRUE)) {
  directCommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        for(i4 in 1:length(bigVec)) {
          big = bigVec[i4]
          
          for(i5 in 1:length(testVec)) {
            test = testVec[i5]
            
            directCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, 
                                          big=big, test=test)
            i = i + 1
          }
        }
      }
    }
  }
  
  save(directCommandArgs, file="directCommandArgs.RData")
}