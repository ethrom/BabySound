#############################
### OUTLINE OF THIS SCRIPT
###
### For each audio property:
### - create distribution plots
### - test assumptions for ANOVA
### - calculate ANOVA
### - calculate follow-up test of potential significant results
###
### Elena Throm, 28/07/2020
#############################


##############
# PREPARATIONS
#############

update.packages(ask = FALSE)

library(psych)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(car)
library(lme4)
library(nlme)
library(ggplot2)
library(shiny)
library(BayesFactor)
library(coda)
library(Matrix)
library(DHARMa) #residuenanalyse glmer
library(multcomp) # post hoc tests
library(lattice)
library(readxl)
library(xlsx)
library(plyr)
library(tidyverse)
library(funModeling)
library(Hmisc)
library(tidyverse)
library(gdata)
library(dplyr)
library(ggcorrplot)
library(GGally)
library(GauPro)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra)
library(dendsort)
library(pheatmap)
library(corrplot)
library(FactoMineR)
library(stringr)
library(ggfortify)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(pracma)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(vioplot)
install.packages("psych")
library(psych)
library("RColorBrewer")

# set working directory
setwd("/Users/cbcd/PhD/BONDS/fMRI/20sounds")
getwd()

# read data
data = read_excel("/Users/cbcd/PhD/BONDS/fMRI/20sounds/SoundProperties20Sounds.xls")
View(data)
colnames(data) <- c("sound", "condition", "mean_harmonic_ratio",	"CV_hr",	"mean_pitch",	"median_pitch","CV_pitch"	,"mean_spectral_decrease",	"median_spectral_decrease",	"std_spectral_decrease"	,"mean_entropy"	,"median_entropy",	"std_entropy",	"mean_flux",	"median_flux","std_flux","mean_centroid","median_centroid",	"std_centroid"	,"mean_kurtosis",	"median_kurtosis"	,"std_kurtosis"	,"mean_spread",	"median_spread"	,"std_spread",	"mean_crest",	"median_crest",	"std_crest"	,"mean_spectralPeak",	"median_spectralPeak"	,"std_spectralPeak",	"mean_flatness"	,"median_flatness","std_flatness", "mean_slope",	"median_slope",	"std_slope",	"mean_skewness","median_skewness","std_skewness","mean_rolloffpoint","median_rolloffpoint","std_rolloffpoint")



#############
### FUNCTIONS
#############


# FUNCTION: vioplot

vio.plot = function(property){
  
  property=as.matrix(property)
  plotname = strcat(c(colnames(property)))
  
  x = as.matrix(data$sound)
  y = gsub("_normalised.wav", "", x)
  
  f = ggplot(data, aes(factor(condition), property, label = y))
  graph = f + geom_violin(trim=FALSE,  fill='#A4A4A4', color="darkred") +  
    scale_x_discrete(breaks=c("1","2","3"),  labels=c("social", "partly social", "nonsocial")) +
    theme(legend.position = "none") +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    labs(title=plotname) +
    geom_boxplot(width = 0.1) +
    geom_jitter(height = 0, width = 0.05) +
    theme(plot.title = element_text(color = "black", size = 20, face = "bold")) +
    geom_text_repel()
  
  # # alternative plot design
  #   plotname = strcat(c(colnames(property)))
  #   vioplot(property ~ condition, data=data, col = "grey",
  #           xlab="condition", ylab= plotname,
  #           names=c("social","partly social","nonsocial"))
  #   text(property ~ condition, labels=y,data=data, cex=.6, font=0.05)
  
  name = colnames(property)
  pdfname = strcat(c("/Users/cbcd/PhD/BONDS/fMRI/20sounds/output/", name, ".png"))
  png(file= pdfname) 
  par(mar=c(5,10,4,1)+.1)
  print(graph)
  dev.off()
  
  return(graph)
}



# FUNCTION: test assumptions for ANOVA

test.assumptions = function(property){
  
  p = unlist(property)
  p = as.numeric(p)
  
  # build the linear model
  model  <- lm(p~ condition, data = data)
  
  # check residuals
  ggqqplot(residuals(model))
  shapiro_test(residuals(model))
  
  # check for normal distribution
  normality_shapiro = shapiro.test(p)
  x = normality_shapiro[["p.value"]]
  ggdensity(p, fill = "lightgray")
  ggqqplot(p)
  
  # check for homogeneity of variances
  c = as.factor(data$condition)
  homogen_levene = leveneTest(p ~ c, data=data)
  y = homogen_levene[1,3]
  
  # check for homogeneity of variances (alternative if normal distribution of residuals not given)
  homogen_fligner = fligner.test(p ~ condition, data)
  z = homogen_fligner[["p.value"]]

  name = colnames(property)
  res = cbind(c(name,x,y,z))
  return(res)
}



# FUNCTION: ANOVA for equally and unequally distributed variances (depending on Levene-test outcome)
calc.anova = function(property){
  
  p = unlist(property)
  p = as.numeric(p)
  
  if (res[3,1] < .05){
    
    res.aov.hetero = oneway.test(p ~ condition, data=data, var.equal=FALSE)
    s = res.aov.hetero[["p.value"]]
  }
    
  else {
    res.aov.homo = aov(p ~ condition, data = data)
    aov = summary(res.aov.homo)
    s = aov[[1]] [1,5]
  }
 
  return(s)
} 



##########
### LOOPS
##########

# LOOP for plotting: for each audio property, create and save violin plot

i = 3
for(i in 3:43) {
  property = data[,i]
  # property = data[,"std_flatness"]
  # property=as.matrix(property)
  vio.plot(property)
  i = i +1
}



# LOOP for ANOVA: for each audio property: calculate ANOVA, write down p values; 
# if significant, calculate post-hoc test

i = 3
for(i in 3:43) {
  property = data[,i] 
  property=as.matrix(property)
  test.assumptions(property)
  P = calc.anova(property)
  # if(P < 0.05){
        write(colnames(property), file = "properties_across_conditions.csv", append = TRUE)
        write(P, file = "properties_across_conditions.csv", append = TRUE)
   # }
    
  if(P < 0.05){
        tt = pairwise.t.test(property, data$condition,
                            p.adjust.method = "bonf")
        t = tt[["p.value"]]
        write(t, file = "properties_across_conditions.csv", append = TRUE)
  }
    
  i = i +1
}



