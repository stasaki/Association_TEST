
library(tidyverse)
library(reshape2)
source("assoc_test.R")

# make synthetic data
# expdata need to have three variables, expression matrix,sample annotation, and gene annotations)
# rownames and colnames should match for all variables.

expdata=list()
expdata$data=rnorm(10*1000)%>%matrix(.,nrow = 1000,ncol = 10)
colnames(expdata$data)=paste0("s",1:10)
rownames(expdata$data)=paste0("g",1:1000)

expdata$sample.ann=data.frame(SAMPLE.ID=paste0("s",1:10),
                              Group=c(rep(1,5),rep(0,5))%>%as.factor(),
                              Cov1=rnorm(10))
rownames(expdata$sample.ann)=expdata$sample.ann$SAMPLE.ID
expdata$gene.ann = data.frame(Name=paste0("g",1:1000))
rownames(expdata$gene.ann)=expdata$gene.ann$Name

# Set parameters for association test 
PARAM=list()
PARAM$PRIMARY_COVS="Group" # Main phenotype you want to investigate
PARAM$ADJUST_COVS="Cov1" # Covariates you want to adjust
PARAM$TEST.METHOD="LIMMA" # Method you want to use LIMMA or LM (just normal lm)
PARAM$spline_func=NULL # set if you want to investigate time-dependent effect
PARAM$INTERACT_COV=NULL # set if you want to investigate interacions
PARAM$BLOCK_COV=NULL # set if you want to do mixed-effect modeling
PARAM$SAVE_ADJUSTED=FALSE # set TRUE if you want to save Covariate-adjusted data
PARAM$OUTLOC="./out_limma/" # set output directory

run_group_test(PARAM,expdata)

PARAM$TEST.METHOD="LM"
PARAM$OUTLOC="./out_lm/" # set output directory
run_group_test(PARAM,expdata)
