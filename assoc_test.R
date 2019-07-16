
# run script
run_group_test<-function(PARAM,expdata){
  library(limma)
  # drop NA
  indx = !is.na(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]])
  expdata$data = expdata$data[,indx]
  expdata$sample.ann = expdata$sample.ann[indx,,drop=FALSE]
  
  # primary covariate has a variation?
  if(length(unique(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]]))<=1){
    print("phenotype has no variation")
    return()
  }
  
  if(is.factor(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]])&&
     min(table(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]]))< 2){
    
    if(is.null(PARAM$ignore_sample_size) || !PARAM$ignore_sample_size){
      print("minimum phenotype group should have 3 samples at least")
      unlink(list.files(PARAM$OUTLOC,full.names = T))
      return()
    }else{
      print("sample size less than 2")
    }
    
  }
  if((expdata$sample.ann%>%nrow) < (2+length(PARAM[["ADJUST_COVS"]]))){
    print("sample size too small")
    return()
  }

  # make interaction covariates and adjusted covariate unique
  PARAM[["ADJUST_COVS"]] = setdiff(PARAM[["ADJUST_COVS"]],PARAM[["INTERACT_COV"]])
  
  if(PARAM$TEST.METHOD=="LIMMA"){
    limma_group_block(expdata,
                      PRIMARY_COV=PARAM[["PRIMARY_COVS"]],
                      ADJUST_COVS=PARAM[["ADJUST_COVS"]],
                      BLOCK_COV=PARAM[["BLOCK_COV"]],
                      INTERACT_COV=PARAM[["INTERACT_COV"]],
                      spline_func=PARAM[["spline_func"]],
                      OUTLOC=PARAM[["OUTLOC"]],
                      eBayes=TRUE,
                      SAVE_ADJUSTED=PARAM[["SAVE_ADJUSTED"]],
                      PARAM=PARAM)
  }else{
    if(!is.null(PARAM[["INTERACT_COV"]])){
      return()
    }
    limma_group_block(expdata,
                      PRIMARY_COV=PARAM[["PRIMARY_COVS"]],
                      ADJUST_COVS=PARAM[["ADJUST_COVS"]],
                      BLOCK_COV=PARAM[["BLOCK_COV"]],
                      INTERACT_COV=PARAM[["INTERACT_COV"]],
                      spline_func=PARAM[["spline_func"]],
                      OUTLOC=PARAM[["OUTLOC"]],
                      eBayes=FALSE,
                      SAVE_ADJUSTED=PARAM[["SAVE_ADJUSTED"]],
                      PARAM=PARAM)
  }
}


limma_group_block <- function(expdata,
                              PRIMARY_COV,ADJUST_COVS,BLOCK_COV,INTERACT_COV=NULL,
                              spline_func=NULL,OUTLOC,
                              B.num=1000,eBayes=TRUE,
                              SAVE_ADJUSTED=FALSE,
                              PARAM,
                              DiffVar=FALSE){
  library(limma)
  library(splines)
  # make output directory
  dir.create(paste0(OUTLOC),recursive = TRUE,showWarnings = FALSE)
  
  # Variable should be there in data
  NEED_FOR_PLOT=c(PRIMARY_COV,ADJUST_COVS,BLOCK_COV,INTERACT_COV)
  
  # check data is ready
  expdata = sanity_check(expdata,PRIMARY_COV,NEED_FOR_PLOT)
  
  # get component
  pheno = expdata$sample.ann
  temp_mat = expdata$data
  gene.ann = expdata$gene.ann
  
  # spline_function
  if(is.null(spline_func)){
    spline_func = PRIMARY_COV
  }else{
    spline_func = paste0(spline_func,"(",PRIMARY_COV,")")
  }
  
  if(!is.null(INTERACT_COV)){
    design_obj=make_design_obj_interaction(pheno,spline_func,PRIMARY_COV,INTERACT_COV,ADJUST_COVS)
    GROUP_COV=INTERACT_COV
  }else{
    design_obj=make_design_obj(pheno,spline_func,PRIMARY_COV,ADJUST_COVS)
    GROUP_COV=PRIMARY_COV
  }
  
  
  if (!is.null(PARAM$VOOM) && PARAM$VOOM){
    print("voom normalization")
    # Normalize using voom, WITH adjusting for known and HV covariates:
    temp_mat = calcNormFactors(DGEList(counts = temp_mat),
                               method=PARAM$calcNormFactors.method)%>%
      voom(., design = design_obj$design,
           normalize.method = PARAM$normalize.method,
           plot=F)
  }
  
  
  # run limma block
  if(DiffVar==TRUE){
    temp_mat_adjusted=diffvar_group_sub(design_obj$indx_coef,design_obj$design,
                                        design_obj$contrast.matrix,design_obj$contrast.group,
                                        NEED_FOR_PLOT,
                                        PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                                        pheno,temp_mat,gene.ann,
                                        OUTLO)    
  }else if(eBayes==TRUE){
    temp_mat_adjusted=limma_group_sub(design_obj$indx_coef,design_obj$design,
                                      design_obj$contrast.matrix,design_obj$contrast.group,
                                      NEED_FOR_PLOT,
                                      PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                                      pheno,temp_mat,gene.ann,
                                      OUTLOC)
  }else{
    temp_mat_adjusted=lm_group_sub(design_obj$indx_coef,design_obj$design,
                                   design_obj$contrast.matrix,design_obj$contrast.group,
                                   spline_func,
                                   NEED_FOR_PLOT,
                                   PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                                   pheno,temp_mat,gene.ann,
                                   OUTLOC)
  }
  
  if (SAVE_ADJUSTED==TRUE){
    saveRDS(temp_mat_adjusted,file=paste0(OUTLOC,"/data_adjusted.rds"))
  }
  return(TRUE)
}


limma_group_sub<-function(indx_coef,design,contrast.matrix,contrast.group,
                          NEED_FOR_PLOT,
                          PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                          pheno,temp_mat,gene.ann,
                          OUTLOC,B.num){
  

  
  # fit with linear model
  if(!is.null(BLOCK_COV)){
    print("mixed effect modeling")
    corfit <- duplicateCorrelation(temp_mat,design,block=pheno[[BLOCK_COV]])
    fit <- lmFit(temp_mat,design,block=pheno[[BLOCK_COV]],correlation=corfit$consensus)
  }else{
    fit <- lmFit(temp_mat,design)
  }
  fit <- eBayes(fit)
  
  
  # get residual
  # back primary variable information
  if(colnames(design)[1]%in%c("(Intercept)","X.Intercept.")){
    back_coef = c(1,indx_coef)
  }else{
    back_coef = indx_coef
  }
  temp_mat_adjusted = residuals(fit,temp_mat)
  temp_mat_adjusted = temp_mat_adjusted+t(design[,back_coef]%*%t(fit$coefficients[,back_coef]))
  #fitting with all
  fitted.value = t(fit$design%*%t(fit$coefficients))
  # fitting with adjusted covs
  adjust.value =  t(fit$design[,-back_coef]%*%t(fit$coefficients[,-back_coef]))
  # adjusting by adjust covs
  fitted.value_adjusted = fitted.value-adjust.value
  
  if (is.null(contrast.matrix)){
    # numeric test
    saveRDS(fit,file=paste(OUTLOC,"/fit_",PRIMARY_COV,".rds",sep=""))
    write.table(design,file=paste(OUTLOC,"/design_",PRIMARY_COV,".txt",sep=""),
                quote=FALSE,sep="\t",row.names=FALSE)
    
    DS=topTable(fit, genelist= rownames(temp_mat),
                number=fit$coefficients%>%nrow(),
                p.value=1,lfc=log2(1),
                coef=indx_coef, adjust="BH")
    suffix = NULL
    p_group = unique(pheno[[GROUP_COV]])
    limma_group_plot(DS,
                     temp_mat,temp_mat_adjusted,
                     fitted.value,fitted.value_adjusted,
                     gene.ann,pheno,
                     PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                     suffix,p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
    
    return(list(data = temp_mat,
                data.adjusted = temp_mat_adjusted,
                fitted = fitted.value,
                fitted.adjusted = fitted.value_adjusted))
  }
  
  fit.const <- contrasts.fit(fit,contrast.matrix)
  fit.const <- eBayes(fit.const) # this is basically ANOVA
  
  # save result
  saveRDS(contrast.matrix,file=paste(OUTLOC,"/contrast.matrix_",PRIMARY_COV,".rds",sep=""))
  saveRDS(fit.const,file=paste(OUTLOC,"/fit_",PRIMARY_COV,".rds",sep=""))
  write.table(design,file=paste(OUTLOC,"/design_",PRIMARY_COV,".txt",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)
  
  
  if (ncol(contrast.matrix)==1){
    print("just two group comparision")
    suffix = contrast.group
    DS=topTable(fit.const, genelist= rownames(temp_mat),
                number=fit.const$coefficients%>%nrow(),
                p.value=1,lfc=log2(1),
                coef=suffix, adjust="BH")
    
    p_group = unique(unlist(strsplit(suffix,"-")))
    limma_group_plot(DS,
                     temp_mat,temp_mat_adjusted,
                     fitted.value,fitted.value_adjusted,
                     gene.ann,pheno,
                     PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                     suffix,p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
  }else{
    print("do anova")
    suffix = "ANOVA"
    DS=topTable(fit.const, genelist= rownames(temp_mat),
                number=fit.const$coefficients%>%nrow(),
                p.value=1,lfc=log2(1),
                coef=colnames(contrast.matrix), adjust="BH")
    p_group = unique(unlist(strsplit(contrast.group,"-")))
    limma_group_plot(DS,
                     temp_mat,temp_mat_adjusted,
                     fitted.value,fitted.value_adjusted,
                     gene.ann,pheno,
                     PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                     suffix,p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
    
    print("do each of comparison")
    
    for (suffix in unique(contrast.group)){
      
      DS=topTable(fit.const, genelist= rownames(temp_mat),
                  number=fit.const$coefficients%>%nrow(),
                  p.value=1,lfc=log2(1),
                  coef=colnames(contrast.matrix)[contrast.group==suffix],
                  adjust="BH")
      p_group = unique(unlist(strsplit(suffix,"-")))
      limma_group_plot(DS,
                       temp_mat,temp_mat_adjusted,
                       fitted.value,fitted.value_adjusted,
                       gene.ann,pheno,
                       PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                       suffix,p_group,
                       NEED_FOR_PLOT,
                       OUTLOC)
    }
  }
  
  return(list(data = temp_mat,
              data.adjusted = temp_mat_adjusted,
              fitted = fitted.value,
              fitted.adjusted = fitted.value_adjusted))
}


limma_group_plot<-function(DS,
                           temp_mat,temp_mat_adjusted,
                           fitted.value,fitted.value_adjusted,
                           gene.ann,pheno,
                           PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                           suffix,p_group,
                           NEED_FOR_PLOT,
                           OUTLOC){
  
  #output
  indx = match(DS[,1],rownames(gene.ann))
  DS = cbind(DS,gene.ann[indx,,drop=FALSE])
  colnames(DS)[1]="Name"
  tryCatch({
    #    write.table(DS,file=paste(OUTLOC,"/stat_",PRIMARY_COV,"_",suffix,".txt",sep=""),
    #                quote=FALSE,sep="\t",row.names=FALSE)
  },error = function(e) {
  })
  saveRDS(DS,file=paste(OUTLOC,"/stat_",PRIMARY_COV,"_",suffix,".rds",sep=""))
  
  tryCatch({
    # make pvalue histgram
    pdf(paste0(OUTLOC,"/pvalue_dist_",PRIMARY_COV,"_",suffix,".pdf"))
    hist(DS$P.Value)
    dev.off()
  }, warning = function(w) {
  }, error = function(e) {
  }, finally = {
  })
  
  # plot top 100
  indx = DS$P.Value<0.05
  indx[is.na(indx)]=FALSE
  SIG_VARS = DS[1:min(6,nrow(DS)),1]
  #SIG_VARS = SIG_VARS[1:6]
  
  
  # plot raw data and adjusted data
  # also spline fit curve 
  
  # take out data for plot
  # unadjusted
  if(class(temp_mat)=="EList"){
    temp_mat=temp_mat$E
  }
  ds = cbind(t(temp_mat[SIG_VARS,,drop=FALSE]),
             pheno[,NEED_FOR_PLOT,drop=FALSE])
  ds = melt(ds,id.vars =NEED_FOR_PLOT)
  
  ds_fitted = cbind(t(fitted.value[SIG_VARS,,drop=FALSE]),
                    pheno[,NEED_FOR_PLOT,drop=FALSE])
  ds_fitted = melt(ds_fitted,id.vars = NEED_FOR_PLOT,
                   value.name = "fit")
  # combine
  ds =merge(ds,ds_fitted,
            by = c(NEED_FOR_PLOT,"variable"))
  ds$Type="raw"
  
  # adjusted
  ds_adjusted = cbind(t(temp_mat_adjusted[SIG_VARS,,drop=FALSE]),
                      pheno[,NEED_FOR_PLOT,drop=FALSE])
  ds_adjusted = melt(ds_adjusted,id.vars = NEED_FOR_PLOT)
  
  ds_fitted = cbind(t(fitted.value_adjusted[SIG_VARS,,drop=FALSE]),
                    pheno[,NEED_FOR_PLOT,drop=FALSE])
  ds_fitted = melt(ds_fitted,id.vars =NEED_FOR_PLOT,
                   value.name = "fit")
  # combine
  ds_adjusted =merge(ds_adjusted,ds_fitted,
                     by = c(NEED_FOR_PLOT,"variable"))
  ds_adjusted$Type = "adjusted"
  
  ds_adjusted$variable = factor(as.character(ds_adjusted$variable),
                                levels=SIG_VARS)
  ds$variable = factor(as.character(ds$variable),
                       levels=SIG_VARS)
  # filter group
  ds = ds[ds[[GROUP_COV]]%in%p_group,]
  ds_adjusted = ds_adjusted[ds_adjusted[[GROUP_COV]]%in%p_group,]
  
  library(ggthemes)
  if(is.numeric(ds[[PRIMARY_COV]])|is.numeric(ds[[GROUP_COV]])){
    g_obj=list()
    g_obj[[1]]=ggplot(ds_adjusted,
                      aes_string(PRIMARY_COV,"value",group=BLOCK_COV,
                                 color=GROUP_COV))+
      facet_wrap(~variable,scales = "free_y",ncol = 3)+
      geom_point()+
      geom_line(aes(y=fit),alpha=.3,size=2)
    g_obj[[2]]=ggplot(ds,
                      aes_string(PRIMARY_COV,"value",group=BLOCK_COV,
                                 color=GROUP_COV))+
      facet_wrap(~variable,scales = "free_y",ncol = 3)+
      geom_point()+
      geom_line(aes(y=fit),alpha=.3,size=2)
    
    plot.nrow = ceiling(length(unique(ds$variable))/3)
    plot.ncol = min(3,length(unique(ds$variable)))
    
    pdf(paste0(OUTLOC,"/plot_",PRIMARY_COV,"_",suffix,".pdf"),
        width = plot.ncol*3+2,height = plot.nrow*2+0.5)
    plot(g_obj[[1]])
    plot(g_obj[[2]])
    dev.off()
    return()
  }
  
  # for classes, plot two version
  g_obj=list()
  g_obj[[1]]=ggplot(ds_adjusted,
                    aes_string(PRIMARY_COV,"value",
                               fill=GROUP_COV))+
    facet_wrap(~variable,scales = "free_y",ncol = 3)+
    geom_boxplot()+theme_bw()+scale_fill_stata()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "grey10"))
  g_obj[[2]]=ggplot(ds,
                    aes_string(PRIMARY_COV,"value",
                               fill=GROUP_COV))+
    facet_wrap(~variable,scales = "free_y",ncol = 3)+
    geom_boxplot()+theme_bw()+scale_fill_stata()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "grey10"))
  
  plot.nrow = ceiling(length(unique(ds$variable))/3)
  plot.ncol = min(3,length(unique(ds$variable)))
  n_groups = length(p_group)
  
  pdf(paste0(OUTLOC,"/plot_",PRIMARY_COV,"_",suffix,".pdf"),
      width = plot.ncol*(1+n_groups*0.5)+2,height = plot.nrow*2+0.5)
  plot(g_obj[[1]])
  plot(g_obj[[2]])
  dev.off()
}


lm_group_sub <- function(indx_coef,design,contrast.matrix,contrast.group,
                         spline_func,
                         NEED_FOR_PLOT,
                         PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                         pheno,temp_mat,gene.ann,
                         OUTLOC,B.num){
  
  library(lme4)
  orig_name = rownames(temp_mat)
  rownames(temp_mat) = paste0("X",1:nrow(temp_mat)) # get fomula can handle
  #fit_ds = cbind(pheno[,c(PRIMARY_COV,ADJUST_COVS,BLOCK_COV),drop=FALSE],t(temp_mat))
  pheno_use = pheno[,c(PRIMARY_COV,ADJUST_COVS,BLOCK_COV),drop=FALSE]
  mat_use = t(temp_mat)
  rownames(temp_mat)= orig_name # back to original ones
  
  DS.ordinary=lapply(1:nrow(temp_mat), function(j){
    #print(j)
    y = paste0("X",j)
    cbind(mat_use[,j,drop=F],pheno_use)%>%
      #dplyr::select_(fit_ds,.dots=c(y,PRIMARY_COV,ADJUST_COVS,BLOCK_COV))%>%
      mutate(r.n. = rownames(.))%>%
      filter(rowSums(is.na(.))==0) -> tmp.data
    
    indx = apply(tmp.data,2,function(x)length(unique(x)))==1
    ADJUST_COVS_rm = colnames(tmp.data)[indx]
    
    
    if(length(unique(tmp.data[!is.na(tmp.data[,1]),2]))==1|
       length(unique(tmp.data[!is.na(tmp.data[,1]),1]))<3){
      return(NULL)
    }
    
    if(is.null(contrast.matrix)){
      # numeric test
      suffix = "quant"
      
      if(is.null(BLOCK_COV)){
        lm.obj = lm(as.formula(paste0(y,
                                      paste(c("~1",spline_func,setdiff(ADJUST_COVS,ADJUST_COVS_rm)),
                                            collapse = "+"))),
                    data = tmp.data, x=TRUE)
        indx_coef = which(regexpr(PRIMARY_COV,names(lm.obj$coefficients),fixed = T)>0)
      }else{
        lm.obj = lmer(as.formula(paste0(y,
                                        paste(c("~1",spline_func,setdiff(ADJUST_COVS,ADJUST_COVS_rm),
                                                paste0("(1|",BLOCK_COV,")")),
                                              collapse = "+"))),
                      data = tmp.data, REML=FALSE)
        indx_coef = which(regexpr(PRIMARY_COV,names(coef(lm.obj)[[1]]),fixed = T)>0)
      }
      DS.ordinary=list()
      DS.ordinary[[suffix]]=myTopTable_coef(lm.obj,indx_coef)
      DS.ordinary[[suffix]]$ID = orig_name[j]
      DS.ordinary[[suffix]]$AveExpr=mean(tmp.data[[y]])
      
      # get residual
      back_coef=c(1,indx_coef)
      fitted = get_fitted_values(lm.obj,back_coef)
      colnames(fitted$raw.value_adjusted)=tmp.data$r.n.
      colnames(fitted$fitted.value)=tmp.data$r.n.
      colnames(fitted$fitted.value_adjusted)=tmp.data$r.n.
      
      
      return(list(DS.ordinary=DS.ordinary,fitted=fitted))
    }
    
    if(is.null(BLOCK_COV)){
      lm.obj = lm(as.formula(paste0(y,
                                    paste(c("~0",spline_func,setdiff(ADJUST_COVS,ADJUST_COVS_rm)),
                                          collapse = "+"))),
                  data = tmp.data, x=TRUE)
      indx_coef = which(regexpr(PRIMARY_COV,names(lm.obj$coefficients),fixed = T)>0)   
    }else{
      lm.obj = lmer(as.formula(paste0(y,
                                      paste(c("~0",spline_func,setdiff(ADJUST_COVS,ADJUST_COVS_rm),
                                              paste0("(1|",BLOCK_COV,")")),
                                            collapse = "+"))),
                    data = tmp.data, REML=FALSE)
      indx_coef = which(regexpr(PRIMARY_COV,names(fixef(lm.obj)),fixed = T)>0)
    }
    
    
    fitted = get_fitted_values(lm.obj,indx_coef)
    colnames(fitted$raw.value_adjusted)=tmp.data$r.n.
    colnames(fitted$fitted.value)=tmp.data$r.n.
    colnames(fitted$fitted.value_adjusted)=tmp.data$r.n.
    
    DS.ordinary=list()
    if (ncol(contrast.matrix)==1){
      suffix = contrast.group
      # single contrast
      DS.ordinary[[suffix]]=myTopTable(lm.obj,contrast.matrix[indx_coef,,drop=F],PRIMARY_COV)
      DS.ordinary[[suffix]]$ID = orig_name[j]
      DS.ordinary[[suffix]]$AveExpr=mean(tmp.data[[y]])
    }else{
      #ANOVA
      suffix="ANOVA"
      DS.ordinary[[suffix]]=myTopTable(lm.obj,contrast.matrix[indx_coef,,drop=F],PRIMARY_COV)
      DS.ordinary[[suffix]]$ID = orig_name[j]
      DS.ordinary[[suffix]]$AveExpr=mean(tmp.data[[y]])
      for (suffix in unique(contrast.group)){
        # each contrast
        DS.ordinary[[suffix]]=myTopTable(lm.obj,
                                         contrast.matrix[indx_coef,colnames(contrast.matrix)[contrast.group==suffix],drop=F],
                                         PRIMARY_COV)
        DS.ordinary[[suffix]]$ID = orig_name[j]
        DS.ordinary[[suffix]]$AveExpr=mean(tmp.data[[y]])
      }
    }
    return(list(DS.ordinary=DS.ordinary,fitted=fitted))
  })
  
  names(DS.ordinary)=orig_name
  indx=!sapply(DS.ordinary, is.null) #remove NA
  DS.ordinary = DS.ordinary[indx]
  temp_mat = temp_mat[indx,]
  
  # recover nessesary information for plotting
  temp_mat_adjusted = lapply(1:length(DS.ordinary), function(x){
    data.frame(SAMPLE.ID = colnames(DS.ordinary[[x]]$fitted$raw.value_adjusted),
               ID=names(DS.ordinary)[x],
               value=t(DS.ordinary[[x]]$fitted$raw.value_adjusted))
  })
  temp_mat_adjusted = do.call(rbind,temp_mat_adjusted)
  group_by(temp_mat_adjusted,SAMPLE.ID,ID)%>%
    spread(SAMPLE.ID,value) -> temp_mat_adjusted
  r.n.=temp_mat_adjusted$ID
  temp_mat_adjusted=as.matrix(temp_mat_adjusted[,-1])
  rownames(temp_mat_adjusted)=r.n.
  indx=match(rownames(temp_mat),rownames(temp_mat_adjusted))
  temp_mat_adjusted=temp_mat_adjusted[indx,]
  indx=match(colnames(temp_mat),colnames(temp_mat_adjusted))
  temp_mat_adjusted=temp_mat_adjusted[,indx]
  
  fitted.value = lapply(1:length(DS.ordinary), function(x){
    data.frame(SAMPLE.ID = colnames(DS.ordinary[[x]]$fitted$fitted.value),
               ID=names(DS.ordinary)[x],
               value=t(DS.ordinary[[x]]$fitted$fitted.value))
  })
  fitted.value = do.call(rbind,fitted.value)
  group_by(fitted.value,SAMPLE.ID,ID)%>%
    spread(SAMPLE.ID,value) -> fitted.value
  r.n.=fitted.value$ID
  fitted.value=as.matrix(fitted.value[,-1])
  rownames(fitted.value)=r.n.
  indx=match(rownames(temp_mat),rownames(fitted.value))
  fitted.value=fitted.value[indx,]
  indx=match(colnames(temp_mat),colnames(fitted.value))
  fitted.value=fitted.value[,indx]
  
  fitted.value_adjusted = lapply(1:length(DS.ordinary), function(x){
    data.frame(SAMPLE.ID = colnames(DS.ordinary[[x]]$fitted$fitted.value_adjusted),
               ID=names(DS.ordinary)[x],
               value=t(DS.ordinary[[x]]$fitted$fitted.value_adjusted))
  })
  fitted.value_adjusted = do.call(rbind,fitted.value_adjusted)
  group_by(fitted.value_adjusted,SAMPLE.ID,ID)%>%
    spread(SAMPLE.ID,value) -> fitted.value_adjusted
  r.n.=fitted.value_adjusted$ID
  fitted.value_adjusted=as.matrix(fitted.value_adjusted[,-1])
  rownames(fitted.value_adjusted)=r.n.
  indx=match(rownames(temp_mat),rownames(fitted.value_adjusted))
  fitted.value_adjusted=fitted.value_adjusted[indx,]
  indx=match(colnames(temp_mat),colnames(fitted.value_adjusted))
  fitted.value_adjusted=fitted.value_adjusted[,indx]
  
  for (suffix in names(DS.ordinary[[1]]$DS.ordinary)){
    DS = lapply(DS.ordinary, function(x){
      return(x$DS.ordinary[[suffix]])
    })
    do.call(rbind,DS)%>%arrange(P.Value)->DS
    
    if(suffix=="ANOVA"){
      p_group = unique(unlist(strsplit(contrast.group,"-")))
    }else if(suffix=="quant"){
      suffix = NULL
      p_group = unique(pheno[[GROUP_COV]])
    }else{
      p_group = unique(unlist(strsplit(suffix,"-"))) 
    }
    limma_group_plot(DS,
                     temp_mat,temp_mat_adjusted,
                     fitted.value,fitted.value_adjusted,
                     gene.ann,pheno,
                     PRIMARY_COV,ADJUST_COVS,BLOCK_COV,GROUP_COV,
                     suffix,p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
  }
  return(list(data = temp_mat,
              data.adjusted = temp_mat_adjusted,
              fitted = fitted.value,
              fitted.adjusted = fitted.value_adjusted))
}




my_mcp<-function (linfct, interaction_average = FALSE, covariate_average = FALSE) 
{
  library(multcomp)
  linfct <- lapply(linfct, function(x) {
    if (is.numeric(x) && !is.matrix(x)) {
      return(matrix(x, nrow = 1))
    }
    else {
      return(x)
    }
  })
  if (is.null(names(linfct))) 
    stop(sQuote("linfct"), " doesn't have a ", sQuote("names"), 
         " attribute")
  classes <- sapply(linfct, function(x) inherits(x, "matrix") || 
                      inherits(x, "character"))
  if (length(linfct) == 1 && linfct[[1]] == "Means") {
    class(linfct) <- "means"
    return(linfct)
  }
  attr(linfct, "interaction_average") <- interaction_average
  attr(linfct, "covariate_average") <- covariate_average
  if (all(classes)) {
    class(linfct) <- "mcp"
    return(linfct)
  }
  stop("Arguments don't consist of either matrices or characters")
}


myTopTable=function(lm.obj,contrast.matrix,PRIMARY_COV){
  linfct =  list(t(contrast.matrix))
  names(linfct)[1]=PRIMARY_COV
  library(multcomp)
  t <- glht(lm.obj, linfct = my_mcp(linfct)) 
  
  if (ncol(contrast.matrix)==1){
    tstat=summary(t,univariate())$test$tstat[1]
    df=summary(t,univariate())$df
    if (df > 0) pvals <- 2 * pt(abs(tstat), df,lower.tail = FALSE) else pvals <- 2 * pnorm(abs(tstat),lower.tail = FALSE)
    return(data.frame(ID=NA,  
                      logFC=summary(t,univariate())$test$coefficients[1],
                      AveExpr=NA,
                      t = summary(t,univariate())$test$tstat[1],
                      P.Value=pvals,
                      row.names = NULL))
  }else{
    # Ftest
    return(data.frame(ID=NA,  
                      t(summary(t,test=Ftest())$test$coefficients),
                      AveExpr=NA,
                      F = summary(t,test=Ftest())$test$fstat[1,1],
                      P.Value= summary(t,test=Ftest())$test$pvalue[1,1],
                      row.names = NULL))
  }
}


myTopTable_coef=function(lm.obj,indx_coef){
  library(multcomp)
  if(class(lm.obj)=="lmerMod"){
    # t <- glht(lm.obj, linfct = paste0(names(coef(lm.obj)[[1]])[indx_coef],"=0"))  
    K=diag(length(coef(lm.obj)[[1]]))[indx_coef,,drop=F]
    rownames(K)=names(coef(lm.obj)[[1]])[indx_coef]
    t <- glht(lm.obj, linfct = K)
  }else{
    #t <- glht(lm.obj, linfct = paste0(names(lm.obj$coefficients)[indx_coef],"=0"))  
    K=diag(length(coef(lm.obj)))[indx_coef,,drop=F]
    rownames(K)=names(coef(lm.obj))[indx_coef]
    t <- glht(lm.obj, linfct = K)
  }
  
  if (length(indx_coef)==1){
    tstat=summary(t,univariate())$test$tstat[1]
    df=summary(t,univariate())$df
    if (df > 0) pvals <- 2 * pt(abs(tstat), df,lower.tail = FALSE) else pvals <- 2 * pnorm(abs(tstat),lower.tail = FALSE)
    return(data.frame(ID=NA,  
                      logFC=summary(t,univariate())$test$coefficients[1],
                      AveExpr=NA,
                      t = summary(t,univariate())$test$tstat[1],
                      P.Value=pvals,
                      row.names = NULL))
  }else{
    # Ftest
    return(data.frame(ID=NA,  
                      t(summary(t,test=Ftest())$test$coefficients),
                      AveExpr=NA,
                      F = summary(t,test=Ftest())$test$fstat[1,1],
                      P.Value= summary(t,test=Ftest())$test$pvalue[1,1],
                      row.names = NULL))
  }
}

get_fitted_values = function(lm.obj,back_coef){
  
  #fitting with all
  fitted.value = t(matrix(fitted(lm.obj)))
  
  if(class(coef(lm.obj))=="coef.mer"){
    #fitted.value_adjusted = rowSums(model.matrix(lm.obj)[,colnames(coef(lm.obj)[[1]])[back_coef]]*
    #                                  coef(lm.obj)[[1]][,back_coef])
    fitted.value_adjusted = t(model.matrix(lm.obj)[,names(fixef(lm.obj))[back_coef]]%*%
                                fixef(lm.obj)[back_coef])
  }else{
    fitted.value_adjusted = t(model.matrix(lm.obj)[,back_coef]%*%coef(lm.obj)[back_coef])
  }
  
  temp_mat_adjusted = fitted.value_adjusted+residuals(lm.obj)
  
  return(list(raw.value_adjusted=temp_mat_adjusted,
              fitted.value=fitted.value,
              fitted.value_adjusted=fitted.value_adjusted))
}




make_interaction_cnst = function(pheno,design,indx_coef,INTERACT_COV){
  
  # make contrast matrix
  ###
  cnst = colnames(design)[indx_coef]
  cnst_group = rep(NA,length(cnst))
  
  i_levels = levels(pheno[[INTERACT_COV]])
  cnst_group = do.call(rbind,strsplit(cnst,":"))[,2]
  cnst_group = paste0(cnst_group,"-",paste0(INTERACT_COV,i_levels[1]))
  
  cnstA=make.names(cnst)
  cnst_groupA=cnst_group
  
  ###
  itable = do.call(rbind,strsplit(colnames(design)[indx_coef],":"))
  
  u_plevel=unique(itable[,1])
  p_position=which(itable[,1]%in%u_plevel[1])
  cnst=apply(t(combn(p_position,2)),1,function(x){
    paste0(make.names(colnames(design)[indx_coef[x]]),collapse = "-")
  })
  cnst=sapply(1:length(u_plevel), function(x){
    gsub(u_plevel[1],u_plevel[x],cnst,fixed = TRUE)
  })
  
  cnst_group=cnst
  for (x in make.names(unique(itable[,1]))){
    cnst_group=gsub(x,"",cnst_group,fixed = TRUE)
  }
  cnst_group=gsub("^\\.","",cnst_group)
  cnst_group=gsub("-\\.","-",cnst_group)
  cnst = c(cnst,cnstA)
  cnst_group = c(cnst_group,cnst_groupA)
  cnst_group = gsub(INTERACT_COV,"",cnst_group)
  return(list(cnst=cnst,cnst_group=cnst_group))
}


check_expdata_easy = function(expdata){
  # check expression data format
  # nessesary
  needed = c("data","gene.ann","sample.ann")
  
  indx = needed%in%names(expdata)
  if (!all(indx)){
    print(paste0("NOT found ",paste0(needed[!indx],collapse=" "),sep=" "))
    stop()
  }
  

  # check class
  if (!(is.matrix(expdata$data)|class(expdata$data)=="EList"|class(expdata$data)=="DGEList")){
    print("expression data should be a matrix or EList")
    stop()
  }
  if (!is.data.frame(expdata$gene.ann)){
    print("gene annotation data should be a data.frame")
    stop()
  }
  if (!is.data.frame(expdata$sample.ann)){
    print("sample annotation data should be a data.frame")
    stop()
  }
  
  # check sample number
  nsample = dim(expdata$data)[2]
  if (dim(expdata$sample.ann)[1]!=nsample){
    print("sample size does not mactch")
    stop()
  }

  
  # check sample name
  sampleid = colnames(expdata$data)
  if (!all(rownames(expdata$sample.ann)==sampleid)){
    print("sample name not mactch")
    stop()
  }

  
  # check gene number
  ngene = dim(expdata$data)[1]
  if (dim(expdata$gene.ann)[1]!=ngene){
    print("gene size does not mactch")
    stop()
  }

  # check gene name
  geneid = rownames(expdata$data)
  if (!all(rownames(expdata$gene.ann)==geneid)){
    print("gene name does not mactch")
    stop()
  }

  
  return(TRUE)
}




sanity_check <- function(expdata,PRIMARY_COV,NEED_FOR_PLOT){
  # check expdata
  if (!check_expdata_easy(expdata)){
    stop()
  }
  
  # checking NA
  if(any(colSums(is.na(expdata$sample.ann[,NEED_FOR_PLOT,drop=FALSE]))>0)){
    print("having NAs")
    stop()
  }
  
  # chekaing pheno type existence
  if(any(!(NEED_FOR_PLOT%in%colnames(expdata$sample.ann)))){
    print("cannot find nessesary phenotypes")
    stop()
  }
  
  # chekaing phenotype class
  indx  = !sapply(expdata$sample.ann[,NEED_FOR_PLOT,drop=FALSE],class)%in%c("factor","numeric","integer")
  if(sum(indx)>0){
    print("factor or numeric")
    stop()
  }
  # redefine factor levels
  for (i in 1:ncol(expdata$sample.ann)){
    if (is.factor(expdata$sample.ann[,i])){
      expdata$sample.ann[,i] = factor(expdata$sample.ann[,i])
    }
  }
  
  if(!is.numeric(expdata$sample.ann[[PRIMARY_COV]])){
    if(all.equal(make.names(expdata$sample.ann[[PRIMARY_COV]]),
                 as.character(expdata$sample.ann[[PRIMARY_COV]]))!=TRUE){
      
      print("variable name is not compatible with limma")
      print("renaming...")
      expdata$sample.ann[[PRIMARY_COV]] = 
        make.names(expdata$sample.ann[[PRIMARY_COV]])%>%as.factor()
      # stop()
    }
  }
  
  
  return(expdata)
}



make_design_obj_interaction<-function(pheno,spline_func,PRIMARY_COV,INTERACT_COV,ADJUST_COVS){
  # make design matrix
  if (is.numeric(pheno[[PRIMARY_COV]])){
    design <- model.matrix(as.formula(paste0("~1+",
                                             paste(c(paste0(spline_func,"*",INTERACT_COV),
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    #colnames(design)
    indx_coef = which(regexpr(spline_func,colnames(design),fixed = T)>0&
                        regexpr(":",colnames(design),fixed = T)>0&
                        regexpr(INTERACT_COV,colnames(design),fixed = T)>0)
    if(is.numeric(pheno[[INTERACT_COV]])|length(indx_coef)==1){
      colnames(design)=make.names(colnames(design))
      contrast.matrix = NULL
      contrast.group = NULL
    }else{
      cnst.obj = make_interaction_cnst(pheno,design,indx_coef,INTERACT_COV)
      colnames(design)=make.names(colnames(design))
      contrast.matrix <- makeContrasts(contrasts=cnst.obj$cnst, levels=design)
      contrast.group = cnst.obj$cnst_group
    }
  }else{
    
    
    p_cof_levels = levels(pheno[[PRIMARY_COV]])
    
    if(length(p_cof_levels)!=2){
      print("now only can handel two classes")
      stop()
    }
    
    design <- model.matrix(as.formula(paste0("~0+",
                                             paste(c(paste0(spline_func,"*",INTERACT_COV),
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef = which(regexpr(spline_func,colnames(design),fixed = T)>0&
                        regexpr(":",colnames(design),fixed = T)>0&
                        regexpr(INTERACT_COV,colnames(design),fixed = T)>0)
    colnames(design) = gsub(PRIMARY_COV,"",colnames(design))
    
    if(is.numeric(pheno[[INTERACT_COV]])){
      colnames(design)=make.names(colnames(design))
      contrast.matrix = NULL
      contrast.group = NULL
    }else{
      
      colnames(design) = gsub(INTERACT_COV,"",colnames(design))
      colnames(design)=make.names(colnames(design))
      
      if(length(indx_coef)==1){
        contrast.matrix = NULL
        contrast.group = NULL
      }else{
        # make constraints
        cnst = colnames(design)[indx_coef]
        for (ii in 1:(length(indx_coef)-1)){
          for (jj in (ii+1):length(indx_coef)){
            cnst = c(cnst,paste(colnames(design)[indx_coef[jj]],
                                "-",colnames(design)[indx_coef[ii]],sep=""))
          }
        }
        
        contrast.matrix <- makeContrasts(contrasts=cnst, levels=design)
        colname_contrast.matrix = gsub(paste0(p_cof_levels[length(p_cof_levels)],"."),"",
                                       colnames(contrast.matrix))
        colname_contrast.matrix[1:length(indx_coef)]=
          paste0(colname_contrast.matrix[1:length(indx_coef)],"-",levels(pheno[[INTERACT_COV]])[1])
        colnames(contrast.matrix) = colname_contrast.matrix
        contrast.group = colnames(contrast.matrix)
      }
    }
  }
  
  colnames(design) = gsub(INTERACT_COV,"",colnames(design))
  colnames(contrast.matrix) = gsub(INTERACT_COV,"",colnames(contrast.matrix))
  rownames(contrast.matrix) = gsub(INTERACT_COV,"",rownames(contrast.matrix))
  #contrast.group = gsub(INTERACT_COV,"",contrast.group)
  return(list(indx_coef=indx_coef,design=design,
              contrast.matrix=contrast.matrix,
              contrast.group=contrast.group))
}

make_design_obj<-function(pheno,spline_func,PRIMARY_COV,ADJUST_COVS){
  
  # make design matrix
  if (is.numeric(pheno[[PRIMARY_COV]])){
    design <- model.matrix(as.formula(paste0("~1+",
                                             paste(c(spline_func,ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef = which(regexpr(PRIMARY_COV,colnames(design),fixed = T)>0)
    
    contrast.matrix=NULL
    contrast.group=NULL
  }else{
    design <- model.matrix(as.formula(paste0("~0+",
                                             paste(c(spline_func,ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef = which(regexpr(PRIMARY_COV,colnames(design),fixed = T)>0)
    colnames(design) = gsub(PRIMARY_COV,"",colnames(design))
    
    # make constraints
    cnst = c()
    for (ii in 1:(length(indx_coef)-1)){
      for (jj in (ii+1):length(indx_coef)){
        cnst = c(cnst,paste(colnames(design)[indx_coef[jj]],
                            "-",colnames(design)[indx_coef[ii]],sep=""))
      }
    }
    
    #print(cnst)
    contrast.matrix <- makeContrasts(contrasts=cnst, levels=design)
    contrast.group = colnames(contrast.matrix)
  }
  return(list(indx_coef=indx_coef,design=design,
              contrast.matrix=contrast.matrix,
              contrast.group=contrast.group))
}


