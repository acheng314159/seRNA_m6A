#Organize RF pvalue plot
RF_P_bar=function(RF_p,Clinic.factor.temp,WER,p_cutoff=0.25){
  k=length(Clinic.factor.temp)+length(WER)
  RF_peak=strsplit2(RF_p[,1],split = "\\[")[,1]
  RF_peak=strsplit2(RF_peak,split = "_|/")
  RF_peak=RF_peak[,(grep("temp",RF_peak[1,])+1):(grep("Construct",RF_peak[1,])-1)]
  RF_peak=paste(RF_peak[,1],":",RF_peak[,2],"-",RF_peak[,3],sep="")
  RF_p[,1]=strsplit2(RF_p[,1],split = "\\[")[,2]
  RF_p[,k]=gsub("\\]","",RF_p[,k])
  RF_p[,k*2]=gsub("\\]","",RF_p[,k*2])
  
  i=1
  n=grep("target.gene",RF_p[,1])
  for(i in n){
    temp=RF_p[i,1:k]
    temp=ordered(temp,levels=temp)
    temp=rep(temp,2)
    temp=temp[order(temp)]
    RF_p[i,]=as.character(temp)
  }
  
  RF_p=t(RF_p)
  RF_p_plot=RF_p[grep("\\)",RF_p[,2]),]
  RF_p_plot=rbind(RF_p_plot,RF_peak)
  
  # q=dim(RF_p_plot)[2]/2
  # n=seq(2,2*q,by=2)
  # both_sig=c()
  # for(i in n){
  #   if(length(which(as.numeric(gsub(")","",RF_p_plot[,i]))<p_cutoff))>1){
  #     both_sig=append(both_sig,RF_p_plot[which(rownames(RF_p_plot)=="RF_peak"),i])
  #   }
  # }
  
  n=seq(1,2*q,by=2)
  q=dim(RF_p_plot)[2]/2
  RF_p_plot_dcast=c()
  for(i in n){
    RF_p_plot_dcast=rbind(RF_p_plot_dcast,cbind(RF_p_plot[1:k,i],RF_p_plot[1:k,i+1],RF_p_plot[k+1,i]))
    
  }
  RF_p_plot_dcast[,2]=gsub("\\)","",RF_p_plot_dcast[,2])
  RF_p_plot_dcast[,1]=gsub(" ","",RF_p_plot_dcast[,1])
  RF_p_plot_dcast[,2]=gsub(" ","",RF_p_plot_dcast[,2])
  RF_p_plot_dcast=data.frame(RF_p_plot_dcast)
  RF_p_plot_dcast$X2=as.numeric(RF_p_plot_dcast$X2)
  RF_p_plot_dcast=data.frame(RF_p_plot_dcast,type="Non")
  RF_p_plot_dcast[which(RF_p_plot_dcast$X2<p_cutoff),"type"]="Sig"
  return(RF_p_plot_dcast)
}


#RBP corr with m6A-seRNAs
gene_corr_peak=function(peak_for_analysis,gene_for_analysis,peak_matrix,gene_matrix){
  corr_peak=peak_matrix[peak_for_analysis,]
  corr_gene=gene_matrix[gene_for_analysis,]
  write.table(corr_peak,"temp/corr_peak.txt",quote = F,row.names = T,col.names = T,sep="\t")
  write.table(corr_gene,"temp/corr_gene.txt",quote = F,row.names = T,col.names = T,sep="\t")
  i=1
  n=nrow(corr_gene)
  script=c()
  label=abs(rnorm(1))
  while(i<=n){
    script=rbind(script,paste(sep="","Rscript correlation.R temp/corr_peak.txt temp/corr_gene.txt ",
                              rownames(corr_gene)[i]," ",label," temp/"))
    i=i+1
  }
  write.table(script,"run.sh",quote = F,row.names = F,col.names = F)
  system("cat run.sh| xargs -iCMD -P50 bash -c CMD")
  system(paste("cat temp/*",label,"*.txt > temp/temp",sep=""))
  result=read.table("temp/temp",sep="\t",header=F)
  system("rm temp/*")
  colnames(result)=c("Pvalue","Corr","gene","Peak.id")
  return(result)
}


#different expression
diff.exp.gene.DESeq2=function(exprSet,group_list){
  exprSet <- round(exprSet,digits=0)
  colData <- data.frame(row.names=colnames(exprSet),
                        group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(exprSet),
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
  res <- results(dds,
                 contrast=c("group_list","KO","Ctrl"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  return(DESeq2_DEG)
}

diff.exp.gene.edgeR=function(exprSet,group_list){
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>0) >=0
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  colnames(fit)
  lrt <- glmLRT(fit,contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  edgeR_DEG =nrDEG 
  return(edgeR_DEG)
}

