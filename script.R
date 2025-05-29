options(stringsAsFactors = F)
options(scipen = 20000)
library(data.table)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(reshape2)
library(ggsci)
library(RColorBrewer)
library(DiffBind)
library(webshot)
library(networkD3)


##merge m6Apeak
system("cat *macs2peak.bed | sortBed -i - |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>4{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' |sort -u >meet5S_macs2peak_allpeak.bed")
system("cat *metpeak.bed | sortBed -i - |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>4{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' |sort -u >meet5S_metpeak_allpeak.bed")
system("cat meet5S_macs2peak_allpeak.bed meet5S_metpeak_allpeak.bed | sortBed -i - |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>0{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' |sort -u >overlappeak.bed")

##remove m6Am 5'UTR peaks 
system("perl m6A_annotate_forGTF_xingyang_v2.pl /Sum/lirui/database/hg19/gencode.v27lift37.annotation.gtf /Sum/lirui/project/98S_callpeak/overlappeak.bed  /Sum/lirui/project/98S_callpeak/annotation/overlappeak")

overlap.anno=read.table("/Sum/lirui/project/98S_callpeak/annotation/overlappeak.anno.txt",header = F,sep="\t")
colnames(overlap.anno)=c("Chr","Start","End","Peak.id","Chr.gene","Start.gene","End.gene","Transcript.id",
                         "Nouse","Strand","Gene","Gene.type","Gene.site","Peak.position","Ensembl.gene.id","Level.1.gene.type",
                         "Level.2.gene.type")
rownames(overlap.anno)=overlap.anno$Peak.id

anno_5UTR=overlap.anno[which(overlap.anno$Gene.site=="5UTR"),]

gtf.temp=fread("/Sum/lirui/database/hg19/gencode.v27lift37.annotation.gtf",sep="\t",skip = 5,data.table = F)
gtf.temp=cbind(gtf.temp,Transcript.id=strsplit2(strsplit2(gtf.temp$V9,split = "transcript_id ")[,2],split = ";")[,1])
gtf.temp$Transcript.id=gsub("\"","",gtf.temp$Transcript.id)
gtf.temp$Gene.site=gtf.temp$V3
gtf.temp$UTR.temp=strsplit2(gtf.temp[,9],split = ";")[,9]
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp==" exon_number 1"),"Gene.site"]="5UTR"
gtf.temp[which(gtf.temp$Gene.site=="UTR" & gtf.temp$UTR.temp!=" exon_number 1"),"Gene.site"]="3UTR"
gtf.temp$UTR.temp=NULL

write.table(strsplit2(anno_5UTR$Peak.id,split=":|-"),"/Sum/lirui/project/98S_callpeak/UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /Sum/lirui/database/hg19/genome.fa -bed /Sum/lirui/project/98S_callpeak/UTR5.peak.bed -fo /Sum/lirui/project/98S_callpeak/UTR5.peak.fa")
system("homer2 find -i /Sum/lirui/project/98S_callpeak/UTR5.peak.fa -m /Sum/lirui/database/BCA.motif -p 50 > /Sum/lirui/project/98S_callpeak/BCA_peak_offset.txt")
BCA_in_5UTR_offset=read.table("/Sum/lirui/project/98S_callpeak/BCA_peak_offset.txt",header=F)

anno_5UTR=overlap.anno[unique(BCA_in_5UTR_offset$V1),]
anno_5UTR=merge(gtf.temp,anno_5UTR,by="Transcript.id",all.y=T)
anno_5UTR=anno_5UTR[which(anno_5UTR$V3=="UTR"),]
anno_5UTR$temp.start=anno_5UTR$V4-anno_5UTR$Start
anno_5UTR$temp.end=anno_5UTR$V5-anno_5UTR$Start
anno_5UTR[which(anno_5UTR$temp.start>0),"temp.start"]=1
anno_5UTR[which(anno_5UTR$temp.start<0),"temp.start"]=(-1)
anno_5UTR[which(anno_5UTR$temp.end>0),"temp.end"]=1
anno_5UTR[which(anno_5UTR$temp.end<0),"temp.end"]=(-1)
anno_5UTR=anno_5UTR[which((anno_5UTR$temp.start*anno_5UTR$temp.end)<=0),]
anno_5UTR.bed=cbind(anno_5UTR$V1,anno_5UTR$V4,anno_5UTR$V5,anno_5UTR$Peak.id,".",anno_5UTR$V7)
write.table(anno_5UTR.bed,"/Sum/lirui/project/98S_callpeak/UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /Sum/lirui/database/hg19/genome.fa -bed /Sum/lirui/project/98S_callpeak/UTR5.peak.bed -s -name -fo /Sum/lirui/project/98S_callpeak/UTR5.peak.fa")

temp.utr5=read.table("/Sum/lirui/project/98S_callpeak/UTR5.peak.fa",sep="\n")
temp.utr5=cbind(temp.utr5,substr(temp.utr5[,1],1,1))

i=2
n=nrow(temp.utr5)
temp.utr5=cbind(temp.utr5,type=NA)
while(i<=n){
  if(temp.utr5[i,2]=="A"){
    temp.utr5[c(i-1,i),"type"]="m6Am"
  }
  i=i+2
}
m6Am=na.omit(temp.utr5)
m6Am=m6Am[grep(">",m6Am$V1),]
m6Am=gsub(">","",m6Am$V1)
m6Am=strsplit2(m6Am,split="[(]")[,1]

overlap.anno=overlap.anno[setdiff(rownames(overlap.anno),m6Am),]

peak <- m6Am
merge.bed <- fread("/Sum/lirui/project/98S_callpeak/overlappeak.bed",
                   data.table = FALSE,sep = "\t",header = FALSE)
filter.bed <- merge.bed[-which(merge.bed$V4%in%peak),]
write.table(filter.bed,"/Sum/lirui/project/98S_callpeak/removem6Am/mergepeak_removem6Am.bed4",row.names = F,col.names = F,quote = F,sep="\t")

system("intersectBed -a /Sum/lirui/project/98S_callpeak/removem6Am/mergepeak_removem6Am.bed4 -b filter_17467seRNA_locihg19.bed4 -f 0.2 -wa -wb | awk -F "\t" '{print $5"\t"$6"\t"$7"\t"$8}' |sort -u > filter_5375m6AseRNA_locihg19.bed4")


####search RRACH motif in seRNA
system("bedtools getfasta -fi /data1/database/human/hg19/genome.fa -bed filter_17467seRNA_locihg19.bed4 -fo filter_17467seRNA_locihg19.fa")
system("seqkit seq filter_17467seRNA_locihg19.fa --dna2rna > filter_17467seRNA_locihg19.RNA.fa")
system("perl runsramp.pl /Sum/lirui/project/motif/filter_17467seRNA_locihg19.fa /Sum/lirui/project/motif/filter_17467seRNA_locihg19_sramp.txt full")

####diff m6AseRNA
m6A.df <- fread("/Sum/lirui/project/98S/normalized_m6Alevel.txt",
                data.table = FALSE,header = FALSE)
Tumor.samp <- fread("/Sum/lirui/project/samplelist/65tumorID.txt",
                    data.table = FALSE,header = FALSE)

Tumor.m6a <- m6A.df[,Tumor.samp$V1]

Normal.m6a <- m6A.df[,-which(colnames(m6A.df)%in%Tumor.samp$V1)]

pair.df <- fread("/Sum/lirui/project/samplelist/33T_33N_pairinformation.txt",
                 data.table = FALSE)
pair.df$Tumor <- paste(pair.df$Tumor,"T",sep = "")
pair.df$Normal <- paste(pair.df$Normal,"N",sep = "")

pair.T.m6a <- Tumor.m6a[,pair.df$Tumor]
pair.N.m6a <- Normal.m6a[,pair.df$Normal]
pair.dm.df <- cbind(pair.N.m6a,pair.T.m6a)
pair.dm.df$eRNA <- rownames(pair.dm.df)

#wilcox test for tumor and normal
row.wilcox.test=function(common_peak,tumor_sample_in_peak,normal_sample_in_peak,paired.set=F){
  i=1
  n=nrow(common_peak)
  test=c()
  while(i<=n){
    pvalue.paired=wilcox.test(as.numeric(common_peak[i,tumor_sample_in_peak]),
                              as.numeric(common_peak[i,normal_sample_in_peak]),paired = paired.set)
    log2FC.paired=log2(mean(as.numeric(common_peak[i,tumor_sample_in_peak]))/
                         mean(as.numeric(common_peak[i,normal_sample_in_peak])))
    test=rbind(test,c(rownames(common_peak)[i],pvalue.paired$p.value,log2FC.paired))
    i=i+1
  }
  test=data.frame(test)
  colnames(test)=c("Peak.id","P.value.paired","log2FC.paired")
  test$P.value.paired=as.numeric(test$P.value.paired)
  test=cbind(test,FDR.paired=p.adjust(test$P.value.paired,method = "fdr"))
  test$log2FC.paired=as.numeric(test$log2FC.paired)
  return(test)
}

diff.common.peak=row.wilcox.test(pair.dm.df,pair.df$Tumor,pair.df$Normal,paired.set = T)
write.table(diff.common.peak,"/Sum/lirui/project/DM/5375m6AseRNA/pairedTvsN_wilcoxtest.txt",row.names = F,
            quote = FALSE,col.names = TRUE,sep = "\t")

rownames(diff.common.peak)=diff.common.peak$Peak.id
significan.diff.common.peak=diff.common.peak[which(diff.common.peak$FDR.paired<0.1),]
significan.diff.common.peak=significan.diff.common.peak[order(significan.diff.common.peak$log2FC.paired,decreasing = T),]
hyper_peaks=significan.diff.common.peak[which(significan.diff.common.peak$log2FC.paired>=log2(3/2)),]
hypo_peaks=significan.diff.common.peak[which(significan.diff.common.peak$log2FC.paired<=log2(2/3)),]

Not.pair.T.dm.df <- Tumor.m6a[,-which(colnames(Tumor.m6a)%in%pair.df$Tumor)]
Not.pair.dm.df <- cbind(pair.N.m6a,Not.pair.T.dm.df)
diff.common.peak.unpaired.normal=row.wilcox.test(Not.pair.dm.df,colnames(Not.pair.T.dm.df),
                                                 pair.df$Normal,paired.set = F)
write.table(diff.common.peak.unpaired.normal,"/Sum/lirui/project/DM/5375m6AseRNA/unpairedTvsN_wilcoxtest.txt",quote = F,row.names = F,sep="\t")
sig.diff.common.peak.unpaired.normal=diff.common.peak.unpaired.normal[which(diff.common.peak.unpaired.normal$FDR.paired<0.1),]

hyper_unpair=sig.diff.common.peak.unpaired.normal[which(sig.diff.common.peak.unpaired.normal$log2FC.paired>=log2(3/2)),]
hypo_unpair=sig.diff.common.peak.unpaired.normal[which(sig.diff.common.peak.unpaired.normal$log2FC.paired<=log2(2/3)),]

hyper.peak.id <- intersect(hyper_unpair$Peak.id,hyper_peaks$Peak.id)
hypo.peak.id <- intersect(hypo_unpair$Peak.id,hypo_peaks$Peak.id)

write.table(data.frame(hyper=hyper.peak.id),
            "/Sum/lirui/project/DM/5375m6AseRNA/pairTvsN_unpairTvsN_wilcoxtest_hyperseRNA.txt",
            row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

write.table(data.frame(hyper=hypo.peak.id),
            "/Sum/lirui/project/DM/5375m6AseRNA/pairTvsN_unpairTvsN_wilcoxtest_hyposeRNA.txt",
            row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")


#run randomforest
run_RF=function(clinic_data.for.analysis,target.peak.matrix,gene_matrix,target.gene,Clinic.factor.temp){
  get_script=build_RF_script(clinic_data.for.analysis=clinic_data.for.analysis,target.peak.matrix=target.peak.matrix,
                             gene_matrix=gene_matrix,target.gene=target.gene,Clinic.factor.temp=Clinic.factor.temp)
  if(length(get_script)<1){return(NULL)}
  write.table(get_script,"/data/lirui/Randomforest/RF_script/combine_temp/get_script.sh",
              quote = F,row.names = F,col.names = F)
  
  #################################################################
  #sh get_script.sh
  #awk puthon3 chr*.py
  #grep "{" for contribution
  #grep "[" for pvalue
  #################################################################
  system("sh /data/lirui/Randomforest/RF_script/combine_temp/get_script.sh")
  system("ls /data/lirui/Randomforest/RF_script/combine_temp/*_Construct_model_and_output.py|awk \'{print \"python3 \" $1\" > \"$1\".result\"}\'|xargs -iCMD -P30 bash -c CMD")
  system("grep {  /data/lirui/Randomforest/RF_script/combine_temp/*_Construct_model_and_output.py.result > /data/lirui/Randomforest/RF_script/combine_temp/RBP.contribution.randomforest.txt")
  system("grep \"^\\[\" /data/lirui/Randomforest/RF_script/combine_temp/*_Construct_model_and_output.py.result > /data/lirui/Randomforest/RF_script/combine_temp/RBP.pvalue.randomforest.txt")
  system("rm -r /data/lirui/Randomforest/RF_script/combine_temp/*_Construct_model_and_output.py*")
  
  #random forest contribution
  contribution=read.table("/data/lirui/Randomforest/RF_script/combine_temp/RBP.contribution.randomforest.txt",sep=",")
  
  contribution.WER_plot=RF_CDF_data(contribution,Clinic.factor.temp = Clinic.factor.temp,target.gene=target.gene)
  
  #random forest pvalue
  RF_p=read.table("/data/lirui/Randomforest/RF_script/combine_temp/RBP.pvalue.randomforest.txt",sep=",",fill=T)
  
  RF_p_plot_dcast=RF_P_bar(RF_p,Clinic.factor.temp,target.gene,p_cutoff=0.1)
  
  result=list(contribution=contribution.WER_plot,RF_p_plot_dcast=RF_p_plot_dcast)
  result=contribution.WER_plot
  return(result)
}

#make randomforest script
build_RF_script=function(gene_matrix,clinic_data.for.analysis,target.peak.matrix,target.gene,Clinic.factor.temp){
  gene_matrix <- gene_matrix  ###nuclear RBPs
  clinic_data.for.analysis <- clinical
  target.peak.matrix <- target.peak.matrix
  target.gene <- target.gene
  Clinic.factor.temp <- Clinic.factor.temp
  j=1
  k=nrow(target.peak.matrix)
  clinic_data.for.analysis=clinic_data.for.analysis[colnames(target.peak.matrix),Clinic.factor.temp]
  clinic_data.for.analysis=data.frame(target.gene=t(gene_matrix[target.gene,colnames(target.peak.matrix)]),clinic_data.for.analysis)
  get_script=c()
  while(j<=k){
    peak.temp=rownames(target.peak.matrix)[j]
    target.temp=target.peak.matrix[j,]
    p=1
    q=nrow(clinic_data.for.analysis)
    data_in_script=c()
    while(p<=q){
      data_in_script=append(data_in_script,paste("[",paste(clinic_data.for.analysis[p,],collapse = ","),"]"))
      p=p+1
    }
    temp.script=paste("sed 's/lr_data/[",paste(data_in_script,collapse = ","),
                      "]/g' /data/lirui/Randomforest/RF_script/Construct_model_and_output.py |sed 's/xingyang_target/[",
                      paste(target.temp,collapse = ","),
                      "]/g'|sed 's/feature_list_lr/[",paste("\"target.gene.",target.gene,"\"",collapse = ",",sep=""),
                      ",\"Gender\",\"Age\",\"grade\",\"Smoking\",\"Drinking\",\"Neural\",\"Vascular\",\"Lymphnode\",\"stage\"]/g'|sed 's/temp_dir_xingyang/",
                      gsub("-","_",gsub(":","_",peak.temp)),"_temp/g' >  /data/lirui/Randomforest/RF_script/combine_temp/",
                      gsub("-","_",gsub(":","_",peak.temp)),"_Construct_model_and_output.py",sep="")
    get_script=append(get_script,temp.script)
    j=j+1
  }
  return(get_script)
}

#Organize RF CDF plot
RF_CDF_data=function(contribution,Clinic.factor.temp,target.gene){
  contribution_peak=strsplit2(contribution[,1],split = "[{]")[,1]
  contribution_peak=strsplit2(contribution_peak,split = "_|/")
  contribution_peak=contribution_peak[,(grep("Construct",contribution_peak[1,])-3):
                                        (grep("Construct",contribution_peak[1,])-1)]
  contribution_peak=paste(contribution_peak[,1],":",contribution_peak[,2],"-",contribution_peak[,3],sep="")
  
  k=length(Clinic.factor.temp)+length(target.gene)
  contribution[,1]=strsplit2(contribution[,1],split = "[{]")[,2]
  contribution[,k]=gsub("}","",contribution[,k])
  contribution=as.character(as.matrix(contribution))
  contribution=cbind(contribution,rep(contribution_peak,k))
  contribution=data.frame(contribution)
  
  contribution.target.gene=contribution[grep("target.gene",contribution[,1]),]
  contribution.target.gene=cbind(contribution.target.gene,strsplit2(contribution.target.gene[,1],split=":"))
  contribution.target.gene[,4]=as.numeric(contribution.target.gene[,4])
  contribution.target.gene[,3]=gsub("target.gene.","",contribution.target.gene[,3])
  contribution.target.gene[,3]=gsub(" ","",contribution.target.gene[,3])
  contribution.target.gene=data.frame(contribution.target.gene)
  contribution.target.gene$contribution=NULL
  contribution.target.gene$X2=as.numeric(contribution.target.gene$X2)
  contribution.target.gene[which(contribution.target.gene$X2>1),"X2"]=1
  contribution.target.gene[which(contribution.target.gene$X2<(-1)),"X2"]=-1
  contribution.target.gene=contribution.target.gene[!duplicated(paste(contribution.target.gene[,1],":",contribution.target.gene[,2],sep="")),]
  contribution.target.gene_plot=contribution.target.gene[,2:3]
  return(contribution.target.gene_plot)
}

write.table(contribution.target.gene,
            "/data/lirui/seRNA.RBP.contriscore.txt",
            row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")

#Organize RF pvalue plot
RF_P_bar=function(RF_p,Clinic.factor.temp,target.gene,p_cutoff=0.25){
  k=length(Clinic.factor.temp)+length(target.gene)
  RF_peak=strsplit2(RF_p[,1],split = "\\[")[,1]
  RF_peak=strsplit2(RF_peak,split = "_|/")
  RF_peak=RF_peak[,(grep("Construct",RF_peak[1,])-3):(grep("Construct",RF_peak[1,])-1)]
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
  
  n=seq(2,2*q,by=2)
  q=dim(RF_p_plot)[2]/2
  both_sig=c()
  for(i in n){
    if(length(which(as.numeric(gsub(")","",RF_p_plot[,i]))<p_cutoff))>1){
      both_sig=append(both_sig,RF_p_plot[26,i])
    }
  }
  
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

write.table(RF_p_plot_dcast,
            "/data/lirui/seRNA.RBP.RF.Pvalue.txt",
            row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")

#WER corr with peaks
gene_corr_peak=function(peak_for_analysis,gene_for_analysis,peak_matrix,gene_matrix){
  # peak_for_analysis <- rownames(target.peak.matrix)
  # peak_matrix <- target.peak.matrix
  # gene_matrix <- target.gene_matrix
  # gene_for_analysis <- rownames(gene_matrix)
  corr_peak=peak_matrix[peak_for_analysis,]
  corr_gene=gene_matrix[gene_for_analysis,]
  write.table(corr_peak,"/data/lirui/Randomforest/RF_script/combine_temp/temp/corr_peak.txt",quote = F,row.names = T,col.names = T,sep="\t")
  write.table(corr_gene,"/data/lirui/Randomforest/RF_script/combine_temp/temp/corr_gene.txt",quote = F,row.names = T,col.names = T,sep="\t")
  i=1
  n=nrow(corr_gene)
  script=c()
  label=abs(rnorm(1))
  while(i<=n){
    script=rbind(script,paste(sep="","Rscript /data/lirui/Randomforest/RF_script/combine_temp/correlation.R /data/lirui/Randomforest/RF_script/combine_temp/temp/corr_peak.txt /data/lirui/Randomforest/RF_script/combine_temp/temp/corr_gene.txt ",
                              rownames(corr_gene)[i]," ",label," /data/lirui/Randomforest/tmp/"))
    i=i+1
  }
  write.table(script,"/data/lirui/Randomforest/RF_script/combine_temp/run.sh",quote = F,row.names = F,col.names = F)
  system("cat /data/lirui/Randomforest/RF_script/combine_temp/run.sh| xargs -iCMD -P50 bash -c CMD")
  system(paste("cat /data/lirui/Randomforest/tmp/*",label,"*.txt > /data/lirui/Randomforest/tmp/temp",sep=""))
  result=read.table("/data/lirui/Randomforest/tmp/temp",sep="\t",header=F)
  system("rm /data/lirui/Randomforest/tmp/*")
  colnames(result)=c("Pvalue","Corr","gene","Peak.id")
  return(result)
}
write.table(result,
            "/data/lirui/seRNA.RBP.corr.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")

RBP_hyper_peak_corr <- fread("/data/lirui/seRNA.RBP.corr.txt",
                             data.table = FALSE)

not_peak_for_analysis <- fread("/data/lirui/normalized_m6Alevel.txt",
                               data.table = FALSE)
rownames(not_peak_for_analysis) <- not_peak_for_analysis$peakid
hyper.peak <- fread("/data/lirui/hyperseRNAlist.txt",
                    data.table = FALSE)
not_peak_for_analysis <- not_peak_for_analysis[-which(rownames(not_peak_for_analysis)%in%hyper.peak$final.hyper),]

#Non diff peaks RBP correlation
RBP_non_peak_corr=gene_corr_peak(peak_for_analysis = rownames(not_peak_for_analysis),
                                 gene_for_analysis = rownames(gene_matrix),
                                 peak_matrix = not_peak_for_analysis,
                                 gene_matrix = gene_matrix)


#fisher test
sig_score_hyper_peak_corr=RBP_hyper_peak_corr
sig_score_non_peak_corr=RBP_non_peak_corr
temp=table(RBP_hyper_peak_corr[which(abs(RBP_hyper_peak_corr$Corr)>0.25 & 
                                       RBP_hyper_peak_corr$Pvalue<0.05),"gene"])
temp=temp[order(temp,decreasing = T)]
sig_RBP=names(temp)
i=1
n=length(sig_RBP)
enrich_test_RBP=c()
while(i<=n){
  p.temp=fisher.test(
    matrix(c(nrow(sig_score_hyper_peak_corr[which(sig_score_hyper_peak_corr$gene==sig_RBP[i]&
                                                    abs(sig_score_hyper_peak_corr$Corr)>0.25&
                                                    sig_score_hyper_peak_corr$Pvalue<0.05),])[1],
             nrow(sig_score_non_peak_corr[which(sig_score_non_peak_corr$gene==sig_RBP[i]&
                                                  abs(sig_score_non_peak_corr$Corr)>0.25&
                                                  sig_score_hyper_peak_corr$Pvalue<0.05),])[1],
             nrow(sig_score_hyper_peak_corr[which(sig_score_hyper_peak_corr$gene==sig_RBP[i]&
                                                    (abs(sig_score_hyper_peak_corr$Corr)<0.25|
                                                       sig_score_hyper_peak_corr$Pvalue>0.05)),])[1],
             nrow(sig_score_non_peak_corr[which(sig_score_non_peak_corr$gene==sig_RBP[i]&
                                                  (abs(sig_score_non_peak_corr$Corr)<0.25|
                                                     sig_score_hyper_peak_corr$Pvalue>0.05)),])[1]),2,2),
    alternative = "greater")
  enrich_test_RBP=rbind(enrich_test_RBP,c(P=p.temp$p.value,OR=p.temp$estimate,RBP=sig_RBP[i]))
  i=i+1
}
enrich_test_RBP=data.frame(enrich_test_RBP)
enrich_test_RBP$P=as.numeric(enrich_test_RBP$P)
enrich_test_RBP$OR.odds.ratio=as.numeric(enrich_test_RBP$OR.odds.ratio)

enrich_test_RBP=enrich_test_RBP[order(enrich_test_RBP$P),]
enrich_test_RBP$RBP=ordered(enrich_test_RBP$RBP,levels=enrich_test_RBP$RBP)
enrich_test_RBP$logP=-log10(enrich_test_RBP$P)

enrich_test_RBP_plot_hyper=enrich_test_RBP
enrich_test_RBP_plot=enrich_test_RBP

####circos plot of correlation analysis between CFL1 and hyperseRNA 
hyper.CFL1.corr <- fread("seRNA.RBP.corr.txt",
                         data.table = FALSE)
CFL1.corr <- hyper.CFL1.corr[which(hyper.CFL1.corr$gene=="CFL1"),]
CFL1.corr <- CFL1.corr[which(CFL1.corr$Pvalue<0.05&CFL1.corr$Corr>0.3),]
CFL1.corr$Corr <- paste("value=",CFL1.corr$Corr,sep="")
CFL1.corr$hs <- sapply(strsplit(CFL1.corr$Peak.id,":"),"[",1)
CFL1.corr$region <- sapply(strsplit(CFL1.corr$Peak.id,":"),"[",2)
CFL1.corr$start <- sapply(strsplit(CFL1.corr$region,"-"),"[",1)
CFL1.corr$end <- sapply(strsplit(CFL1.corr$region,"-"),"[",2)
CFL1.corr$hs <- gsub("chr","hs",CFL1.corr$hs)
ciscos.df <- data.frame(cfl1.hs="hs11",cfl1.start=65590493,cfl1.end=65629497,
                        peak.hs=CFL1.corr$hs,peak.start=CFL1.corr$start,peak.end=CFL1.corr$end,
                        R=CFL1.corr$Corr)
write.table(ciscos.df,
            "/Sum/lirui/project/circos/CFL1.txt",
            row.names = FALSE,col.names = FALSE,quote = FALSE,sep=" ")

system("perl /Sum/lirui/soft/circos-0.69-9/bin/circos -outputdir /Sum/lirui/project/circos/CFL1/test/  -outputfile  CFL1.png -conf circos_run.conf")


###CFL1 diff expression
allfiles=list.files(path = "/Sum/lirui/project/Tumor",full.names = TRUE)

files1 <-list.files(path = "/Sum/lirui/project/Tumor") 

fpkmFileLabels <- gsub(".genes.results","",files1)


datanew <- fread(allfiles[1],data.table = FALSE,header = TRUE)
datanew <- datanew[,c(1,5)]

colnames(datanew) <- c("gene_id",fpkmFileLabels[1])


for (i in 2:length(allfiles)) {
  data1 <- fread(allfiles[i],data.table = FALSE,header = TRUE)
  
  data1 <- data1[,c(1,5)]
  colnames(data1) <- c("gene_id",fpkmFileLabels[i])
  datanew <- merge(x=datanew,y=data1,by="gene_id")
  
}

Tumor.df <- datanew

allfiles=list.files(path = "/Sum/lirui/project/Normal",full.names = TRUE)

files1 <-list.files(path = "/Sum/lirui/project/Normal") 

fpkmFileLabels <- gsub(".genes.results","",files1)


datanew.N <- fread(allfiles[1],data.table = FALSE,header = TRUE)
datanew.N <- datanew.N[,c(1,5)]

colnames(datanew.N) <- c("gene_id",fpkmFileLabels[1])


for (i in 2:length(allfiles)) {
  data1 <- fread(allfiles[i],data.table = FALSE,header = TRUE)
  
  data1 <- data1[,c(1,5)]
  colnames(data1) <- c("gene_id",fpkmFileLabels[i])
  datanew.N <- merge(x=datanew.N,y=data1,by="gene_id")
  
}


Normal.df <- datanew.N

#####挑出不同的几种样本###

pair.df <- fread("/Sum/lirui/project/samplelist/33T_33N_pairinformation.txt",
                 data.table = FALSE)
pair.df$Tumor <- paste(pair.df$Tumor,"T",sep = "")
pair.df$Normal <- paste(pair.df$Normal,"N",sep = "")
All.t <- fread("/Sum/lirui/project/samplelist/65tumorID.txt",
               data.table = FALSE,header = FALSE)
rownames(Normal.df) <- Normal.df$gene_id
rownames(Tumor.df) <- Tumor.df$gene_id
pair.N.exp <- Normal.df[,pair.df$Normal]
pair.T.exp <- Tumor.df[,pair.df$Tumor]
all.T.exp <- Tumor.df[,All.t$V1]
single.T.exp <- all.T.exp[,-which(colnames(all.T.exp)%in%pair.df$Tumor)]
merge.df <- cbind(pair.N.exp,pair.T.exp,single.T.exp,all.T.exp)
merge.df$gene <- sapply(strsplit(rownames(merge.df),"_"),"[",1)
merge.df$symbol <- sapply(strsplit(rownames(merge.df),"_"),"[",3)
CFL1.fpkm <- merge.df[which(merge.df$symbol%in%c("CFL1")),]

rownames(CFL1.fpkm) <- CFL1.fpkm$symbol
CFL1.fpkm.163S <- CFL1.fpkm[,1:163]
CFL1.fpkm.163S <- as.data.frame(t(CFL1.fpkm.163S))
CFL1.fpkm.163S$sampleid <- rownames(CFL1.fpkm.163S)
CFL1.fpkm.163S$group <- c(rep("paired_normal",time=33),rep("pairer_tumor",time=33),rep("unpair_tumor",time=32),
                        rep("all_tumor",time=65))
write.table(CFL1.fpkm.163S,
            "/Sum/lirui/project/PAAD_pairNT_unpairT_AllSam_CFL1_count.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")

###CFL1/METTL3KO seRNA diff methylation
seRNA.m6A.df <- fread("/Sum/lirui/project/WT_KO_normalized.count.txt",
                     data.table = FALSE,header = TRUE)
seRNA.m6A.df$WT1 <- as.numeric(as.character(seRNA.m6A.df$WT1_IP))/as.numeric(as.character(seRNA.m6A.df$WT1_input))
seRNA.m6A.df$WT2 <- as.numeric(as.character(seRNA.m6A.df$WT2_IP))/as.numeric(as.character(seRNA.m6A.df$WT2_input))

seRNA.m6A.df$KO1 <- as.numeric(as.character(allseRNA.df$KO1_IP))/as.numeric(as.character(seRNA.m6A.df$KO1_input))
seRNA.m6A.df$KO2 <- as.numeric(as.character(allseRNA.df$KO2_IP))/as.numeric(as.character(seRNA.m6A.df$KO2_input))

seRNA.m6A.df[which(is.na(seRNA.m6A.df$WT1)),9] <- 0
seRNA.m6A.df[which(is.na(seRNA.m6A.df$WT2)),10] <- 0
seRNA.m6A.df[which(is.na(seRNA.m6A.df$KO1)),11] <- 0
seRNA.m6A.df[which(is.na(seRNA.m6A.df$KO2)),12] <- 0

seRNA.m6A.df$KO1vsWT1 <- as.numeric(as.character(allseRNA.df$KO1))/as.numeric(as.character(seRNA.m6A.df$WT1))
seRNA.m6A.df$KO2vsWT2 <- as.numeric(as.character(allseRNA.df$KO2))/as.numeric(as.character(seRNA.m6A.df$WT2))


seRNA.m6A.df$type="unchanged"
seRNA.m6A.df[which(seRNA.m6A.df$KO1vsWT1>1.2 & seRNA.m6A.df$KO2vsWT2>1.2),15] <- "hyper"
seRNA.m6A.df[which(seRNA.m6A.df$KO1vsWT1<(5/6) & seRNA.m6A.df$KO2vsWT2<(5/6)),15] <- "hypo"


# compare histone levels of high/low m6A-seRNA local chromatin
setwd("/Sum/lirui/project/HighLow_m6A_histone/")
system("computeMatrix reference-point --referencePoint center -p 6 -R PANC1_highm6AseRNA.bed PANC1_lowm6AseRNA.bed -S PANC1_H3K27ac.bw -o m6AseRNA_H3K27ac.10kb.gz --outFileSortedRegions m6AseRNA_H3K27ac.heatmap.10kb.bed  -a 5000 -b 5000")
system("plotProfile -m m6AseRNA_H3K27ac.10kb.gz -out m6AseRNA_H3K27ac.10kb.pdf --plotType fill --colors '#75C3F8' '#2B4071' --plotHeight 7 --plotWidth 8")

system("computeMatrix reference-point --referencePoint center -p 6 -R PANC1_highm6AseRNA.bed PANC1_lowm6AseRNA.bed -S PANC1_H3K4me3.bw -o m6AseRNA_H3K4me3.10kb.gz --outFileSortedRegions m6AseRNA_H3K4me3.10kb.heatmap.bed  -a 5000 -b 5000")
system("plotProfile -m m6AseRNA_H3K4me3.10kb.gz -out m6AseRNA_H3K4me3.10kb.pdf --plotType fill --colors '#B6CDBB' '#4F5743' --plotHeight 7 --plotWidth 8")

system("computeMatrix reference-point --referencePoint center -p 6 -R PANC1_highm6AseRNA.bed PANC1_lowm6AseRNA.bed -S PANC1_H3K4me1.bw -o m6AseRNA_H3K4me1.10kb.gz --outFileSortedRegions m6AseRNA_H3K4me1.10kb.heatmap.bed  -a 5000 -b 5000")
system("plotProfile -m m6AseRNA_H3K4me1.10kb.gz -out m6AseRNA_H3K4me1.10kb.pdf --plotType fill --colors '#D3C3C3' '#5D3D37' --plotHeight 7 --plotWidth 8")


##CFL1KO H3K4me3 diff 
setwd("/Sum/lirui/project/CFL1_histone/")
dbObj <- dba(sampleSheet=paste("NC12_SI12.csv",sep=""))
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)  
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1) ###NC vs CFL1KO
out.edgeR <- as.data.frame(comp1.edgeR)
edgeR.loss <- out.edgeR[which(out.edgeR$p.value<0.05 & out.edgeR$Fold >0.5),]
edgeR.gain <- out.edgeR[which(out.edgeR$p.value<0.05 & out.edgeR$Fold < (-0.5)),]
edgeR.loss$region <- paste(edgeR.loss$seqnames,":",edgeR.loss$start,"-",edgeR.loss$end,sep = "")
edgeR.gain$region <- paste(edgeR.gain$seqnames,":",edgeR.gain$start,"-",edgeR.gain$end,sep = "")
write.table(edgeR.loss,
            "/Sum/lirui/project/Histone/siCFL1_vs_WT_H3K4me3_edgeR_loss.txt",
            row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
write.table(edgeR.gain,
            "/Sum/lirui/project/Histone/siCFL1_vs_WT_H3K4me3_edgeR_gain.txt",
            row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
out.edgeR$log10pvalue <- log10(out.edgeR$p.value)*(-1)
out.edgeR$Signif <- "Not Signif."
out.edgeR[which(out.edgeR$p.value <0.05 & out.edgeR$Fold>= 0.5),13] <- "Loss H3K4me3"
out.edgeR[which(out.edgeR$p.value <0.05 & out.edgeR$Fold<=(-0.5)),13] <- "Gain H3K4me3"
plot.df <- out.edgeR[which(out.edgeR$Signif%in%c("Loss H3K4me3","Gain H3K4me3")),]
plot.df$LogFC <- as.numeric(as.character(plot.df$Fold))*(-1)  ##Fold represents NC versus CFL1KO 
pdf("/Sum/lirui/project/plot/CFL1KOvsWT_H3K4me3_edgeR_lossgain_histogram.pdf",
    width = 5,height = 5.5)
p <- ggplot(plot.df, aes(x=LogFC, fill=Signif))+  
  geom_histogram(binwidth = 0.05,size=0.2)+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  scale_fill_manual(values=c("#CE524E","#2F6782"))+
  xlab("Log2FC of H3K4me3 enrichment [CFL1 KO/KO-Control]")+
  ylab("Counts")+
  theme(legend.position = "top", 
        legend.text = element_text(size = 14, face = "plain", colour="black"))+
  theme(axis.text.x = element_text(size = 14, face = "plain", angle = 0, colour="black"),
        axis.text.y = element_text(size = 14, face = "plain", hjust=0, colour="black"),
        axis.title = element_text(size = 14, face = "plain", colour="black"))+
  scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,2))+
  scale_y_continuous(limits = c(0,250),expand = c(0,0))+
  annotate("text",x=-3,y=190,label="Loss 2382",size=5,colour="#2F6782")+
  annotate("text",x=3,y=190,label="Gain 262",size=5,colour="#CE524E")+
  theme(plot.margin = margin(1,1,1,1, "cm"))
print(p)
dev.off()

gain.macs2.homer <- fread("/Sum/lirui/project/Histone/siCFL1_vs_WT_H3K4me3_edgeR_gainanno.txt",
                              data.table = FALSE)
gain.macs2.homer$Annotation <- as.character(gain.macs2.homer$Annotation)
gain.macs2.homer$Annotation2 <- gain.macs2.homer$Annotation
gain.macs2.homer$Annotation2 <- sapply(strsplit(gain.macs2.homer$Annotation," \\("),"[",1)
gain.homer.statistic <- as.data.frame(table(gain.macs2.homer$Annotation2))
gain.homer.statistic$percent <- (gain.homer.statistic$Freq/sum(gain.homer.statistic$Freq))*100
gain.homer.statistic$percent <- round(gain.homer.statistic$percent,2)


loss.macs2.homer <- fread("/Sum/lirui/project/Histone/siCFL1_vs_WT_H3K4me3_edgeR_loss.anno",
                              data.table = FALSE)
loss.macs2.homer$Annotation <- as.character(loss.macs2.homer$Annotation)
loss.macs2.homer$Annotation2 <- loss.macs2.homer$Annotation
loss.macs2.homer$Annotation2 <- sapply(strsplit(loss.macs2.homer$Annotation," \\("),"[",1)
loss.homer.statistic <- as.data.frame(table(loss.macs2.homer$Annotation2))
loss.homer.statistic$percent <- (loss.homer.statistic$Freq/sum(loss.homer.statistic$Freq))*100
loss.homer.statistic$percent <- round(loss.homer.statistic$percent,2)
gain.homer.statistic$type <- "gain"
loss.homer.statistic$type <- "loss"
anno.df <- rbind(loss.homer.statistic,gain.homer.statistic)
anno.df$Var1 <- factor(anno.df$Var1,levels = rev(c("intron","Intergenic","promoter-TSS","exon","TTS","3' UTR","non-coding","5' UTR")))
anno.df$type <- factor(anno.df$type,levels = c("loss","gain"))

pdf("/Sum/lirui/project/plot/H3K4me3Cuttag_CFL1KOvsWT_lossgainpeakanno_percentagebarplot.pdf",width = 7,height = 5)
p <- ggplot(data=anno.df,aes(type,Freq,fill=Var1))+
  geom_bar(stat="identity", position="fill",color="black", width=0.8,size=0.25)+
  theme_test()+
  scale_fill_manual(values=c('#E8D1CA','#BEBFC9','#B2B8A9','#CDDDE3','#8E7DA0','#D68670','#82545F','#294E5F'))+
  xlab("")+ylab("Percentage")+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.title=element_text(size=14,face="plain",color="black"),
    legend.position = "right"
  )+coord_flip()
print(p)
dev.off()


##YTHDC2KO diff H3K4me3
setwd("/Sum/lirui/project/DC2_histone/")
dbObj <- dba(sampleSheet=paste("NC12_SI12.csv",sep=""))
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)  
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)  ##NC vs siMLL1
out.edgeR <- as.data.frame(comp1.edgeR)
edgeR.loss <- out.edgeR[which(out.edgeR$p.value<0.05 & out.edgeR$Fold >0.5),]
edgeR.gain <- out.edgeR[which(out.edgeR$p.value<0.05 & out.edgeR$Fold < (-0.5)),]
edgeR.loss$region <- paste(edgeR.loss$seqnames,":",edgeR.loss$start,"-",edgeR.loss$end,sep = "")
edgeR.gain$region <- paste(edgeR.gain$seqnames,":",edgeR.gain$start,"-",edgeR.gain$end,sep = "")
write.table(edgeR.loss,
            "siYTHDC2_vs_WT_H3K4me3_edgeR_loss.txt",
            row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
write.table(edgeR.gain,
            "siYTHDC2_vs_WT_H3K4me3_edgeR_gain.txt",
            row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
out.edgeR$log10pvalue <- log10(out.edgeR$p.value)*(-1)
out.edgeR$Signif <- "Not Signif."
out.edgeR[which(out.edgeR$p.value <0.05 & out.edgeR$Fold>= 0.5),13] <- "Loss H3K4me3"
out.edgeR[which(out.edgeR$p.value <0.05 & out.edgeR$Fold<=(-0.5)),13] <- "Gain H3K4me3"
plot.df <- out.edgeR[which(out.edgeR$Signif%in%c("Loss H3K4me3","Gain H3K4me3")),]
plot.df$LogFC <- as.numeric(as.character(plot.df$Fold))*(-1)  ##Fold represents NC versus DC2KO
pdf("/Sum/lirui/project/plot/DC2KOvsWT_H3K4me3_edgeR_lossgain_histogram.pdf",
    width = 5,height = 5.5)
p <- ggplot(plot.df, aes(x=LogFC, fill=Signif))+  
  geom_histogram(binwidth = 0.05,size=0.2)+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  scale_fill_manual(values=c("#CE524E","#2F6782"))+
  xlab("Log2FC of H3K4me3 enrichment [YTHDC2 KO/KO-Control]")+
  ylab("Counts")+
  theme(legend.position = "top", 
        legend.text = element_text(size = 14, face = "plain", colour="black"))+
  theme(axis.text.x = element_text(size = 14, face = "plain", angle = 0, colour="black"),
        axis.text.y = element_text(size = 14, face = "plain", hjust=0, colour="black"),
        axis.title = element_text(size = 14, face = "plain", colour="black"))+
  scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,2))+
  scale_y_continuous(limits = c(0,250),expand = c(0,0))+
  annotate("text",x=-3,y=190,label="Loss 2229",size=5,colour="#2F6782")+
  annotate("text",x=3,y=190,label="Gain 340",size=5,colour="#CE524E")+
  theme(plot.margin = margin(1,1,1,1, "cm"))
print(p)
dev.off()

###m6A peak distribution in mRNA
m6A_anno=fread("/Sum/lirui/project/CFL1KO_m6A/annotation/control12_merged2sam.anno.txt",header=F,sep="\t",data.table = FALSE)
WT_peak_freq=c(as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="5UTR"),14])),
               as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="CDS"),14]))+100,
               as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="3UTR"),14]))+200)
plot(table(WT_peak_freq)/length(WT_peak_freq),xlim=c(0,300),ylim=c(0,0.020),lwd=2,ylab="m6A coding peak density",
     xlab=NA,main=NA,col="#00468BFF",xaxt="n",type="l") #xingyang
abline(v=100,lty=2,col="black",lwd=1)
abline(v=200,lty=2,col="black",lwd=1)
text(x=c(80,180,280),y=-0.005,pos=2,labels=c("5' UTR","CDS","3' UTR"),xpd=TRUE,font=2)
par(new=TRUE)  ####make two lines in the same plot

KO_m6A_anno=fread("/Sum/lirui/project/CFL1KO_m6A/annotation/KOCFL1_merged2sam.anno.txt",
                  header=F,sep="\t",quote="",data.table = FALSE)
KO_peak_freq=c(as.numeric(as.vector(KO_m6A_anno[which(KO_m6A_anno[,13]=="5UTR"),14])),
               as.numeric(as.vector(KO_m6A_anno[which(KO_m6A_anno[,13]=="CDS"),14]))+100,
               as.numeric(as.vector(KO_m6A_anno[which(KO_m6A_anno[,13]=="3UTR"),14]))+200)
plot(table(KO_peak_freq)/length(KO_peak_freq),xlim=c(0,300),ylim=c(0,0.020),lwd=2,
     main=NA,col="#ED0000FF",xaxt="n",type="l",axes = FALSE,xlab = "", ylab = "")


##ATAC diff analysis
setwd("/Sum/lirui/project/ATAC/")
dbObj <- dba(sampleSheet="finalCFL1.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)  
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj) 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)   ###NC vs CFL1
out.Deseq2 <- as.data.frame(comp1.deseq)  
siCFL1vsWT.loss <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold >0.5),]
siCFL1vsWT.gain <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold < (-0.5)),]
write.table(out.Deseq2,
            "NC_vs_siCFL1_Deseq2_diffbind_result.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siCFL1vsWT.loss,
            "Deseq2_diffbind_siCFL1vsCtrl_loss.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siCFL1vsWT.gain,
            "Deseq2_diffbind_siCFL1vsCtrl_gain.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")

dbObj <- dba(sampleSheet="finalMLL1.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)  
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj) 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1) ###NC vs siMLL1
out.Deseq2 <- as.data.frame(comp1.deseq)
siMLL1vsWT.loss <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold >0.5),]
siMLL1vsWT.gain <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold < (-0.5)),]
write.table(out.Deseq2,
            "NC_vs_siMLL1_Deseq2_diffbind_result.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siMLL1vsWT.loss,
            "Deseq2_diffbind_siMLL1vsCtrl_loss.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siMLL1vsWT.gain,
            "Deseq2_diffbind_siMLL1vsCtrl_gain.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")


dbObj <- dba(sampleSheet="finalYTHDC2.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)  
plot(dbObj) 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1) ###NC vs siDC2
out.Deseq2 <- as.data.frame(comp1.deseq)
siDC2vsWT.loss <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold >0.5),]
siDC2vsWT.gain <- out.Deseq2[which(out.Deseq2$FDR<0.05 & out.Deseq2$Fold < (-0.5)),]
write.table(out.Deseq2,
            "NC_vs_siYTHDC2_Deseq2_diffbind_result.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siDC2vsWT.loss,
            "Deseq2_diffbind_siYTHDC2vsCtrl_loss.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
write.table(siDC2vsWT.gain,
            "Deseq2_diffbind_siYTHDC2vsCtrl_gain.txt",
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")


#####overlap target gene hallmarks pathway
Hallmark.df <- fread("/Sum/lirui/project/targetgene/Hallmark.txt",
                     data.table = FALSE)
Hallmark.df <- Hallmark.df[which(as.numeric(as.character(Hallmark.df$`P-value`))<0.05),]
Hallmark.df <- Hallmark.df[order(Hallmark.df$`P-value`,decreasing = FALSE),]
Hallmark.df$gene_number <- sapply(strsplit(Hallmark.df$Overlap,"/"),"[",1)
colnames(Hallmark.df) <- c("pathway","P","Combined_Score","Genes","Gene_number")
Hallmark.df <- Hallmark.df[1:11,]
data_long <- Hallmark.df
data_long$source <- "539 target gene"
data_long <- data_long[,c(6,1,5)]
colnames(data_long) <- c("source","target","value")
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1
ColourScal ='d3.scaleOrdinal() .range(["#d9534f","#B589BC","#77A1EC","#A3C4A5","#EED785","#F48F89","#7FC28E","#CD7560","#6276B2","#EDB869","#85C1D7"])'
# Make the Network
P <- sankeyNetwork(Links = data_long, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=15, nodePadding=20)
saveNetwork(P,"/Sum/lirui/project/targetgene/sankeyNetworkplot.html")
if(!is_phantomjs_installed()){
  install_phantomjs()
}
is_phantomjs_installed()
webshot("/Sum/lirui/project/targetgene/sankeyNetworkplot.html" , 
        "/Sum/lirui/project/targetgene/sankeyNetworkplot.pdf")