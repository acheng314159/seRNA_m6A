#Rscript correlation.R <peak matrix> <gene matrix> <gene> <type_lable> <output_prefix>
Args=commandArgs(trailingOnly = TRUE)

data_matrix=read.table(Args[1],row.names = 1,header = TRUE,sep="\t")
gene_exp_data_m=read.table(Args[2],row.names = 1,header = TRUE,sep="\t")

colnames(data_matrix)=gsub(".ip","",colnames(data_matrix))
colnames(gene_exp_data_m)=gsub(".ip","",colnames(gene_exp_data_m))

gene_exp_data_m=gene_exp_data_m[,intersect(colnames(gene_exp_data_m),colnames(data_matrix))]
data_matrix=data_matrix[,intersect(colnames(gene_exp_data_m),colnames(data_matrix))]
#rownames(gene_exp_data_m)=toupper(rownames(gene_exp_data_m))
#dim(gene_exp_data_m)
#dim(data_matrix)
temp1=as.numeric(gene_exp_data_m[Args[3],])
#if(sd(temp1)<1){
#	q()
#}
j=1
k=nrow(data_matrix)
cor_data=c()
while(j<=k){
  temp=cor.test(temp1,as.numeric(data_matrix[j,]),method="spearman")
  #if(abs(temp$estimate)>0.25){
  cor_data=rbind(cor_data,c(round(temp$p.value,5),round(temp$estimate,5),Args[3],rownames(data_matrix)[j]))
  #}
  j=j+1
  #	print(j)
  #print(cor_data)
}
cor_data=data.frame(cor_data)

#colnames(cor_data)=c("pvalue","cor","factor","target")
#print(paste(Args[6],"_",Args[4],".cor_data.txt",sep=""))
write.table(cor_data,paste(Args[5],Args[3],"_",Args[4],".cor_data.txt",sep=""),col.names = F,row.names = FALSE,quote = F,sep="\t")
