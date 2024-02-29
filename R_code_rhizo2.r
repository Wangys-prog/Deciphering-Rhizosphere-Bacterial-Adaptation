library(seqinr)
install.packages("seqinr")
library(seqinr)
all_fasta <- read.fasta('bacteria_dna-sequences.fasta')
head(names(all_fasta))
####目标names######
target_id<- read.csv("ASV_ID.csv")
target_id[1,]
target_id[,1]
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]
names(sub_fasta)
# 写出文件
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_target_occ15.fasta')
target_id[,2]
#这一步把名字换成想要的样子
names(sub_fasta) <- gsub(names(target_id[,2]))
#这一步把名字换成想要的样子
names(sub_fasta) <- names(target_id[,2])
names(sub_fasta)
#这一步把名字换成想要的样子
names(sub_fasta) <- gsub(":.*", "", names(target_id[,2]))
names(sub_fasta)
#这一步把名字换成想要的样子
names(sub_fasta) <- gsub(target_id[,2])
names(sub_fasta) <- gsub(":.*", "", target_id[,2])
names(sub_fasta)
####目标names######
target_id<- read.csv("ASV_ID.csv")
target_id
View(target_id)
####目标names######
target_id<- read.csv("ASV_ID.csv",header =0)
View(target_id)
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]
names(sub_fasta)
#这一步把名字换成想要的样子
names(sub_fasta) <- gsub(target_id[,2])
names(sub_fasta) <- gsub(":.*", "", target_id[,2])
names(sub_fasta)
# 写出文件
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_target_occ15.fasta')
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_target_occ15.fasta')
####目标names######
target_id<- read.csv("ASV_ID.csv",header =0)
View(target_id)
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,3]]
names(sub_fasta)
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]
names(sub_fasta)
gsub("ASV_.*?(_.*$)","\\1",target_id$V2)
names(result1)=paste0(names(result1),gsub("ASV_.*?(_.*$)","\\1",target_id$V2))
#这一步把名字换成想要的样子
result1 =sub_fasta
names(result1)=paste0(names(result1),gsub("ASV_.*?(_.*$)","\\1",target_id$V2))
names(result1)
loc = match(target_id[,1],names(sub_fasta))
loc
names(sub_fasta)[loc,]
names(sub_fasta)[loc]
names(sub_fasta)[,loc]
sub_fasta[loc,]
sub_fasta
sub_fasta.items()
fasta.items()
fasta.items()
install.packages("propr")
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
setwd("../16s_function_prediction")
install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)
library(Tax4Fun2)
buildReferenceData(path_to_working_directory = ".", use_force = FALSE, install_suggested_packages = TRUE)
buildReferenceData(path_to_working_directory = '.', use_force = FALSE, install_suggested_packages = TRUE)
buildReferenceData(path_to_working_directory = '.', use_force = T, install_suggested_packages = TRUE)
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
buildDependencies(path_to_reference_data = "./Tax4Fun2_ReferenceData_v2")
getwd()
setwd("../core_csvs")
setwd("../core_asvs")
fungi_asv<- read.csv("fungi_asv_select_rhizo2.csv",header=T,row.names=1)
#去掉全为0的列
colnames(fungi_asv)
ncol(fungi_asv)
nrow(fungi_asv)
fungi_asv2<-fungi_asv[which(rowSums(fungi_asv) !=0),]
fungi_asv3<-fungi_asv2[,which(colSums(fungi_asv2) != 0)]
A<-fungi_asv3
B<-A/sum(colSums(A)) #计算col 列
C=t(B)
rowSums(C)
colSums(C)
write.csv(C,"fungi_asv_abundance.csv")
#例如只保留相对丰度总和高于 0.01% 的CSV
fungi_asv_filter<- C[,which(colSums(C)>= 0.0001)]
ncol(fungi_asv_filter)
nrow(fungi_asv_filter)
write.csv(fungi_asv_filter,"fungi_asv_filter.csv")
fungi_asv_filter1 <- fungi_asv_filter
#大于0的赋值为1
fungi_asv_filter1[fungi_asv_filter1>0] <- 1
write.csv(fungi_asv_filter1,"fungi_asv_filter1.csv")
fungi_core_abundance<-read.csv("./fungi_core_abundance.csv",header = T,row.names = 1)
View(fungi_core_abundance)
library(ggplot2)
names(fungi_core_abundance)
p3<- ggplot(fungi_core_abundance, aes(Relative.abundance...,Occurrence.frequency...,colour=Occurrence.frequency...)) +
geom_point(size = 3) +
xlab("Relative abundance (%)") +
ylab("Occurrence frequency (%)")+
geom_hline(yintercept=15, linetype="dotted",size=1.2)+
geom_vline(xintercept=0.01, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
ggsave(p3, file='fungi_core_point.pdf', width=7, height=8)
ggplot(fungi_core_abundance, aes(Relative.abundance...,Occurrence.frequency...,colour=Occurrence.frequency...)) +
geom_point(size = 3) +
xlab("Relative abundance (%)") +
ylab("Occurrence frequency (%)")+
geom_hline(yintercept=15, linetype="dotted",size=1.2)+
geom_vline(xintercept=0.01, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
ggplot(fungi_core_abundance, aes(Relative.abundance...,Occurrence.frequency...,colour=Occurrence.frequency...)) +
geom_point(size = 3) +
xlab("Relative abundance (%)") +
ylab("Occurrence frequency (%)")+
geom_hline(yintercept=10, linetype="dotted",size=1.2)+
geom_vline(xintercept=0.01, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
p3<-ggplot(fungi_core_abundance, aes(Relative.abundance...,Occurrence.frequency...,colour=Occurrence.frequency...)) +
geom_point(size = 3) +
xlab("Relative abundance (%)") +
ylab("Occurrence frequency (%)")+
geom_hline(yintercept=10, linetype="dotted",size=1.2)+
geom_vline(xintercept=0.01, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
ggsave(p3, file='fungi_core_point.pdf', width=7, height=8)
###vegan####
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
fungi_filter_group4<-cbind(fungi_asv_filter,fungi_group4)
colnames(fungi_filter_group4)
View(fungi_filter_group4)
fungi_filter_group4_sum<- t(aggregate(. ~ group4, data = fungi_filter_group4,sum))
fungi_filter_group4_sum<- fungi_filter_group4_sum[c(-1),]
View(fungi_filter_group4_sum)
fungi_filter_group4<-cbind(fungi_asv_filter,fungi_group4)
View(fungi_filter_group4)
fungi_filter_group4_sum<- t(aggregate(. ~ group4, data = fungi_filter_group4,sum))
View(fungi_filter_group4_sum)
###vegan####
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
fungi_filter_group4<-cbind(fungi_asv_filter,fungi_group4)
fungi_filter_group4_sum<- t(aggregate(. ~ group4, data = fungi_filter_group4,sum))
View(fungi_filter_group4_sum)
fungi_filter_group4_sum<- fungi_filter_group4_sum[c(-1),]
ss <- t(fungi_filter_group4_sum)
install.packages("ggvenn")
install.packages("ggvenn")
library(ggvenn) #加载ggvenn包
gene_set1 <-names(ss[1,])[ss[1,]>0]
gene_set2 <- names(ss[2,])[ss[2,]>0]
base::intersect(gene_set1, gene_set2) #获取两个基因集的交集
union(gene_set1, gene_set2) #获取两个基因集的并集
setdiff(gene_set1, gene_set2) #获取gene_set1减去gene_set2的差集
#[1] "BIRC3"   "MAP3K14" "MAP3K8"  "CREB5"   "CCL2"    "CCL20"
setdiff(gene_set2, gene_set1) #获取gene_set2减去gene_set1的差集
#[1] "FOSB"    "JUND"    "FOSL1"   "TNFAIP3" "CXCL8"   "CXCL10"
a <- list(`Bacteria_Bulk` = gene_set1,
`Bacteria_Rhizosphere` = gene_set2) #将基因集变成列表变量
p2 <- ggvenn(a, show_elements = FALSE,fill_color = c("green", "red"),
label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
text_size = 4)
p2
ss <- t(fungi_filter_group4_sum)
install.packages("ggvenn")
library(ggvenn) #加载ggvenn包
gene_set1 <-names(ss[1,])[ss[1,]>0]
gene_set2 <- names(ss[2,])[ss[2,]>0]
base::intersect(gene_set1, gene_set2) #获取两个基因集的交集
union(gene_set1, gene_set2) #获取两个基因集的并集
setdiff(gene_set1, gene_set2) #获取gene_set1减去gene_set2的差集
#[1] "BIRC3"   "MAP3K14" "MAP3K8"  "CREB5"   "CCL2"    "CCL20"
setdiff(gene_set2, gene_set1) #获取gene_set2减去gene_set1的差集
#[1] "FOSB"    "JUND"    "FOSL1"   "TNFAIP3" "CXCL8"   "CXCL10"
a <- list(`Bacteria_Bulk` = gene_set1,
`Bacteria_Rhizosphere` = gene_set2) #将基因集变成列表变量
#[1] "FOSB"    "JUND"    "FOSL1"   "TNFAIP3" "CXCL8"   "CXCL10"
a <- list(`Fungi_Bulk` = gene_set1,
`Fungi_Rhizosphere` = gene_set2) #将基因集变成列表变量
p2 <- ggvenn(a, show_elements = FALSE,fill_color = c("green", "red"),
label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
text_size = 4)
p2
ggsave(p2, file='fungi_core_venn.pdf', width=9, height=6)
View(asv_occ15)
fungi_taxa_rhizo2<- read.csv("./fungi_taxa_rhizo2.csv",header = T,row.names = 1)
fungi_taxa_rhizo2<- read.csv("./fungi_taxa_rhizo2.csv",header = T,row.names = 1)
row.names(fungi_taxa_rhizo2)<-fungi_taxa_rhizo2[,1]
fungi_taxa_rhizo2<-fungi_taxa_rhizo2[,-1]
asv_occ10<- read.csv("clipboard")
View(asv_occ10)
asv_occ10_taxa<-fungi_taxa_rhizo2[asv_occ10$ASV_ID,]
fungi_asv_select_rhizo2<- read.csv("fungi_asv_select_rhizo2.csv",header = T,row.names = 1)
asv_occ10_table<-fungi_asv_select_rhizo2[asv_occ10$ASV_ID,]
asv_occ10_table2<-asv_occ10_table[which(rowSums(asv_occ10_table) !=0),]
asv_occ10_table3<-asv_occ10_table2[,which(colSums(asv_occ10_table2) != 0)]
View(asv_occ10_table3)
ncol(asv_occ10_table3)
nrow(asv_occ10_table)
View(asv_occ10_table)
write.csv(asv_occ10_table,"fungi_asv_occ10_table.csv")
c<- asv_occ10_taxa$c
c
o<- asv_occ10_taxa$o
f<- asv_occ10_taxa$f
g<- asv_occ10_taxa$g
asv_occ10_table_c<- cbind(asv_occ10_table3,c)
asv_occ10_table_o<- cbind(asv_occ10_table3,o)
asv_occ10_table_f<- cbind(asv_occ10_table3,f)
asv_occ10_table_g<- cbind(asv_occ10_table3,g)
occ10_c_sum<-aggregate(. ~c, data =asv_occ10_table_c,sum) #忽略参数不对提示
occ10_o_sum<-aggregate(. ~o, data =asv_occ10_table_o,sum)
occ10_f_sum<-aggregate(. ~f, data =asv_occ10_table_f,sum)
occ10_g_sum<-aggregate(. ~g, data =asv_occ10_table_g,sum)
occ10_g_sum[1,1]<-"unassigned"
rownames(occ10_g_sum)<-occ10_g_sum[,1]
occ10_g_sum1<-occ10_g_sum[,-1]
g_occ10_table_abundance<-t(t(occ10_g_sum1)/colSums(occ10_g_sum1))
View(g_occ10_table_abundance)
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
group <-fungi_group4[colnames(g_occ10_table_abundance),]
occ10_g_sum2<- data.frame(cbind(t(g_occ10_table_abundance),group))
View(occ10_g_sum2)
write.csv(occ10_g_sum2,"fungi_occ10_g_sum2.csv")
occ10_g_sum2 <- read.csv("fungi_occ10_g_sum2.csv",header = T,row.names = 1)
occ10_g_sum3 <- aggregate(. ~group, data = as.data.frame(occ10_g_sum2), mean)
View(occ10_g_sum3)
write.csv(occ10_g_sum3,"fungi_occ10_g_sum3.csv")
occ10_f_sum[1,1]<-"unassigned"
rownames(occ10_f_sum)<-occ10_f_sum[,1]
occ10_f_sum1<-occ10_f_sum[,-1]
f_occ10_table_abundance<-t(t(occ10_f_sum1)/colSums(occ10_f_sum1))
View(f_occ10_table_abundance)
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
group <-fungi_group4[colnames(f_occ10_table_abundance),]
occ10_f_sum2<- data.frame(cbind(t(f_occ10_table_abundance),group))
View(occ10_f_sum2)
write.csv(occ10_f_sum2,"fungi_occ10_f_sum2.csv")
occ10_f_sum2 <- read.csv("fungi_occ10_f_sum2.csv",header = T,row.names = 1)
occ10_f_sum3 <- aggregate(. ~group, data = as.data.frame(occ10_f_sum2), mean)
View(occ10_f_sum3)
write.csv(occ10_f_sum3,"fungi_occ10_f_sum3.csv")
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
getwd()
setwd("../core_asvs")
occ10_o_sum
occ10_o_sum[1,1]<-"unassigned"
rownames(occ10_o_sum)<-occ10_o_sum[,1]
occ10_o_sum1<-occ10_o_sum[,-1]
o_occ10_table_abundance<-t(t(occ10_o_sum1)/colSums(occ10_o_sum1))
View(o_occ10_table_abundance)
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
group <-fungi_group4[colnames(o_occ10_table_abundance),]
occ10_o_sum2<- data.frame(cbind(t(o_occ10_table_abundance),group))
View(occ10_o_sum2)
write.csv(occ10_o_sum2,"fungi_occ10_o_sum2.csv")
occ10_o_sum2 <- read.csv("fungi_occ10_o_sum2.csv",header = T,row.names = 1)
occ10_o_sum3 <- aggregate(. ~group, data = as.data.frame(occ10_o_sum2), mean)
View(occ10_o_sum3)
write.csv(occ10_o_sum3,"fungi_occ10_o_sum3.csv")
occ10_c_sum[1,1]<-"unassigned"
rownames(occ10_c_sum)<-occ10_c_sum[,1]
occ10_c_sum1<-occ10_c_sum[,-1]
c_occ10_table_abundance<-t(t(occ10_c_sum1)/colSums(occ10_c_sum1))
View(c_occ10_table_abundance)
fungi_group4<- read.csv("fungi_group4.csv",header = T,row.names = 1)
group <-fungi_group4[colnames(c_occ10_table_abundance),]
occ10_c_sum2<- data.frame(cbind(t(c_occ10_table_abundance),group))
View(occ10_c_sum2)
write.csv(occ10_c_sum2,"fungi_occ10_c_sum2.csv")
occ10_c_sum2 <- read.csv("fungi_occ10_c_sum2.csv",header = T,row.names = 1)
occ10_c_sum3 <- aggregate(. ~group, data = as.data.frame(occ10_c_sum2), mean)
View(occ10_c_sum3)
write.csv(occ10_c_sum3,"fungi_occ10_c_sum3.csv")
##转置求相对丰度
####bar_plot########
library(ggplot2)
library(RColorBrewer)
library(hrbrthemes)
library(ggsci)
fungi_core_barplot<- read.csv("fungi_core_barplot_data.csv",header = T)
name2<- factor(fungi_core_barplot$taxa_1,ordered=TRUE,levels=c("Class_Bulk","Class_Rhizosphere",
"Order_Bulk","Order_Rhizosphere",
"Family_Bulk","Family_Rhizosphere",
"Genus_Bulk","Genus_Rhizosphere"))
p_bar<-ggplot(data=fungi_core_barplot,mapping=aes(x=Relative.abundance..., y=name2,fill=Taxa))+
geom_bar(stat="identity",position=position_stack(0.90),color="black", width=.7)+
labs(x = "Relative abundance(%)", y = "") +
theme_set(theme_bw())+
geom_text(aes(label=group1),size=4,hjust = 0, vjust =0,angle=45)+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15),
axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=0,size=15),
axis.title.x = element_text(face = "bold",colour="black",size=15),
axis.title.y = element_text(face = "bold",colour="black",size=15),
legend.position='none')
ggplot(data=fungi_core_barplot,mapping=aes(x=Relative.abundance..., y=name2,fill=Taxa))+
geom_bar(stat="identity",position=position_stack(0.90),color="black", width=.7)+
labs(x = "Relative abundance(%)", y = "") +
theme_set(theme_bw())+
geom_text(aes(label=group1),size=4,hjust = 0, vjust =0,angle=45)+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15),
axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=0,size=15),
axis.title.x = element_text(face = "bold",colour="black",size=15),
axis.title.y = element_text(face = "bold",colour="black",size=15),
legend.position='none')
fungi_core_barplot<- read.csv("fungi_core_barplot_data.csv",header = T)
ggplot(data=fungi_core_barplot,mapping=aes(x=Relative.abundance..., y=name2,fill=Taxa))+
geom_bar(stat="identity",position=position_stack(0.90),color="black", width=.7)+
labs(x = "Relative abundance(%)", y = "") +
theme_set(theme_bw())+
geom_text(aes(label=group1),size=4,hjust = 0, vjust =0,angle=45)+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15),
axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=0,size=15),
axis.title.x = element_text(face = "bold",colour="black",size=15),
axis.title.y = element_text(face = "bold",colour="black",size=15),
legend.position='none')
ggsave(p_bar, file='fungi_core_barplot.pdf', width=10, height=6)
p_bar<-ggplot(data=fungi_core_barplot,mapping=aes(x=Relative.abundance..., y=name2,fill=Taxa))+
geom_bar(stat="identity",position=position_stack(0.90),color="black", width=.7)+
labs(x = "Relative abundance(%)", y = "") +
theme_set(theme_bw())+
geom_text(aes(label=group1),size=4,hjust = 0, vjust =0,angle=45)+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15),
axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=0,size=15),
axis.title.x = element_text(face = "bold",colour="black",size=15),
axis.title.y = element_text(face = "bold",colour="black",size=15),
legend.position='none')
ggsave(p_bar, file='fungi_core_barplot.pdf', width=10, height=6)
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
getwd()
setwd("../16s_function_prediction")
getwd()
setwd("./bacteria_core_blast_Ref100NR")
dir()
######Taxa4Fun############################
bacteria_function_prediction<- read.csv("./bacteria_functional_prediction.csv",header = T,row.names = 1)
bacteria_group4
View(bacteria_group4)
bacteria_function_prediction_t<- cbind(t(bacteria_function_prediction),bacteria_group4)
bacteria_function_prediction_t<- cbind(t(bacteria_function_prediction[,-1]),bacteria_group4)
names(bacteria_function_prediction_t)
rownames(bacteria_function_prediction_t)
bacteria_function_prediction_t<- t(bacteria_function_prediction[,-1])
######Taxa4Fun############################
bacteria_ko_prediction<- read.csv("./bacteria_functional_prediction.csv",header = T,row.names = 1)
bacteria_ko_prediction_t<- t(bacteria_ko_prediction[,-1])
bacteria_Ko_dist <- vegdist(bacteria_class_select_t, "bray")
######Taxa4Fun############################
libaray(vegan)
######Taxa4Fun############################
library(vegan)
bacteria_ko_prediction<- read.csv("./bacteria_functional_prediction.csv",header = T,row.names = 1)
bacteria_Ko_dist <- vegdist(bacteria_class_select_t, "bray")
bacteria_ko_prediction_t<- t(bacteria_ko_prediction[,-1])
bacteria_Ko_dist <- vegdist(bacteria_class_select_t, "bray")
bacteria_Ko_dist <- vegdist(bacteria_ko_select_t, "bray")
bacteria_Ko_dist <- vegdist(bacteria_ko_prediction_t, "bray")
names(bacteria_ko_prediction[,-1])
View(bacteria_ko_prediction_t)
bacteria_Ko_dist <- vegdist(bacteria_ko_prediction_t, "bray")
bacteria_ko_prediction<- read.csv("./bacteria_functional_prediction2.csv",header = T,row.names = 1)
bacteria_ko_prediction_t<- t(bacteria_ko_prediction[,-1])
bacteria_Ko_dist <- vegdist(bacteria_ko_prediction_t, "bray")
bacteria_ko_prediction_t2<-bacteria_ko_prediction_t[which(rowSums(bacteria_ko_prediction_t) !=0),]
bacteria_ko_prediction_t3<-bacteria_ko_prediction_t2[,which(colSums(bacteria_ko_prediction_t2) != 0)]
bacteria_Ko_dist <- vegdist(bacteria_ko_prediction_t3, "bray")
nmds_dis <- metaMDS(bacteria_Ko_dist, k = 2)
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
write.csv(nmds_dis_site,"bacteria_ko_nmds_dis_sit.csv")
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("bacteria_group4.csv",header = T,row.names=1)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
View(map)
ggplot(merged, aes(MDS1, MDS2)) +
geom_point(aes(colour = factor(group1),shape=factor(group4)), size = 3) +
xlab("MDS1") +
ylab("MDS2")+
geom_hline(yintercept = 0, linetype='dotted', size = 1.2)+
geom_vline(xintercept=0, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
library(ggplot2)
ggplot(merged, aes(MDS1, MDS2)) +
geom_point(aes(colour = factor(group4),shape=factor(group4)), size = 3) +
xlab("MDS1") +
ylab("MDS2")+
geom_hline(yintercept = 0, linetype='dotted', size = 1.2)+
geom_vline(xintercept=0, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
?anosim
View(merged)
####################
dune.ano <- with(merged, anosim(nmds_dis, merged$group4))
####################
dune.ano <- with(merged, anosim(nmds_dis, group4))
####################
dune.ano <- with(merged, anosim(bacteria_Ko_dist, group4))
summary(dune.ano)
ggsave(p1, file='bacteria_ko_nmds.pdf', width=9, height=6)
p1<-ggplot(merged, aes(MDS1, MDS2)) +
geom_point(aes(colour = factor(group4),shape=factor(group4)), size = 3) +
xlab("MDS1") +
ylab("MDS2")+
geom_hline(yintercept = 0, linetype='dotted', size = 1.2)+
geom_vline(xintercept=0, linetype="dotted",size=1.2)+
theme_set(theme_bw())+
theme(axis.line=element_blank(),
panel.background = element_rect(fill = "white", color="black",size=2),
panel.grid = element_blank(),
text=element_text(face = "bold",colour="black",size=15),
legend.key=element_rect(fill='transparent', color='transparent'),
axis.text = element_text(face = "bold",colour="black",size=15))
ggsave(p1, file='bacteria_ko_nmds.pdf', width=9, height=6)
source("F:/wangys_uestc_data/rhizosphere_microbe/core_asvs/R_code_core_asvs.R")
library(propr)
install.packages("propr")
library(propr)
counts <- matrix(rpois(20*50, 100), 20, 50)
View(counts)
group <- sample(c("A", "B"), size = 20, replace = TRUE)
group
devtools::install_github("tpq/propr")
library(propr)
pr <- propr(counts, # rows as samples, like it should be
metric = "rho", # or "phi", "phs", "cor", "vlr"
ivar = "clr", # or can use "iqlr" instead
alpha = NA, # use to handle zeros
p = 100) # used by updateCutoffs
pr
updateCutoffs(pr,
cutoff = seq(0, 1, .05), # cutoffs at which to estimate FDR
ncores = 1) # parallelize here
?propr
show(object)

############差异性物种#######
# https://zhuanlan.zhihu.com/p/547246100
library(tidyverse)
library(rstatix)
bacteria_asv_genus<- read.csv("bacteria_asv_genus_top50.csv")
dir()
bacteria_group4<- read.csv("bacteria_group4.csv")
# 表达量+分组 合并，整理
df = bacteria_asv_genus %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(genus,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(genus) %>%
  t_test(value ~ group4,var.equal=T) # t_test 方法源于rstatix包; var.equal=T等方差
dfP
# 对P值FDR校正
dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  # p.adjust对P值FDR校正,算法选择BH
dfP_FDR

#整理合并表格
dfdata = dfP %>%
  as_tibble() %>%
  left_join(dfFC) %>%
  left_join(dfP_FDR)
dfdata
write.csv(dfdata,"bacteria_differential_species.csv")
###绘制火山图#####
library(ggplot2)
dataset<- read.csv("bacteria_differential_species.csv",header = T,row.names = 1)
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dataset$change = ifelse(dataset$FDR < cut_off_pvalue & abs(dataset$FC) >= cut_off_logFC, 
                        ifelse(dataset$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)
p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","red"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
    axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p
ggsave(p, file='bacteria_differential_volcano.pdf', width=8, height=6)
install.packages("ggrepel")
library(ggrepel)
p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","red"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  geom_label_repel(aes(label =genus),data=dataset, nudge_y = 2,  alpha = 0.6  )+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p
ggsave(p, file='bacteria_differential_volcano2.pdf', width=8, height=6)
#########fungi###############
library(tidyverse)
library(rstatix)
fungi_asv_genus<- read.csv("fungi_asv_genus_top60.csv")
dir()
fungi_group4<- read.csv("fungi_group4.csv")
# 表达量+分组 合并，整理
df = fungi_asv_genus %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(fungi_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(genus,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(genus) %>%
  t_test(value ~ group4,var.equal=T) # t_test 方法源于rstatix包; var.equal=T等方差
dfP
# 对P值FDR校正
dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  # p.adjust对P值FDR校正,算法选择BH
dfP_FDR

#整理合并表格
dfdata = dfP %>%
  as_tibble() %>%
  left_join(dfFC) %>%
  left_join(dfP_FDR)
dfdata
write.csv(dfdata,"fungi_differential_species.csv")
###绘制火山图#####
library(ggplot2)
dataset<- read.csv("fungi_differential_species.csv",header = T,row.names = 1)
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dataset$change = ifelse(dataset$FDR < cut_off_pvalue & abs(dataset$FC) >= cut_off_logFC, 
                        ifelse(dataset$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)
p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","red"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-5,22))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p
ggsave(p, file='fungi_differential_volcano.pdf', width=8, height=6)
install.packages("ggrepel")
library(ggrepel)
p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","red"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-5,22))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  geom_label_repel(aes(label =genus),data=dataset, nudge_y = 2,  alpha = 0.6  )+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p
ggsave(p, file='fungi_differential_volcano2.pdf', width=8, height=6)

######计算遗传聚类######
library(ape)
library(tidyverse)

path = "……/R语言计算DNA序列的遗传距离"

setwd(path)

fas.F = read.FASTA("seqdumpB.fas")

mat = dist.dna(fas.F,as.matrix = T)

write.csv(mat,file = "genetic distance.csv")

#####筛选差异属######
##bacteria### "g_Dyadobacter","g_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
#"g_Flavobacterium","g_Mucilaginibacter","g_Luteolibacter","g_Streptomyces","g_Caulobacter",
#"g_Pseudomonas","g_Burkholderia-Caballeronia-Paraburkholderia","g_Pedobacter","g_Massilia"
#"g_S0134_terrestrial_group","g_Chryseobacterium","g_Novosphingobium"
bacteria_asv_select_rhizo<- read.csv("bacteria_asv_select_rhizo.csv",header = T,row.names = 1)
g_Dyadobacter <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Dyadobacter"))
g_Rhizobium <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"))
g_Flavobacterium <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Flavobacterium"))
g_Mucilaginibacter <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Mucilaginibacter"))
g_Luteolibacter <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Luteolibacter"))
g_Streptomyces <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Streptomyces"))
g_Caulobacter <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Caulobacter"))
g_Pseudomonas <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Pseudomonas"))
g_Paraburkholderia <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Burkholderia-Caballeronia-Paraburkholderia"))
g_Pedobacter <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Pedobacter"))
g_Massilia <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Massilia"))
g_S0134 <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_S0134_terrestrial_group"))
g_Chryseobacterium <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Chryseobacterium"))
g_Novosphingobium <- bacteria_asv_select_rhizo %>% filter(str_detect(g,"g_Novosphingobium"))

write.csv(g_Dyadobacter,"bacteria_g_Dyadobacter.csv")
write.csv(g_Rhizobium,"bacteria_g_Rhizobium.csv")
write.csv(g_Flavobacterium,"bacteria_g_Flavobacterium.csv")
write.csv(g_Mucilaginibacter,"bacteria_g_Mucilaginibacter.csv")
write.csv(g_Luteolibacter,"bacteria_g_Luteolibacter.csv")
write.csv(g_Streptomyces,"bacteria_g_Streptomyces.csv")
write.csv(g_Caulobacter,"bacteria_g_Caulobacter.csv")
write.csv(g_Pseudomonas,"bacteria_g_Pseudomonas.csv")
write.csv(g_Paraburkholderia,"bacteria_g_Paraburkholderia.csv")
write.csv(g_Pedobacter,"bacteria_g_Pedobacter.csv")
write.csv(g_Massilia,"bacteria_g_Massilia.csv")
write.csv(g_S0134,"bacteria_g_S0134.csv")
write.csv(g_Chryseobacterium,"bacteria_g_Chryseobacterium.csv")
write.csv(g_Novosphingobium,"bacteria_g_Novosphingobium.csv")





####根据序列号抽提fasta###########

library(seqinr)
all_fasta <- read.fasta('bacteria_dna-sequences.fasta')
head(names(all_fasta))
####目标names######
target_id<-g_Novosphingobium
View(target_id)
rownames(target_id)
target_id$asv_id
target_id[,1]
#选取
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(target_id)]
names(sub_fasta)
# 写出文件
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Novosphingobium.fasta')

##seqkit###
## seqkit replace --ignore-case --kv-file rename.txt --pattern "^(\w+)" --replacement "{kv}" genome.fa -o genome.new.fa

# seqkit replace --ignore-case --kv-file bacteria_g_Caulobacter.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Caulobacter.fasta -o bacteria_g_Caulobacter_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Rhizobium.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Rhizobium.fasta -o bacteria_g_Rhizobium_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Flavobacterium.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Flavobacterium.fasta -o bacteria_g_Flavobacterium_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Mucilaginibacter.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Mucilaginibacter.fasta -o bacteria_g_Mucilaginibacter_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Luteolibacter.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Luteolibacter.fasta -o bacteria_g_Luteolibacter_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Streptomyces.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Streptomyces.fasta -o bacteria_g_Streptomyces_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Dyadobacter.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Dyadobacter.fasta -o bacteria_g_Dyadobacter_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Pseudomonas.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Pseudomonas.fasta -o bacteria_g_Pseudomonas_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Paraburkholderia.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Paraburkholderia.fasta -o bacteria_g_Paraburkholderia_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Pedobacter.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Pedobacter.fasta -o bacteria_g_Pedobacter_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Massilia.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Massilia.fasta -o bacteria_g_Massilia_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_S0134.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_S0134.fasta -o bacteria_g_S0134_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Chryseobacterium.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Chryseobacterium.fasta -o bacteria_g_Chryseobacterium_new.fasta
# seqkit replace --ignore-case --kv-file bacteria_g_Novosphingobium.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_g_Novosphingobium.fasta -o bacteria_g_Novosphingobium_new.fasta


######mafft##########
# mafft --auto bacteria_g_Caulobacter_new.fasta > bacteria_g_Caulobacter_new_mafft.fasta
# mafft --auto bacteria_g_Rhizobium_new.fasta > bacteria_g_Rhizobium_new_mafft.fasta
# mafft --auto bacteria_g_Flavobacterium_new.fasta > bacteria_g_Flavobacterium_new_mafft.fasta
# mafft --auto bacteria_g_Mucilaginibacter_new.fasta > bacteria_g_Mucilaginibacter_new_mafft.fasta
# mafft --auto bacteria_g_Luteolibacter_new.fasta > bacteria_g_Luteolibacter_new_mafft.fasta
# mafft --auto bacteria_g_Streptomyces_new.fasta > bacteria_g_Streptomyces_new_mafft.fasta
# mafft --auto bacteria_g_Dyadobacter_new.fasta > bacteria_g_Dyadobacter_new_mafft.fasta
# mafft --auto bacteria_g_Pseudomonas_new.fasta > bacteria_g_Pseudomonas_new_mafft.fasta
# mafft --auto bacteria_g_Paraburkholderia_new.fasta > bacteria_g_Paraburkholderia_new_mafft.fasta
# mafft --auto bacteria_g_Pedobacter_new.fasta > bacteria_g_Pedobacter_new_mafft.fasta
# mafft --auto bacteria_g_Massilia_new.fasta > bacteria_g_Massilia_new_mafft.fasta
# mafft --auto bacteria_g_S0134_new.fasta > bacteria_g_S0134_new_mafft.fasta
# mafft --auto bacteria_g_Chryseobacterium_new.fasta > bacteria_g_Chryseobacterium_new_mafft.fasta
# mafft --auto bacteria_g_Novosphingobium_new.fasta > bacteria_g_Novosphingobium_new_mafft.fasta

####
library(ggplot2)
library(ape)
library(dplyr)
setwd("../Niche_differentiation/bacteria_niche")
ex.dna4 <- read.dna("bacteria_g_Caulobacter_new_mafft.fasta", format = "fasta")

#计算遗传距离
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance
image(ex.dna4)
image(ex.dna4, "-", "blue") #只显示缺失值为蓝色
image(ex.dna4, c("g", "c"), "green") # 显示G+C为绿色
image(ex.dna4, "g", "yellow", "black") # 显示G为换色，背景为黑色
dev.off ()
g_Caulobacter2<-select(g_Caulobacter,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Caulobacter2)<-g_Caulobacter2[,1]
g_Caulobacter3<-g_Caulobacter2[,-1]
g_Caulobacter4<-g_Caulobacter3[which(rowSums(g_Caulobacter3) !=0),]
g_Caulobacter5<-g_Caulobacter4[,which(colSums(g_Caulobacter4) != 0)]

###############################
counts <- matrix(rpois(20*50, 100), 20, 50)
group <- sample(c("A", "B"), size = 20, replace = TRUE)
devtools::install_github("tpq/propr")
library(propr)
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
pr2<-updateCutoffs(pr,
              cutoff = seq(0.05, 0.05, .05), # cutoffs at which to estimate FDR
              ncores = 1) # parallelize here
pd <- propd(counts,
            group, # a vector of 2 or more groups
            alpha = NA, # whether to handle zeros
            weighted = TRUE, # whether to weigh log-ratios
            p = 100) # used by updateCutoffs
setDisjointed(pd)
setEmergent(pd)
pd <- updateF(pd,
              moderated = FALSE, # moderate stats with limma-voom
              ivar = "clr") # used for moderation
?getResults # get results in long-format
?getMatrix # get results as a square matrix
?getAdj # get an adjacency matrix

getMatrix(pr)
getResults(pr2)
getAdj(pr) 
getMatrix(pd)
getResults(pd)
getAdj(pd)
write.csv(getResults(pr),"bacteria_g_Caulobacter_rho.csv")
###################
library(ape)
library(ggplot2)
library(propr)
counts <- as.matrix(t(g_Caulobacter5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
ex.dna4 <- read.dna("bacteria_g_Caulobacter_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4)
#ggsave(p, file='bacteria_g_Caulobacter_dnadist.pdf', width=7, height=8)
dev.off()
#计算遗传距离
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
#然后将上三角阵转换为NA，再使用reshape2::melt函数即可：
## 替换上三角阵为NA
# ex.dna4.distance[upper.tri(ex.dna4.distance)] = NA
# ex.dna4.distance2<-reshape2::melt(ex.dna4.distance, na.rm=T)

write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Caulobacter_dnadist.csv")

ex.gene4.distance<- dist.gene(ex.dna4, method = "pairwise")
ex.gene4.distance
write.csv(as.matrix(ex.gene4.distance),"bacteria_g_Caulobacter_genedist.csv")
################将对称矩阵转换成两两对应关系#################
mtrx2cols = function(m1,m2,val1,val2){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt], val2= m2[lt]) #按列依次获取矩阵下半角的元素
  names(res)[3:4] = c(val1,val2) #对后两列重命名，支持多个矩阵合并
  return(res)
}
mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Caulobacter_dnadist.csv",header = T,row.names = 1)
View(ex_dna4_distance)
res = mtrx2cols(ex_dna4_distance,'val')
write.csv(res,"bacteria_g_Caulobacter_dnadist2.csv")
ex_gene4_distance<- read.csv("bacteria_g_Caulobacter_genedist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_gene4_distance,'val')
write.csv(res,"bacteria_g_Caulobacter_genedist2.csv")
rho1<- read.csv("bacteria_g_Caulobacter_rho_merge.csv",header = T,row.names = 1)
library(splines)
library(ggpmisc)
set.seed(1234)
my.formula <-rho1$lrv~ rho1$nucl_dist
p1<- ggplot(rho1, aes(nucl_dist, lrv)) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula = y ~ x ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label..,sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.1, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Caulobacter") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(family='', hjust=0.2, vjust=1,face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p1, file='bacteria_g_Caulobacter_geomline.pdf', width=9, height=7)
p1
#################bacteria_g_Chryseobacterium##################
########读取文件###########
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
g_Chryseobacterium<-read.csv("bacteria_g_Chryseobacterium.csv",header = T,row.names = 1)
g_Chryseobacterium2<-select(g_Chryseobacterium,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Chryseobacterium2)<-g_Chryseobacterium2[,1]
g_Chryseobacterium3<-g_Chryseobacterium2[,-1]
g_Chryseobacterium4<-g_Chryseobacterium3[which(rowSums(g_Chryseobacterium3) !=0),]
g_Chryseobacterium5<-g_Chryseobacterium4[,which(colSums(g_Chryseobacterium4) != 0)]
write.csv(g_Chryseobacterium5,"bacteria_g_Chryseobacterium2.csv")
counts <- as.matrix(t(g_Chryseobacterium5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Chryseobacterium_rho.csv")
ex.dna4 <- read.dna("bacteria_g_Chryseobacterium_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Chryseobacterium_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Chryseobacterium_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Chryseobacterium_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Chryseobacterium_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Chryseobacterium_rho.csv",header = T,row.names = 1)

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
          |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Chryseobacterium_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)

rho2<- read.csv("bacteria_g_Chryseobacterium_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p2<- ggplot(rho2, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.1, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Chryseobacterium") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.1, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p2, file='bacteria_g_Chryseobacterium_geomline.pdf', width=9, height=7)

###############bacteria_g_Dyadobacter######################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
g_Dyadobacter<-read.csv("bacteria_g_Dyadobacter.csv",header = T,row.names = 1)
g_Dyadobacter2<-select(g_Dyadobacter,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Dyadobacter2)<-g_Dyadobacter2[,1]
g_Dyadobacter3<-g_Dyadobacter2[,-1]
g_Dyadobacter4<-g_Dyadobacter3[which(rowSums(g_Dyadobacter3) !=0),]
g_Dyadobacter5<-g_Dyadobacter4[,which(colSums(g_Dyadobacter4) != 0)]
write.csv(g_Dyadobacter5,"bacteria_g_Dyadobacter2.csv")
counts <- as.matrix(t(g_Dyadobacter5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Dyadobacter_rho.csv")
ex.dna4 <- read.dna("bacteria_g_Dyadobacter_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Dyadobacter_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Dyadobacter_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Dyadobacter_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Dyadobacter_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Dyadobacter_rho.csv",header = T,row.names = 1)

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
          |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Dyadobacter_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho3<- read.csv("bacteria_g_Dyadobacter_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p3<- ggplot(rho3, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Dyadobacter") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.1, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p3, file='bacteria_g_Dyadobacter_geomline.pdf', width=9, height=7)
p3

######################Flavobacterium########################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Flavobacterium<-read.csv("bacteria_g_Flavobacterium.csv",header = T,row.names = 1)
g_Flavobacterium2<-select(g_Flavobacterium,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Flavobacterium2)<-g_Flavobacterium2[,1]
g_Flavobacterium3<-g_Flavobacterium2[,-1]
g_Flavobacterium4<-g_Flavobacterium3[which(rowSums(g_Flavobacterium3) !=0),]
g_Flavobacterium5<-g_Flavobacterium4[,which(colSums(g_Flavobacterium4) != 0)]
write.csv(g_Flavobacterium5,"bacteria_g_Flavobacterium2.csv")
counts <- as.matrix(t(g_Flavobacterium5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Flavobacterium_rho.csv")
all_fasta<- read.fasta('bacteria_g_Flavobacterium_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Flavobacterium5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Flavobacterium_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Flavobacterium_new2_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Flavobacterium_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Flavobacterium_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Flavobacterium_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Flavobacterium_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Flavobacterium_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Flavobacterium_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho4<- read.csv("bacteria_g_Flavobacterium_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p4<- ggplot(rho4, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Flavobacterium") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.5, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p4, file='bacteria_g_Flavobacterium_geomline.pdf', width=9, height=7)
p4

#####################bacteria_g_Luteolibacter########################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Luteolibacter<-read.csv("bacteria_g_Luteolibacter.csv",header = T,row.names = 1)
g_Luteolibacter2<-select(g_Luteolibacter,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Luteolibacter2)<-g_Luteolibacter2[,1]
g_Luteolibacter3<-g_Luteolibacter2[,-1]
g_Luteolibacter4<-g_Luteolibacter3[which(rowSums(g_Luteolibacter3) !=0),]
g_Luteolibacter5<-g_Luteolibacter4[,which(colSums(g_Luteolibacter4) != 0)]
write.csv(g_Luteolibacter5,"bacteria_g_Luteolibacter2.csv")
counts <- as.matrix(t(g_Luteolibacter5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Luteolibacter_rho.csv")
all_fasta<- read.fasta('bacteria_g_Luteolibacter_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Luteolibacter5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Luteolibacter_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Luteolibacter_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Luteolibacter_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Luteolibacter_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Luteolibacter_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Luteolibacter_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Luteolibacter_rho.csv",header = T,row.names = 1)
sum(is.na(res))
sum(is.na(rho))
res[is.na(res)] = 0
rho[is.na(rho)] = 0
sum(is.na(rho))

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
         ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
                           |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Luteolibacter_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho5<- read.csv("bacteria_g_Luteolibacter_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p5<- ggplot(rho5, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.8, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Luteolibacter") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.1, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p5, file='bacteria_g_Luteolibacter_geomline.pdf', width=9, height=7)
p5
###################g_Massilia################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Massilia<-read.csv("bacteria_g_Massilia.csv",header = T,row.names = 1)
g_Massilia2<-select(g_Massilia,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Massilia2)<-g_Massilia2[,1]
g_Massilia3<-g_Massilia2[,-1]
g_Massilia4<-g_Massilia3[which(rowSums(g_Massilia3) !=0),]
g_Massilia5<-g_Massilia4[,which(colSums(g_Massilia4) != 0)]
write.csv(g_Massilia5,"bacteria_g_Massilia2.csv")
counts <- as.matrix(t(g_Massilia5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Massilia_rho.csv")
all_fasta<- read.fasta('bacteria_g_Massilia_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Massilia5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Massilia_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Massilia_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Massilia_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Massilia_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Massilia_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Massilia_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Massilia_rho.csv",header = T,row.names = 1)

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Massilia_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho6<- read.csv("bacteria_g_Massilia_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p6<- ggplot(rho6, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Massilia") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.2, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p6, file='bacteria_g_Massilia_geomline.pdf', width=9, height=7)
p6
#############bacteria_g_Mucilaginibacter##################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Mucilaginibacter<-read.csv("bacteria_g_Mucilaginibacter.csv",header = T,row.names = 1)
g_Mucilaginibacter2<-select(g_Mucilaginibacter,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Mucilaginibacter2)<-g_Mucilaginibacter2[,1]
g_Mucilaginibacter3<-g_Mucilaginibacter2[,-1]
g_Mucilaginibacter4<-g_Mucilaginibacter3[which(rowSums(g_Mucilaginibacter3) !=0),]
g_Mucilaginibacter5<-g_Mucilaginibacter4[,which(colSums(g_Mucilaginibacter4) != 0)]
write.csv(g_Mucilaginibacter5,"bacteria_g_Mucilaginibacter2.csv")
counts <- as.matrix(t(g_Mucilaginibacter5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Mucilaginibacter_rho.csv")
all_fasta<- read.fasta('bacteria_g_Mucilaginibacter_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Mucilaginibacter5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Mucilaginibacter_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Mucilaginibacter_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Mucilaginibacter_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Mucilaginibacter_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Mucilaginibacter_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Mucilaginibacter_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Mucilaginibacter_rho.csv",header = T,row.names = 1)

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Mucilaginibacter_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho7<- read.csv("bacteria_g_Mucilaginibacter_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p7<- ggplot(rho7, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Mucilaginibacter") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.05, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p7, file='bacteria_g_Mucilaginibacter_geomline.pdf', width=9, height=7)
p7
######################Novosphingobium######################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Novosphingobium<-read.csv("bacteria_g_Novosphingobium.csv",header = T,row.names = 1)
g_Novosphingobium2<-select(g_Novosphingobium,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Novosphingobium2)<-g_Novosphingobium2[,1]
g_Novosphingobium3<-g_Novosphingobium2[,-1]
g_Novosphingobium4<-g_Novosphingobium3[which(rowSums(g_Novosphingobium3) !=0),]
g_Novosphingobium5<-g_Novosphingobium4[,which(colSums(g_Novosphingobium4) != 0)]
write.csv(g_Novosphingobium5,"bacteria_g_Novosphingobium2.csv")
counts <- as.matrix(t(g_Novosphingobium5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Novosphingobium_rho.csv")
all_fasta<- read.fasta('bacteria_g_Novosphingobium_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Novosphingobium5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Novosphingobium_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Novosphingobium_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Novosphingobium_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Novosphingobium_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Novosphingobium_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Novosphingobium_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Novosphingobium_rho.csv",header = T,row.names = 1)

merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
    }
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Novosphingobium_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho8<- read.csv("bacteria_g_Novosphingobium_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p8<- ggplot(rho8, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Novosphingobium") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.05, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p8, file='bacteria_g_Novosphingobium_geomline.pdf', width=9, height=7)
p8
##############Paraburkholderia####################

library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Paraburkholderia<-read.csv("bacteria_g_Paraburkholderia.csv",header = T,row.names = 1)
g_Paraburkholderia2<-select(g_Paraburkholderia,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Paraburkholderia2)<-g_Paraburkholderia2[,1]
g_Paraburkholderia3<-g_Paraburkholderia2[,-1]
g_Paraburkholderia4<-g_Paraburkholderia3[which(rowSums(g_Paraburkholderia3) !=0),]
g_Paraburkholderia5<-g_Paraburkholderia4[,which(colSums(g_Paraburkholderia4) != 0)]
write.csv(g_Paraburkholderia5,"bacteria_g_Paraburkholderia2.csv")
counts <- as.matrix(t(g_Paraburkholderia5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Paraburkholderia_rho.csv")
all_fasta<- read.fasta('bacteria_g_Paraburkholderia_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Paraburkholderia5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Paraburkholderia_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Paraburkholderia_new2_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Paraburkholderia_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}

ex_dna4_distance<- read.csv("bacteria_g_Paraburkholderia_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Paraburkholderia_dnadist2.csv")


####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Paraburkholderia_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Paraburkholderia_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Paraburkholderia_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho9<- read.csv("bacteria_g_Paraburkholderia_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p9<- ggplot(rho9, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Paraburkholderia") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.02, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p9, file='bacteria_g_Paraburkholderia_geomline.pdf', width=9, height=7)
p9
#################bacteria_g_Pedobacter######################
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Pedobacter<-read.csv("bacteria_g_Pedobacter.csv",header = T,row.names = 1)
g_Pedobacter2<-select(g_Pedobacter,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Pedobacter2)<-g_Pedobacter2[,1]
g_Pedobacter3<-g_Pedobacter2[,-1]
g_Pedobacter4<-g_Pedobacter3[which(rowSums(g_Pedobacter3) !=0),]
g_Pedobacter5<-g_Pedobacter4[,which(colSums(g_Pedobacter4) != 0)]
write.csv(g_Pedobacter5,"bacteria_g_Pedobacter2.csv")
counts <- as.matrix(t(g_Pedobacter5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Pedobacter_rho.csv")
all_fasta<- read.fasta('bacteria_g_Pedobacter_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Pedobacter5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Pedobacter_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Pedobacter_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Pedobacter_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Pedobacter_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Pedobacter_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Pedobacter_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Pedobacter_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Pedobacter_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho10<- read.csv("bacteria_g_Pedobacter_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p10<- ggplot(rho10, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.95, label.y.npc = 0.5,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Pedobacter") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.1, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p10, file='bacteria_g_Pedobacter_geomline.pdf', width=9, height=7)
p10
#####################bacteria_g_Pseudomonas##############
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Pseudomonas<-read.csv("bacteria_g_Pseudomonas.csv",header = T,row.names = 1)
g_Pseudomonas2<-select(g_Pseudomonas,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Pseudomonas2)<-g_Pseudomonas2[,1]
g_Pseudomonas3<-g_Pseudomonas2[,-1]
g_Pseudomonas4<-g_Pseudomonas3[which(rowSums(g_Pseudomonas3) !=0),]
g_Pseudomonas5<-g_Pseudomonas4[,which(colSums(g_Pseudomonas4) != 0)]
write.csv(g_Pseudomonas5,"bacteria_g_Pseudomonas2.csv")
counts <- as.matrix(t(g_Pseudomonas5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Pseudomonas_rho.csv")
all_fasta<- read.fasta('bacteria_g_Pseudomonas_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Pseudomonas5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Pseudomonas_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Pseudomonas_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Pseudomonas_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Pseudomonas_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Pseudomonas_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Pseudomonas_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Pseudomonas_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Pseudomonas_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho11<- read.csv("bacteria_g_Pseudomonas_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p11<- ggplot(rho11, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.95, label.y.npc = 0.9,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Pseudomonas") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.1, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p11, file='bacteria_g_Pseudomonas_geomline.pdf', width=9, height=7)
p11
############bacteria_g_Rhizobium#############
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Rhizobium<-read.csv("bacteria_g_Rhizobium.csv",header = T,row.names = 1)
g_Rhizobium2<-select(g_Rhizobium,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Rhizobium2)<-g_Rhizobium2[,1]
g_Rhizobium3<-g_Rhizobium2[,-1]
g_Rhizobium4<-g_Rhizobium3[which(rowSums(g_Rhizobium3) !=0),]
g_Rhizobium5<-g_Rhizobium4[,which(colSums(g_Rhizobium4) != 0)]
write.csv(g_Rhizobium5,"bacteria_g_Rhizobium2.csv")
counts <- as.matrix(t(g_Rhizobium5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Rhizobium_rho.csv")
all_fasta<- read.fasta('bacteria_g_Rhizobium_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Rhizobium5)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Rhizobium_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Rhizobium_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Rhizobium_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Rhizobium_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Rhizobium_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Rhizobium_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Rhizobium_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Rhizobium_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho12<- read.csv("bacteria_g_Rhizobium_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p12<- ggplot(rho12, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.95, label.y.npc = 0.95,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Rhizobium") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.05, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p12, file='bacteria_g_Rhizobium_geomline.pdf', width=9, height=7)
p12
###############bacteria_g_S0134#############
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_S0134<-read.csv("bacteria_g_S0134.csv",header = T,row.names = 1)
g_S01342<-select(g_S0134,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_S01342)<-g_S01342[,1]
g_S01343<-g_S01342[,-1]
g_S01344<-g_S01343[which(rowSums(g_S01343) !=0),]
g_S01345<-g_S01344[,which(colSums(g_S01344) != 0)]
write.csv(g_S01345,"bacteria_g_S01342.csv")
counts <- as.matrix(t(g_S01345))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_S0134_rho.csv")
all_fasta<- read.fasta('bacteria_g_S0134_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_S01345)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_S0134_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_S0134_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_S0134_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_S0134_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_S0134_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_S0134_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_S0134_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_S0134_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho13<- read.csv("bacteria_g_S0134_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p13<- ggplot(rho13, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.95, label.y.npc = 0.95,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("S0134") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.08, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p13, file='bacteria_g_S0134_geomline.pdf', width=9, height=7)
p13
#########Streptomyces############
library(ape)
library(ggplot2)
library(dplyr)
library(propr)
library(splines)
library(ggpmisc)
library(seqinr)
g_Streptomyces<-read.csv("bacteria_g_Streptomyces.csv",header = T,row.names = 1)
g_Streptomyces2<-select(g_Streptomyces,-c(Feature.ID,ASV_ID,d,p,c,o,f,g,s))
rownames(g_Streptomyces2)<-Streptomyces2[,1]
g_Streptomyces3<-g_Streptomyces2[,-1]
g_Streptomyces4<-g_Streptomyces3[which(rowSums(g_Streptomyces3) !=0),]
g_Streptomyces5<-g_Streptomyces4[,which(colSums(g_Streptomyces4) != 0)]
write.csv(g_Streptomyces5,"bacteria_g_Streptomyces2.csv")
counts <- as.matrix(t(g_Streptomyces5))
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
getMatrix(pr)
getResults(pr)
write.csv(getResults(pr),"bacteria_g_Streptomyces_rho.csv")
all_fasta<- read.fasta('bacteria_g_Streptomyces_new.fasta')
sub_fasta <- all_fasta[names(all_fasta) %in% rownames(g_Streptomyces)]
names(sub_fasta)
write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_g_Streptomyces_new2.fasta')
ex.dna4 <- read.dna("bacteria_g_Streptomyces_new_mafft.fasta", format = "fasta")
checkAlignment(ex.dna4) # width=14.02, height=8.63
ex.dna4.distance <- dist.dna(ex.dna4, model = "K80")
ex.dna4.distance #1-ex.dna4.distance
write.csv(as.matrix(ex.dna4.distance),"bacteria_g_Streptomyces_dnadist.csv")

###将dna_dist转化成两两对应关系#######

mtrx2cols = function(m1,m2,val1){
  lt = lower.tri(m1)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(m1,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(m1,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = m1[lt]) #按列依次获取矩阵下半角的元素 #对后两列重命名，支持多个矩阵合并
  return(res)
}
ex_dna4_distance<- read.csv("bacteria_g_Streptomyces_dnadist.csv",header = T,row.names = 1)
res = mtrx2cols(ex_dna4_distance,'val')
names(res)<- c("Partner","Pair","val1") 
write.csv(res,"bacteria_g_Streptomyces_dnadist2.csv")

####将矩阵对应关系填入rho文件###############
setwd("../Niche_differentiation/bacteria_niche")
res<- read.csv("bacteria_g_Streptomyces_dnadist2.csv",header = T,row.names = 1)
rho<- read.csv("bacteria_g_Streptomyces_rho.csv",header = T,row.names = 1)
merge_data= function(res,rho){
  data3=data.frame()
  for (i in 1:nrow(res)){
    for (j in 1:nrow(res)){
      if(!is.na(rho$Partner[i])&
         !is.na(res$Partner[j])&
         !is.na(res$Pair[j])&
         !is.na(rho$Pair[j])
      ){if (rho$Partner[i]==res$Partner[j]&rho$Pair[i]==res$Pair[j] 
            |rho$Partner[i]==res$Pair[j]&rho$Pair[i]==res$Partner[j])
      {data3 = rbind(data3,cbind(rho[i,],res[j,]))}
      }
      
    }
    
  }
  return(data3)
}

merge_data2<-merge_data(res,rho) 
merge_data2
write.csv(merge_data2,"bacteria_g_Streptomyces_rho_dnadist_merge.csv")
library(ggplot2)
library(splines)
library(ggpmisc)
#####删除大于distance0.5#########################
rho14<- read.csv("bacteria_g_Streptomyces_rho_dnadist_merge.csv",header = T,row.names = 1)
set.seed(1234)

p14<- ggplot(rho14, aes(as.numeric(val1), as.numeric(lrv))) + 
  geom_point(colour="black", size = 4) + 
  geom_smooth(method = "lm", formula =y ~ ns(x, 1) ,colour="black",lwd=1.5, lty=1)+
  stat_poly_eq(formula = y ~ ns(x, 1), 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~"),size=2.5),
               size = 8,label.x.npc = 0.9, label.y.npc = 0.95,
               parse = TRUE) +
  xlab("Nucleotide pairwise distances") + 
  ylab("Rho proportionality")+
  ggtitle("Streptomyces") +
  theme_set(theme_bw())+
  theme(plot.title=element_text(hjust=0.07, vjust=1,family='', face='bold', colour='black', size=20, margin=margin(t=40,b=-30)),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))
ggsave(p14, file='bacteria_g_Streptomyces_geomline.pdf', width=9, height=7)
p14
####图片组合#######https://zhuanlan.zhihu.com/p/384189537
p15 <- cowplot::plot_grid(p4, p9, p11, p12,p14, nrow = 3, labels = letters[1:5],label_size = 25)
p15
ggsave(p15, file='all_bacteria_g_remarkable_geomline.pdf', width=22, height=20)
p16 <- ggpubr::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, labels = c('a', 'b', 'c', 'd',"e"),size=20,font.label = list(color = 'black'))#将p1-p4四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D，颜色为红色(通过font.label = list()修改)，无法通过label.color = 'red'或其他方式修改。
p16
ggsave(p16, file='all_bacteria_g_remarkable_geomline.pdf', width=22, height=20)
?plot_grid
p17 <- cowplot::plot_grid(p1,p2,p3,p5,p6,p7,p8,p10,p13, nrow = 3,ncol=3,labels = letters[1:9],label_size = 25)
p17
ggsave(p17, file='all_bacteria_g_unremarkable_geomline.pdf', width=25, height=20)


