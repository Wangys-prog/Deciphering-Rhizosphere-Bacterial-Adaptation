bacteria_asv_rhizo<- read.csv("./bacteria_asv_select_rhizo.csv",header = T,row.names = 1)
bacteria_asv_select_t <- t(bacteria_asv_rhizo)

set.seed(123)
select_train <- sample(120, 120*0.7)
ASV_train <- bacteria_asv_select_t[-select_train, ]
ASV_test <- bacteria_asv_select_t[select_train, ]
ASV_train.forest <- randomForest(as.factor(group4) ~ ., data = ASV_train, importance = TRUE) 
#训练集自身测试
train_predict <- predict(ASV_train.forest, ASV_train)
compare_train <- table(train_predict,  ASV_train$group4) 
sum(diag(compare_train)/sum(compare_train))  
compare_train
#使用测试集评估 
test_predict <- predict(ASV_train.forest, ASV_test) 
compare_test <- table(ASV_test$group4, test_predict, dnn = c('Actual', 'Predicted'))
compare_test
#或者使用函数 importance()
importance_ASV <- data.frame(importance(ASV_train.forest))
head(importance_ASV)

#可以根据某种重要性的高低排个序，例如根据“Mean Decrease Accuracy”指标
importance_ASV <- importance_ASV[order(importance_ASV$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_ASV)

#输出表格
write.table(importance_ASV, 'bacteria_importance_ASV.txt', sep = '\t', col.names = NA, quote = FALSE)

#读取表格

data=read.csv("OTU.csv",row.names = 1,header=T)

#ASVs with fewer than 10 reads were removed.

data[data < 10] <- 0 

#过滤方法

#install.packages('Hmisc')

library(Hmisc)

#处理的属水平，读取属水平genus物种丰度表格

genus <- read.csv('otu_genus.csv', row.names = 1,header=T)

genus <- data 

#完成由绝对丰度转化为相对丰度

A=genus  #行=样本 列=OTU

C=A/rowSums(A)

genus=t(C)

#可选事先过滤一些低丰度或低频的类群，行=OTU

#例如只保留相对丰度总和高于 0.01% 的属

genus <- genus[which(rowSums(genus) >= 0.0001), ] 

genus1 <- genus

#大于0的赋值为1

genus1[genus1>0] <- 1

#上一步大于0变成1后，下一步合算1多少个，保留下来

genus <- genus[which(rowSums(genus1) >=52 ), ]  

#例如只保留在 52个及以上样本中出现的属 按20%计算

write.table(genus,"otu_filter_0.0001.csv",sep=",",row.names=TRUE,col.names=TRUE) 



bacteria_asv<- read.csv("bacteria_asv_select_rhizo2.csv",header=T,row.names=1)

#去掉全为0的列
colnames(bacteria_asv)
ncol(bacteria_asv)
View(rowSums(bacteria_asv))
nrow(bacteria_asv)
bacteria_asv2<-bacteria_asv[which(rowSums(bacteria_asv) !=0),]
bacteria_asv3<-bacteria_asv2[,which(colSums(bacteria_asv2) != 0)]
View(bacteria_asv3)
A<-bacteria_asv3
B<-A/sum(colSums(A)) #计算col 列
C=t(B)
rowSums(C)
colSums(C)
View(C)
write.csv(C,"bacteria_asv_abundance.csv")
#例如只保留相对丰度总和高于 0.01% 的CSV
bacteria_asv_filter<- C[,which(colSums(C)>= 0.0001)]
ncol(bacteria_asv_filter)
nrow(bacteria_asv_filter)
write.csv(bacteria_asv_filter,"bacteria_asv_filter.csv")
bacteria_asv_filter1 <- bacteria_asv_filter
#大于0的赋值为1
bacteria_asv_filter1[bacteria_asv_filter1>0] <- 1
write.csv(bacteria_asv_filter1,"bacteria_asv_filter1.csv")
#上一步大于0变成1后，下一步合算1多少个，保留下来  97%
bacteria_asv_filter_broad<- bacteria_asv_filter[,which(colSums(bacteria_asv_filter1)>=4149)]
ncol(bacteria_asv_filter_broad)
nrow(bacteria_asv_filter)
#例如只保留在 52个及以上样本中出现的属 按20%计算
write.table(genus,"otu_filter_0.0001.csv",sep=",",row.names=TRUE,col.names=TRUE) 

###vegan####
bacteria_group4<- read.csv("bacteria_group4.csv",header = T,row.names = 1)

bacteria_filter_group4<-cbind(bacteria_asv_filter,bacteria_group4)
colnames(bacteria_filter_group4)
View(bacteria_filter_group4)
bacteria_filter_group4_sum<- t(aggregate(. ~ group4, data = bacteria_filter_group4,sum))
View(bacteria_filter_group4_sum)
colnames(bacteria_filter_group4_sum)<-bacteria_filter_group4_sum[1,]
bacteria_filter_group4_sum<- bacteria_filter_group4_sum[c(-1),]
ss <- t(bacteria_filter_group4_sum)
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
p1 <- ggvenn(a, show_elements = TRUE,fill_color = c("white", "white"),
             label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
             text_size = 3)
p1


p2 <- ggvenn(a, show_elements = FALSE,fill_color = c("green", "red"),
             label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
             text_size = 4)
p2
ggsave(p2, file='bacteria_core_venn.pdf', width=9, height=6) 


fill_color：表示填充颜色，默认值是c("blue","yellow", "green", "red")；


bacteria_core_abundance<-read.csv("./bacteria_core_abundance.csv",header = T,row.names = 1)
View(bacteria_core_abundance)
library(ggplot2)
p3<- ggplot(bacteria_core_abundance, aes(Relative.abundance...,Occurrence.frequency...,colour=Occurrence.frequency...)) + 
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
ggsave(p3, file='bacteria_core_point.pdf', width=7, height=8)

bacteria_taxa_rhizo2<- read.csv("./bacteria_taxa_rhizo2.csv",header = T,row.names = 1)
row.names(bacteria_taxa_rhizo2)<-bacteria_taxa_rhizo2[,1]
bacteria_taxa_rhizo2<-bacteria_taxa_rhizo2[,-1]
asv_occ15<- read.csv("clipboard")
View(asv_occ15)
asv_occ15_taxa<-bacteria_taxa_rhizo2[asv_occ15$ASV_ID,] 
bacteria_asv_select_rhizo2<- read.csv("bacteria_asv_select_rhizo2.csv",header = T,row.names = 1)
asv_occ15_table<-bacteria_asv_select_rhizo2[asv_occ15$ASV_ID,]
asv_occ15_table2<-asv_occ15_table[which(rowSums(asv_occ15_table) !=0),]
asv_occ15_table3<-asv_occ15_table2[,which(colSums(asv_occ15_table2) != 0)]
View(asv_occ15_table3)
ncol(asv_occ15_table3)
nrow(asv_occ15_table)
View(asv_occ15_table)
write.csv(asv_occ15_table,"bacteria_asv_occ15_table.csv")

c<- asv_occ15_taxa$c
o<- asv_occ15_taxa$o
f<- asv_occ15_taxa$f
g<- asv_occ15_taxa$g
asv_occ15_table_c<- cbind(asv_occ15_table3,c)
asv_occ15_table_o<- cbind(asv_occ15_table3,o)
asv_occ15_table_f<- cbind(asv_occ15_table3,f)
asv_occ15_table_g<- cbind(asv_occ15_table3,g)


occ15_c_sum<-aggregate(. ~c, data =asv_occ15_table_c,sum) #忽略参数不对提示
occ15_o_sum<-aggregate(. ~o, data =asv_occ15_table_o,sum)
occ15_f_sum<-aggregate(. ~f, data =asv_occ15_table_f,sum)
occ15_g_sum<-aggregate(. ~g, data =asv_occ15_table_g,sum)

occ15_g_sum[1,1]<-"unassigned"
rownames(occ15_g_sum)<-occ15_g_sum[,1]
occ15_g_sum1<-occ15_g_sum[,-1]
g_occ15_table_abundance<-t(t(occ15_g_sum1)/colSums(occ15_g_sum1))
View(g_occ15_table_abundance)

bacteria_group4<- read.csv("bacteria_group4.csv",header = T,row.names = 1)
group <-bacteria_group4[colnames(g_occ15_table_abundance),]
occ15_g_sum2<- data.frame(cbind(t(g_occ15_table_abundance),group))
View(occ15_g_sum2)
write.csv(occ15_g_sum2,"bacteria_occ15_g_sum2.csv")
occ15_g_sum2 <- read.csv("bacteria_occ15_g_sum2.csv",header = T,row.names = 1)
occ15_g_sum3 <- aggregate(. ~group, data = as.data.frame(occ15_g_sum2), mean)
View(occ15_g_sum3)
write.csv(occ15_g_sum3,"bacteria_occ15_g_sum3.csv")
##转置求相对丰度
####bar_plot########
library(ggplot2)
library(RColorBrewer)
library(hrbrthemes)
library(ggsci)
bacteria_core_barplot<- read.csv("bacteria_core_barplot_data.csv",header = T)
name2<- factor(bacteria_core_barplot$taxa_1,ordered=TRUE,levels=c("Class_Bulk","Class_Rhizosphere",
                                                   "Order_Bulk","Order_Rhizosphere",
                                                   "Family_Bulk","Family_Rhizosphere",
                                                   "Genus_Bulk","Genus_Rhizosphere"))

p_bar<-ggplot(data=bacteria_core_barplot,mapping=aes(x=Relative.abundance..., y=name2,fill=Taxa))+
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

ggsave(p_bar, file='bacteria_core_barplot.pdf', width=10, height=6)

######提取序列######

library(seqinr)
all_fasta <- read.fasta('bacteria_dna-sequences.fasta')

#fasta names######

names(all_fasta) 

####目标names######使用mgsub函数
target_id<- read.csv("ASV_ID.csv",header =0)

#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]

names(sub_fasta)

write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'bacteria_target_occ15.fasta')

#这一步把名字换成想要的样子
# seqkit replace --ignore-case --kv-file ASV_ID.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_target_occ15.fasta -o bacteria_target_occ15_new.fasta

##rename.txt 就是改名列表，第一列是原ID，第二列是新ID，中间用tab隔开。 genome.fa 是需要改ID的文件名，genome.new.fa 是新生成的改ID后的文件名。特别要注意的是列表中一定要包含所有的ID，不然他会将列表中不包含的ID改成空白
# picrust2_pipeline.py -s otu.fasta -i otu_table.txt -o picrust2_result -p 4

# https://www.jianshu.com/p/06d3a96c5717

# https://blog.csdn.net/woodcorpse/article/details/106554521

########真菌core##############################

fungi_asv<- read.csv("fungi_asv_select_rhizo2.csv",header=T,row.names=1)

#去掉全为0的列
colnames(fungi_asv)
ncol(fungi_asv)
View(rowSums(fungi_asv))
nrow(fungi_asv)
fungi_asv2<-fungi_asv[which(rowSums(fungi_asv) !=0),]
fungi_asv3<-fungi_asv2[,which(colSums(fungi_asv2) != 0)]
View(fungi_asv3)
A<-fungi_asv3
B<-A/sum(colSums(A)) #计算col 列
C=t(B)
rowSums(C)
colSums(C)
View(C)
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
View(fungi_filter_group4_sum)
colnames(fungi_filter_group4_sum)<-fungi_filter_group4_sum[1,]
fungi_filter_group4_sum<- fungi_filter_group4_sum[c(-1),]
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
a <- list(`Fungi_Bulk` = gene_set1,
          `Fungi_Rhizosphere` = gene_set2) #将基因集变成列表变量


p2 <- ggvenn(a, show_elements = FALSE,fill_color = c("green", "red"),
             label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
             text_size = 4)
p2
ggsave(p2, file='fungi_core_venn.pdf', width=9, height=6) 

######################################################

fungi_taxa_rhizo2<- read.csv("./fungi_taxa_rhizo2.csv",header = T,row.names = 1)
row.names(fungi_taxa_rhizo2)<-fungi_taxa_rhizo2[,1]
fungi_taxa_rhizo2<-fungi_taxa_rhizo2[,-1]
asv_occ10<- read.csv("clipboard")
asv_occ10<- read.csv("fungi_occ_10_id.csv",header =T)
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
##转置求相对丰度
####bar_plot########
library(ggplot2)
library(RColorBrewer)
library(hrbrthemes)
library(ggsci)
fungi_core_barplot<- read.csv("fungi_core_barplot_data2.csv",header = T)
name2<- factor(fungi_core_barplot$Taxa1,ordered=TRUE,levels=c("Class_Bulk","Class_Rhizosphere",
                                                                  "Order_Bulk","Order_Rhizosphere",
                                                                  "Family_Bulk","Family_Rhizosphere",
                                                                  "Genus_Bulk","Genus_Rhizosphere"))

p_bar<-ggplot(data=fungi_core_barplot,mapping=aes(x=Relative.abundance...., y=name2,fill=Taxa))+
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


######提取序列######

library(seqinr)
all_fasta <- read.fasta('fungi_dna-sequences.fasta')

#fasta names######

names(all_fasta) 

####目标names######使用mgsub函数
target_id<- read.csv("fungi_occ_10_id.csv",header =0)

#选取
sub_fasta <- all_fasta[names(all_fasta) %in% target_id[,1]]

names(sub_fasta)

write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'fungi_target_occ10.fasta')

#这一步把名字换成想要的样子
# seqkit replace --ignore-case --kv-file ASV_ID.txt --pattern "^(\w+)" --replacement "{kv}" bacteria_target_occ15.fasta -o bacteria_target_occ15_new.fasta

seqkit replace --ignore-case --kv-file fungi_occ_10_id.txt --pattern "^(\w+)" --replacement "{kv}" fungi_target_occ10.fasta -o fungi_target_occ10_new.fasta

##rename.txt 就是改名列表，第一列是原ID，第二列是新ID，中间用tab隔开。 genome.fa 是需要改ID的文件名，genome.new.fa 是新生成的改ID后的文件名。特别要注意的是列表中一定要包含所有的ID，不然他会将列表中不包含的ID改成空白
# picrust2_pipeline.py -s otu.fasta -i otu_table.txt -o picrust2_result -p 4

# https://www.jianshu.com/p/06d3a96c5717

# https://blog.csdn.net/woodcorpse/article/details/106554521

# bioawk -c fastx '{print}' bacteria_target_occ15_new.fasta| sort -k1,1V | awk '{print ">"$1;print $2}' >bacteria_target_occ15_new_sort.fasta

# picrust2_pipeline.py -s otu.fasta -i otu_table.txt -o picrust2_result -p 4

# picrust2_pipeline.py -s bacteria_target_occ15_new_sort.fasta -i bacteria_asv_occ15_table.txt -o bacteria_picrust2_result -p 4

# picrust2_pipeline.py -s fungi_target_occ10_new_sort.fasta -i fungi_asv_occ10_table.txt -o fungi_picrust2_result -p 4

######Taxa4Fun############################
library(vegan)
bacteria_ko_prediction<- read.csv("./bacteria_functional_prediction2.csv",header = T,row.names = 1)
bacteria_ko_prediction_t<- t(bacteria_ko_prediction[,-1])
bacteria_ko_prediction_t2<-bacteria_ko_prediction_t[which(rowSums(bacteria_ko_prediction_t) !=0),]
bacteria_ko_prediction_t3<-bacteria_ko_prediction_t2[,which(colSums(bacteria_ko_prediction_t2) != 0)]

names(bacteria_ko_prediction_t)
rownames(bacteria_ko_prediction_t)

bacteria_Ko_dist <- vegdist(bacteria_ko_prediction_t3, "bray")
nmds_dis <- metaMDS(bacteria_Ko_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_ko_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("bacteria_group4.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")
library(ggplot2)
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
####################
dune.ano <- with(merged, anosim(bacteria_Ko_dist, group4))
summary(dune.ano)


