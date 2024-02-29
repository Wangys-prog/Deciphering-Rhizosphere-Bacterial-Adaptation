library(vegan)
library(ape)
library(tidyverse)
library(ggplot2)

require(data.table)
data<-as.data.frame(fread("bacteria_feature_table100_7000.tsv"))
colnames(data)
rownames(data)
rownames(data)=data[,1]
rownames(data)
data[2,1]
##删去第一列
data<-data[,-1]
data[2,1]
write.csv(data,"bacteria_table100.csv")


###bacteria_class
bacteria_class<-read.csv("bacteria_class_sum.csv",header=T,row.names = 1)
View(bacteria_class)
colnames(bacteria_class)
rownames(bacteria_class)
ncol(bacteria_class)


new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_class[ ,colnames(bacteria_class) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)

write.csv(newdata2,"bacteria_class_select100_7000.csv")

#fungi_feature_table_sel<- read.csv("fungi_feature_table_select.csv",header = T,row.names = 1)
#fungi_feature_table_sel <- select(fungi_feature_table_sel,-c(1))
##取出第一列
#rownames(newdata2)=newdata2[,1]
##将第一列删除
#newdata2=newdata2[,-1]

rownames(newdata2)
colnames(newdata2)

ncol(newdata2)
nrow(newdata2)

bacteria_class_select<-data.frame(newdata2)

#删除小于sum2的行
# fungi_feature_table_sel$Part_sum <- rowSums(fungi_feature_table_sel[c(1:75303),])
# View(fungi_feature_table_sel$Part_sum)
# nrow(fungi_feature_table_sel)
# fungi_feature_table_sel2<- fungi_feature_table_sel[-which(fungi_feature_table_sel$Part_sum<2),]
# nrow(fungi_feature_table_sel2)
# ncol(fungi_feature_table_sel2)

# write.csv(fungi_feature_table_sel2,"fungi_feature_table_sel2.csv")

#删除part_sum列
# fungi_feature_table_sel3<-select(fungi_feature_table_sel,-c(Part_sum))

# nrow(fungi_feature_table_sel3)
# ncol(fungi_feature_table_sel3)

# write.csv(fungi_feature_table_sel3,"fungi_feature_table_sel3.csv")

# fungi_feature_table_sel3<- read.csv("fungi_feature_table_sel3.csv",header = T,row.names = 1)

bacteria_class_select_t<-t(bacteria_class_select)
# View(fungi_feature_table_sel3_t)
nrow(bacteria_class_select_t)
bacteria_class_dist <- vegdist(bacteria_class_select_t, "bray")
bacteria_map<-  read.csv("bacteria_group2.csv",header = T,row.names=1)
# res <- pcoa(feature_dist)
# head(feature_dist)
# View(feature_dist)

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(bacteria_class_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_class_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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

p <- ggplot(data = merged, aes(MDS1, MDS2)) +
  geom_point(size=3,aes(colour = factor(merged$group1)),shape = group1)) +
  stat_ellipse(aes(fill = group1), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +   #添加置信椭圆，注意不是聚类
  scale_color_manual(values =color[1:length(unique(map$group1))]) +
  scale_fill_manual(values = color[1:length(unique(map$group1))]) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
  #, legend.position = 'none'
  geom_vline(xintercept = 0, color = 'gray', size = 1) +
  geom_hline(yintercept = 0, color = 'gray', size = 1)+
  #geom_text(data =merged, aes(label = Row.names,x =MDS1, y = MDS2), size=4, check_overlap = TRUE)+
  xlab("MDS1") + ylab("MDS2")+
  theme_bw()+ 
  theme(rect = element_rect(colour = "black", size = 2, linetype = 1),
        panel.grid.major=element_line(colour=NA),
        text=element_text(face = "bold",colour="black",size = 20))
p

##bacteria_phylum
bacteria_phylum<-read.csv("bacteria_phylum_sum.csv",header=T,row.names = 1)
View(bacteria_phylum)
colnames(bacteria_phylum)
rownames(bacteria_phylum)
ncol(bacteria_phylum)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_phylum[ ,colnames(bacteria_phylum) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_phylum_select100_7000.csv")
bacteria_phylum_select<-data.frame(newdata2)
bacteria_phylum_select_t<-t(bacteria_phylum_select)
nrow(bacteria_phylum_select_t)
bacteria_phylum_dist <- vegdist(bacteria_phylum_select_t, "bray")
# res <- pcoa(bacteria_phylum_dist)
nmds_dis <- metaMDS(bacteria_phylum_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_phylum_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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


##bacteria_family
bacteria_family<-read.csv("bacteria_family_sum.csv",header=T,row.names = 1)
View(bacteria_family)
colnames(bacteria_family)
rownames(bacteria_family)
ncol(bacteria_family)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_family[ ,colnames(bacteria_family) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_family_select100_7000.csv")
bacteria_family_select<-data.frame(newdata2)
bacteria_family_select_t<-t(bacteria_family_select)
nrow(bacteria_family_select_t)
bacteria_family_dist <- vegdist(bacteria_family_select_t, "bray")
# res <- pcoa(bacteria_family_dist)
nmds_dis <- metaMDS(bacteria_family_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_family_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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

##bacteria_genus
bacteria_genus<-read.csv("bacteria_genus_sum.csv",header=T,row.names = 1)
View(bacteria_genus)
colnames(bacteria_genus)
rownames(bacteria_genus)
ncol(bacteria_genus)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_genus[ ,colnames(bacteria_genus) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_genus_select100_7000.csv")
bacteria_genus_select<-data.frame(newdata2)
bacteria_genus_select_t<-t(bacteria_genus_select)
nrow(bacteria_genus_select_t)
bacteria_genus_dist <- vegdist(bacteria_genus_select_t, "bray")
# res <- pcoa(bacteria_genus_dist)
nmds_dis <- metaMDS(bacteria_genus_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_genus_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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

##bacteria_order
bacteria_order<-read.csv("bacteria_order_sum.csv",header=T,row.names = 1)
View(bacteria_order)
colnames(bacteria_order)
rownames(bacteria_order)
ncol(bacteria_order)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_order[ ,colnames(bacteria_order) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_order_select100_7000.csv")
bacteria_order_select<-data.frame(newdata2)
bacteria_order_select_t<-t(bacteria_order_select)
nrow(bacteria_order_select_t)
bacteria_order_dist <- vegdist(bacteria_order_select_t, "bray")
# res <- pcoa(bacteria_genus_dist)
nmds_dis <- metaMDS(bacteria_order_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_order_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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

##bacteria_species
bacteria_species<-read.csv("bacteria_species_sum.csv",header=T,row.names = 1)
View(bacteria_species)
colnames(bacteria_species)
rownames(bacteria_species)
ncol(bacteria_species)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_species[ ,colnames(bacteria_species) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_species_select100_7000.csv")
bacteria_species_select<-data.frame(newdata2)
bacteria_species_select_t<-t(bacteria_species_select)
nrow(bacteria_species_select_t)
bacteria_species_dist <- vegdist(bacteria_species_select_t, "bray")
# res <- pcoa(bacteria_genus_dist)
nmds_dis <- metaMDS(bacteria_species_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_species_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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

##bacteria_asv
bacteria_asv<-read.csv("bacteria_table100.csv",header=T,row.names = 1)
View(bacteria_asv)
colnames(bacteria_asv)
rownames(bacteria_asv)
ncol(bacteria_asv)
new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria_asv[ ,colnames(bacteria_asv) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
write.csv(newdata2,"bacteria_asv_select100_7000.csv")
bacteria_asv_select<-data.frame(newdata2)
bacteria_asv_select_t<-t(bacteria_asv_select)
nrow(bacteria_asv_select_t)
bacteria_asv_dist <- vegdist(bacteria_asv_select_t, "bray")
# res <- pcoa(bacteria_genus_dist)
nmds_dis <- metaMDS(bacteria_asv_dist, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"bacteria_asv_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
map <- read.csv("bacteria_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
write.csv(merged,"fungi_class_merged.csv")
View(merged)
p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(group1),shape=factor(group1)), size = 3) + 
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
ggsave (p,filename = "bacteria_asv_nmds.pdf",width = 9,height = 6)

# R ASV绝对值转化成相对丰度，即每一列除以对应的列之和
bacteria_asv_abo<-read.csv('bacteria_asv_select100_7000.csv',header = T,row.names=1)
bacteria_asv_abo<-t(bacteria_asv_abo)
bacteria_asv_abo<-bacteria_asv_abo/rowSums(bacteria_asv_abo)


#species_network
bacteria_species_abo<-read.csv('bacteria_species_select100_7000.csv',header = T,row.names=1)
bacteria_species_abo<-t(bacteria_species_abo)

bacteria_species_abundance<-bacteria_species_abo/rowSums(bacteria_species_abo)
ncol(bacteria_species_abundance)
nrow(bacteria_species_abundance)
# bacteria_species_abundance<-t(bacteria_species_abundance)
write.csv(bacteria_species_abundance,file='bacteria_species_abundance.csv')
ncol(bacteria_species_abundance)
nrow(bacteria_species_abundance)
newdata3<-bacteria_species_abundance[which(rowSums(bacteria_species_abundance)>0),]
ncol(newdata3)
nrow(newdata3)
bacteria_species_abundance2<-newdata3[,which(colSums(newdata3)>0)]
ncol(bacteria_species_abundance2)
nrow(bacteria_species_abundance2)

write.csv(bacteria_species_abundance2,file='bacteria_species_abundance2.csv')

# #读取属水平的物种丰度表，

# genus <- genus[which(rowSums(genus) >= 0.005), ]    #只保留相对丰度总和高于 0.005 的属

#例如只保留在 5 个及以上样本中出现的属
# genus1 <- genus
# genus1[genus1>0] <- 1
# genus <- genus[which(rowSums(genus1) >= 5), ]

##genus_network

bacteria_genus_abo<-read.csv('bacteria_genus_select100_7000.csv',header = T,row.names=1)
bacteria_genus_abo<-t(bacteria_genus_abo)

bacteria_genus_abundance<-bacteria_genus_abo/rowSums(bacteria_genus_abo)
ncol(bacteria_genus_abundance)
nrow(bacteria_genus_abundance)
# bacteria_genus_abundance<-t(bacteria_genus_abundance)
write.csv(bacteria_genus_abundance,file='bacteria_genus_abundance.csv')
ncol(bacteria_genus_abundance)
nrow(bacteria_genus_abundance)
newdata3<-bacteria_genus_abundance[which(rowSums(bacteria_genus_abundance)>0),]
ncol(newdata3)
nrow(newdata3)
bacteria_genus_abundance2<-newdata3[,which(colSums(newdata3)>0)]
ncol(bacteria_genus_abundance2)
nrow(bacteria_genus_abundance2)

write.csv(bacteria_genus_abundance2,file='bacteria_genus_abundance2.csv')

##phylum_network

bacteria_phylum_abo<-read.csv('bacteria_phylum_select100_7000.csv',header = T,row.names=1)
bacteria_phylum_abo<-t(bacteria_phylum_abo)
bacteria_phylum_abundance<-bacteria_phylum_abo/rowSums(bacteria_phylum_abo)
ncol(bacteria_phylum_abundance)
nrow(bacteria_phylum_abundance)
# bacteria_phylum_abundance<-t(bacteria_phylum_abundance)
write.csv(bacteria_phylum_abundance,file='bacteria_phylum_abundance.csv')
ncol(bacteria_phylum_abundance)
nrow(bacteria_phylum_abundance)
newdata3<-bacteria_phylum_abundance[which(rowSums(bacteria_phylum_abundance)>0),]
ncol(newdata3)
nrow(newdata3)
bacteria_phylum_abundance2<-newdata3[,which(colSums(newdata3)>0)]
ncol(bacteria_phylum_abundance2)
nrow(bacteria_phylum_abundance2)
write.csv(bacteria_phylum_abundance2,file='bacteria_phylum_abundance2.csv')
 
                             
##asv_network
bacteria_asv_abo<-read.csv('bacteria_asv_select100_7000.csv',header = T,row.names=1)
bacteria_asv_abo<-t(bacteria_asv_abo)
bacteria_asv_abundance<-bacteria_asv_abo/rowSums(bacteria_asv_abo)
ncol(bacteria_asv_abundance)
nrow(bacteria_asv_abundance)
# bacteria_asv_abundance<-t(bacteria_asv_abundance)
write.csv(bacteria_asv_abundance,file='bacteria_asv_abundance.csv')
ncol(bacteria_asv_abundance)
nrow(bacteria_asv_abundance)
newdata3<-bacteria_asv_abundance[which(rowSums(bacteria_asv_abundance)>0),]
ncol(newdata3)
nrow(newdata3)
bacteria_asv_abundance2<-newdata3[,which(colSums(newdata3)>0)]
ncol(bacteria_asv_abundance2)
nrow(bacteria_asv_abundance2)
write.csv(bacteria_asv_abundance2,file='bacteria_asv_abundance2.csv')

rownames(bacteria_asv_abundance2)
rownames(bacteria_asv_abundance2)
colnames(bacteria_asv_abundance2)
write.csv(rownames(bacteria_asv_abundance2),"bacteria_asv_taxa.csv")
write.csv(colnames(bacteria_asv_abundance2),"bacteria_asv_sample.csv")

bacteria_asv_abundance3<- read.csv("bacteria_asv_abundance2.csv")

