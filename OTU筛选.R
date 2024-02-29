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

###fungi_class
fungi_class<-read.csv("fungi_class_sum.csv",header=T,row.names = 1)
View(fungi_class)
colnames(fungi_class)
rownames(fungi_class)
ncol(fungi_class)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_class[ ,colnames(fungi_class) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_class_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_class_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_class_select<-data.frame(newdata4)

fungi_class_select_t<-t(fungi_class_select)

View(fungi_class_select)
nrow(fungi_class_select_t)
fungi_class_dist <- vegdist(fungi_class_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_class_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_class_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_class_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_class_nmds.pdf",width = 9,height = 6)

##fungi_species

fungi_species<-read.csv("fungi_species_sum.csv",header=T,row.names = 1)
View(fungi_species)
colnames(fungi_species)
rownames(fungi_species)
ncol(fungi_species)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_species[ ,colnames(fungi_species) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_species_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_species_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_species_select<-data.frame(newdata4)

fungi_species_select_t<-t(fungi_species_select)

View(fungi_species_select)
nrow(fungi_species_select_t)
fungi_species_dist <- vegdist(fungi_species_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_species_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_species_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_species_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_species_nmds.pdf",width = 9,height = 6)

##fungi_phylum
fungi_phylum<-read.csv("fungi_phylum_sum.csv",header=T,row.names = 1)
View(fungi_phylum)
colnames(fungi_phylum)
rownames(fungi_phylum)
ncol(fungi_phylum)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_phylum[ ,colnames(fungi_phylum) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_phylum_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_phylum_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_phylum_select<-data.frame(newdata4)

fungi_phylum_select_t<-t(fungi_phylum_select)

View(fungi_phylum_select)
nrow(fungi_phylum_select_t)
fungi_phylum_dist <- vegdist(fungi_phylum_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_phylum_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_phylum_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_phylum_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_phylum_nmds.pdf",width = 9,height = 6)


##fungi_family

fungi_family<-read.csv("fungi_family_sum.csv",header=T,row.names = 1)
View(fungi_family)
colnames(fungi_family)
rownames(fungi_family)
ncol(fungi_family)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_family[ ,colnames(fungi_family) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_family_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_family_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_family_select<-data.frame(newdata4)

fungi_family_select_t<-t(fungi_family_select)

View(fungi_family_select)
nrow(fungi_family_select_t)
fungi_family_dist <- vegdist(fungi_family_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_family_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_family_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_family_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_family_nmds.pdf",width = 9,height = 6)



##fungi_genus

fungi_genus<-read.csv("fungi_genus_sum.csv",header=T,row.names = 1)
View(fungi_genus)
colnames(fungi_genus)
rownames(fungi_genus)
ncol(fungi_genus)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_genus[ ,colnames(fungi_genus) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_genus_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_genus_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_genus_select<-data.frame(newdata4)

fungi_genus_select_t<-t(fungi_genus_select)

View(fungi_genus_select)
nrow(fungi_genus_select_t)
fungi_genus_dist <- vegdist(fungi_genus_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_genus_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_genus_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_genus_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_genus_nmds.pdf",width = 9,height = 6)

##fungi_order

fungi_order<-read.csv("fungi_order_sum.csv",header=T,row.names = 1)
View(fungi_order)
colnames(fungi_order)
rownames(fungi_order)
ncol(fungi_order)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_order[ ,colnames(fungi_order) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 

write.csv(newdata4,"fungi_order_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_order_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_order_select<-data.frame(newdata4)

fungi_order_select_t<-t(fungi_order_select)

View(fungi_order_select)
nrow(fungi_order_select_t)
fungi_order_dist <- vegdist(fungi_order_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_order_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_order_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_order_group2.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_order_nmds.pdf",width = 9,height = 6)


##fungi_asv
fungi_asv<-read.csv("fungi_new_feature_table_select.csv",header=T,row.names = 1)
View(fungi_asv)
colnames(fungi_asv)
rownames(fungi_asv)
ncol(fungi_asv)

new_d2 <- scan("fungi_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-fungi_asv[ ,colnames(fungi_asv) %in% new_d2]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
# remove rows 
df <- df[rowSums(df[-(1:7)]) !=0, ]
# remove columns    
df <- df[c(1:7,7 + which(colSums(df[-(1:7)]) !=0))]

write.csv(newdata4,"fungi_asv_select100_7000.csv")

write.csv(colnames(newdata4),"fungi_asv_group.csv")
rownames(newdata4)
colnames(newdata4)

ncol(newdata4)
nrow(newdata4)

fungi_asv_select<-data.frame(newdata4)

fungi_asv_select_t<-t(fungi_asv_select)

View(fungi_asv_select)
nrow(fungi_asv_select_t)
fungi_asv_dist <- vegdist(fungi_asv_select_t, "bray")

#计算mnds矩阵距离
#https://www.jianshu.com/p/c80a1cb9d5e9

nmds_dis <- metaMDS(fungi_asv_dist, k = 2)
#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191
#样方得分
nmds_dis_site <- data.frame(nmds_dis$points)
View(nmds_dis_site)
write.csv(nmds_dis_site,"fungi_asv_nmds_dis_sit.csv")
View(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
View(nmds_dis_site)
# nmds_dis_site<- read.csv("fungi_nmds_dis_sit_2.csv",header = T,row.names = 1)
map <- read.csv("fungi_asv_group.csv",header = T,row.names=1)
View(map)
#nmds_dis_site$group <- map$group
merged=merge(nmds_dis_site,map,by="row.names",all.x=TRUE)
View(merged)
color1=c("#00A087B2","#3C5488B2","#F39B7FB2")

p<-ggplot(merged, aes(MDS1, MDS2)) + 
  geom_point(aes(colour = factor(Fungi_GroupB),shape=factor(Fungi_GroupB)), size = 3) + 
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
ggsave (p,filename = "fungi_asv_nmds.pdf",width = 9,height = 6)


