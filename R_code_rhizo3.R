##########core_ASV################
# https://zhuanlan.zhihu.com/p/547246100
setwd("../core_asvs/")
library(tidyverse)
library(rstatix)
bacteria_asv_occ15<- read.csv("bacteria_asv_occ15_table.csv")
dir()
bacteria_group4<- read.csv("bacteria_group4.csv")
# 表达量+分组 合并，整理
df = bacteria_asv_occ15 %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(ASV_ID,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(ASV_ID) %>%
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
write.csv(dfdata,"bacteria_differential_core_asvs.csv")
###绘制火山图#####
library(ggplot2)
dataset<- read.csv("bacteria_differential_core_asvs.csv",header = T,row.names = 1)
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
ggsave(p, file='bacteria_differential_core_asv_volcano.pdf', width=8, height=6)
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
  geom_label_repel(aes(label =ASV_ID),data=dataset, nudge_y = 2,  alpha = 0.6  )+
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
ggsave(p, file='bacteria_differential_core_asv_volcano2.pdf', width=8, height=6)

###################功能预测差异分析##########################
setwd("../16s_function_prediction/bacteria_core_blast_Ref100NR/")

bacteria_function<- read.csv("bacteria_functional_prediction2.csv")
bacteria_group4<- read.csv("bacteria_group4.csv")
# 表达量+分组 合并，整理
df = bacteria_function %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(KO,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(KO) %>%
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
write.csv(dfdata,"bacteria_differential_KO.csv")

###绘制火山图#####
library(ggplot2)
dataset<- read.csv("bacteria_differential_KO.csv",header = T,row.names = 1)
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dataset$change = ifelse(dataset$FDR < cut_off_pvalue & abs(dataset$FC) >= cut_off_logFC, 
                        ifelse(dataset$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)
write.csv(dataset,"bacteria_differential_KO.csv")

p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("blue","red","gray"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
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
ggsave(p, file='bacteria_differential_ko_volcano.pdf', width=8, height=6)
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
  geom_label_repel(aes(label =KO),data=dataset, nudge_y = 2,  alpha = 0.6  )+
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
ggsave(p, file='bacteria_differential_core_asv_volcano2.pdf', width=8, height=6)

################################################################
bacteria_differential_KO_up <- read.csv("bacteria_functional_prediction_up.csv",header = T,row.names = 1)
bacteria_differential_KO_up3<-bacteria_differential_KO_up[,-1]
bacteria_differential_KO_up4<-bacteria_differential_KO_up3[which(rowSums(bacteria_differential_KO_up3) !=0),]
bacteria_differential_KO_up5<-bacteria_differential_KO_up4[,which(colSums(bacteria_differential_KO_up4) != 0)]
ncol(bacteria_differential_KO_up)
ncol(bacteria_differential_KO_up5)
nrow(bacteria_differential_KO_up)
nrow(bacteria_differential_KO_up5)

write.csv(bacteria_differential_KO_up5,"bacteria_differential_KO_up5.csv")
bacteria_differential_KO_up5$description<- bacteria_differential_KO_up[,1]
bacteria_group4_2<- read.csv("bacteria_group4_2.csv")
df = bacteria_differential_KO_up5 %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4_2,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(description,group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列
df
write.csv(df,"bacteria_differential_KO_up5_sum.csv")

##bacteria_differential_KO2_sum <- aggregate(. ~description, data = bacteria_differential_KO_up5,sum)

####根据group求总和######

library(dplyr)
dft_c <- dft %>% 
  group_by(group) %>% 
  arrange(group,dat) %>% 
  mutate(cum_num = cumsum(num)) %>% 
  ungroup()
dft_c

##########################
##
setwd("../16s_function_prediction/bacteria_core_blast_Ref100NR/")
bacteria_ko_top<- read.csv("bacteria_differential_KO4.csv",header = T)

p<- ggplot(bacteria_ko_top, aes(x=description2, y=value)) + 
  geom_bar(stat = "identity",aes(fill =group),position = position_stack(reverse = TRUE))+
  labs(x="",
     y="Count")+
  theme_set(theme_bw())+
  theme(legend.position="top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p1<- p+coord_flip()

ggsave(p1, file='bacteria_differential_ko_geombar.pdf', width=15, height=7)


##############pathway########################
#############################################
library(dplyr)
library(tidyverse)
library(rstatix)
bacteria_pathway<- read.csv("bacteria_pathway_prediction2.csv")
bacteria_group4<- read.csv("bacteria_group4.csv")
# 表达量+分组 合并，整理
df = bacteria_pathway %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(pathway,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk,.groups = "drop_last")  # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(pathway) %>%
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
write.csv(dfdata,"bacteria_differential_pathway.csv")

###绘制火山图#####
library(ggplot2)
dataset<- read.csv("bacteria_differential_pathway.csv",header = T,row.names = 1)
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dataset$change = ifelse(dataset$FDR < cut_off_pvalue & abs(dataset$FC) >= cut_off_logFC, 
                        ifelse(dataset$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)
write.csv(dataset,"bacteria_differential_pathway.csv")

p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","red","gray"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
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
ggsave(p, file='bacteria_differential_pathway_volcano.pdf', width=8, height=6)
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
  geom_label_repel(aes(label =pathway),data=dataset, nudge_y = 2,  alpha = 0.6  )+
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
ggsave(p, file='bacteria_differential_core_pathway_volcano2.pdf', width=8, height=6)

################################################################
bacteria_differential_KO_up <- read.csv("bacteria_functional_prediction_up.csv",header = T,row.names = 1)
bacteria_differential_KO_up3<-bacteria_differential_KO_up[,-1]
bacteria_differential_KO_up4<-bacteria_differential_KO_up3[which(rowSums(bacteria_differential_KO_up3) !=0),]
bacteria_differential_KO_up5<-bacteria_differential_KO_up4[,which(colSums(bacteria_differential_KO_up4) != 0)]
ncol(bacteria_differential_KO_up)
ncol(bacteria_differential_KO_up5)
nrow(bacteria_differential_KO_up)
nrow(bacteria_differential_KO_up5)

write.csv(bacteria_differential_KO_up5,"bacteria_differential_KO_up5.csv")
bacteria_differential_KO_up5$description<- bacteria_differential_KO_up[,1]
bacteria_group4_2<- read.csv("bacteria_group4_2.csv")
df = bacteria_differential_KO_up5 %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4_2,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(description,group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列
df
write.csv(df,"bacteria_differential_KO_up5_sum.csv")

##bacteria_differential_KO2_sum <- aggregate(. ~description, data = bacteria_differential_KO_up5,sum)

####根据group求总和######

library(dplyr)
dft_c <- dft %>% 
  group_by(group) %>% 
  arrange(group,dat) %>% 
  mutate(cum_num = cumsum(num)) %>% 
  ungroup()
dft_c

##########################
##
setwd("../16s_function_prediction/bacteria_core_blast_Ref100NR/")
bacteria_ko_top<- read.csv("bacteria_differential_KO4.csv",header = T)

p<- ggplot(bacteria_ko_top, aes(x=description2, y=value)) + 
  geom_bar(stat = "identity",aes(fill =group),position = position_stack(reverse = TRUE))+
  labs(x="",
       y="Count")+
  theme_set(theme_bw())+
  theme(legend.position="top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p1<- p+coord_flip()

ggsave(p1, file='bacteria_differential_ko_geombar.pdf', width=15, height=7)

###############功能冗余指数#######################
####https://www.sci666.com.cn/57925.html#######
bacteria_aFRI <- read.table('absolute_functional_redundancy.txt',sep = '\t',header = TRUE,row.names = 1)
#names(bacteria_aFRI)[4150]
#bacteria_aFRI2 <- select(bacteria_aFRI,-description)
# ncol(bacteria_aFRI2 )
####计算分组平均值
df = bacteria_aFRI %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(row.names(bacteria_aFRI),group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列
df

write.csv(df,"bacteria_aFRI_avg.csv")

#####差异性aFRI########
bacteria_aFRI <- read.csv('absolute_functional_redundancy.csv')
bacteria_aFRI2 <- select(bacteria_aFRI,-1)
bacteria_group4<- read.csv("bacteria_group4.csv")
df = bacteria_aFRI %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(KO,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk,.groups = "drop_last")  # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(KO) %>%
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

cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dfdata$change = ifelse(dfdata$FDR < cut_off_pvalue & abs(dfdata$FC) >= cut_off_logFC, 
                        ifelse(dfdata$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
write.csv(dfdata,"bacteria_differential_aFRI.csv")

####计算分组平均值
bacteria_aFRI <- read.table('absolute_functional_redundancy.txt',sep = '\t',header = TRUE,row.names = 1)
bacteria_aFRI2 <- select(bacteria_aFRI,-description)

df = bacteria_aFRI %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(description,group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列
df
write.csv(df,"bacteria_aFRI_avg.csv")

df = bacteria_aFRI %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(description,group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列

#############################################
###############功能冗余指数#######################
####https://www.sci666.com.cn/57925.html#######
setwd("../16s_function_prediction/bacteria_core_blast_Ref100NR/")
bacteria_aFRI_avg2<- read.csv("bacteria_aFRI_avg2.csv",header = T,row.names=1)
#添加分组信息#####
# x <- ifelse(str_detect(group$title,"Control"),"control","treat")
#group$group <- x
bacteria_aFRI_avg2[which(bacteria_aFRI_avg2$log <= 0),"group"] <- "red"
bacteria_aFRI_avg2[which(bacteria_aFRI_avg2$log > 0),"group"] <- "blue"
View(bacteria_aFRI_avg2)
library(ggplot2)
p<- ggplot(bacteria_aFRI_avg2, aes(x=num, y=log,color = group)) + 
  geom_point(size=2)+
  scale_color_manual(values = c("#f88421","#006b7b")) +
  labs(x="",
       y="log ratio(FRIRhizosphere/FRIBulk)")+
  theme_set(theme_bw())+
  geom_hline(yintercept=c(0),lty=4,col="black",lwd=1.5) +
  theme(legend.position="top",
        legend.text=element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",linewidth=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

ggsave(p, file='bacteria_differential_FRI_log.pdf', width=8, height=7)


#####################################################
library(dplyr)
library(tidyr)
library(rstatix)
bacteria_function<- read.csv("bacteria_functional_prediction3.csv")
bacteria_group4<- read.csv("bacteria_group4.csv")
# 表达量+分组 合并，整理
df = bacteria_function %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(description,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(description) %>%
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
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dfdata$change = ifelse(dfdata$FDR < cut_off_pvalue & abs(dfdata$FC) >= cut_off_logFC, 
                        ifelse(dfdata$FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dfdata)
write.csv(dfdata,"bacteria_differential_description.csv")

###绘制火山图#####
library(ggplot2)
dataset<- read.csv("bacteria_differential_description.csv",header = T,row.names = 1)

p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("blue","red","gray"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
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
ggsave(p, file='bacteria_differential_description_volcano.pdf', width=8, height=6)

library(ggplot2)
dataset<- read.csv("bacteria_differential_description2.csv",header = T,row.names = 1)
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
  geom_label_repel(aes(label =description),data=dataset, nudge_y = 2,  alpha = 0.6)+
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
ggsave(p, file='bacteria_differential_core_pathway_volcano2.pdf', width=8, height=6)



##########
#####根据上调description提取列表####
description<- read.csv("bacteria_differential_description2.csv",header = T,row.names = 1)
a<- description$description
bacteria_function$description
###筛选######
bacteria_function_up<-bacteria_function %>% filter(description %in% a)
nrow(description)
nrow(bacteria_function_up)
View(bacteria_function_up)
write.csv(bacteria_function_up,"bacteria_function_description_table_up.csv")
bacteria_differential_description_table_up<- read.csv("bacteria_function_description_table_up.csv",header = T,row.names = 1)
nrow(bacteria_differential_description_table_up)
bacteria_differential_description_table_up3<-bacteria_differential_description_table_up[,-1]
bacteria_differential_description_table_up4<-bacteria_differential_description_table_up3[which(rowSums(bacteria_differential_description_table_up3) !=0),]
bacteria_differential_description_table_up5<-bacteria_differential_description_table_up4[,which(colSums(bacteria_differential_description_table_up4) != 0)]
ncol(bacteria_differential_description_table_up5)
nrow(bacteria_differential_description_table_up)
nrow(bacteria_differential_description_table_up5)
write.csv(bacteria_differential_description_table_up5,"bacteria_differential_description_table_up5.csv")
bacteria_differential_description_table_up5$description<- bacteria_differential_description_table_up[,1]
b<-names(bacteria_differential_description_table_up5)
bacteria_group4<- read.csv("bacteria_group4.csv")
names(bacteria_group4)
bacteria_group4$sample
bacteria_group5<- bacteria_group4 %>% filter(sample %in% b)
View(bacteria_group5)
names(bacteria_group5)

df = bacteria_differential_description_table_up5 %>%
  as_tibble() %>%
  pivot_longer(-description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(bacteria_group5,by=c("sample" = "sample"))%>%# 与分组数据合并
  group_by(description,group4)%>%
  summarise(value_mean = mean(value),.groups = "drop")%>% #分组计算平均值
  arrange(desc(value_mean)) #倒序排列
df
write.csv(df,"bacteria_differential_description_table_up5_avg.csv")

################################
################################
bacteria_differential_description2<- read.csv("bacteria_differential_description2.csv",header = T,row.names = 1)
c<-bacteria_differential_description2$description[c(1:30)]
c
bacteria_differential_description_table_up5_avg<- read.csv("bacteria_differential_description_table_up5_avg.csv",header = T,row.names = 1)
bacteria_differential_description_table_up5_avg$description

bacteria_differential_description_table_up5_avg2<- bacteria_differential_description_table_up5_avg %>% filter(description %in% c)
bacteria_differential_description_table_up5_avg2
write.csv(bacteria_differential_description_table_up5_avg2,"bacteria_differential_description_table_up5_avg2.csv")

bacteria_differential_description3<- merge(bacteria_differential_description2,bacteria_differential_description_table_up5_avg2, by = "description")
View(bacteria_differential_description3)
write.csv(bacteria_differential_description3,"bacteria_differential_description3.csv")
p<- ggplot(bacteria_differential_description3, aes(x=description, y=value_mean)) + 
  geom_bar(stat = "identity",aes(fill =group4),position = position_stack(reverse = TRUE))+
  labs(x="",
       y="Count")+
  theme_set(theme_bw())+
  theme(legend.position="top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))

p1<- p+coord_flip()
p1
ggsave(p1, file='bacteria_differential_description_geombar.pdf', width=15, height=7)


#######################fungi_function#######################
library(dplyr)
library(tidyr)
library(rstatix)
setwd("../18s_ITS_function_prediction/picrust2Dypu577a4kKTvhWdbhUkFo_result")
fungi_function<- read.csv("Sample_EC_abundance3.csv")
fungi_group4<- read.csv("fungi_group4.csv")
# 表达量+分组 合并，整理
df = fungi_function %>%
  as_tibble() %>%
  pivot_longer(-Description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(fungi_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(Description,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(Description) %>%
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
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dfdata$change = ifelse(dfdata$FDR < cut_off_pvalue & abs(dfdata$FC) >= cut_off_logFC, 
                       ifelse(dfdata$FC> cut_off_logFC ,'Up','Down'),
                       'Stable')
head(dfdata)
write.csv(dfdata,"fungi_differential_EC_description.csv")

###绘制火山图#####
library(ggplot2)
dataset<- read.csv("fungi_differential_EC_description.csv",header = T,row.names = 1)

p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("blue","red","gray"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",linewidth=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size =15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))
p
ggsave(p, file='fungi_differential_EC_description.csv_volcano.pdf', width=8, height=6)

#################################

setwd("../18s_ITS_function_prediction/picrust2Dypu577a4kKTvhWdbhUkFo_result")
fungi_function<- read.csv("MetaCyc_pathway_abundance2.csv")
fungi_group4<- read.csv("fungi_group4.csv")
# 表达量+分组 合并，整理
df = fungi_function %>%
  as_tibble() %>%
  pivot_longer(-Description,names_to = "sample",values_to = "value") %>%  # 转换成长数据
  left_join(fungi_group4,by=c("sample" = "sample")) # 与分组数据合并
df
#计算FC
dfFC = df %>%
  group_by(Description,group4) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = group4,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = Rhizosphere/Bulk)                           # 实验组/对照组 计算差异倍数FC
dfFC
# t_test 计算P 值
dfP = df %>%
  group_by(Description) %>%
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
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 1 
dfdata$change = ifelse(dfdata$FDR < cut_off_pvalue & abs(dfdata$FC) >= cut_off_logFC, 
                       ifelse(dfdata$FC> cut_off_logFC ,'Up','Down'),
                       'Stable')
head(dfdata)
write.csv(dfdata,"fungi_differential_pathway_description.csv")

###绘制火山图#####
library(ggplot2)
dataset<- read.csv("fungi_differential_pathway_description.csv",header = T,row.names = 1)

p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = FC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("blue","red","gray"))+
  #横坐标刻度修改
  scale_x_continuous(limits = c(-2.5,10))+
  # 辅助线
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",linewidth=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size =15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))
p
ggsave(p, file='fungi_differential_pathway_description.csv_volcano.pdf', width=8, height=6)

#####ggclusternet#################
setwd("../core_network")
library(dplyr)
bacteria_core_otu<- read.csv("bacteria_asv_occ15_table2.csv",row.names = 1,header = T)
fungi_core_otu<- read.csv("fungi_asv_occ10_table2.csv",row.names = 1,header = T)
View(bacteria_core_otu)
View(fungi_core_otu)
a<- names(fungi_core_otu)
bacteria_core_otu_new<- bacteria_core_otu[ ,colnames(bacteria_core_otu) %in% a]
View(bacteria_core_otu_new)
ncol(bacteria_core_otu_new)
nrow(bacteria_core_otu_new)
ncol(fungi_core_otu)
nrow(fungi_core_otu)
b<-names(bacteria_core_otu_new)
fungi_core_otu_new<- fungi_core_otu[ ,colnames(fungi_core_otu) %in% b]
ncol(fungi_core_otu_new)
nrow(fungi_core_otu_new)

bacteria_fungi_core_otu_new<- rbind(bacteria_core_otu_new,fungi_core_otu_new)
View(bacteria_fungi_core_otu_new)

bacteria_fungi_core_otu_new2<-bacteria_fungi_core_otu_new[which(rowSums(bacteria_fungi_core_otu_new) !=0),]
bacteria_fungi_core_otu_new3<-bacteria_fungi_core_otu_new2[,which(colSums(bacteria_fungi_core_otu_new2) != 0)]
ncol(bacteria_fungi_core_otu_new3)
nrow(bacteria_fungi_core_otu_new3)
write.csv(bacteria_fungi_core_otu_new3,"bacteria_fungi_core_asv_new3.csv")
bacteria_fungi_core_otu_new3<- read.csv("bacteria_fungi_core_asv_new3.csv",header = T,row.names = 1)
bacteria_fungi_core_otu_new3_group<- read.csv("bacteria_fungi_core_asv_new3_group.csv",header = T)
bacteria_fungi_core_otu_new3_t<- t(bacteria_fungi_core_otu_new3)
bacteria_fungi_core_otu_new3_t <- transform(bacteria_fungi_core_otu_new3_t,group =bacteria_fungi_core_otu_new3_group$group )

bacteria_fungi_core_otu_new3_t_mean <- aggregate(. ~group, data = bacteria_fungi_core_otu_new3_t,mean)
write.csv(bacteria_fungi_core_otu_new3_t_mean,"bacteria_fungi_core_otu_new3_t_mean.csv")
bacteria_fungi_core_otu_new3_t_mean_t<-t(bacteria_fungi_core_otu_new3_t_mean)
write.csv(bacteria_fungi_core_otu_new3_t_mean_t,"bacteria_fungi_core_otu_new3_t_mean_t.csv")
bacteria_fungi_core_otu_new3_t_mean_t<- read.csv("bacteria_fungi_core_otu_new3_t_mean_t.csv",header = T,row.names = 1)
bacteria_fungi_core_otu_new3_t_mean_t_relative <- t(bacteria_fungi_core_otu_new3_t_mean_t)/rowSums(t(bacteria_fungi_core_otu_new3_t_mean_t))
write.csv(bacteria_fungi_core_otu_new3_t_mean_t_relative,"bacteria_fungi_core_otu_new3_t_mean_t_relative.csv")

bacteria_asv_id<- read.csv("bacteria_ASV_ID.csv",header = 0)
View(bacteria_asv_id)
bacteria_tax_select_sort<- read.csv("bacteria_tax_select_sort.csv")
View(bacteria_tax_select_sort)
c<- bacteria_asv_id$V1
d<-bacteria_tax_select_sort$X
bacteria_tax_select_sort2<- bacteria_tax_select_sort[d%in%c,]
View(bacteria_tax_select_sort2)
write.csv(bacteria_tax_select_sort2,"bacteria_tax_select_sort2.csv")

###转化为相对丰度#############
library(BiodiversityR)
library(ggplot2)
asv_relative <- t(bacteria_fungi_core_otu_new3)/rowSums(t(bacteria_fungi_core_otu_new3))  #转化为相对丰度
rowSums(asv_relative)
write.csv(asv_relative,"asv_relative.csv")
ggclusternet
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq") # 需要30多分钟，安装时间比较长
install.packages('igraph')
install.packages('network')
install.packages('sna')
# 以下需要20多分钟，安装时间比较长
install.packages('tidyverse')  
# 以下需要20分钟左右，目的是安装ggClusterNet
install.packages('devtools') 
install.packages('ggalluvial')
devtools::install_github("taowenmicro/ggClusterNet") 
######
###https://mp.weixin.qq.com/s/xlay-DL6YU8rM51eaQ02HA
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(tidyfst)
metadata <- read.delim("metadata_bulk.txt",row.names =1)
metadata <- read.csv("metadata_bulk.csv",row.names = 1)
setwd("../core_network/16S_ITS_network_bulk/")
otutab <- read.delim("asv_relative_bulk4.txt", row.names=1)
taxonomy  = read.delim("bacteria_fungi_taxa_zipi.txt", row.names=1)
ps = phyloseq(sample_data(metadata), otu_table(as.matrix(otutab), taxa_are_rows=TRUE), tax_table(as.matrix(taxonomy)))
ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))
ps = phyloseq(sample_data(metadata), otu_table(as.matrix(otutab), taxa_are_rows=TRUE), tax_table(as.matrix(taxonomy)))

#-提取丰度最高的指定数量的otu进行构建网络
#----------计算相关#----
result = corMicro(ps = ps,
                  N = 155,
                  method.scale = "TMM",
                  r.threshold=0.6,
                  p.threshold=0.05,
                  method = "spearman")
#--提取相关矩阵
cor = result[[1]]
cor
#-网络中包含的OTU的phyloseq文件提取
ps_net = result[[3]]
ps_net 
# 以下为输出信息
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 150 taxa and 18 samples ]
#sample_data() Sample Data:       [ 18 samples by 10 sample variables ]
#tax_table()   Taxonomy Table:    [ 150 taxa by 7 taxonomic ranks ]

###这是网络布局的基础，无论是什么聚类布局，都需要制作一个分组文件，这个文件有两列，一列是节点，一列是分组信息，这个分组信息名称为：group。这个文件信息就是用于对节点进行分组，然后按照分组对节点归类，使用包中可视化函数计算节点位置。
# 注意分组文件的格式，分为两列，第一列是网络中包含的OTU的名字，第二列是分组信息，同样的分组标记同样的字符。
#--人工构造分组信息：将网络中全部OTU分为五个部分，等分
netClu = data.frame(ID = row.names(otu_table),group =rep(1,length(row.names(otu_table)))[1:length(row.names(otu_table))] )
netClu$group = as.factor(netClu$group)
head(netClu)
##PolygonClusterG 根据分组，计算布局位置坐标
#--------计算布局#---------
result2 = PolygonClusterG (cor = cor,nodeGroup =netClu  )
node = result2[[1]]
head(node)
tax_table = ps_net %>%
  vegan_tax() %>%
  as.data.frame()
#---node节点注释#-----------
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
head(edge)
p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                             data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p
####按照微生物分类不同设定分组
#----------计算相关#----
result = corMicro (ps = ps,
                   N = 155,
                   method.scale = "TMM",
                   r.threshold=0.6,
                   p.threshold=0.05,
                   method = "spearman"
)
#--提取相关矩阵
cor = result[[1]]
head(cor)
#-网络中包含的OTU的phyloseq文件提取
ps_net = result[[3]]
#-导出otu表格
otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()
tax = ps_net %>% vegan_tax() %>%
  as.data.frame()
tax$filed = tax$Phylum
group2 <- data.frame(ID = row.names(tax),group = tax$Phylum)
group2 <- data.frame(ID = row.names(tax),group = tax$Kingdom)
group2$group  =as.factor(group2$group)

result2 = PolygonClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]
# ---node节点注释#-----------
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

##################################
##################################
###按照网络模块分析定义分组
library(ggraph)
netClu  = modulGroup( cor = cor,cut = NULL,method = "cluster_fast_greedy" )

result2 = model_maptree_group(cor = cor,
                              nodeGroup = group2,
)

# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]

# ---node节点注释#-----------
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)

nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Model_",nodes2$group,sep = "")

#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

##############################
##############################

result2 = model_maptree2(cor = cor,
                         method = "cluster_fast_greedy"
)

# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]

# ---node节点注释#-----------
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)

nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Model_",nodes2$group,sep = "")

#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

#####################
#######model_igraph布局
result = cor_Big_micro(ps = ps,
                       N = 1000,
                       r.threshold=0.6,
                       p.threshold=0.05,
                       method = "spearman"
)

#--提取相关矩阵
cor = result[[1]]
dim(cor)

result2 <- model_igraph(cor = cor,
                        method = "cluster_fast_greedy",
                        seed = 12
)
node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod

#---node节点注释#-----------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)

tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
  select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

library(ggnewscale)

p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                              data = edge2, size = 1) +
  scale_colour_manual(values = col0) 

ggsave("./cs1.pdf",p1,width = 16,height = 14)
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
  scale_colour_manual(values = col) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2
ggsave("./network_bulk_model.pdf",p2,width = 16,height = 14)

############节点模块化可视化
result4 = nodeEdge(cor = cor)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
p <- res[[1]]
p



################################
################################
###https://blog.csdn.net/qazplm12_3/article/details/125967064
result = cor_Big_micro(ps = ps,
                       N = 155,
                       r.threshold=0.6,
                       p.threshold=0.05,
                       method = "spearman"
)

#--提取相关矩阵
cor = result[[1]]
dim(cor)

result2 <- model_igraph(cor = cor,
                        method = "cluster_fast_greedy",
                        seed = 12
)
node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod

#---node节点注释#-----------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)

tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% dplyr::select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color
write.csv(edge2,"network_color.csv")
edge2<- read.csv("network_color.csv",header = T,row.names = 1)
View(edge2)
write.csv(dat,"dat.csv")
dat<- read.csv("dat.csv",header = T,row.names = 1)
View(dat)
library(ggnewscale)
p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                              data = edge2, size = 1) +
  scale_colour_manual(values = col0) 

# ggsave("./cs1.pdf",p1,width = 16,height = 14)
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =group), data = dat,size = 4) +
  #scale_colour_manual(values = col) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2
# ggsave("./cs2.pdf",p2,width = 16,height = 14)
result4 = nodeEdge(cor = cor)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
p <- res[[1]]

####################################
##igraph 包计算网络模块library(igraph)#输入数据示例，邻接矩阵#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adjacency_unweight <- read.delim('./16S_ITS_network/adjacency_matrix/0.6adjacency_matrx.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight 
head(adjacency_unweight)[1:6]    
#邻接矩阵类型的网络文件#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络igraph igraph    
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph    #igraph 的邻接列表
#计算节点度
V(igraph)$degree <- igraph::degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  "nodes_id" = V(igraph)$name, 
  "degree" = V(igraph)$degree, 
  "modularity" = as.numeric(V(igraph)$modularity))
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, './16S_ITS_network/meta-network/network_topological_features/nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度(Zi)和模块间连通度(Pi)
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight <- read.delim('./16S_ITS_network/adjacency_matrix/0.6adjacency_matrx.txt', row.names = 1, sep = '\t', check.names = FALSE)

adjacency_unweight <- read.delim('./adjacency_matrix/0.6adjacency_matrx.txt', row.names = 1, sep = '\t', check.names = FALSE)
#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('./meta-network/network_topological_features/real_meta_network_NodeLevel_top2.txt', row.names = 1, sep = '\t', check.names = FALSE)

nodes_list <- read.delim('./meta-network/network_topological_features/real_meta_network_NodeLevel_top.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度(Zi)和模块间连通度(Pi)
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
#计算模块内连通度（Zi）和模块间连通度（Pi）

#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'Modularity Class')

head(zi_pi) 

write.table(zi_pi, './meta-network/network_topological_features/zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)
##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
zi_pi <- na.omit(zi_pi)  #NA 值最好去掉，不要当 0 处理

zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'

zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs' 

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),    
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black'),    
        panel.background = element_blank(), 
        legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +geom_hline(yintercept = 2.5)


ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  
  geom_point(aes(color = type), alpha = 0.5, size = 3) +
  
  scale_color_manual(values = c('gray','red','blue','purple'),
                     
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        
        panel.background = element_blank(), legend.key = element_blank()) +
  
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  
  geom_vline(xintercept = 0.62) +geom_hline(yintercept = 2.5)
  

##################################
#####################################
##zipi#####
library(ggnewscale)
library(ggplot2)
setwd("../16S_ITS_network_bulk/")
zi_pi <- read.delim('./meta-network/network_topological_features/zi_pi_result.txt', row.names = 1, sep = '\t', check.names = FALSE)
View(zipi_table)
names(zipi_table)

zi_pi <- na.omit(zi_pi)  #NA 值最好去掉，不要当 0 处理

zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'

zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs' 

p3<- ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 1, size = 4) +
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(-2,3.5))+
  scale_color_manual(values = c('#6E8B3E','red','blue','purple'),    
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))+
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept=0.62,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2.5,lty=4,col="black",lwd=0.8)+
  geom_text(aes(0.625,-1,label="B_ASV_3679"),hjust=-0.03,vjust=1.5)
p3
ggsave("./Zi_Pi_bulk.pdf",p3,width = 10,height = 7)  

########
node_feature<- read.delim("./meta-network/network_topological_features/zi_pi_result.txt",
                          row.names = 1, sep = '\t', check.names = FALSE)

View(node_feature)

p4<- ggplot(node_feature, aes(degree,betweenesscentrality)) +
  geom_point(aes(color = group), alpha = 1, size = 4) +
  scale_color_manual(values = c('red','blue'),    
                     limits = c('Bacterial_ASV', 'Fungal_ASV'))+
  scale_x_continuous(limits = c(-7,30))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=17),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=17))+
  labs(x = 'Degree', y = 'Betweeness centrality', color = '') +
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)+
  geom_text(aes(8,3294.644733,label="B_ASV_3853"),hjust=-0.05,vjust=-0.5)+
  geom_text(aes(4,2256,label="F_ASV_13948"),hjust=0.1,vjust=-0.8)+
  geom_text(aes(2,2185,label="F_ASV_23484"),hjust=0.9,vjust=-0.7)+
  geom_text(aes(4,2122.221731,label="B_ASV_1355"),hjust=-0.06,vjust=-0.2)+
  geom_text(aes(2,2112,label="F_ASV_19554"),hjust=1.1,vjust=0.6)+
  geom_text(aes(4,2100.499509,label="B_ASV_3679"),hjust=-0.05,vjust=1)+
  geom_text(aes(3,2037,label="F_ASV_13511"),hjust=-0.03,vjust=1.5)

p4
ggsave("./degree_Betweeness_bulk.pdf",p4,width = 10,height = 7) 

##################

p5<- ggplot(node_feature, aes(degree,closnesscentrality)) +
  geom_point(aes(color = group), alpha = 1, size = 4) +
  scale_color_manual(values = c('#BA55D3','#48D1CC'),    
                     limits = c('Bacterial_ASV', 'Fungal_ASV'))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=17),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=17))+
  labs(x = 'Degree', y = 'Closness centrality', color = '') +
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)
  #geom_text(aes(17,0.01154769,label="B_ASV_15566"),hjust=1,vjust=-0.5)+
  #geom_text(aes(5,0.011530398,label="F_ASV_19554"),hjust=0.1,vjust=-0.7)+
  #geom_text(aes(2,0.011529535,label="F_ASV_23484"),hjust=-0.1,vjust=0.7)+
  #geom_text(aes(10,0.011526946,label="F_ASV_8189"),hjust=0.1,vjust=-0.7)

p5
ggsave("./degree_closnesscentrality_bulk.pdf",p5,width = 10,height = 7) 

####################zipi_rhizosphere###############
setwd("../16S_ITS_network_rhizosphere/")

##计算模块内连通度(Zi)和模块间连通度(Pi)
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight <- read.delim('./adjacency_matrix/0.6adjacency_matrx.txt', row.names = 1, sep = '\t', check.names = FALSE)
#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('./meta-network/network_topological_features/real_meta_network_NodeLevel_top.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度(Zi)和模块间连通度(Pi)
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
#计算模块内连通度（Zi）和模块间连通度（Pi）

#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'Modularity Class')

head(zi_pi) 

write.table(zi_pi, './meta-network/network_topological_features/zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)
##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
zi_pi <- na.omit(zi_pi)  #NA 值最好去掉，不要当 0 处理

zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'

zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'

zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs' 

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),    
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black'),    
        panel.background = element_blank(), 
        legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +geom_hline(yintercept = 2.5)


ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  
  geom_point(aes(color = type), alpha = 0.5, size = 3) +
  
  scale_color_manual(values = c('gray','red','blue','purple'),
                     
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        
        panel.background = element_blank(), legend.key = element_blank()) +
  
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  
  geom_vline(xintercept = 0.62) +geom_hline(yintercept = 2.5)


##################################
#####################################
##zipi#####
library(ggnewscale)
library(ggplot2)

zi_pi <- read.delim('./meta-network/network_topological_features/zi_pi_result.txt', row.names = 1, sep = '\t', check.names = FALSE)
View(zipi_table)
names(zipi_table)
node_feature_rhizosphere<- read.delim("./meta-network/network_topological_features/zi_pi_result.txt",
                                      row.names = 1, sep = '\t', check.names = FALSE)
node_feature_rhizosphere <- na.omit(node_feature_rhizosphere)  #NA 值最好去掉，不要当 0 处理

node_feature_rhizosphere[which(node_feature_rhizosphere$within_module_connectivities < 2.5 & node_feature_rhizosphere$among_module_connectivities < 0.62),'type'] <- 'Peripherals'

node_feature_rhizosphere[which(node_feature_rhizosphere$within_module_connectivities <= 2.5 & node_feature_rhizosphere$among_module_connectivities > 0.62),'type'] <- 'Connectors'

node_feature_rhizosphere[which(node_feature_rhizosphere$within_module_connectivities > 2.5 & node_feature_rhizosphere$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'

node_feature_rhizosphere[which(node_feature_rhizosphere$within_module_connectivities > 2.5 & node_feature_rhizosphere$among_module_connectivities > 0.62),'type'] <- 'Network hubs' 

p7<- ggplot(node_feature_rhizosphere, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 1, size = 4) +
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(-2,3.5))+
  scale_color_manual(values = c('#6E8B3E','red','blue','purple'),    
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))+
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept=0.62,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2.5,lty=4,col="black",lwd=0.8)+
  geom_text(aes(0.625,-1.42313077572628,label="B_ASV_12040"),hjust=-0.03,vjust=1.5)
p7
ggsave("./Zi_Pi_rhizosphere.pdf",p7,width = 10,height = 7)  

########


View(node_feature)

p8<- ggplot(node_feature_rhizosphere, aes(degree,betweenesscentrality)) +
  geom_point(aes(color = group), alpha = 1, size = 4) +
  scale_color_manual(values = c('red','blue'),    
                     limits = c('Bacterial_ASV', 'Fungal_ASV'))+
  scale_x_continuous(limits = c(-7,30))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=17),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=17))+
  labs(x = 'Degree', y = 'Betweeness centrality', color = '') +
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)+
  geom_text(aes(23,1468.645126,label="F_ASV_5263"),hjust=-0.05,vjust=-0.5)+
  geom_text(aes(24,1192.057171,label="B_ASV_15299"),hjust=0.2,vjust=-0.9)+
  geom_text(aes(13,1176.80578,label="F_ASV_8189"),hjust=-0.02,vjust=-0.7)+
  geom_text(aes(10,1135.923433,label="F_ASV_3909"),hjust=-0.1,vjust=0.4)
  

p8
ggsave("./degree_Betweeness_rhizosphere.pdf",p8,width = 10,height = 7) 

##################
library(ggplot2)
p9<- ggplot(node_feature_rhizosphere, aes(degree,closnesscentrality)) +
  geom_point(aes(color = group), alpha = 1, size = 4) +
  scale_color_manual(values = c('#BA55D3','#48D1CC'),    
                     limits = c('Bacterial_ASV', 'Fungal_ASV'))+
  theme(legend.position="right",
        legend.title = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=17),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=17))+
  labs(x = 'Degree', y = 'Closness centrality', color = '') +
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)
#geom_text(aes(17,0.01154769,label="B_ASV_15566"),hjust=1,vjust=-0.5)+
#geom_text(aes(5,0.011530398,label="F_ASV_19554"),hjust=0.1,vjust=-0.7)+
#geom_text(aes(2,0.011529535,label="F_ASV_23484"),hjust=-0.1,vjust=0.7)+
#geom_text(aes(10,0.011526946,label="F_ASV_8189"),hjust=0.1,vjust=-0.7)

p9
ggsave("./degree_closnesscentrality_rhizosphere.pdf",p9,width = 10,height = 7) 

p10<- cowplot::plot_grid(p3, p4, p5, p7,p8,p9, ncol=2, labels = letters[1:6])
p10



setwd("../species_annotation")
fungi_asv<- read.csv("fungi_asv_select_rhizo.csv",header = T,row.names = 1)
fungi_asv_distance<- vegdist(t(fungi_asv), method = 'bray')
fungi_nmds_before_group<- read.csv("fungi_nmds_before_group.csv",header = T,row.names = 1)
dune.ano <- with(fungi_nmds_before_group, anosim(fungi_asv_distance, predict_group))

df_nmds <- metaMDS(otu.distance, k = 2)
summary(df_nmds)
df_nmds_stress <- df_nmds$stress
df_nmds_stress

bacteria_asv<- read.csv("bacteria_asv_select_rhizo.csv",header = T,row.names = 1)
bacteria_asv_distance<- vegdist(t(bacteria_asv), method = 'bray')
bacteria_nmds_before_group<- read.csv("bacteria_nmds_before_group.csv",header = T,row.names = 1)
bac.dune.ano <- with(bacteria_nmds_before_group, anosim(bacteria_asv_distance, predict_group))


setwd("../alpha_diversity")
library(tidyverse)
library(rstatix)
library(ggpubr)
set.seed(123)
data("anxiety", package = "datarium")
anxiety %>% sample_n_by(group, size = 1)
anxiety <- anxiety %>%
  gather(key = "time", value = "score", t1, t2, t3) %>%
  convert_as_factor(id, time)
set.seed(123)
anxiety %>% sample_n_by(group, time, size = 1)
stat.test <- anxiety %>%
  group_by(group) %>%
  pairwise_t_test(
    score ~ time, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test
# 绘制图形
bxp <- ggboxplot(
  anxiety, x = "group", y = "score",
  color = "time", palette = "jco"
)
# 添加显著性标记
stat.test <- stat.test %>% add_xy_position(x = "group")
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE, tip.length = 0
)


bacteria_alpha_diversity<- read.csv("bacteria_alpha_diversity4.csv",header = T,row.names = 1)
pairwise.wilcox.test(bacteria_alpha_diversity$shannon_entropy, bacteria_alpha_diversity$group6, p.adjust.method = "bonf")
pairwise.wilcox.test(bacteria_alpha_diversity$shannon_entropy, bacteria_alpha_diversity$group4, p.adjust.method = "bonf")

fungi_alpha_diversity<- read.csv("fungi_alpha_diversity3.csv",header = T,row.names = 1)
pairwise.wilcox.test(fungi_alpha_diversity$shannon_entropy, fungi_alpha_diversity$group6, p.adjust.method = "bonf")
pairwise.wilcox.test(fungi_alpha_diversity$shannon_entropy, fungi_alpha_diversity$group4, p.adjust.method = "bonf")


stat.test <- bacteria_alpha_diversity %>%
  group_by(group6) %>%
  pairwise_t_test(
    shannon_entropy ~ group6, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details

stat.test
names(bacteria_alpha_diversity)

setwd("./16S_ITS_network_bulk")
asv_relative_bulk<-read.csv("asv_relative_bulk.csv",header = T,row.names = 1)
asv_relative_bulk2<-asv_relative_bulk[which(rowSums(asv_relative_bulk) !=0),]
asv_relative_bulk3<-asv_relative_bulk2[,which(colSums(asv_relative_bulk2) != 0)]
write.csv(asv_relative_bulk3,"asv_relative_bulk3.csv")

setwd("../16S_ITS_network_rhizosphere")
asv_relative_rhizosphere<-read.csv("asv_relative_rhizosphere.csv",header = T,row.names = 1)
asv_relative_rhizosphere2<-asv_relative_rhizosphere[which(rowSums(asv_relative_rhizosphere) !=0),]
asv_relative_rhizosphere3<-asv_relative_rhizosphere2[,which(colSums(asv_relative_rhizosphere2) != 0)]
write.csv(asv_relative_rhizosphere3,"asv_relative_rhizosphere3.csv")


#####ggplot_bar###############
library(ggplot2)
topp<- read.table("clipboard",header = T)
View(topp)
p1<- ggplot(data=topp,mapping=aes(x=name, y =value,fill=group))+
  geom_bar(stat="identity",width=0.8,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Average_degree") +
  theme_set(theme_bw())+
  theme(legend.position = "top",
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=45,size=20),
        axis.title.x = element_text(face = "bold",colour="black",size=20),
        axis.title.y = element_text(face = "bold",colour="black",size=20))
p1
ggsave(p1, file='sample_degree_path.pdf', width=7, height=9)

topp2<- read.table("clipboard",header = T,row.names = 1)

p3<- ggplot(data=topp2,mapping=aes(x=name, y =value,fill=group))+
  geom_bar(stat="identity",width=0.8,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "value") +
  theme_set(theme_bw())+
  theme(legend.position = "top",
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=45,size=20),
        axis.title.x = element_text(face = "bold",colour="black",size=20),
        axis.title.y = element_text(face = "bold",colour="black",size=20))
p3
setwd("../core_network/")
ggsave(p3, file='meta_degree_path.pdf', width=7, height=9)
