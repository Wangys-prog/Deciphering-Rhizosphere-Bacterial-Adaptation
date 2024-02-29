ASV_train.cv <- replicate(5, rfcv(ASV_train[-ncol(ASV_train)], ASV_train$group4, cv.fold = 10,step = 1.5), simplify = FALSE)
ASV_train.cv <- replicate(5, rfcv(ASV_train[-ncol(ASV_train)], as.factor(ASV_train$group4), cv.fold = 10,step = 1.5), simplify = FALSE)
bacteria_genus_sum_t<- read.csv("bacteria_genus_sum_t.csv",header = T,row.names = 1)
getwd()
setwd("../species_annotation")
bacteria_genus_sum_t<- read.csv("bacteria_genus_sum_t.csv",header = T,row.names = 1)
#将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
select_train <- sample(120, 120*0.7)
ASV_train <- bacteria_genus_sum_t[-select_train, ]
ASV_test <- bacteria_genus_sum_t[select_train, ]
ASV_train.forest <- randomForest(as.factor(group4) ~ ., data = ASV_train, importance = TRUE)
library(randomForest)
library(ggplot2)
bacteria_asv_taxa_genus_sum <- aggregate(. ~ g, data = bacteria_asv_taxa,sum)
ASV_train.forest <- randomForest(as.factor(group4) ~ ., data = ASV_train, importance = TRUE)
ASV_train.forest
plot(randomForest::margin(ASV_train.forest), main = '观测值被判断正确的概率图')
# 原来的做图代码如下：出现Error in unit(c(t, r, b, l), unit) :  'list' object cannot be coerced to type 'double'
# plot(margin(otu_train.forest, otu_traingroups), main = '观测值被判断正确的概率图')
#训练集自身测试
train_predict <- predict(ASV_train.forest, ASV_train)
compare_train <- table(train_predict,  ASV_train$group4)
sum(diag(compare_train)/sum(compare_train))
compare_train
#使用测试集评估
test_predict <- predict(ASV_train.forest, ASV_test)
compare_test <- table(ASV_test$group4, test_predict, dnn = c('Actual', 'Predicted'))
compare_test
##关键 ASVs 识别
#查看表示每个变量（ASVs）重要性的得分
#summary(otu_train.forest)
importance_ASV <- ASV_train.forest$importance
head(importance_ASV)
#或者使用函数 importance()
importance_ASV <- data.frame(importance(ASV_train.forest))
head(importance_ASV)
#可以根据某种重要性的高低排个序，例如根据“Mean Decrease Accuracy”指标
importance_ASV <- importance_ASV[order(importance_ASV$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_ASV)
#输出表格
write.table(importance_ASV, 'bacteria_importance_ASV.txt', sep = '\t', col.names = NA, quote = FALSE)
#作图展示 top30 重要的 OTUs
varImpPlot(ASV_train.forest, n.var = min(30, nrow(ASV_train.forest$importance)), main = 'Top 30 - variable importance')
##交叉验证帮助选择特定数量的 OTUs
#5 次重复十折交叉验证
set.seed(123)
#作图展示 top30 重要的 OTUs
varImpPlot(ASV_train.forest, n.var = min(30, nrow(ASV_train.forest$importance)), main = 'Top 30 - variable importance')
##交叉验证帮助选择特定数量的 OTUs
#5 次重复十折交叉验证
set.seed(123)
ASV_train.cv <- replicate(5, rfcv(ASV_train[-ncol(ASV_train)], as.factor(ASV_train$group4), cv.fold = 10,step = 1.5), simplify = FALSE)
library(dplyr)
fungi_genus<- read.csv("./fungi_asv_taxa_genus_sum.csv",header = T,row.names = 1)
fungi_importance<- read.csv("./fungi_importance_genus.csv",header = T,row.names = 1)
fungi_top60<- fungi_importance[c(1:60),]
fungi_asv_genus_top60<- fungi_genus[rownames(fungi_top60),]
write.csv(fungi_asv_genus_top60,"fungi_asv_genus_top60.csv")
fungi_asv_genus_top60<- read.csv("./fungi_asv_genus_top60.csv",header = T,row.names = 1)
setwd("../")
setwd("./species_annotation")
fungi_genus<- read.csv("./fungi_asv_taxa_genus_sum.csv",header = T,row.names = 1)
fungi_importance<- read.csv("./fungi_importance_genus.csv",header = T,row.names = 1)
fungi_top60<- fungi_importance[c(1:60),]
fungi_asv_genus_top60<- fungi_genus[rownames(fungi_top60),]
write.csv(fungi_asv_genus_top60,"fungi_asv_genus_top60.csv")
fungi_asv_genus_top60<- read.csv("./fungi_asv_genus_top60.csv",header = T,row.names = 1)
fungi_asv_genus_top60_t<- t(fungi_asv_genus_top60)
group <- read.csv("./fungi_heatmap_group.csv",header = T,row.names = 1)
fungi_asv_genus_top60_t_group<- cbind(fungi_asv_genus_top60_t,group)
datamean=group_by(fungi_asv_genus_top60_t_group, group) %>% summarize_each(funs(mean))
View(fungi_asv_genus_top60_t_group)
ncol(fungi_asv_genus_top60_t_group)
ncol(fungi_asv_genus_top60_t)
datamean=group_by(fungi_asv_genus_top60_t_group, group) %>% summarize_each(funs(mean))
write.csv(datamean,"fungi_datamean.csv")
fungi_datamean<- read.csv("fungi_datamean.csv",header = T,row.names = 1)
data_matrix <- as.matrix(select(fungi_datamean,-1))
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-2)))
data_matrix3<-scale(t(data_matrix2),center = T)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
library(gplots)
library(RColorBrewer)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
col<- colorRampPalette(c("red","white","blue"))(256)
group<- fungi_datamean$group2
plot_color = c('chocolate','cyan4')[as.factor(group)]
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
data_matrix <- as.matrix(select(fungi_datamean,-1))
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-2)))
data_matrix3<-scale(t(data_matrix2),center = T)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
plot_color = c('chocolate','cyan4')[as.factor(group)]
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
group<- fungi_datamean$group2
plot_color = c('chocolate','cyan4')[as.factor(group)]
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
col<- colorRampPalette(c("red","yellow","blue"))(256)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6))
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
data_matrix <- as.matrix(select(fungi_datamean,-1))
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
colnames(data_matrix2)
data_matrix3<-scale(t(data_matrix2),center = T)
col<- colorRampPalette(c("red","yellow","blue"))(256)
group<- fungi_datamean$group2
plot_color = c('chocolate','cyan4')[as.factor(group)]
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
col<- colorRampPalette(c("red","white","blue"))(256)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
colnames(fungi_datamean)
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
colnames(data_matrix2)
data_matrix3<-scale(t(data_matrix2),center = T)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6)
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
ColSideColors =plot_color,
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=plot_color,
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=fungi_datamean$group2,
)
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
data_matrix <- as.matrix(select(fungi_datamean,-1))
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
data_matrix3<-scale(t(data_matrix2),center = T)
col<- colorRampPalette(c("red","white","blue"))(256)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=fungi_datamean$group2,
)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
)
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
data_matrix3<-scale(t(data_matrix2),center = T)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
)
?pheatmap
ann_colors=list(class=c(L='#009933',R='#CC33CC',F='#FDDCA9'))
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_colors=ann_colors)
group<- fungi_datamean$group2
group
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_row=class)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_row=group)
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
data_matrix3<-scale(t(data_matrix2),center = T)
group<- fungi_datamean$group2
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_row=group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_col=group)
ncol(data_matrix3)
ncol(group)
group<- data.frame(fungi_datamean$group2)
View(group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_col=group)
ncol(group)
nrow(group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_row ==group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_row =group)
row.names(data_matrix3)
colnames(data_matrix3)
group<- read.csv("fungi_heatmap_group2.csv",header=T,row.names=1)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,annotation_col=group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_row=group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_row=as.data.frame(group))
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_row=as.data.frame(fungi_datamean[,"group2"], row.names=rownames(fungi_datamean)))
nrow(group)
rownames(group) <- colnames(data_matrix3)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_row=group)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group)
ann_colors=list(group=c(Bulk='chocolate',Rhizophere='cyan4'))
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
ann_colors=list(group=c(Bulk=#CC6666',Rhizophere='cyan4'))
ann_colors=list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))
ann_colors<-list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
group<- read.csv("fungi_heatmap_group2.csv",header=T,row.names=1)
ann_colors<-list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
rownames(group) <- colnames(data_matrix3)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
ann_colors <- list(group=Sexcolor)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
Sexcolor <- c("red","#016D06")
ann_colors <- list(group=Sexcolor)
pheatmap::pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
library(pheatmap)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors=ann_colors)
ann_colors<-list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))
ann_colors
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
anno_color<-list(Treatment=c(BL="blue",GL="green",RL="red",WL="grey"))
anno_color
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_color)
anno_color<-list(Treatment=c(BL="blue",GL="green",RL="red",WL="grey"))
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_color)
ann_colors<-list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_color)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
bacteria_datamean<- read.csv("bacteria_datamean2.csv",header = T,row.names = 1)
data_matrix2 <- as.matrix(select(bacteria_datamean,c(-1,-2)))
data_matrix3<-scale(t(data_matrix2),center = T)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
group<- bacteria_datamean$group2
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
rownames(group) <- colnames(data_matrix3)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
group<- as.data.frame(bacteria_datamean$group2)
rownames(group) <- colnames(data_matrix3)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)
getwd()
setwd("../")
setwd("./species_annotation")
library(pheatmap)
library(dplyr)
library(gplots)
library(RColorBrewer)
bacteria_datamean<- read.csv("bacteria_datamean2.csv",header = T,row.names = 1)
data_matrix2 <- as.matrix(select(bacteria_datamean,c(-1,-2)))
data_matrix3<-scale(t(data_matrix2),center = T)
heatmap.2(data_matrix3,
trace="none",#不显示trace
ColSideColors =plot_color,
col=col,#修改热图颜色
density.info = "none",#图例取消density
key.xlab ='Correlation',
key.title = "",
cexRow = 1,cexCol = 1,#修改横纵坐标字体
Rowv = F,Colv = F, #去除聚类
margins = c(6, 6))
pheatmap(data_matrix3,cluster_row = FALSE,cluster_col=FALSE,annotation_col=group,annotation_colors = ann_colors)
bacteria_datamean<- read.csv("bacteria_datamean2.csv",header = T,row.names = 1)
data_matrix2 <- as.matrix(select(bacteria_datamean,c(-1,-2)))
data_matrix3<-scale(t(data_matrix2),center = T)
heatmap.2(data_matrix3,trace="none",#不显示traceColSideColors =plot_color,col=col,#修改热图颜色density.info = "none",#图例取消densitykey.xlab ='Correlation',
key.title = "",cexRow = 1,cexCol = 1,#修改横纵坐标字体Rowv = F,Colv = F, #去除聚类margins = c(6, 6))
pheatmap(data_matrix3,
cluster_row = FALSE,
cluster_col=FALSE,
annotation_col=group,
annotation_colors = ann_colors)

setwd("./species_annotation")
