require(data.table)
##select_sample
fungi_class<-read.csv("fungi_asv_select10_4000.csv",header=T,row.names = 1)
View(fungi_class)
##select_sample
fungi_asv<-read.csv("fungi_asv_select10_4000.csv",header=T,row.names = 1)
colnames(fungi_asv)
rownames(fungi_asv)
ncol(fungi_asv)
sample_id <- scan("fungi_group.csv",what = "character")
View(sample_id)
##选择列
newdata2<-fungi_asv[ ,colnames(fungi_asv) %in% sample_id]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
sample_id <- scan("fungi_group2.csv",what = "character")
##选择列
newdata2<-fungi_asv[ ,colnames(fungi_asv) %in% sample_id]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
rownames(fungi_asv)
colnames(fungi_asv)
ncol(fungi_asv)
colnames(fungi_asv)
colnames(fungi_asv) %in% sample_id
sample_id <- scan("fungi_group2.csv",what = "character")
##选择列
newdata2<-fungi_asv[ ,colnames(fungi_asv) %in% sample_id]
nrow(newdata2)
ncol(newdata2)
#去掉全为0的行和列
colnames(newdata2)
View(rowSums(newdata2))
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
write.csv(newdata4,"fungi_asv_select_rhizo.csv")
write.csv(colnames(newdata4),"fungi_group3.csv")
rownames(newdata4)
colnames(newdata4)
ncol(newdata4)
nrow(newdata4)
fungi_asv_select<-data.frame(newdata4)
fungi_asv_select_t<-t(fungi_asv_select)
write.csv(colnames(newdata4),"fungi_group3.csv")

#####pcoa####
library(vegan)
library(ape)
library(ggplot2)
fungi_asv_select<-data.frame(newdata4)
fungi_asv_select <- read.csv("fungi_asv_select_rhizo.csv",header=T,row.names=1)
fungi_asv_select_t<-t(fungi_asv_select)
fungi_dm<- vegdist(fungi_asv_select_t, "bray")
res2 <- pcoa(fungi_dm)
pcoa <- cmdscale (fungi_dm,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

write.csv(pc12,"fungi_bray_pcoa_values.csv")
write.csv(res2$vectors,"fungi_bray_pcoa_vectors.csv")

varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)


###########bacteria##########
##select_sample
bacteria_asv<-read.csv("bacteria_asv_original0.05.csv",header=T,row.names = 1)
sample_id <- scan("bacteria_group2.csv",what = "character")
colnames(bacteria_asv)
rownames(bacteria_asv)
##选择列
newdata2<-bacteria_asv[ ,colnames(bacteria_asv) %in% sample_id]
nrow(newdata2)
ncol(newdata2)
View(newdata2)
#去掉全为0的行和列
colnames(newdata2)
newdata3<-newdata2[which(rowSums(newdata2) !=0),]
newdata4<-newdata2[,which(colSums(newdata3) != 0)]
write.csv(colnames(newdata4),"bacteria_group3.csv")
rownames(newdata4)
colnames(newdata4)
write.csv(newdata4,"bacteria_asv_select_rhizo.csv")
rownames(newdata4)
ncol(newdata4)
nrow(newdata4)

#####pcoa####
library(vegan)
library(ape)
library(ggplot2)
bacteria_asv_select<-data.frame(newdata4)
bacteria_asv_select <- read.csv("bacteria_asv_select_rhizo.csv",header=T,row.names=1)
bacteria_asv_select_t<-t(bacteria_asv_select)
bacteria_dm<- vegdist(bacteria_asv_select_t, "bray")
res2 <- pcoa(bacteria_dm)
write.csv(res2$values,"bacteria_bray_pcoa_values.csv")
write.csv(res2$vectors,"bcteria_bray_pcoa_vectors.csv")

pcoa <- cmdscale (bacteria_dm,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

write.csv(pc12,"bacteria_bray_pcoa_values.csv")

library(ggplot2)

pcoa_vector<- read.csv("bacteria_bray_pcoa_values.csv",header = T,row.names = 1)

# df <- merge(pc12,group,by="samples")
# color=c("#1597A5","#FFC24B","#FEB3AE")

color=c("#1597A5","#FFC24B")
p1<-ggplot(data=pcoa_vector,aes(x=V1,y=V2,colour=group6))+
  guides(fill=guide_legend(title=NULL))+
  geom_point(size=2.5)+
  geom_vline(xintercept = 0,lty="dashed",size=1.2)+
  geom_hline(yintercept = 0,lty="dashed",size=1.2)+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  xlab("PC1 6.78%")+
  ylab("PC2 3.45%")+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20),
        axis.title.x = element_text(face = "bold",colour="black",size=20),
        axis.title.y = element_text(face = "bold",colour="black",size=20))

p1


pcoa_vector<- read.csv("fungi_bray_pcoa_values.csv",header = T,row.names = 1)
View(pcoa_vector)
# df <- merge(pc12,group,by="samples")
# color=c("#1597A5","#FFC24B","#FEB3AE")

color=c("#1597A5","#FFC24B")
p1<-ggplot(data=pcoa_vector,aes(x=V1,y=V2,
                                color=group6))+
  guides(fill=guide_legend(title=NULL))+
  geom_point(size=2.5)+
  geom_vline(xintercept = 0,lty="dashed",size=1.2)+
  geom_hline(yintercept = 0,lty="dashed",size=1.2)+
  xlab("PC1 6.19%")+
  ylab("PC2 4.12%")+
  scale_fill_manual(values = c("#1597A5","#FFC24B"))+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20),
        axis.title.x = element_text(face = "bold",colour="black",size=20),
        axis.title.y = element_text(face = "bold",colour="black",size=20))

p1

#进行Anosim分析
library(vegan)
bacteria_pcoa_vector<- read.csv("bacteria_bray_pcoa_values.csv",header = T,row.names = 1)
anosim=anosim(bacteria_asv_select_t, bacteria_pcoa_vector$group4, permutations=999)
summary(anosim)

#进行Anosim分析
library(vegan)
fungi_pcoa_vector<- read.csv("fungi_bray_pcoa_values.csv",header = T,row.names = 1)
anosim=anosim(fungi_asv_select_t, fungi_pcoa_vector$group4, permutations=999)
summary(anosim)


library(ggplot2)
bacteria_pcoa <-
g <- ggplot(pCoaVecsTaxa,aes(x=PCo1,y=PCo2,color=P.copri)) + 
  geom_point(size=1.25) + 
  theme_classic() + 
  scale_color_gradientn(colours = rainbow(5)) + 
  xlab(paste0('PCo1 ',round(xVar,2),' %')) + 
  ylab(paste0('PCo2 ',round(yVar,2),' %')) + 
  theme(text = element_text(size = 16))
print(g)


pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
{library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])}
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

pairwise.adonis(t(otu), group$site, sim.method="bray", p.adjust.m= "bonferroni")

pairwise.adonis(bacteria_asv_select_t, bacteria_pcoa_vector$group5, sim.method="bray", p.adjust.m= "bonferroni")

fungi_pcoa_vector<- read.csv("fungi_bray_pcoa_values.csv",header = T,row.names = 1)
pairwise.adonis(fungi_asv_select_t, fungi_pcoa_vector$group5, sim.method="bray", p.adjust.m= "bonferroni")


bacteria_alpha<- read.csv("bacteria_alpha_diversity4.csv",header = T,row.names = 1)
label<- factor(bacteria_alpha$group6,ordered=TRUE,levels=c("Bulk","Apiaceae",
                                                   "Burseraceae","Cruciferae",
                                                   "Cucurbitaceae","Gramineae",
                                                   "Lamiaceae","Leguminosae",
                                                   "Musaceae","Solanaceae",
                                                   "Taxodiaceae"))
library(RColorBrewer)
mycol= brewer.pal(n =11, name = "Set3")
#View(bacteria_alpha)
ggplot(data=bacteria_alpha, mapping = aes(x=label, y=fisher_alpha, group=label,fill=label)) +
  geom_boxplot()+
  scale_fill_manual(values = mycol) +
  #geom_errorbar(aes(ymin=chao1+SD,ymax=chao1-SD),width=0.2) +
  xlab("") +
  ylab("Fisher_alpha")+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=30,size=15),
        axis.text = element_text(face = "bold",colour="black",size=15))

fungi_alpha2<- read.csv("fungi_alpha_diversity3.csv",header = T,row.names = 1)
label2<- factor(fungi_alpha2$group6,ordered=TRUE,levels=c("Bulk","Burseraceae","Cuculidae",
                                                          "Cruciferae","Gramineae", "Lamiaceae",
                                                          "Leguminosae","Malvaceae","Musaceae",
                                                          "SolaBulkceae","Taxodiaceae","Vitaceae"))
library(RColorBrewer)
mycol= brewer.pal(n =11, name = "Set3")
#View(bacteria_alpha)
ggplot(data=fungi_alpha2, mapping = aes(x=label2, y=simpson, group=label2,fill=label2)) +
  geom_boxplot()+
  scale_fill_manual(values = mycol) +
  #geom_errorbar(aes(ymin=chao1+SD,ymax=chao1-SD),width=0.2) +
  xlab("") +
  ylab("Simpson")+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=30,size=15),
        axis.text = element_text(face = "bold",colour="black",size=15))

sample_id <- scan("bacteria_taxa_rhizo.csv",what = "character")
##选择列
bacteria_taxa_select_sort<- read.csv("./bacteria_tax_select_sort.csv",header = T,row.names = 1)
bacteria_rhizo<-bacteria_taxa_select_sort[rownames(bacteria_taxa_select_sort) %in% sample_id,]
View(bacteria_rhizo)


sample_id <- scan("fungi_taxa_rhizo.csv",what = "character")
##选择列
fungi_taxa_select_sort<- read.csv("./fungi_tax_select_sort_new.csv",header = T,row.names = 1)
fungi_rhizo<-fungi_taxa_select_sort[rownames(fungi_taxa_select_sort) %in% sample_id,]
View(fungi_rhizo)
write.csv(fungi_rhizo,"fungi_taxa_rhizo2.csv")

install.packages("randomForest")
library(randomForest)

############### 随机森林 #####################
#随机森林
# install.packages("randomForest")
# install.packages("ROCR")
# install.packages("Hmisc")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)

#制作性状表格

trait=as.data.frame(row.names(signres_norm_matrix_forest))
trait=cbind(trait,substr(trait$`row.names(signres_norm_matrix_forest)`,1,12))
colnames(trait)[2]="barcode"
good_trait=cbind(good_OS_patients[,1],rep("good",length(good_OS_patients[,1])))
poor_trait=cbind(poor_OS_patients[,1],rep("poor",length(poor_OS_patients[,1])))
merge_trait=as.data.frame(rbind(good_trait,poor_trait))
colnames(merge_trait)[1]="barcode"
trait=merge(trait,merge_trait,by="barcode")
trait=trait[,-1]
colnames(trait)[1]="barcode"
# barcode=as.character(trait$barcode)
# status=trait$V2
# trait=data.frame(barcode,status)
forest_trait=cbind(row.names(signres_norm_matrix_forest),signres_norm_matrix_forest)
colnames(forest_trait)[1]="barcode"
forest_trait=merge(trait,forest_trait,by="barcode")
colnames(forest_trait)[2]="trait"
row.names(forest_trait)=forest_trait[,1]
forest_trait=forest_trait[,-1]
#将数据集分为2/3训练集和1/3测试集,并查看数据集基本属性。
ind=sample(2,nrow(forest_trait),replace=TRUE, prob=c(0.67,0.33))
train=forest_trait[ind==1,]
test=forest_trait[ind==2,]
#选取mtry节点值,选择误差最小时对应的i值
#每次都set.seed一下，保证数据的可重复性
errset=1:(length(names(train))-1)
for (i in 1:length(names(train))-1){
  set.seed(350)
  mtry_fit=randomForest(trait ~ .,train,mtry=i)
  err=mean(mtry_fit$err.rate)
  print(err)
  errset[i]=err
}
best_mtry=which(errset==min(errset),arr.ind=TRUE)
#选择ntree，选择图像趋于稳定时的树值
set.seed(350)
model <- randomForest(trait~., train, mtry=best_mtry[1], ntree = 1000)
plot(model)
#大约400棵树比较稳定了

#训练集进行随机森林
set.seed(350)
rf <- randomForest(trait~., train, mtry=best_mtry[1], ntree = 400, importance = TRUE)

#查看变量的重要性
importance <- importance(x=rf)

#绘制变量的重要性图
varImpPlot(rf)

#最后验证并预测

pred1<-predict(rf,data=train)
Freq1<-table(pred1,train$trait)
sum(diag(Freq1))/sum(Freq1)
plot(margin(rf,test$trait))
#全体数据用于训练集
set.seed(350)
rf_output=randomForest(trait ~ ., forest_trait, importance = TRUE, mtry=best_mtry[1], ntree = 400, proximity=TRUE)

#查看变量的重要性
importance <- importance(x=rf_output)

#绘制变量的重要性图
varImpPlot(rf_output)


#####随机森林差异物种###################
###https://cloud.tencent.com/developer/article/1870681

install.packages("randomForest")
library(randomForest)
library(ggplot2)

bacteria_asv_rhizo<- read.csv("./bacteria_asv_select_rhizo.csv",header = T,row.names = 1)
fungi_asv_rhizo<- read.csv("./fungi_asv_select_rhizo.csv",header = T,row.names = 1)

bacteria_asv_select_t <- t(bacteria_asv_rhizo)
fungi_asv_select_t <- t(fungi_asv_rhizo)

taxa<- read.csv("./bacteria_taxa_rhizo2.csv",header = T,row.names = 1)
taxa<- read.csv("./fungi_taxa_rhizo2.csv",header = T,row.names = 1)

g<- taxa$g

bacteria_asv_taxa<- cbind(bacteria_asv_select,g)
fungi_asv_taxa<- cbind(fungi_asv_select,g)

bacteria_asv_taxa_df<- data.frame(bacteria_asv_taxa)

fungi_asv_taxa_df<- data.frame(fungi_asv_taxa)

bacteria_asv_taxa_genus_sum <- aggregate(. ~ g, data = bacteria_asv_taxa,sum)

fungi_asv_taxa_genus_sum <- aggregate(. ~ g, data = fungi_asv_taxa,sum)

write.csv(bacteria_asv_taxa_genus_sum,"bacteria_asv_taxa_genus_sum.csv")

write.csv(fungi_asv_taxa_genus_sum,"fungi_asv_taxa_genus_sum.csv")

bacteria_asv_taxa_genus_sum <- read.csv("bacteria_asv_taxa_genus_sum.csv",header=T,row.names = 1)

fungi_asv_taxa_genus_sum <- read.csv("fungi_asv_taxa_genus_sum.csv",header=T,row.names = 1)

bacteria_genus_sum_t<- t(bacteria_asv_taxa_genus_sum)

fungi_genus_sum_t<- t(fungi_asv_taxa_genus_sum)

write.csv(bacteria_genus_sum_t,"bacteria_genus_sum_t.csv")

write.csv(fungi_genus_sum_t,"fungi_genus_sum_t.csv")

bacteria_genus_sum_t<- read.csv("bacteria_genus_sum_t.csv",header = T,row.names = 1)

fungi_genus_sum_t<- read.csv("fungi_genus_sum_t.csv",header = T,row.names = 1)

#将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
select_train <- sample(120, 120*0.7)
ASV_train <- bacteria_genus_sum_t[-select_train, ]
ASV_test <- bacteria_genus_sum_t[select_train, ]
ASV_train.forest <- randomForest(as.factor(group4) ~ ., data = ASV_train, importance = TRUE) 
plot(randomForest::margin(ASV_train.forest), main = '观测值被判断正确的概率图') 

set.seed(123)
select_train <- sample(120, 120*0.7)
ASV_train <- fungi_genus_sum_t[-select_train, ]
ASV_test <- fungi_genus_sum_t[select_train, ]
ASV_train.forest <- randomForest(as.factor(group4) ~ ., data = ASV_train, importance = TRUE) 
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

write.table(importance_ASV, 'fungi_importance_genus.txt', sep = '\t', col.names = NA, quote = FALSE)

###此处“Mean Decrease Accuracy”和“Mean Decrease Gini”为随机森林模型中的两个重要指标。其中，“mean decrease accuracy”表示随机森林预测准确性的降低程度，该值越大表示该变量的重要性越大；“mean decrease gini”计算每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性。该值越大表示该变量的重要性越大。

##到这一步，可从中筛选一些关键OTUs作为代表物种，作为有效区分两种环境的生物标志物。

###该图展示了其中top30关键的OTUs，将它们划分为“关键OTUs”的依据为模型中的两个重要指标（两个指标下各自包含30个OTUs，默认由高往低排）。

#作图展示 top30 重要的 OTUs
varImpPlot(ASV_train.forest, n.var = min(60, nrow(ASV_train.forest$importance)), main = 'Top 60 - variable importance')

##交叉验证帮助选择特定数量的 OTUs
#5 次重复十折交叉验证
set.seed(123)
ASV_train.cv <- replicate(5, rfcv(ASV_train[-ncol(ASV_train)], as.factor(ASV_train$group4), cv.fold = 10,step = 1.5), simplify = FALSE)
ASV_train.cv

#提取验证结果绘图
ASV_train.cv <- data.frame(sapply(ASV_train.cv, '[[', 'error.cv'))
ASV_train.cv$asvs <- rownames(ASV_train.cv)
ASV_train.cv <- reshape2::melt(ASV_train.cv, id = 'asvs')
ASV_train.cv$asvs <- as.numeric(as.character(ASV_train.cv$asvs))

#拟合线图
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线
pdf("bacteria_asv_cross.pdf",width=9,height=6)
pdf("fungi_asv_cross.pdf",width=9,height=6)

p <- ggplot(ASV_train.cv, aes(asvs, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p
#交叉验证曲线展示了模型误差与用于拟合的OTUs数量之间的关系。误差首先会随OTUs数量的增加而减少，
#开始时下降非常明显，但到了特定范围处，下降幅度将不再有显著变化，甚至有所增加。再根据简约性原则，
#大致选择重要性排名前30的OTUs就可以了。
#大约提取前 30 个重要的 OTUs
p + geom_vline(xintercept = 60)
dev.off()
#根据 OTUs 重要性排序后选择，例如根据“Mean Decrease Accuracy”指标
importance_ASV <- importance_ASV[order(importance_ASV$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_ASV)

#输出表格
write.table(importance_ASV[1:60, ], 'fungi_importance_asv_top60.txt', sep = '\t', col.names = NA, quote = FALSE)

##简约分类器
#选择 top30 重要的 OTUs，例如上述已经根据“Mean Decrease Accuracy”排名获得
ASV_select <- rownames(importance_ASV)[1:50]
ASV_select <- rownames(importance_ASV)[1:60]

#数据子集的训练集和测试集
ASV_train_top50 <- ASV_train[ ,c(ASV_select, 'group4')]
ASV_test_top50 <- ASV_test[ ,c(ASV_select, 'group4')]

ASV_train_top60 <- ASV_train[ ,c(ASV_select, 'group4')]
ASV_test_top60 <- ASV_test[ ,c(ASV_select, 'group4')]

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
ASV_train.forest_50 <- randomForest(as.factor(group4) ~ ., data = ASV_train_top50, importance = TRUE)
ASV_train.forest_50
plot(margin(ASV_train.forest_50, ASV_test_top50$group4), main = '观测值被判断正确的概率图')

set.seed(123)
ASV_train.forest_60 <- randomForest(as.factor(group4) ~ ., data = ASV_train_top60, importance = TRUE)
ASV_train.forest_60
plot(margin(ASV_train.forest_60, ASV_test_top60$group4), main = '观测值被判断正确的概率图')


#训练集自身测试
train_predict <- predict(ASV_train.forest_50, ASV_train_top50)
compare_train <- table(train_predict, ASV_train_top50$group4)
compare_train

train_predict <- predict(ASV_train.forest_60, ASV_train_top60)
compare_train <- table(train_predict, ASV_train_top60$group4)
compare_train

#使用测试集评估
test_predict <- predict(ASV_train.forest_50, ASV_test_top50)
compare_test <- table(ASV_test_top50$group4, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

test_predict <- predict(ASV_train.forest_60, ASV_test_top60)
compare_test <- table(ASV_test_top60$group4, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

##NMDS 排序图中展示分类
#NMDS 降维
nmds <- vegan::metaMDS(bacteria_genus_sum_t[,!names(bacteria_genus_sum_t) %in% c("group4")], distance = 'bray')
nmds <- vegan::metaMDS(bacteria_genus_sum_t[,rownames(importance_ASV)[1:50]], distance = 'bray')

nmds <- vegan::metaMDS(fungi_genus_sum_t[,!names(fungi_genus_sum_t) %in% c("group4")], distance = 'bray')

result <- nmds$points
result <- as.data.frame(cbind(result, rownames(result)))
write.csv(result,"bacteria_nmds_before.csv")
write.csv(result,"fungi_nmds_before.csv")

#获得上述的分类预测结果
predict_group <- c(train_predict, test_predict)
predict_group <- as.character(predict_group[rownames(result)])

#作图
colnames(result)[1:3] <- c('NMDS1', 'NMDS2', 'samples')
result$NMDS1 <- as.numeric(as.character(result$NMDS1))
result$NMDS2 <- as.numeric(as.character(result$NMDS2))
result$samples <- as.character(result$samples)
result <- cbind(result, predict_group)
head(result)
write.csv(result,"bacteria_nmds_before_group.csv")
write.csv(result,"fungi_nmds_before_group.csv")
result <- read.csv("bacteria_nmds_before_group2.csv",header = T,row.names = 1)
group6<- factor(result$group6,ordered=TRUE,levels=c("Bulk","Apiaceae",
                                                           "Burseraceae","Cruciferae",
                                                           "Cucurbitaceae","Gramineae",
                                                           "Lamiaceae","Leguminosae",
                                                           "Musaceae","Solanaceae"))


result<- read.csv("fungi_nmds_before_group2.csv",header = T,row.names = 1)
label2<- factor(result$group6,ordered=TRUE,levels=c("Bulk","Burseraceae","Cuculidae",
                                                          "Cruciferae","Gramineae", "Lamiaceae",
                                                          "Leguminosae","Malvaceae","Musaceae",
                                                          "SolaBulkceae","Taxodiaceae","Vitaceae"))
pdf("bacteria_nmds_before_group6.pdf",width=9,height=6)
pdf("fungi_nmds_before_group6.pdf",width=9,height=6)

ggplot(result, aes(NMDS1, NMDS2, color = predict_group)) +
  geom_polygon(data = plyr::ddply(result, 'predict_group', function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA) +
  geom_point()+
  xlab("NMDS1") +
  ylab("NMDS2")+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=0,size=15),
        axis.text = element_text(face = "bold",colour="black",size=15))
dev.off()

library(dplyr)
bacteria_asv_genus_top50<- read.csv("./bacteria_asv_genus_top50.csv",header = T,row.names = 1)
bacteria_asv_genus_top50_t<- t(bacteria_asv_genus_top50)
View(bacteria_asv_genus_top50_t)
ncol(bacteria_asv_genus_top50_t)
group <- read.csv("./bacteria_heatmap_group.csv",header = T,row.names = 1)
View(group)
bacteria_asv_genus_top50_t_group<- cbind(bacteria_asv_genus_top50_t,group)
View(bacteria_asv_genus_top50_t_group)
ncol(bacteria_asv_genus_top50_t_group)
datamean=group_by(bacteria_asv_genus_top50_t_group, group) %>% summarize_each(funs(mean))
datamean
write.csv(datamean,"bacteria_datamean.csv")
######FAST#### https://blog.csdn.net/woodcorpse/article/details/1065
library(ggplot2)
bacteria_importance<- read.csv("./bacteria_importance_genus.csv",header = T,row.names = 1)
bacteria_top50<- bacteria_importance[c(1:50),]
bacteria_top50_2<- bacteria_top50[order(bacteria_top50$MeanDecreaseAccuracy,decreasing = TRUE),]
x_label<- factor(rownames(bacteria_top50_2), levels=c("Unassigned","g__Rhizobium","g__Streptomyces","g__SC.I.84", 
                                                      "g__Vicinamibacteraceae","g__Nitrospira","g__Flavobacterium",
                                                      "g__Asticcacaulis","g__Dyadobacter","g__Pseudomonas","g__Luteolibacter",
                                                      "g__Burkholderia.Caballeronia.Paraburkholderia","g__MND1","g__TRA3.20",                                   
                                                      "g__Mucilaginibacter" ,"g__Dongia","g__Bacillus","g__RB41" ,"g__Gaiella","g__Gemmatimonas",                              
                                                      "g__Subgroup_7","g__Pedobacter","g__67.14","g__Massilia","g__WD2101_soil_group",
                                                      "g__Bradyrhizobium", "g__Nocardioides", "g__Haliangium", "g__Pseudolabrys","g__MB.A2.108" ,"g__Nitrosospira",
                                                      "g__Bryobacter", "g__Caulobacter", "g__Sphingomonas","g__Rokubacteriales", 
                                                      "g__S0134_terrestrial_group","g__KF.JG30.B3","g__Rubrobacter","g__CCD24",  "g__Candidatus_Solibacter",
                                                      "g__KD4.96","g__Subgroup_5", "g__BIrii41" ,"g__Reyranella","g__Skermanella","g__Chryseobacterium", 
                                                      "g__Acidibacter","g__Halomonas", "g__Nitrososphaeraceae","g__Novosphingobium"))
View(x_label)
library(forcats)
fct_rev(x_label)

library(rlist)
pdf("bacteria_top50_taxa.pdf",width=6,height=9)
ggplot(bacteria_top50_2,aes(fct_rev(x_label),MeanDecreaseAccuracy)) +
  geom_point(shape=1)+
  xlab("MeanDecreaseAccuracy") +
  ylab("")+
  coord_flip()+
  theme_set(theme_bw())+
  theme(legend.position = "none")+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=12),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=0,size=12),
        axis.text.y = element_text(vjust=1.1,hjust=1,angle=0,size=10),
        axis.text = element_text(face = "bold",colour="black",size=12))
dev.off()


fungi_importance<- read.csv("./fungi_importance_genus.csv",header = T,row.names = 1)
fungi_top60<- fungi_importance[c(1:60),]
fungi_top60_2<- fungi_top60[order(fungi_top60$MeanDecreaseAccuracy,decreasing = TRUE),]
x_label<- factor(rownames(fungi_top60_2), levels=c("g__Mortierella","g__Talaromyces","g__Podospora","g__Westerdykella","Unassigned",           
                                                  "g__Rhodosporidium","g__Neurospora","g__Gibberella","g__Mollisia","g__Ustilago",          
                                                   "g__Hypocrea","g__Ophiocordyceps","g__Nectria","g__Aspergillus","g__Leptodontidium",    
                                                   "g__Malassezia","g__Teberdinia","g__Dokmaia","g__Capronia","g__Lecythophora",      
                                                  "g__Eupenicillium","g__Sakaguchia","g__Waitea","g__Davidiella","g__Cryptococcus",      
                                                   "g__Toxicocladosporium","g__Archaeorhizomyces","g__Rhizopus","g__Glomerella","g__Stachybotrys",      
                                                   "g__Myceliophthora","g__Sarocladium","g__Microdochium", "g__Zopfiella","g__Clavulinopsis",     
                                                   "g__Oidiodendron","g__Umbelopsis","g__Emericellopsis","g__Gibellulopsis","g__Cephalotrichum",    
                                                   "g__Ceratobasidium","g__Phoma","g__Chalara","g__Khuskia","g__Phaeosphaeriopsis", 
                                                   "g__Didymella","g__Hygrocybe","g__Acremonium","g__Cyphellophora","g__Cortinarius",       
                                                   "g__Paecilomyces","g__Arthrobotrys","g__Humicola","g__Eurotium","g__Harzia",            
                                                   "g__Actinomucor","g__Pseudogymnoascus","g__Penicillium","g__Neonectria","g__Chaetomium"))

pdf("fungi_top60_taxa.pdf",width=6,height=9)
ggplot(fungi_top60_2, aes(fct_rev(x_label),MeanDecreaseAccuracy)) +
  geom_point(shape=1)+
  xlab("MeanDecreaseAccuracy") +
  ylab("")+
  coord_flip()+
  theme_set(theme_bw())+
  theme(legend.position = "none")+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=12),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=0,size=12),
        axis.text.y = element_text(vjust=1.1,hjust=1,angle=0,size=10),
        axis.text = element_text(face = "bold",colour="black",size=12))
dev.off()

#####heatmap######

bacteria_genus<- read.csv("./bacteria_asv_taxa_genus_sum.csv",header = T,row.names = 1)
View(bacteria_genus)
bacteria_importance<- read.csv("./bacteria_importance_genus.csv",header = T,row.names = 1)
bacteria_top50<- bacteria_importance[c(1:50),]
View(bacteria_top50)
bacteria_Asv_genus_top50<- bacteria_genus[rownames(data.frame(bacteria_top50)),]
View(bacteria_Asv_genus_top50)
setdiff(rownames(bacteria_top50),rownames(bacteria_Asv_genus_top50))
bacteria_Asv_genus_top50<- read.csv("bacteria_Asv_genus_top50.csv",header = T,row.names = 1)
View(bacteria_Asv_genus_top50)
bacteria_importance<- read.csv("./bacteria_importance_genus2.csv",header = T)
View(bacteria_importance)

pdf("bacteria_top50_taxa_barplot.pdf",width=7,height=9)
ggplot(bacteria_importance, aes(x=value,y=species)) +
  geom_bar(aes(fill = group),stat = "identity",
           position ="dodge",width=0.9)+
  xlab("Value") +
  ylab("")+
  theme_set(theme_bw())+
  theme(legend.position = c(0.8,0.9))+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=12),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text.x = element_text(vjust=0.9,hjust=0.8,angle=0,size=12),
        axis.text.y = element_text(vjust=1.1,hjust=1,angle=0,size=10),
        axis.text = element_text(face = "bold",colour="black",size=12))
dev.off()



library(dplyr)
bacteria_datamean<- read.csv("bacteria_datamean2.csv",header = T,row.names = 1)
View(bacteria_datamean)
data_matrix <- as.matrix(select(bacteria_datamean,-1))
data_matrix2 <- as.matrix(select(bacteria_datamean,c(-1,-2)))
df %>%
  group_by(team) %>%
  mutate(percent = points/sum(points))

data_matrix2 %>%
  mutate(percent = points/sum(points))

My_percnet <- function(x){
   for (i in 1:ncol(x)){
     print(mean(i))
     res <- rbind(x,round(x[,i]/sum(x[,i])*100,2))
   } 
  return(res)
}

data_matrix2_pre<- My_percnet(t(data_matrix2))
View(data.frame(data_matrix2_pre))

require(tidyverse)
df <- read.table("~/ATGC.txt", header = T)
df %>%
  mutate(Percentage=paste0(round(Num/sum(Num)*100,2),"%"))
My_percnet <- function(x){
  return(paste0(round(Num/sum(Num)*100,2),"%"))
}


ncol(bacteria_datamean)
ncol(data_matrix)
View(data_matrix)
group<- bacteria_datamean$group2
library(gplots)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(7,"BrBG"))(10)#换个好看的颜色
hM <- format(round(data_matrix, 2))#对数据保留2位小数

data_matrix3<-scale(t(data_matrix2),center = T)
View(data_matrix3)

heatmap(df, scale = "none")
col<- colorRampPalette(c("red","white","blue"))(256)
col <- colorRampPalette(brewer.pal(8,"Set1"))(20)

heatmap(df, scale = "none", col =  col)
group<- as.data.frame(bacteria_datamean$group2)
rownames(group) <- colnames(data_matrix3)

plot_color = c('chocolate','cyan4')[as.factor(group)]
plot_color
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

pheatmap(data_matrix3, 
         cluster_row = FALSE,
         cluster_col=FALSE,
         annotation_col=group,
         annotation_colors = ann_colors)


######fungi heatmap######

library(dplyr)
library(gplots)
library(RColorBrewer)
fungi_genus<- read.csv("./fungi_asv_taxa_genus_sum.csv",header = T,row.names = 1)
View(fungi_genus)
fungi_importance<- read.csv("./fungi_importance_genus.csv",header = T,row.names = 1)
fungi_top60<- fungi_importance[c(1:60),]
View(fungi_top60)
fungi_asv_genus_top60<- fungi_genus[rownames(fungi_top60),]
View(fungi_asv_genus_top60)
write.csv(fungi_asv_genus_top60,"fungi_asv_genus_top60.csv")
fungi_asv_genus_top60<- read.csv("./fungi_asv_genus_top60.csv",header = T,row.names = 1)
fungi_asv_genus_top60_t<- t(fungi_asv_genus_top60)
View(fungi_asv_genus_top60_t)
ncol(fungi_asv_genus_top60_t)
group <- read.csv("./fungi_heatmap_group.csv",header = T,row.names = 1)
View(group)
fungi_asv_genus_top60_t_group<- cbind(fungi_asv_genus_top60_t,group)
View(fungi_asv_genus_top60_t_group)
ncol(fungi_asv_genus_top60_t_group)
datamean=group_by(fungi_asv_genus_top60_t_group, group) %>% summarize_each(funs(mean))
datamean
write.csv(datamean,"fungi_datamean.csv")
fungi_datamean<- read.csv("fungi_datamean2.csv",header = T,row.names = 1)
View(fungi_datamean)

data_matrix <- as.matrix(select(fungi_datamean,-1))
data_matrix2 <- as.matrix(select(fungi_datamean,c(-1,-6)))
data_matrix3<-scale(t(data_matrix2),center = T)
View(data_matrix3)

col<- colorRampPalette(c("red","white","blue"))(256)
col <- colorRampPalette(brewer.pal(8,"Set1"))(20)
group<- fungi_datamean$group2
plot_color = c('chocolate','cyan4')[as.factor(group)]
plot_color
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
library(pheatmap)

group<- read.csv("fungi_heatmap_group2.csv",header=T,row.names=1)
ann_colors<-list(group=c(Bulk='#CC6666',Rhizophere='cyan4'))

rownames(group) <- colnames(data_matrix3)
View(group)
pheatmap(data_matrix3, 
                   cluster_row = FALSE,
                   cluster_col=FALSE,
                   annotation_col=group,
         annotation_colors = ann_colors)

anno_color<-list(Treatment=c(BL="blue",GL="green",RL="red",WL="grey"))

library(dendextend)
library(pheatmap)
Cor  = cor(dat,method = "pearson")

dend = hclust(dist(Cor))
dend = color_branches(dend, k = 2)


library(gplots)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#换个好看的颜色
hM <- format(round(Ca, 2))#对数据保留2位小数

heatmap.2(data_matrix3,
          trace="none",#不显示trace
          col=coul,#修改热图颜色
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv = F, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)






