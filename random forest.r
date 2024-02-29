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

#作图展示 top60 重要的 OTUs
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
library(splines) 
pdf("bacteria_asv_cross.pdf",width=9,height=6)

p <- ggplot(ASV_train.cv, aes(asvs, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p

#输出表格
write.table(importance_ASV, 'bacteria_importance_ASV.txt', sep = '\t', col.names = NA, quote = FALSE)
#读取表格


