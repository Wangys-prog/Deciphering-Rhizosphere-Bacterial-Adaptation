library(vegan)
library(ape)
library(tidyverse)
library(ggplot2)

###bacteria
bacteria<-read.csv("bacteria.csv",header=T,row.names = 1)
View(bacteria)
colnames(bacteria)
rownames(bacteria)
ncol(bacteria)

new_d2 <- scan("bacteria_id_select.csv",what = "character")
View(new_d2)
##选择列
newdata2<-bacteria[ ,colnames(bacteria) %in% new_d2]


rownames(newdata2)
colnames(newdata2)

ncol(newdata2)
nrow(newdata2)

bacteria_select<-data.frame(newdata2)

bacteria_select_t<-t(bacteria_select)
nrow(bacteria_select_t)
bacteria_dist <- vegdist(bacteria_select_t, "bray")
bacteria_map<-  read.csv("bacteria_group2.csv",header = T,row.names=1)

#计算mnds矩阵距离

nmds_dis <- metaMDS(bacteria_dist, k = 2)

#应力函数值，一般不大于 0.2 为合理
nmds_dis$stress
# nmds_dis$stress
# 0.001747191

#样方得分

nmds_dis_site <- data.frame(nmds_dis$points)

nmds_dis_site$name <- rownames(nmds_dis_site)

map <- read.csv("bacteria_group2.csv",header = T,row.names=1)

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
		

#进行Anosim分析
library(vegan)
anosim=anosim(bacteria_select_t, group4, permutations=999)
summary(anosim)

####adonis分析####
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
