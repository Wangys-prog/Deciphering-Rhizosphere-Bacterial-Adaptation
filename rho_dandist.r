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
####图片组合#######
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
