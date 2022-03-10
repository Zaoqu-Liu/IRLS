rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(GSVA)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
load('/Users/CRC-data-lnc/TCGA.rda')
rm(CRC.lncRNA)
load('/Users/Tool and Data/cellMarker_ssGSEA.Rdata')

# -------------------------------------------------------------------------

gsva_data <- gsva(as.matrix(CRC.mRNA),cellMarker, method = "ssgsea")
ss <- gsva_data

# -------------------------------------------------------------------------

dir.create('ConsensusCluster/')
## 一致性聚类
results = ConsensusClusterPlus(as.matrix(ss),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               tmyPal = c('navy','darkred'),
                               title='ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")

icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'pdf')

## PAC = Proportion of ambiguous clustering 模糊聚类比例
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK

PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=4,shape=21,color='darkred',fill='orange')+
  ggtitle('Proportion of ambiguous clustering')+
  xlab('Cluster number K')+ylab(NULL)+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=13))
ggsave(filename = 'ConsensusCluster/PAC.pdf',width = 3.8,height = 4)

## 保存分型信息
clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)

head(sub)

my <- results[[2]][["ml"]]
library(pheatmap)
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
pheatmap(1-my,show_colnames = F,show_rownames = F,
         treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'complete',
         color = colorRampPalette(c("white","#C75D30"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_col = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_colors = list(Cluster=c('C2'='#B5739D','C1'='#4E8279')))
library(export)
graph2pdf(file='cluster2.pdf',width=5.5,height=4.5)

# -------------------------------------------------------------------------

ss2 <- merge(sub,t(ss),by.x=1,by.y=0)
ss2 <- pivot_longer(ss2,3:30,names_to = 'cell',values_to = 'value')

ggplot(ss2,aes(cell,value,fill=Cluster))+
  geom_boxplot(outlier.colour = NA)+
  stat_compare_means(label = 'p.signif')+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))

# -------------------------------------------------------------------------

TCGA_clin <- merge(sub,TCGA_clin,by=1)
my <- TCGA_clin[,c(1,2,7,8,12,13)]%>%column_to_rownames('Sample')
my$Age <- ifelse(my$Age>65,'>65','≤65')
table(my$Age)
my$Gender <- Hmisc::capitalize(my$Gender)
my <- my[order(my$Cluster,my$Age,my$Stage,my$Status,my$Gender),]

ee <- t(scale(t(ee)))
table(my$Cluster)

# -------------------------------------------------------------------------

my[is.na(my)] <- 'NA'
my$Age <- factor(my$Age,levels = c('≤65','>65'))
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Stage <- factor(my$Stage,levels = c('I','II','III','IV','NA'))
my$Cluster <- factor(my$Cluster)
my$Status <- factor(my$Status)

# -------------------------------------------------------------------------

Cluster <- c('#4E8279','#B5739D')
names(Cluster) <- levels(my$Cluster)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
table(my$Gender)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
table(my$Stage)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick','White')
names(Stage) <- levels(my$Stage)
table(my$Recurrence)
table(my$Status)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)

# -------------------------------------------------------------------------

Top = HeatmapAnnotation(Cluster=my$Cluster,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col=list(Cluster = Cluster,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))

Heatmap(ee,name='Z-score',
        top_annotation = Top,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c('#21b6af','white','#eeba4d')),#49b0d9
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,316),rep(2,268)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold")),
        column_gap = unit(2,'mm')
) 
library(export)
graph2pdf(file='cell-heatmap.pdf',width=9,height=5.5)

# -------------------------------------------------------------------------


ee <- as.data.frame(ss)[,rownames(my)]
tt <- cbind(Cluster=as.character(my$Cluster),t(ee))
tt <- as.data.frame(tt)
tt2 <- pivot_longer(tt,cols = 2:29,names_to = 'cell',values_to = 'value')
tt2$value <- as.numeric(tt2$value)

# -------------------------------------------------------------------------

source("GeomSplitViolin.R")

ggplot(tt2, aes(cell,value, fill = Cluster)) + 
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75), #画4分位线
                    trim = T, #是否修剪小提琴图的密度曲线
                    linetype = "solid", #周围线的轮廓
                    color = "black", 
                    size = 0.35,
                    na.rm = T,
                    position ="identity")+ #周围线粗细
  ylab("Relative Infiltration") + xlab(NULL) +
  scale_fill_manual(values = c('#21b6af','#eeba4d'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,size=10),
        legend.position=c(0.05,0.003),legend.justification = c(0,0),
        axis.title.y = element_text(size=12))
ggsave(filename = 'cell-boxplot.pdf',width = 8,height=4.3)

# -------------------------------------------------------------------------

save(TCGA_clin,ss,file = 'Cluster+ssgsea.Rda')


