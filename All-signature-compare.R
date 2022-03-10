rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(compareC)
library(ggsci)
load('all-signature.rda')

uni <- data.frame()
for (i in names(tmp2)) {
  dd <- tmp2[[i]]
  for(j in colnames(dd)[4:ncol(dd)]){
    scox <- summary(coxph(Surv(OS.time,OS)~get(j),dd))
    p <- compareC(dd$OS.time,dd$OS,dd$IRLS,dd[,j])$pval
    uni <- rbind(uni,data.frame(ID=i,A=j,
                                HR=scox$conf.int[,1],
                                HR.95L=scox$conf.int[,3],
                                HR.95R=scox$conf.int[,4],
                                pvalue=scox$coefficients[,5],
                                cindex=scox$concordance[1]%>%as.numeric(),
                                cse=scox$concordance[2]%>%as.numeric(),
                                CP=p))
  }
}

unique(uni$ID)
dd <- uni[uni$ID=='TCGA-CRC',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[1])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="TCGA-CRC")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.89,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.7,0.9),limits = c(0.4,0.94))

# -------------------------------------------------------------------------

unique(uni$ID)
dd <- uni[uni$ID=='GSE17536',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[2])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE17536")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.71,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.6,0.7))

# -------------------------------------------------------------------------

unique(uni$ID)
dd <- uni[uni$ID=='GSE17537',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL


ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[3])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE17537")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.85,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.4,0.6,0.8),limits = c(0.38,0.88))

# -------------------------------------------------------------------------

unique(uni$ID)
dd <- uni[uni$ID=='GSE29621',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL


ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[4])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE29621")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.79,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.4,0.6,0.8),limits = c(0.36,0.82))


unique(uni$ID)
dd <- uni[uni$ID=='GSE38832',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL
ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[5])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE38832")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.78,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.4,0.6,0.8))

unique(uni$ID)
dd <- uni[uni$ID=='GSE39582',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[6])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE39582")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.69,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.6,0.7))


unique(uni$ID)
dd <- uni[uni$ID=='GSE72970',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[7])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE72970")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.69,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.6,0.7))

unique(uni$ID)
dd <- uni[uni$ID=='Meta-Cohort',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))
rownames(dd) <- NULL

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(10)[8])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="Meta-Cohort")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.69,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.6,0.7))

# -------------------------------------------------------------------------

uni$x <- ifelse(uni$pvalue<0.05&uni$HR>1,'Risky',
                ifelse(uni$pvalue<0.05&uni$HR<1,'Protective','P > 0.05'))
xx <- pivot_wider(uni[,c(1,2,10)],names_from = 'ID',values_from = 'x')
xx <- column_to_rownames(xx,'A')

library(ComplexHeatmap)
library(circlize)

ss <- read.csv('signature.csv')[,3:4]
ss <- distinct(ss,Author,.keep_all = T)
ss <- ss[ss$Author%in%rownames(xx),]
ss[107,] <- c('lncRNA','IRLS')
rownames(ss) <- ss$Author
ss <- ss[order(ss$Type,ss$Author),]
rownames(ss) <- NULL
xx <- xx[ss$Author,]

right = HeatmapAnnotation(Signature=ss$Type,
                          annotation_legend_param=list(labels_gp = gpar(fontsize = 12),border = T,
                                                       title_gp = gpar(fontsize = 12,fontface = "bold"),
                                                       ncol=1),
                          border = T,which = 'column',
                          col=list(Signature = c('lncRNA'=pal_npg(alpha = 0.8)(10)[3],'mRNA'=pal_nejm(alpha = 0.8)(8)[3])),
                          show_annotation_name = F,
                          annotation_name_gp = gpar(fontsize = 12))
Heatmap(t(xx),col = c('#F2F2F2',pal_npg()(2)[2],pal_npg()(2)[1]),name = 'Type',
        rect_gp = gpar(col='grey70'),
        top_annotation = right,
        row_names_side = 'left',
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 12), border = T,
                                  title_gp = gpar(fontsize = 12, fontface = "bold")))



