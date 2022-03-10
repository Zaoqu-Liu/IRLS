rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(export)
library(ggsci)

load('sampleclin.rda')

# -------------------------------------------------------------------------

kmp <- function(data,legend,main){
  mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                             font.x = c(14,"plain", "black"),
                             font.y = c(14,"plain", "black"),
                             legend = "top")
  cut <- surv_cutpoint(data,'OS.time','OS','RS',minprop = 0.1)
  cat <- surv_categorize(cut)
  fit <- survfit(Surv(OS.time,OS)~RS,cat)
  pp <- ggsurvplot(fit,data = cat,
                   palette= pal_nejm()(2),
                   conf.int=FALSE,size=1.3,
                   pval=T,pval.method = T,
                   legend.labs=c('High-risk','Low-risk'), 
                   legend.title="",
                   legend=legend,
                   xlab="Time(years)",
                   ylab='Overall survival',
                   ggtheme = mytheme)
  return(ggpar(pp,main = main))
}

table(my$IC)
kmp(my,legend = c(0.84,0.97),main = NULL)
graph2pdf(file='All-patients.pdf',width=4.2,height=4)

# -------------------------------------------------------------------------

kmp <- function(data,legend,main){
  mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                             font.x = c(14,"plain", "black"),
                             font.y = c(14,"plain", "black"),
                             legend = "top")
  cut <- surv_cutpoint(data,'RFS.time','RFS','RS',minprop = 0.1)
  cat <- surv_categorize(cut)
  fit <- survfit(Surv(RFS.time,RFS)~RS,cat)
  pp <- ggsurvplot(fit,data = cat,
                   palette= pal_nejm()(2),
                   conf.int=FALSE,size=1.3,
                   pval=T,pval.method = T,
                   legend.labs=c('High-risk','Low-risk'), 
                   legend.title="",
                   legend=legend,
                   xlab="Time(years)",
                   ylab='Recurrence-free survival',
                   ggtheme = mytheme)
  return(ggpar(pp,main = main))
}

kmp(my,legend = c(0.84,0.97),main = NULL)
graph2pdf(file='All-patients-rfs.pdf',width=4.2,height=4)

#kmp(my[my$IC==1,],legend = c(0.84,0.97),main = 'Patients with ICI (n = 65)')
#graph2pdf(file='Patients with immunotherapy-rfs.pdf',width=4.2,height=4)

#kmp(my[my$IC==0,],legend = c(0.84,0.97),main = 'Patients without ICI (n = 167)')
#graph2pdf(file='Patients without immunotherapy-rfs.pdf',width=4.2,height=4)

# -------------------------------------------------------------------------

x <- summary(coxph(Surv(OS.time,OS)~.,my[,c(1,2,5:13)]))
y <- data.frame(id=rownames(x$coefficients),
                HR=x$coefficients[,2],
                HR.95L=x$conf.int[,"lower .95"],
                HR.95H=x$conf.int[,'upper .95'],
                pvalue=x$coefficients[,"Pr(>|z|)"])

y <- y[c(1:6,10,7:9),]
y[,-1] <- apply(y[,-1],2,as.numeric)

rt <- y
gene <- rt$id
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr," (",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

if(T){
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  #绘制森林图左边的基因信息
  xlim = c(0,2.6)
  par(mar=c(4,2.5,2,0))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=1
  text(0.5,n:1,gene,adj=0,cex=text.cex)
  text(1.65,n:1,pVal,adj=1,cex=text.cex);
  text(1.5+0.2,n+1,'P-value',cex=text.cex,font=2,adj=1)
  text(2.6,n:1,Hazard.ratio,adj=1,cex=text.cex)
  text(2.45,n+1,'HR (95% CI)',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0,6)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=0,code=3,length=0.05,col="black",lwd=2.5)
  abline(v=1,col="gray",lty=2,lwd=1.5)
  boxcolor = '#67AB9F'
  points(as.numeric(hr), n:1, pch = 15, cex=2,col = boxcolor)
  axis(1)
}
library(export)
graph2pdf(file="inhouse-multicox.pdf", width = 10,height = 4.7)

# -------------------------------------------------------------------------

x <- summary(coxph(Surv(RFS.time,RFS)~.,my[,c(3,4,5:13)]))
y <- data.frame(id=rownames(x$coefficients),
                HR=x$coefficients[,2],
                HR.95L=x$conf.int[,"lower .95"],
                HR.95H=x$conf.int[,'upper .95'],
                pvalue=x$coefficients[,"Pr(>|z|)"])
y <- y[c(1:6,10,7:9),]

rt <- y
gene <- rt$id
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr," (",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

if(T){
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  #绘制森林图左边的基因信息
  xlim = c(0,2.6)
  par(mar=c(4,2.5,2,0))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=1
  text(0.5,n:1,gene,adj=0,cex=text.cex)
  text(1.65,n:1,pVal,adj=1,cex=text.cex);
  text(1.5+0.2,n+1,'P-value',cex=text.cex,font=2,adj=1)
  text(2.6,n:1,Hazard.ratio,adj=1,cex=text.cex)
  text(2.45,n+1,'HR (95% CI)',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0,5)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=0,code=3,length=0.05,col="black",lwd=2.5)
  abline(v=1,col="gray",lty=2,lwd=1.5)
  boxcolor = '#67AB9F'
  points(as.numeric(hr), n:1, pch = 15, cex=2,col = boxcolor)
  axis(1)
}
graph2pdf(file="inhouse-multicox-rfs.pdf", width = 10,height = 4.7)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

tt <- timeROC(my$OS.time,my$OS,my$RS,cause = 1,weighting = 'marginal',times = c(1,3,5),ROC = T)
tp <- tt$TP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'tp')
fp <- tt$FP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'fp')

dd <- tp
dd$fp <- fp$fp
dd$time <- ifelse(dd$time=='t=1',"1-Year = 0.840",
                  ifelse(dd$time=='t=3','3-Year = 0.776','5-Year = 0.818'))

ggplot(dd,aes(fp,tp,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_nejm()+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
ggsave(filename = 'inhouse-timeROC.pdf',width = 4.3,height = 4)

# -------------------------------------------------------------------------

fit <- roc(x$RR,x$RS,auc=T)

dd <- data.frame(x=1-fit$specificities,y=fit$sensitivities)
dd <- dd[order(dd$y,dd$x),]
dd$AUC <- 'IRLS'
tmp <- dd
ggplot(tmp,aes(x,y))+
  geom_line(aes(group=AUC))

tmp$AUC2 <- ifelse(tmp$AUC=='IRLS','IRLS= 0.897',ifelse(tmp$AUC=='PD-L1','PD-L1 = 0.686 ***','CD8A= 0.725 **'))

tmp$AUC2 <- factor(tmp$AUC2,levels = unique(tmp$AUC2))
ggplot(tmp,aes(x,y,color=AUC2))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_nejm()+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
ggsave(filename = 'IC-inhouse-timeROC.pdf',width = 4.3,height = 4)



if(F){
  plot.roc(fit,print.thres = T,print.auc = T)
  pdf(file = 'IC-roc.pdf',width = 4,height = 4)
  plot.roc(fit,
           axes=T, ## 是否显示xy轴
           legacy.axes=T,
           main=NULL, ## Title
           col= "steelblue", ## 曲线颜色
           lty=1, ## 曲线形状
           lwd=3, ## 曲线粗细
           identity=T, ## 是否显示对角线
           identity.col="grey60", ## 对角线颜色
           identity.lty=2, ## 对角线形状
           identity.lwd=2, ## 对角线粗细
           print.thres=F, ## 是否输出cut-off值
           print.thres.pch=20, ## cut-off点的形状
           print.thres.col="red", ## cut-off点和文本的颜色
           print.thres.cex=1.2, 
           print.auc=T, ## 是否显示AUC
           print.auc.pattern="AUC = 0.896", ## 展示AUC的格式
           auc.polygon.border='darkred',
           print.auc.x=0.48, ## AUC值的X位置
           print.auc.y=0.13, ## AUC值的Y位置
           print.auc.cex=1.2, ## AUC值的放大倍数
           print.auc.col='black', ## ACU值的颜色
           auc.polygon=TRUE, ## 是否将ROC下面积上色
           auc.polygon.col='skyblue', 
           max.auc.polygon=TRUE,
           max.auc.polygon.col='WhiteSmoke',
           max.auc.polygon.lty=1
  )
  dev.off()
}
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

library(plyr)
ggdata <- table(x$gg,x$RR)%>%as.data.frame()
ggdata2 <- ddply(ggdata,'Var1',transform,percent_Freq=Freq/sum(Freq)*100)
ggdata2 <- ddply(ggdata2,'Var1',transform,lable=cumsum(percent_Freq)-0.5*percent_Freq)
ggdata2$ll <- paste0(round(ggdata2$percent_Freq/100,2)*100,'%')
ggdata2$Var2 <- factor(ifelse(ggdata2$Var2==0,'NR','R'),levels = c('R','NR'))

ggplot(ggdata2,aes(Var1,percent_Freq,fill=Var2))+
  geom_bar(stat = 'identity',width = 0.85)+
  xlab(NULL)+ylab('Fraction (%)')+
  geom_text(aes(label=ll),size=3.8,
            position = position_stack(vjust = 0.5),
            color='white')+
  scale_fill_manual(values = pal_npg(alpha = 0.9)(2)[1:2])+
  theme_classic()+
  scale_x_discrete(expand = c(0.3,0.2))+
  scale_y_continuous(expand = c(0.01,0.0))+
  ggtitle('****')+
  theme(legend.position = 'right',
        plot.title = element_text(hjust=0.5,face = 'plain'),
        axis.title = element_text(size=13),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=12))+
  labs(fill=NULL)
ggsave(filename = 'IC-bar.pdf',height = 4,width = 3.2)

chisq.test(table(x$gg,x$RR))

# -------------------------------------------------------------------------

ggplot(x2,aes(reorder(ID,RS),RS,fill=Res))+
  geom_bar(stat = 'identity',width = 0.7,color='grey30',position = position_dodge2(width = 0.9))+
  scale_fill_nejm()+labs(y='IRLS score')+
  theme_classic(base_rect_size = 2)+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.01,1),
        legend.justification = c(0,1))

# -------------------------------------------------------------------------

library(compareC)
tt <- my
dd <- data.frame()
for (i in colnames(my)[5:13]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  p <- compareC(tt$OS.time,tt$OS,tt$RS,tt[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}

ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('C-index (Compared with IRLS)')+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.85))+
  geom_text(aes(y=0.80,label=ll),size=5)












