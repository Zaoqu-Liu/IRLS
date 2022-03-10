rm(list = ls())
library(estimate)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ImmuLncRNA)
load('/Users/CRC-data-lnc/TCGA.rda')
load('/Users/CRC-data-lnc/GEO-mRNA-data/GEO-Rdata/TCGA.rda')
rm(TCGA_clin)
load('Purity.Rdata')
CRC.lncRNA <- CRC.lncRNA[apply(CRC.lncRNA,1,function(x){sum(x>0)>ncol(CRC.lncRNA)/3}),]

# -------------------------------------------------------------------------

set.seed(123456)
system.time(
  result <- immu.LncRNA(CRC.mRNA,CRC.lncRNA,adjusted=TRUE,Tumour_purity,pathways,)
)

sig <- subset(result$fgseaRes_all,padj < 0.05 & sigValue > 0.995)
sig_immune_related_lncRNA <- unique(sig$lncRNA)
save(result,sig,sig_immune_related_lncRNA,file = 'immune_related_result.Rdata')

ggplot(bardata,aes(log2(Freq+1),reorder(Var1,Freq),fill=Var1))+
  geom_bar(stat = 'identity',width = 0.7)+
  theme_classic(base_rect_size = 1.5)+
  ylab('')+xlab('log2(Number of immune-related lncRNAs)')+
  ggtitle('ImmuLnc analysis')+
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        axis.title.x = element_text(size = 10),
        legend.position = 'none')+
  scale_y_discrete(expand = c(0.04,0),)+
  scale_x_continuous(expand = c(0,0.1))+
  scale_fill_manual(values = c(pal_npg(alpha = 0.8)(10),pal_d3(alpha = 0.8)(10))[-c(17:20)])
ggsave(filename = 'Immune-lncRNA-type.pdf',width = 5,height = 5)


