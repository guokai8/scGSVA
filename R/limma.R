limmat<-function(obj,group=NULL){
    library(limma) ##limma做差异分析
    exprSet <- kegg
    group <- factor(meta,levels = c("3","6"),ordered = F)## 分组变成向量，并且限定leves的顺序， levels里面，把对照组放在前面
    design <- model.matrix(~group)# 构建比较矩阵
    colnames(design) <- levels(group)
    fit <- lmFit(exprSet,design)#线性模型拟合
    fit1 <- eBayes(fit)#贝叶斯检验
    allDiff=topTable(fit1,adjust='fdr',coef=2,number=Inf)
    write.table(allDiff,"kegg.txt",col.names=T,row.names=T,sep="\t")
    #根据allDiff找出你自己感兴趣的通路,例如
    up <- c("REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",'REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE','REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE','REACTOME_NEUROTOXICITY_OF_CLOSTRIDIUM_TOXINS','REACTOME_GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATIO')
    down <- c('REACTOME_CREB1_PHOSPHORYLATION_THROUGH_NMDA_RECEPTOR_MEDIATED_ACTIVATION_OF_RAS_SIGNALING','REACTOME_RAS_ACTIVATION_UPON_CA2_INFLUX_THROUGH_NMDA_RECEPTOR','REACTOME_LONG_TERM_POTENTIATION','REACTOME_UNBLOCKING_OF_NMDA_RECEPTORS_GLUTAMATE_BINDING_AND_ACTIVATION','REACTOME_CLEC7A_DECTIN_1_INDUCES_NFAT_ACTIVATIO')
    TEST <- c(up,down)
    p <- allDiff
    p$ID <- rownames(p)
    q <- p[TEST,]
    group1 <- c(rep("6",5),rep("3",5))
    df <- data.frame(ID = q$ID, score = q$t,group=group1 )
    # 按照score排序
    sortdf <- df[order(df$score),]
    sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列
    head(sortdf)
    p <- ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity') +
        coord_flip() +
        theme_bw() + #去除背景色
        theme(panel.grid =element_blank())+
        theme(panel.border = element_rect(size = 0.6))
    pdf("gsva_dif.pdf",width=12,height=8)
}

