setwd("C:/Users/shi lab/Desktop/")
vec<-read.table("snps.all_filtered_plink.eigenvec",header=TRUE,sep=" ")
pdf("pca.pdf",width=4,height = 4)
ggplot(pca, aes(x=PC1, y=PC2, shape=population, colour=species))+
# 点的大小
    geom_point(size=2)+
# 点的形状，根据群体数目选择要用多少个、哪些点
    scale_shape_manual(values=c(1:15))+
# 点的颜色
    scale_colour_manual(values=c("#47B8E0","#FF7473","#FFC952"))+
# 不要背景的灰色
    theme_bw()+
# 标题
    ggtitle("pca")+
# panel.grid.major,主网格线格式，类似于excel里面那些主要的分隔
# panel.grid.minor,副网格线格式
# 可以通过panel.grid.(major/minor).(x/y)分别控制（x/y）轴方向（主/副）网格线格式
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0.5,size=rel(1.5),family="Times",face="bold.italic",color = "red"))
dev.off()


# 也可以通过ggthemeassist包来调整
library(ggThemeAssist)
plot<-ggplot(pca, aes(x=PC1, y=PC2, shape=population, colour=species))+geom_point(size=5)+scale_shape_manual(values=c(1:15))+
    +     scale_colour_manual(values=c("#47B8E0","#FF7473","#FFC952"))
ggThemeAssistGadget(plot)