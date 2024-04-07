# Editable figures 
# Updated on 25 Feb
# Yanmei Zhang
library('RamEx')
library('RamanEx')
library_all()
options(warn=1)
library('magrittr')

setwd('D:\\Scientific research\\合作项目\\张宇星\\文章撰写')
# ---------------- Fig 2 --------------------
# Fig 2.C
path <- 'D:\\Scientific research\\合作项目\\张宇星\\文章撰写\\原Fig 2D'
data_2 <- read.spec(path, group.index = 1, cutoff = c(500,3150))
data_2 %<>%  pre.smooth(.,m = 0, p = 5, w = 11, delta.wav = 2)  %>% BG(., cell.index = 2, cal_mean = FALSE)
data_2 %<>%  pre.baseline %>% pre.normalize(., method = 'CH') %>% BG(., cell.index = 2, cal_mean = TRUE)

name_group <- c('type','cell')
data_2 <- folder_reader(path, name_group)
data_2 <- S_G_smooth(data_2)
data_2_sub <- tapply(data_2, data_2$type, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
data_2_sub <- do.call(rbind, list.map(data_2_sub,as.matrix(unlist(.))))
data_2$spc <- data_2_sub
means <- aggregate(data_2$spc, by=list(data_2$type, data_2$cell), mean)  # 五合一
colnames(means) <- c('type', 'cell', colnames(means)[-c(1:2)])
data_2 <- new("hyperSpec",data = means[,1:2],spc = means[,-c(1:2)],wavelength = data_2@wavelength)
data_2 <- data_2[data_2$cell!= 'BG']  # 继续预处理
data_2 <- Raman_baseline(data_2)
data_2 <- Normalization(data_2,from = 2500,to = 3000,nor_cal = "max")
data_2$type %<>% gsub('PANC1','P0',.) %>% gsub('P2G8','P1',.) %>% gsub('P2D6','P2',.) %>% gsub('P1G8','P3',.) %>% gsub('P1E6','P4',.) %>% as.factor()
p <- mean.spec(data_2$spc, group = data_2$type, gap=0.3)
p 
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 5, height = 4.5, append = F, title = 'Fig 2.C') 
# Fig 2.D
library(ggplot2)
library(ggtext)
library('ggsci')
which(data_2@wavelength>500 & data_2@wavelength <1800)
set.seed(1334)
data.red <- data.frame(Rtsne::Rtsne(data_2$spc[,1:633],dims = 2, perplexity = 5, theta = 0.5,verbose = FALSE,max_iter = 10000)$Y)
colnames(data.red) <- c('tSNE 1', 'tSNE 2')
p <- ggplot(data.red, aes(`tSNE 1`, `tSNE 2`, color = as.factor(data_2$type))) +
  geom_point() +
  theme_bw() +
  stat_ellipse(level = 0.8) +
  labs(x = 'tSNE 1', y = 'tSNE 1') +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 15, family = "myFont"),
    legend.background = element_blank(), text = element_text(color = "black"),
    axis.text.x = element_text(size = 15, angle = 0, family = "myFont"),
    axis.text.y = element_text(size = 15, family = "myFont"),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.4, "lines"),
    axis.title = element_text(family = "myFont", size = 15)
  )
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 2.D') 

# ---------------- Fig 3 --------------------
# Fig 3.B ROC
library(multiROC)
load("D:/Scientific research/合作项目/张宇星/8.16分析/RandomForest_train_dataset.Rdata")
load("D:/Scientific research/合作项目/张宇星/8.16分析/Randomforest_3groups_0815.Rdata")
label_test %<>% gsub('ASPC1','high',.) %>% gsub('BXPC3','medium',.) %>% gsub('PANC1','medium',.) %>% gsub('MIAPACA2','low',.) %>% factor(.,levels=c('high','medium','low'))
pred_test <- predict(model.rf, data_test, type='prob') %>% data.frame
colnames(pred_test) <- paste(colnames(pred_test),'_pred')
encoder_test <- model.matrix(~ type -1,data=label_test)
colnames(encoder_test) <- c('high_true','medium_true','low_true')
colnames(pred_test) <- c('high_pred_rf','low_pred_rf','medium_pred_rf')
roc <- multi_roc(cbind(encoder_test, pred_test))
plot_roc <- plot_roc_data(roc)
plot_roc <- plot_roc[plot_roc$Group!='Macro',]
plot_roc$type <- 'single'
plot_roc$type[plot_roc$Group=='Micro'] <- 'Micro'
plot_roc$Group %<>% factor(.,levels=c('high','medium','low','Micro'))
cols <- c('#F8766D','#00BA38','#619CFF','grey40')

p <- ggplot(plot_roc, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), linewidth=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linewidth=1,
               colour='grey', linetype = 'dashed') +
  scale_color_manual(values = cols)+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        panel.grid = element_blank(),
        legend.background = element_rect(fill=NULL, linewidth=0.5, 
                                         linetype="solid",colour ="black")
  )
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 3.B') 


# Fig 3.D
which(test_ASPC1@wavelength > 500 & test_ASPC1@wavelength < 1800)
pred.1 <- predict(model.rf, test_ASPC1$spc[,1:628])
test_ASPC1$strain %<>% gsub('ASPC1-','ASPC1-hnRNPAB-',.) %>% gsub('PANC1-','PANC1-hnRNPLL-',.)
pred.1 %<>% factor(.,levels = c('high','medium','low'))
p <- pre.plot( test_ASPC1$strain, pred.1)
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 3.D') 
# Fig 3.E
pred.2 <- predict(model.rf, test_PANC1$spc[,1:628])
test_PANC1$strain %<>% gsub('P2G8','P1',.) %>% gsub('P2D6','P2',.) %>% gsub('P1G8','P3',.) %>% gsub('P1E6','P4',.) %>% as.factor()
pred.2 %<>% factor(.,levels = c('high','medium','low'))
p <- pre.plot( test_PANC1$strain, pred.2)
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 3.E') 

# Fig ROC Curve
library('multiROC')
test_label <- 
true_label <- dummies::dummy(test_df$Species, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)

roc_res <- multi_roc(final_df, force_diag=F)
pr_res <- multi_pr(final_df, force_diag=F)

plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)
# ---------------- Fig 4 --------------------
# Fig 4.A

# Fig 4.B

# Fig 4.D

# ---------------- Fig 5 --------------------

# Fig 5.C
path <- 'D:\\Scientific research\\合作项目\\张宇星\\10.01分析\\ASPC1_NC_cellculture_mousemodel\\Mouse_tissue'
# Old pipeline
name_group <- c('type','group','cell')
test_mouse <- folder_reader(path, name_group)
test_mouse$filename <- list.files(path, pattern = "*.txt", full.names = TRUE, include.dirs = T, recursive = T)
test_mouse$mouse <- str_split_i(test_mouse$filename, '/', 2)
test_mouse <- S_G_smooth(test_mouse)
test_mouse_sub <- tapply(test_mouse, test_mouse$mouse, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
test_mouse_sub <- do.call(rbind, list.map(test_mouse_sub,as.matrix(unlist(.))))
test_mouse$spc <- test_mouse_sub
means <- aggregate(test_mouse$spc, by=list(test_mouse$mouse,test_mouse$group, test_mouse$cell), mean)  # 五合一
colnames(means) <- c('mouse','group', 'cell', colnames(means)[-c(1:3)])
test_mouse <- new("hyperSpec",data = means[,1:3],spc = means[,-c(1:3)],wavelength = test_mouse@wavelength)
test_mouse <- test_mouse[test_mouse$cell!= 'BG']  # 继续预处理
test_mouse <- Raman_baseline(test_mouse)
test_mouse <- Normalization(test_mouse,from = 2500,to = 3000,nor_cal = "max")
data_mouse <- test_mouse$spc[,6:633]
colnames(data_mouse) <- colnames(data_train)

pred.1 <- predict(model.rf, data_mouse)
pred.1 %<>% factor(.,levels = c('high','medium','low'))
p <- pre.plot(test_mouse$mouse, pred.1)
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 5.C') 

# Fig S.1
pred.2 <- predict(model.rf, data_mouse)
pred.2 %<>% factor(.,levels = c('high','medium','low'))
p <- pre.plot(test_set$mouse, pred.2)
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 8, append = T, title = 'Fig S.1') 

# Fig 5.C Mean
data_mouse <- aggregate(data_mouse, by=list(test_set$mouse, test_set$group), mean)
pred.1 <- predict(model.rf, data_mouse[,-c(1:2)])
pred.1 %<>% factor(.,levels = c('high','medium','low'))
p <- pre.plot(data_mouse[,2], pred.1)
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 4.5, height = 4, append = T, title = 'Fig 5.C.mean') 

# ---------------- Fig S --------------------
# Fig S.2
load("D:/Scientific research/合作项目/张宇星/8.16分析/Randomforest_3groups_0815.Rdata")
importance <- data.frame(wave=as.numeric(colnames(data_train)), intensity=colMeans(data_train), gini=model.rf$importance)
names(importance) <- c('wave','intensity','gini')
data_box <- data.frame(data_train[, colnames(data_train) %in% importance$wave[importance$gini>0.7]])
data_box$label <- label_train
data_box$label %<>% factor(.,levels=c('high','medium','low'))
data_box <- setNames(melt(data_box),c('label','wave','values'))

p <- ggplot(data_box, aes(x = label, y = values, group = label, fill = label)) +  #[I_CH$color != 'batch(1-2)',]
  theme_bw() + geom_point(aes(group = label, color = label), position = position_jitter(0.2)) +
  geom_boxplot(aes(group = label), alpha = 0.2, outlier.alpha = 0) + 
  # geom_signif(comparisons=list(c('ASPC1','BXPC3'),c('BXPC3','PANC1'),c('PANC1','MIAPACA2')),map_signif_level = T,y_position = c(1.05,1.03,1))+
  # geom_hline(aes(yintercept = 0), colour = "red", linetype = "dashed", size = 1) +
  facet_wrap(~wave, scales = 'free')+
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_markdown(size = 15), 
        legend.text = element_markdown(size = 15), 
        legend.position = "none",
        legend.background = element_blank(), text = element_text(color = "black"),
        axis.title.y = element_text(size = 20, angle = 90), 
        axis.text.x = element_text(size = 8, angle = 0), 
        axis.text.y = element_text(size = 8), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"), 
        axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 20)) + 
  xlab("Group") + ylab("Raman Intensity (2937)")
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 8, height = 7, append = T, title = 'Fig S.2') 

# ---------------- Fig 6 --------------------
# Fig 6.A
importance <- data.frame(wave=as.numeric(colnames(data_train)), intensity=colMeans(data_train), gini=model.rf$importance)
names(importance) <- c('wave','intensity','gini')
p <- ggplot(importance, aes(wave,gini))+
  geom_line(aes(color=gini),linewidth=1) + 
  geom_hline(yintercept =0.75, linetype='dashed', color='red', linewidth=1)+
  # scale_color_gradient2(low = "blue", mid = "lightblue", midpoint = 0)+
  labs(y = "Normalized Intensity (a.u.)") +
  xlab(expression(paste("Wavenumber (cm"^{ -1 }, ")"))) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(600,800, 1000,1200,1400,1600,1800)
  )+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 6, height = 4, append = T, title = 'Fig 6.A') 

# Fig 6.2

# Fig 6.C
library(readxl)
library(stringr)
library(reshape2)
lipids <- read_excel('D:\\Scientific research\\合作项目\\张宇星\\Lipid analysis\\张宇星-人细胞-BQ-ZYX20231215-LC-RXB.xlsx')[,c(1,2,6:17)]
lipids$C_chain <- as.character(str_split(lipids$CompoundName, pattern = ' ', simplify = T)[,2]) 
lipids$C_chain <- as.character(str_split(lipids$C_chain, pattern = '-FA', simplify = T)[,1]) 
aa <-str_split(lipids$C_chain, pattern = '-|:', simplify = T)
aa[aa==''] <- 0
lipids$tol_C <- as.numeric(aa[,1])+as.numeric(aa[,3])
lipids$tol_d <- as.numeric(aa[,2])+as.numeric(aa[,4])

lipids_1 <- lipids
lipids_1[,3:14] <- t(t(lipids_1[,3:14]) / colSums(lipids_1[,3:14]))
lipids_1 <- reshape2::melt(lipids_1, id.vars = c('Class', 'CompoundName','C_chain','tol_C','tol_d'))
lipids_1$group <- lipids_1$variable %>% gsub('_1','',.) %>% gsub('_2','',.) %>% gsub('_3','',.) 
lipids_1$group %<>% gsub('ASPC1','High',.) %>% gsub('BXPC3','Medium',.) %>% gsub('PANC1','Medium',.) %>% gsub('MIAPACA2','Low',.) %>% factor(., levels = c('High','Medium','Low'))

# BMP, FFA, DAG, TAG, SM, PE-O, PE-P, PE, PC, 'Cer','DCER','LPC','LPG','LPI',''CE','LCER'
lipids_1$Class <- str_split_i(lipids_1$CompoundName, pattern = ' ',1)
l <- 'HexCer'
# S3
lipids_temp <- lipids_1[lipids_1$Class==l,]
lipids_temp$CompoundName <- str_split_i(lipids_temp$CompoundName, pattern = ' ',2)
p <- ggplot(lipids_temp, aes(x=CompoundName, y=value, group=group,fill=group))+ theme_bw()+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun.data='mean_sd',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+
  facet_wrap(~Class,ncol=1,scales = 'free')+
  # ylim(0,0.012)+
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_markdown(size = 15), 
        legend.text = element_markdown(size = 15), 
        legend.background = element_blank(), text = element_text(color = "black"),
        axis.title.y = element_text(size = 20, angle = 90), 
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 15), 
        strip.background = element_rect(fill = "white"), 
        axis.ticks = element_line(linewidth = 1), 
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15))
p
ggsave(str_c(l,'.png'), width = 5, height = 5,limitsize = FALSE)
eoffice::topptx(figure = p, filename = "Lipids_figures.pptx", width = 8, height = 5, append = T, title = str_c('Fig 6.C.1-',l)) 

# 6.C
library(ggpubr)
lipids_temp <- lipids_1[lipids_1$Class==l,]
lipids_temp <- aggregate(lipids_temp$value, by=list(lipids_temp$tol_d, lipids_temp$variable, lipids_temp$group), sum)
colnames(lipids_temp) <- c('tol_d','variable','group','value')
lipids_temp$tol_d %<>% as.factor()
p <- ggplot(lipids_temp, aes(x=tol_d, y=value, group=group,fill=group))+ theme_bw()+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun.data='mean_sd',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_markdown(size = 15), 
        legend.text = element_markdown(size = 15), 
        legend.background = element_blank(), text = element_text(color = "black"),
        axis.title.y = element_text(size = 20, angle = 90), 
        axis.text.x = element_text(size = 15, angle = 90), 
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"), 
        axis.ticks = element_line(linewidth = 1), 
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15))
p
ggsave(str_c(l,'_degrees.png'), width = 6, height = 6,limitsize = FALSE)
eoffice::topptx(figure = p, filename = "Lipids_figures.pptx", width = 5, height = 5, append = T, title = str_c('Fig 6.C.1-',l)) 

p <- ggplot(lipids_temp, aes(x=group, y=value, group=group,fill=group))+ theme_bw()+
  stat_compare_means(method = "t.test",aes(group=tol_d,label =..p.signif..),comparisons=list(c('High','Medium'),c('Medium','Low'),c('High','Low')))+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun.data='mean_sd',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+
  facet_wrap(~tol_d, ncol=12)+
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_markdown(size = 15), 
        legend.text = element_markdown(size = 15), 
        legend.background = element_blank(), text = element_text(color = "black"),
        axis.title.y = element_text(size = 20, angle = 90), 
        axis.text.x = element_text(size = 15, angle = 90), 
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"), 
        axis.ticks = element_line(linewidth = 1), 
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15))
p
ggsave(str_c(l,'_degrees_2.png'), width = 6, height = 7,limitsize = FALSE)
eoffice::topptx(figure = p, filename = "Lipids_figures.pptx", width = 5, height = 6, append = T, title = str_c('Fig 6.C.1-',l)) 



# Fig S.1



# Fig 6.C
select_lipid <- c('BMP','CER','DCER','LPC','LPG','LPI','PC','PE','CE','LCER')
lipids_1 <- lipids[,c(1,3:14)]
lipids_1[,2:13] <- t(t(lipids_1[,2:13]) / colSums(lipids_1[,2:13]))
lipids_1 <- aggregate(lipids_1[,2:13], by=list(lipids_1$Class),sum)
lipids_1 <- melt(lipids_1, id.vars = c('Group.1'))
lipids_1$group <- lipids_1$variable %>% gsub('_1','',.) %>% gsub('_2','',.) %>% gsub('_3','',.) 
lipids_1$group %<>% gsub('ASPC1','High',.) %>% gsub('BXPC3','Medium',.) %>% gsub('PANC1','Medium',.) %>% gsub('MIAPACA2','Low',.) %>% factor(., levels = c('High','Medium','Low'))
names(lipids_1) <- c('Class','variable','value','group')

p <- ggplot(lipids_1[lipids_1$Class %in% select_lipid,], aes(x=group, y=value, group=group,fill=group))+ theme_bw()+
  stat_compare_means(method = "t.test",aes(group=Class,label =..p.signif..),comparisons=list(c('High','Medium'),c('Medium','Low'),c('High','Low')))+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun.data='mean_sd',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+
  facet_wrap(~Class, ncol=5, scale='free')+
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_markdown(size = 15), 
        legend.text = element_markdown(size = 15), 
        legend.background = element_blank(), text = element_text(color = "black"),
        axis.title.y = element_text(size = 20, angle = 90), 
        axis.text.x = element_text(size = 15, angle = 90), 
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"), 
        axis.ticks = element_line(linewidth = 1), 
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15))
p
eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 12, height = 8, append = T, title = 'Fig 6.C') 


# Fig 6.B
name_group <- c('type','group','cell')
data_path_1 <- 'D:\\Scientific research\\合作项目\\张宇星\\10.01分析\\ASPC1_NC_cellculture_mousemodel\\Cell_culture'
mouse_1 <- folder_reader(data_path_1, name_group)
mouse_1$filename <- list.files(data_path_1, pattern = "*.txt", full.names = TRUE, include.dirs = T, recursive = T)
mouse_1$mouse <- str_split_i(mouse_1$filename, '/', 2)
mouse_1 <- S_G_smooth(mouse_1)
test_set_sub <- tapply(mouse_1, mouse_1$mouse, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
test_set_sub <- do.call(rbind, list.map(test_set_sub,as.matrix(unlist(.))))
mouse_1$spc <- test_set_sub
means <- aggregate(mouse_1$spc, by=list(mouse_1$mouse,mouse_1$group, mouse_1$cell), mean)  # 五合一
colnames(means) <- c('mouse','group', 'cell', colnames(means)[-c(1:3)])
mouse_1 <- new("hyperSpec",data = means[,1:3],spc = means[,-c(1:3)],wavelength = mouse_1@wavelength)
mouse_1 <- mouse_1[mouse_1$cell!= 'BG'] 
mouse_1$type <- 'cell'
data_path_2 <- 'D:\\Scientific research\\合作项目\\张宇星\\10.01分析\\ASPC1_NC_cellculture_mousemodel\\Mouse_tissue'
mouse_2 <- folder_reader(data_path_2, name_group)
mouse_2$filename <- list.files(data_path_2, pattern = "*.txt", full.names = TRUE, include.dirs = T, recursive = T)
mouse_2$mouse <- str_split_i(mouse_2$filename, '/', 2)
mouse_2 <- S_G_smooth(mouse_2)
test_set_sub <- tapply(mouse_2, mouse_2$mouse, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
test_set_sub <- do.call(rbind, list.map(test_set_sub,as.matrix(unlist(.))))
mouse_2$spc <- test_set_sub
means <- aggregate(mouse_2$spc, by=list(mouse_2$mouse,mouse_2$group, mouse_2$cell), mean)  # 五合一
colnames(means) <- c('mouse','group', 'cell', colnames(means)[-c(1:3)])
mouse_2 <- new("hyperSpec",data = means[,1:3],spc = means[,-c(1:3)],wavelength = mouse_2@wavelength)
mouse_2 <- mouse_2[mouse_2$cell!= 'BG'] 
mouse_2$type <- 'tissue'

mouse_all <- rbind2(mouse_1, mouse_2[mouse_2$group=='NC',])
mouse_all <- Raman_baseline(mouse_all)
mouse_all <- Normalization(mouse_all,from = 2500,to = 3000,nor_cal = "max")

spec_mean <- data.frame(wave=mouse_all@wavelength, cell=colMeans(mouse_all$spc[mouse_all$type=='cell',]), tissue=colMeans(mouse_all$spc[mouse_all$type=='tissue',]))
spec_mean$diff <- (spec_mean$cell - spec_mean$tissue)
spec_mean <- melt(spec_mean, id.vars = 'wave')
setnames(spec_mean,c('wave','variable','spec'))
spec_mean$diff <- spec_mean$spec
spec_mean$diff[spec_mean$variable %in% c('cell','tissue')] <- NA
spec_mean$spec[spec_mean$variable == 'diff'] <- NA
p <- ggplot(spec_mean, aes(color=variable, fill=variable))+
  geom_line(aes(x=wave,y=spec),linewidth=1) + 
  geom_line(aes(x= wave,y=diff*3-0.5),linewidth=1) + 
  geom_hline(yintercept =-0.5, linetype='dashed', color='red', linewidth=1)+
  scale_y_continuous(name='Normalized Intensity (a.u.)', breaks=c(0,0.5,1),
                     sec.axis = sec_axis(~(.+0.5)/3,name='Intensity Difference (a.u.)',breaks = c(-0.05,0,0.05,0.1)))+
  xlab(expression(paste("Wavenumber (cm"^{ -1 }, ")"))) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(500,1000,1500,2000,2500, 3000)
  )+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill='white'))
p

eoffice::topptx(figure = p, filename = "All_figures.pptx", width = 6, height = 4.5, append = T, title = 'Fig 6.C') 

# ---------------- Fig Lipids --------------------
# ------------- Correlations between lipids ----------------------
library(readxl)
library(stringr)
library(reshape2)
lipids <- read_excel('D:\\Scientific research\\合作项目\\张宇星\\Lipid analysis\\张宇星-人细胞-BQ-ZYX20231215-LC-RXB.xlsx')[,c(1,2,6:17)]
lipids$C_chain <- as.character(str_split(lipids$CompoundName, pattern = ' ', simplify = T)[,2]) 
lipids$C_chain <- as.character(str_split(lipids$C_chain, pattern = '-FA', simplify = T)[,1]) 
aa <-str_split(lipids$C_chain, pattern = '-|:', simplify = T)
aa[aa==''] <- 0
lipids$tol_C <- as.numeric(aa[,1])+as.numeric(aa[,3])
lipids$tol_d <- as.numeric(aa[,2])+as.numeric(aa[,4])

lipids_2 <- lipids[,c(1,3:14)]
lipids_2[,-1] <- t(t(lipids_2[,-1]) / colSums(lipids_2[,-1]))
lipids_2 <- aggregate(lipids_2[,-1], by=list(lipids_2$Class), sum)
row.names(lipids_2) <- lipids_2$Group.1
lipids_2 <- cor(t(lipids_2[,-1]))
library(pheatmap)
heat_map <- pheatmap(lipids_2,
                     color=c(colorRampPalette(colors = c("blue","white"))(20),colorRampPalette(colors = c("white","red"))(20)), 
                     cluster_cols = T, cluster_rows = T, legend=TRUE, clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",
                     # cellwidth = 30, cellheight = 40,
                     scale = 'none', trace='none',
                     angle_col='45',show_rownames = T,
                     # gaps_row = c(5,14),gaps_col = c(5,14),
                     cutree_rows = 3,cutree_cols=3,
                     density.info = "none",
                     # annotation_col = band,
                     # annotation_colors=col
                     )
heat_map

# ------------------------------------- New ---------------------------------
# Fig 6A
library('RamanEx')
library_all()
load("D:/Scientific research/合作项目/张宇星/8.16分析/RandomForest_train_dataset.Rdata")
load("D:/Scientific research/合作项目/张宇星/8.16分析/Randomforest_3groups_0815.Rdata")
importance <- data.frame(wave=as.numeric(colnames(data_train)), intensity=colMeans(data_train), gini=model.rf$importance)
data_box <- data.frame(data_train[, colnames(data_train) %in% importance$wave[importance$MeanDecreaseGini>0.75]])
colnames(data_box) <- round(as.numeric(gsub('X','',colnames(data_box))),0)
data_box <- data_box[,!colnames(data_box) %in% c('1800','1434','1461')]
data_box$label <- label_train
data_box$label %<>% factor(.,levels=c('high','medium','low'))
cal_info <- function(y,x){
  lm_info <- summary(lm(y ~ as.numeric(x)))
  aa <- c(lm_info$r.squared, lm_info$coefficients[2,1])
  return(aa)
}
# data_mean <- reshape2::melt(aggregate(data_box[,1:12], by=list(data_box$label),FUN=mean))

data_plot <- t(apply(data_box[,1:12],2,function(x)cal_info(x, data_box$label)))
data_plot <- cbind(data_plot,row.names(data_plot))
colnames(data_plot) <- c('r2','coeff','variable')
# order <- data_plot[order(colMeans(data_box[,-13]), decreasing = T),3]
order <- c("1449",  "1302",  "1435", "1129", "1078", "1065",  "854", "1559","1699", "1352", "1668", "1465" )
data_box <- reshape2::melt(data_box)
data_plot <- merge(data_plot, data_box, by='variable')
data_plot$value[data_plot$coeff>0] <- data_plot$value[data_plot$coeff>0]*-1
data_plot$variable %<>% factor(., levels=order)
p <- ggplot(data_plot,aes(variable,value,group=label,fill=label))+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun.data='mean_cl_normal',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+ #mean_sdl
  # stat_compare_means(method = "t.test",aes(group=tol_d,label =..p.signif..),comparisons=list(c('high','medium'),c('medium','low'),c('high','low')))+
  scale_y_continuous(breaks = seq(-0.35, 0.35, 0.1), labels = as.character(abs(seq(-0.35, 0.35, 0.1))), limits = c(-0.35, 0.35))+  
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(),axis.text = element_text(color="black",size=12),axis.title = element_text(color = "black",size=15)) +geom_hline(yintercept = 0, size = 0.4)
p
eoffice::topptx(figure = p, filename = "Mod_figures.pptx", width = 12, height = 4.5, append = T, title = 'Fig 6.B') 

# 
# plot <- ggplot(data_box, aes(x = label, y = values, group = label, fill = label)) +  
#   theme_bw() + geom_point(aes(group = label, color = label), position = position_jitter(0.2)) +
#   geom_boxplot(aes(group = label), alpha = 0.2, outlier.alpha = 0) + 
#   # geom_signif(comparisons=list(c('ASPC1','BXPC3'),c('BXPC3','PANC1'),c('PANC1','MIAPACA2')),map_signif_level = T,y_position = c(1.05,1.03,1))+
#   # geom_hline(aes(yintercept = 0), colour = "red", linetype = "dashed", size = 1) +
#   facet_wrap(~wave, scales = 'free')+
#   theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
#         legend.title = element_markdown(size = 15), 
#         legend.text = element_markdown(size = 15), 
#         legend.position = "none",
#         legend.background = element_blank(), text = element_text(color = "black"),
#         axis.title.y = element_text(size = 20, angle = 90), 
#         axis.text.x = element_text(size = 8, angle = 0), 
#         axis.text.y = element_text(size = 8), 
#         strip.text = element_text(size = 20), 
#         strip.background = element_rect(fill = "white"), 
#         axis.ticks = element_line(size = 1), 
#         axis.ticks.length = unit(0.4, "lines"),
#         axis.title = element_text(size = 20)) + 
#   xlab("Group") + ylab("Raman Intensity (2937)")
# plot

# Fig 6B
lipids_order <- c('DCER','LPC','BMP','LPI','CER','PC','LPG','PE','PI','PS','SM','MGDG','PG','HCER','TAG','LPE','DAG','FFA','LCER','CE')
lipids_1 <- lipids[,c(1,3:14)]
lipids_1[,2:13] <- t(t(lipids_1[,2:13]) / colSums(lipids_1[,2:13]))
lipids_1 <- aggregate(lipids_1[,2:13], by=list(lipids_1$Class),sum)
lipids_1 <- reshape2::melt(lipids_1, id.vars = c('Group.1'))
lipids_1$group <- lipids_1$variable %>% gsub('_1','',.) %>% gsub('_2','',.) %>% gsub('_3','',.) 
lipids_1$group %<>% gsub('ASPC1','High',.) %>% gsub('BXPC3','Medium',.) %>% gsub('PANC1','Medium',.) %>% gsub('MIAPACA2','Low',.) %>% factor(., levels = c('High','Medium','Low'))
names(lipids_1) <- c('Class','variable','value','group')
lipids_1 <- reshape2::dcast(lipids_1[,-2], Class ~ group, fun.aggregate = mean)
lipids_2 <- lipids_1
lipids_2[,2:4] <- t(apply(lipids_2[,2:4],1, scale))
lipids_2 <- reshape2::melt(lipids_2)

# lipids_1 <- aggregate(lipids_1$value, by=list(lipids_1$Class, lipids_1$group), mean)
names(lipids_2) <- c('Class','group','mean')
lipids_2$Class %<>% factor(., levels = lipids_order)
lipids_2$group %<>% factor(., levels = c('Low','Medium','High'))

p <- ggplot(lipids_2,aes(x=Class,y=group)) + 
  geom_tile(aes(fill=mean), colour='grey')+
  # scale_fill_gradient(low="white",high = "#191970") +
  scale_fill_material("blue") +
  # scale_fill_material("cyan") + 
  # scale_fill_material('deep-orange') + 
  # scale_fill_material("light-blue") + 
  # scale_fill_material('teal') +
  # scale_color_npg() +
  # scale_fill_distiller(palette = "Blues", direction = 1)+
  # scale_fill_gsea(alpha = 1) +
  theme_bw() +
  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5) )
  
p

eoffice::topptx(figure = p, filename = "Mod_figures.pptx", width = 10, height = 2.5, append = T, title = 'Fig 6.C') 

lipids_high <- data.frame(Class= lipids_1$Class,prop = lipids_1$High, res = 1-lipids_1$High)

lipids_high$x <- order(lipids_order)
lipids_high$y <- rep(1,20)

library(scatterpie)
p <- ggplot()+
  geom_scatterpie(data=lipids_high,
                  aes(x,y,group=y,r=0.45),
                  cols = c("prop","res"))+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#808080","#FFFFFF"))
p
eoffice::topptx(figure = p, filename = "Mod_figures.pptx", width = 9.2, height = 2, append = T, title = 'Fig 6.C') 

# pie plot of each class
library(RColorBrewer)
lipids_3 <- reshape2::melt(lipids_1)
colnames(lipids_3) <- c('Class','group','prop')
lipids_3$res <- 1-lipids_3$prop
lipids_3$x <- c(rep(1,20), rep(2,20), rep(3,20))
lipids_3$y <- rep(1,60)

lipids_3 <- data.frame(t(lipids_1[,-1]))
colnames(lipids_3) <- lipids_1[,1]
lipids_3$x <- c(1,2,3)
lipids_3$y <- c(1,1,1)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
cols <- sample(col_vector, 20)
pie(rep(1,20),col=cols)

p <- ggplot()+
  geom_scatterpie(data=lipids_3,
                  aes(x,y,group=x,r=0.4),
                  cols = lipids_order)+
  coord_equal()+
  scale_fill_manual(values = cols)+
  theme_void()+
  theme(legend.position = "none")
p
eoffice::topptx(figure = p, filename = "Mod_figures.pptx", width = 10, height = 10, append = T, title = 'Fig S') 


# ---------------- Migration score -----------------
transwell <- read_excel('Cell_metastasis.xlsx')
transwell <- reshape2::melt(transwell, id.var=c('Sample','Type'))
colnames(transwell) <- c('group','Type','variable','score')
weight <- matrix(data=c(1,1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1),ncol=1)

# High, medium, low
pal_npg()(3)
data_box <- train_set$spc[, train_set@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(group=train_set$strain, score=score, batch=train_set$batch)
mig_score$batch %<>% str_split_i(., pattern = '_', 3)
mig_score$group %<>% gsub('ASPC1','High',.) %>% gsub('BXPC3','Medium',.) %>% gsub('PANC1','Medium',.) %>% gsub('MIAPACA2','Low',.) %>% factor(.,levels = c('High','Medium','Low'))

p <- ggplot(mig_score, aes(group,score,group=group))+ #, color=group, fill=group
  # geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  # geom_point(aes(group = group, color = group), position = position_jitter(0.2)) +
  # geom_boxplot(aes(group = group), alpha = 0.2, outlier.alpha = 0) +
  geom_violin(aes(fill = group), trim = FALSE)+
  geom_boxplot(aes(group = group), width=0.12) +
  # stat_summary(fun = mean,color='black',
  #              geom = "errorbar",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  # geom_signif(comparisons=list(c('High','Medium'),c('Medium','Low'),c('High','Low')),map_signif_level = T, test='t.test')+
  scale_y_continuous(name='Raman Migration Score', breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2))+
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF", "#00A087FF"))+
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF", "#00A087FF"))+
  theme_bw() + 
  # coord_cartesian(ylim = c(0.3,1.1))+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Train.pptx", width = 4.5, height = 4, append = T) 
# Calculate sample size
set.seed(42)
sample_size_high <- data.frame(Group='High' ,cal_Minimum_Samplesize_Realtime(mig_score$score[mig_score$group=='High'], initial_number = 1, seq = 1, Max_sample_number = 50))
sample_size_high$Number_sample_max_Pvalue <- max(sample_size_high$Sample_number[sample_size_high$Pvalue<0.95])+1
sample_size_medium <- data.frame(Group='Medium' ,cal_Minimum_Samplesize_Realtime(mig_score$score[mig_score$group=='Medium' ], initial_number = 1, seq = 1, Max_sample_number = 50))
sample_size_medium$Number_sample_max_Pvalue <- max(sample_size_medium$Sample_number[sample_size_medium$Pvalue<0.95])+1
sample_size_low <- data.frame(Group='Low' ,cal_Minimum_Samplesize_Realtime(mig_score$score[mig_score$group=='Low' ], initial_number = 1, seq = 1, Max_sample_number = 51))
sample_size_low$Number_sample_max_Pvalue <- max(sample_size_low$Sample_number[sample_size_low$Pvalue<0.95])+1
sample_size <- rbind(sample_size_high, sample_size_medium, sample_size_low)
sample_size$N <- 1
colnames(sample_size) <- c('Group','Number_sample','Pvalue','Number_sample_max_Pvalue','N')
sample_size$Group %<>% factor(.,levels=c('High','Medium','Low'))
plot <- ggplot(sample_size, aes(x = Number_sample, y = Pvalue,  color = N)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), legend.title = element_blank(), 
        legend.text = element_markdown(size = 15, family = "myFont"), 
        legend.position = "none", legend.background = element_blank(), 
        text = element_text(color = "black"), 
        axis.title.y = element_text(size = 20, angle = 90, family = "myFont"),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 15, angle = 45, family = "myFont", hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15, family = "myFont"), 
        strip.background = element_rect(fill = "white"),
        strip.text = element_markdown(size = 20, family = "myFont"), 
        axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(family = "myFont", size = 20)) + 
  geom_line(aes(group = N, color = N), linewidth=1) + 
  facet_grid(~Group, scales = "free_x") + 
  geom_hline(aes(yintercept = 0.95), colour = "red", linetype = "dashed", linewidth=1) +
  geom_vline(aes(xintercept = Number_sample_max_Pvalue), colour = "red", linetype = "dashed", linewidth=1) + 
  geom_text(aes(x = Number_sample_max_Pvalue+10, y = 0.85, label = paste(Number_sample_max_Pvalue , "cells", sep = " ")), colour = "#990000", size = 5) + 
  xlab("Sample size") + 
  scale_x_continuous(breaks = c(10,20,30,40,50))+
  ylab(expression(paste("Probability (RE"[mean],"<0.95)", sep = "")))
plot
eoffice::topptx(figure = plot, filename = "SI.pptx", width = 8, height = 3, append = T) 

# IRCA Analysis
library(Cairo)
Cairo::CairoPNG( filename = "All_train_local_IRCA.png", width = 7, height = 7,units = "in",dpi = 300)
draw.global.IRCA(train_set$spc, group = paste(train_set$strain, train_set$batch))
dev.off() 
for(group in unique(paste(train_set$strain, train_set$batch))){
  Cairo::CairoPNG( filename = str_c('Global_IRCA_', group, '.png', sep=''), width = 7, height = 7,units = "in",dpi = 300)
  draw.global.IRCA(train_set$spc[paste(train_set$strain, train_set$batch)==group,], group = group)
  dev.off() 
}


irca_corr <- sapply(unique(test_PANC1$batch), function(x){cor_matrix <- Hmisc::rcorr(as.matrix(test_PANC1$spc[test_PANC1$batch==x,]), type = 'pearson')
  cor_matrix$r[cor_matrix$P > 0.05] <- 0
  return(as.vector(cor_matrix$r[upper.tri(cor_matrix$r)]))})
irca_corr[irca_corr > -0.6] <- 0
sample_corr <- cor(irca_corr, method='pearson')
sample_label <- str_split_i(colnames(sample_corr),pattern = 'H2O',1)
sample_corr <- sapply(unique(sample_label), function(x){cor_matrix <- sample_corr[sample_label==x,sample_label==x]
return(mean(cor_matrix[upper.tri(cor_matrix)]))})
# P1E6, P1G8, P2D6, P2G8
pal_npg()(4)
data_box <- test_PANC1$spc[, test_PANC1@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(group=test_PANC1$strain, score=score)

p <- ggplot(mig_score, aes(group,score,group=group))+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  # geom_point(aes(group = group, color = group), position = position_jitter(0.2)) +
  # geom_boxplot(aes(group = group), alpha = 0.2, outlier.alpha = 0) + 
  # geom_violin(aes(fill = group), trim = FALSE)+
  # geom_boxplot(aes(group = group), width=0.12) + 
  stat_summary(fun = mean,color='black',
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  geom_signif(comparisons=list(c('P1E6','P1G8'),c('P2D6','P1G8'),c('P1G8','P2G8')),map_signif_level = T, test='t.test',y_position = c(0.95,0.95,1.0))+
  scale_y_continuous(name='Raman Migration Score', breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1.0))+
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF", "#00A087FF", "#3C5488FF"))+
  scale_color_manual(values = c("#4DBBD5FF","#E64B35FF", "#00A087FF", "#3C5488FF"))+
  theme_bw() + 
  coord_cartesian(ylim = c(0.3,1.1))+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 4.5, height = 4, append = T) 


# ASPC1+PANC1
data_box <- test_ASPC1$spc[, test_ASPC1@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(group=test_ASPC1$strain, score=score)
mig_score$Type <- 'Raman'
mig_score <- rbind(mig_score, transwell[,-3])
mig_score$group %<>% gsub('ASPC1-KO','ASPC1-hnRNPAB-KO',.) %>% gsub('ASPC1-NC','ASPC1-hnRNPAB-NC',.) %>% gsub('PANC1-KO','PANC1-hnRNPLL-KO',.) %>% gsub('PANC1-NC','PANC1-hnRNPLL-NC',.)

mig_temp <- mig_score[mig_score$group %in% c('ASPC1-hnRNPAB-KO','ASPC1-hnRNPAB-NC'),]  # "PANC1-hnRNPLL-NC","PANC1-hnRNPLL-KO"
mig_temp$score[mig_temp$Type=='Raman'] <- (mig_temp$score[mig_temp$Type=='Raman']-0.5)*500
mig_temp$group %<>% factor(., levels=c("PANC1-hnRNPLL-NC","PANC1-hnRNPLL-KO"))
mig_temp$group %<>% factor(., levels=c('ASPC1-hnRNPAB-NC','ASPC1-hnRNPAB-KO'))

p <- ggplot(mig_temp, aes(Type,score,group=group,color=group, fill=group))+
  geom_bar(stat="summary",fun=mean, position="dodge",size=0.5,color="black") +
  stat_summary(fun = mean,color='black',
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  # geom_signif(comparisons=list(c('ASPC1-hnRNPAB-KO','ASPC1-hnRNPAB-NC')),map_signif_level = T, test = 't.test')+
  # stat_summary(fun.data='mean_sd',geom="errorbar",position = position_dodge (.9),width=0.15,size=0.5)+
  scale_y_continuous(name='Cells per field', breaks=c(0,50,100,150),
                     sec.axis = sec_axis(~./500+0.5,name='Raman Migration Score',breaks = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0)))+
  # facet_wrap(~Type)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 4.5, height = 4, append = T) 

# ASPC 1 - RFP
path <- 'D:\\Scientific research\\合作项目\\张宇星\\原始数据\\ASPC1_RFP_NC_KO'
name_group <- c('group','type','cell')
data_4 <- folder_reader(path, name_group)
data_4 <- S_G_smooth(data_4)
data_4$batch <- c(rep(1,130),rep(2,130),rep(3,130),rep(4,130))
data_4_sub <- tapply(data_4, data_4$batch, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
data_4_sub <- do.call(rbind, list.map(data_4_sub,as.matrix(unlist(.))))
data_4$spc <- data_4_sub
means <- aggregate(data_4$spc, by=list(data_4$type, data_4$batch, data_4$cell), mean)  # 五合一
colnames(means) <- c('type','batch', 'cell', colnames(means)[-c(1:3)])
data_4 <- new("hyperSpec",data = means[,c(1,3)],spc = means[,-c(1:3)],wavelength = data_4@wavelength)
data_4 <- data_4[data_4$cell!= 'BG']  # 继续预处理
data_4 <- Raman_baseline(data_4)
data_4 <- Normalization(data_4,from = 2500,to = 3000,nor_cal = "max")
p <- mean.spec(data_4$spc, group = factor(data_4$type), gap=0.3)
p 
pred.1 <- predict(model.rf, data_4$spc[,6:633])
pred.1 %<>% factor(.,levels = c('high','medium','low'))
pre.plot(data_4$type, pred.1)

data_box <- data_4$spc[, data_4@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(group=data_4$type, score=score)
mig_score$Type <- 'Raman'
mig_score <- rbind(mig_score, transwell[,-3])
mig_score$group %<>% gsub('ASPC1-BLVRB-KO','KO2',.) %>% gsub('ASPC1-BLVRB-NC','NC',.) %>% factor(.,levels=c('NC','KO2'))

mig_temp <- mig_score[mig_score$group %in% c('KO2','NC'),]  
mig_temp$score[mig_temp$Type=='Raman'] <- (mig_temp$score[mig_temp$Type=='Raman']-0.6)*300

p <- ggplot(mig_temp, aes(group,score,group=group,color=group, fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  stat_summary(fun = mean,color='black',
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  geom_signif(comparisons=list(c('KO2','NC')),map_signif_level = T, test = 't.test')+
  scale_y_continuous(name='Cells per field', breaks=c(0,50,100,150,200,250,300),
                     sec.axis = sec_axis(~./300+0.6,name='Raman Migration Score',breaks = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0)))+
  facet_grid(~Type)+
  xlab('Type') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 4.5, height = 4, append = T) 
1043
# Mouse
pal_npg()(2)
data_box <- test_mouse$spc[, test_mouse@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(mouse=test_mouse$mouse,group = test_mouse$group, score=score)
mig_score$mouse %<>% gsub('ASPC1_','',.)
mig_score$mouse %<>% factor(., levels = c("KO_1","KO_2","KO_3","KO_4","KO_6","KO_7","KO_8","KO_10",
                                          "NC_1","NC_2","NC_3","NC_4","NC_6","NC_7","NC_8" ))

p <- ggplot(mig_score, aes(mouse,score,group=mouse))+
  # geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  geom_point(aes(group = mouse, color = group), position = position_jitter(0.2)) +
  geom_boxplot(aes(group = mouse, fill=group), alpha = 0.2, outlier.alpha = 0) +
  # geom_violin(aes(fill = group), trim = FALSE)+
  # stat_summary(fun = mean,color='black',
  #              geom = "errorbar",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  # geom_signif(comparisons=list(c('P1E6','P1G8'),c('P2D6','P1G8'),c('P2D6','P2G8')),map_signif_level = T)+
  scale_y_continuous(name='Raman Migration Score', breaks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2))+
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF"))+
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF"))+
  theme_bw() + 
  coord_cartesian(ylim = c(0.3,1.2))+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 10, height = 4, append = T) 

p <- ggplot(mig_score, aes(group,score,group=group))+ #, fill=group, color=group
  # geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  # geom_point(aes(group = group, color = group), position = position_jitter(0.2)) +
  geom_violin(aes(fill = group), trim = FALSE)+
  geom_boxplot(aes(group = group, fill=group), alpha = 0.2, outlier.alpha = 0, width=0.4) +
  # stat_summary(fun = mean,color='black',
  #              geom = "errorbar",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  geom_signif(comparisons=list(c('KO','NC')),map_signif_level = T, test='t.test',y_position = 1.4)+
  scale_y_continuous(name='Raman Migration Score', breaks=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5))+
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF"))+
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF"))+
  theme_bw() + 
  coord_cartesian(ylim = c(0.2,1.5))+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 3, height = 4, append = T) 

# MUM2B  Siha
path <- 'D:\\Scientific research\\合作项目\\张宇星\\原始数据\\MUM2B_Siha'
name_group <- c('group','type','cell')
data_5 <- folder_reader(path, name_group)
data_5 <- S_G_smooth(data_5)
data_5$batch <- c(rep(1,130),rep(2,130),rep(3,130),rep(4,130),rep(5,130),rep(6,130),rep(7,130),rep(8,130))
data_5_sub <- tapply(data_5, data_5$batch, function(data){
  t(t(data$spc) - colMeans(data$spc[data$cell=='BG',])) }, simplify = FALSE, default = 0)
data_5_sub <- do.call(rbind, list.map(data_5_sub,as.matrix(unlist(.))))
data_5$spc <- data_5_sub
means <- aggregate(data_5$spc, by=list(data_5$group,data_5$type, data_5$batch, data_5$cell), mean)  # 五合一
colnames(means) <- c('group','type','batch', 'cell', colnames(means)[-c(1:4)])
data_5 <- new("hyperSpec",data = means[,c(1,4)],spc = means[,-c(1:4)],wavelength = data_5@wavelength)
data_5$type <- means[,2]
data_5 <- data_5[data_5$cell!= 'BG']  # 继续预处理
data_5 <- Raman_baseline(data_5)
data_5 <- Normalization(data_5,from = 2500,to = 3000,nor_cal = "max")

pred.1 <- predict(model.rf, data_5$spc[,6:633])
pred.1 %<>% factor(.,levels = c('high','medium','low'))
pre.plot(paste(data_5$group,data_5$type,sep=' '), pred.1)

data_box <- data_5$spc[, data_5@wavelength %in% importance$wave[importance$gini>0.75]]
score <- data_box %*% weight
mig_score <- data.frame(group=paste(data_5$group,data_5$type,sep=' '), score=score)
mig_score$group %<>% factor(.,levels=c( "MUM2B NC","MUM2B A", "Siha NC","Siha OE" ))

p <- ggplot(mig_score, aes(group,score,group=group))+#,color=group, fill=group
  # geom_bar(stat="summary",fun=mean,position="dodge",size=0.5,color="black") +
  # geom_point(aes(group = group, color = group), position = position_jitter(0.2)) +
  # geom_boxplot(aes(group = group), alpha = 0.2, outlier.alpha = 0) +
  geom_violin(aes(fill = group), trim = FALSE)+
  geom_boxplot(aes(group = group), width=0.2) +
  # stat_summary(fun = mean,color='black',
  #              geom = "errorbar",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x),position = position_dodge (.9),width=0.15,size=0.5)+
  geom_signif(comparisons=list(c("MUM2B NC","MUM2B A"),c("Siha NC","Siha OE")),map_signif_level = T, test = 't.test', y_position = 1.0)+
  scale_y_continuous(name='Raman Migration Score', breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1.0))+
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF", "#00A087FF", "#3C5488FF"))+
  scale_color_manual(values = c("#4DBBD5FF","#E64B35FF", "#00A087FF", "#3C5488FF"))+
  coord_cartesian(ylim = c(0.35,1.05))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'None',
        strip.background = element_rect(fill='white'))
p
eoffice::topptx(figure = p, filename = "Combined data.pptx", width = 4.5, height = 4, append = T) 

# Retrain on mouse tissue
library(randomForest)
for (i in 1:100){
  print(i)
  aa <- c(sample(c('ASPC1_KO_1','ASPC1_KO_2','ASPC1_KO_7','ASPC1_KO_8','ASPC1_KO_6'),2),
        sample(c('ASPC1_KO_3','ASPC1_KO_4','ASPC1_KO_10'),1),
        sample(c('ASPC1_NC_1','ASPC1_NC_3','ASPC1_NC_4','ASPC1_NC_7','ASPC1_NC_8'),2),
        sample(c('ASPC1_NC_2','ASPC1_NC_6'),1))
  mouse_train <- test_mouse$spc[test_mouse$mouse %in% aa,6:633]
  label_mouse_train <- test_mouse$mouse[test_mouse$mouse %in% aa]
  mouse_test <- test_mouse$spc[!test_mouse$mouse %in% aa,6:633]
  label_mouse_test <- test_mouse$mouse[!test_mouse$mouse %in% aa]
  label_mouse_train %<>% gsub('KO_3|KO_4|KO_10|NC_2|NC_6','N',.) %>% 
    gsub(c('KO_1|KO_2|KO_6|KO_7|KO_8|NC_1|NC_3|NC_4|NC_7|NC_8'),'Y',.) %>% factor()
  label_mouse_test %<>% gsub('KO_3|KO_4|KO_10|NC_2|NC_6','N',.) %>% 
    gsub(c('KO_1|KO_2|KO_6|KO_7|KO_8|NC_1|NC_3|NC_4|NC_7|NC_8'),'Y',.)%>% factor()
  model.rf.mouse <- randomForest(mouse_train, as.factor(label_mouse_train), ntree = 100, mtry = 2, replace = TRUE)
  # table(label_mouse_train, predict(model.rf.mouse , mouse_train))
  data_pre <- predict(model.rf.mouse , mouse_test)
  pred_matrix <- as.data.frame(table(label_mouse_test, data_pre))
  base::colnames(pred_matrix) <- c('true_labels', 'pred_labels', 'Freq')
  if(cal_accuracy(pred_matrix)>0.7) cat(cal_accuracy(pred_matrix),aa,'\n')
  pre.plot(label_mouse_test,data_pre)
}

aa <- c('ASPC1_KO_1', 'ASPC1_KO_7', 'ASPC1_KO_4', 'ASPC1_NC_8', 'ASPC1_NC_7', 'ASPC1_NC_6')
cal_accuracy <- function(pred_matrix){
  acc_matrix <- pred_matrix
  acc <- sum(pred_matrix[as.character(pred_matrix[,1])==as.character(pred_matrix[,2]),3])/sum(pred_matrix[,3])
  return(acc)
}



# AUC
cal_acc <- function(x,y){
  conf_mat <- confusionMatrix(factor(x, levels = c(0,1)), factor(y), positive = "1")
  return(conf_mat$byClass[1:2])
}
bb <- cbind(sort(aa$mean)[1:14], sort(aa$mean)[-1])
bb <- rowMeans(bb)

s_s <- lapply(seq(min(aa$mean-0.1),max(aa$mean+0.1),0.001),function(x){
  tt <- rep(0, 15)
  tt[aa$mean>x] <- 1
  return(cal_acc(tt,aa$trans))
})

s_s <- do.call(rbind,s_s)
plot(s_s[,1], 1-s_s[,2], type='l')
lines(0:1,0:1, col='red')
