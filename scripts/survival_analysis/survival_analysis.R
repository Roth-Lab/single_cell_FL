###### survival analysis ######
library(survminer)
library(tidyverse)
library(survival)
library(survMisc)

# fl_rb_fin is a dataframe with proportions of interest and patient meta data (like survival time and CODE_OS)

# FL-RB or FL-RCVP
fl_rb_fin=read.csv('~/Desktop/BRClininfo_w_proportions_CD20-_FL-RB.csv',row.names = 1, check.names = FALSE)
#fl_rb_fin=read.csv('/path/BRClininfo_w_proportions_CD20-_FL-RCVP.csv',row.names = 1, check.names = FALSE)

# KM survival plot no group
p=ggsurvplot(survfit(Surv(fl_rb_fin$`Disease specific survival (y)`, fl_rb_fin$CODE_DSS) ~ 1, data = fl_rb_fin), xlab='Years', ylab='Overall survival probability')
ggsave(plot=p$plot,height=5,width=5,dpi=200, filename="~/Desktop/km_survival_plot_FL-RB.pdf", useDingbats=FALSE)

survdiff(Surv(fl_rb_fin$`Progression free survival (y)`, fl_rb_fin$CODE_PFS) ~ fl_rb_fin$`CD8+LAG3+_max`, data = fl_rb_fin)


# survival metrics
pair1="Disease specific survival (y)"
pair2="CODE_DSS"


# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(fl_rb_fin, time = pair1, event = pair2,
                         variables = c('CD8+LAG3+PD1+_prop_mean', "CD8+LAG3+_prop_mean"))
res.cut

fl_rb_fin$status=ifelse(fl_rb_fin$`CD8+LAG3+_prop_mean`>=0.018,'>=1.8%','<1.8%')

# for foxp3 pattern
fl_rb_fin=fl_rb_fin[fl_rb_fin$`FOXP3 Immunoarchitectural Pattern`!='Uninterpretable',]
fl_rb_fin$status=ifelse(fl_rb_fin$`FOXP3 Immunoarchitectural Pattern` %in% c('Diffuse','Diffuse?'),'Diffuse','Follicular')

fl_rb_fin[,pair2]=as.numeric(fl_rb_fin[,pair2])
fit=survfit(Surv(fl_rb_fin[,pair1], fl_rb_fin[,pair2]) ~ status, data = fl_rb_fin)
ggsurv <- ggsurvplot(fit,
                     pval = TRUE,
                     size = 1,
                     risk.table = TRUE, # Add risk table
                     #risk.table.col = i, # Change risk table color by groups
                     risk.table.y.text.col = TRUE,
                     risk.table.height = 0.25,
                     legend.title = "",
                     ggtheme=theme_bw(base_size = 10),
                     risk.table.fontsize = 2.5,
                     # surv.median.line = "hv", # Specify median survival
                     #ggtheme = theme_bw(), # Change ggplot2 theme
                     xlab = "Time (years)",
                     ylab = paste(pair1,'vs', pair2, sep = " "),
                     palette = c("#666666","#FF0033"),
                     tables.theme = clean_theme())

#ggsave(plot=ggsurv,height=5,width=5,dpi=E00, filename="~/Desktop/km_survival_plot_OS_vs_CODE_OS_FL-RB.pdf", useDingbats=FALSE)
pdf(paste0("/path/km_foxp3_survival_plot_",pair1,"_vs_",pair2,"_FL-RCVP.pdf"),width = 5,height = 6)
print(ggsurv)
dev.off()

# multivariate analysis
tmp=fl_rb_fin[fl_rb_fin$FLIPI != -1,]
#tmp$FLIPIgrp=factor(tmp$FLIPIgrp,levels = c('LOW','INTERMED','HIGH'))
tmp$FLIPIgrp=ifelse(tmp$FLIPIgrp=='HIGH',tmp$FLIPIgrp,'LOW-INTERMED')
tmp$FLIPIgrp=factor(tmp$FLIPIgrp,levels = c('LOW-INTERMED','HIGH'))

tmp$DIAGgrp=ifelse(tmp$DIAG=='FOLL3A','FOLL3A','FOLL1_2')
tmp$DIAGgrp=factor(tmp$DIAGgrp,levels = c('FOLL1_2','FOLL3A'))

#tmp$Proportion=tmp$`CD8+LAG3+_prop_mean`
tmp$Proportion=ifelse(tmp$`CD8+LAG3+_prop_mean`>=0.011,'CD8LAG3_High',tmp$`CD8+LAG3+_prop_mean`)
tmp$Proportion=ifelse(tmp$Proportion<0.011 & tmp$Proportion>=0.004,'CD8LAG3_Mid-High',tmp$Proportion)
tmp$Proportion=ifelse(tmp$Proportion<0.004 & tmp$Proportion>=0.0012,'CD8LAG3_Mid-Low',tmp$Proportion)
tmp$Proportion=ifelse(tmp$Proportion<0.0012,'CD8LAG3_Low',tmp$Proportion)
tmp$Proportion=factor(tmp$Proportion,levels = c('CD8LAG3_Low','CD8LAG3_Mid-Low','CD8LAG3_Mid-High','CD8LAG3_High'))

fit.coxph2 <- coxph(Surv(tmp[,pair1], tmp[,pair2]) ~  status + DIAGgrp + FLIPIgrp,  data = tmp)
summary(fit.coxph2)

p=ggforest(fit.coxph2, data = tmp, main = paste0(pair2, ', ', pair1)) 
pdf(paste0("/path/km_forest_plot_",pair1,"_vs_",pair2,"_FL-RCVP_binary.pdf"),width = 7,height = 6)
print(p)
dev.off()
