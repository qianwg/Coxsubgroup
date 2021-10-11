library(survival);library(dplyr)
lung %>% 
  mutate(status = as.integer(status == 1),
         sex = factor(sex),
         kk = factor(as.integer(pat.karno >= 70)),
         kk1 = factor(as.integer(pat.karno >= 60))) -> lung
table_forest <- TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), 
                                      data=lung, decimal.hr = 3, time_eventrate = 100, line = TRUE)

table_forest <- TableSubgroupMultiCox(Surv(Re_time, Recurrence) ~ AIP2, var_subgroups = c("Diabetes"), 
                                      data=data_unimp, decimal.hr = 3,  line = TRUE)

table_forest_drop <- table_forest[-c(2,6),]
table_forest2 <- table_forest_drop %>% tibble() %>% 
  mutate_if(is.character, as.numeric) 
library(forestplot)
np <- ifelse(!is.na(table_forest_drop$Count), paste(table_forest_drop$Count,' (',table_forest_drop$Percent,')', sep = ''), NA)
HR <- ifelse(!is.na(table_forest_drop$`Point Estimate`), paste0(table_forest_drop$`Point Estimate`, ' (', table_forest_drop$Lower, '-', table_forest_drop$Upper, ')'),NA)
tabletext <- cbind(c("Subgroup", table_forest_drop$Variable),
                   c("No. of Patients (%)",  np),
                   c("Hazard ratio (95% CI)",HR),
                   c("P value",table_forest_drop$`P value`),
                   c("P for interaction",  table_forest_drop$`P for interaction`))


forestplot(labeltext = tabletext,
           mean = c(NA, table_forest2$`Point Estimate`),
           lower = c(NA, table_forest2$Lower),
           upper = c(NA, table_forest2$Upper),
           zero=1,            
           boxsize=0.2,      
           graph.pos= 3,#图放在第四列
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "9"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           xlab="Hazard ratio (95% CI)",
           # xticks=c(0,1,3) ,
           #----------------#字体
           is.summary=c(T, F,T,F,F,T,F,F),#T=粗体
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.1), 
                          xlab=gpar(cex=1), 
                          title=gpar(cex=2)),
           #----------------#线条粗细（x轴、置信区间）
           lwd.zero=1.5,
           lwd.ci=1.5,
           lwd.xaxis=1.5, 
           lty.ci=1,
           ci.vertices =T,
           ci.vertices.height=0.1, 
           clip=c(0,3),
           #----------------#行间距、字间距/box形状                 
           ineheight=unit(6, 'mm'), 
           line.margin=unit(8, 'mm'),
           colgap=unit(6, 'mm'),
           col=fpColors(zero = "#e22e2a",
                        box = '#048410', 
                        lines = 'black'),
           fn.ci_norm="fpDrawCircleCI", 
           title="亚组分析森林图")