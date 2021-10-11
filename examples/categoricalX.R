library(survival);library(dplyr)
lung %>%
  mutate(status = as.integer(status == 1),
         sex = factor(sex),
         kk = factor(as.integer(pat.karno >= 70)),
         kk1 = factor(as.integer(pat.karno >= 60))) -> lung

lung %>%
  mutate(age2 = cut(age, c(min(age), quantile(lung$age)[2], quantile(lung$age)[4], max(age)))) %>%
  mutate(age2 = factor(age2)) ->
  lung

library(jstable)


SubCoxCatX(Surv(time, status) ~ sex, var_subgroup = 'kk', data = lung)

SubMultiCoxCatX(Surv(time, status) ~ sex, var_subgroups = c("kk",'kk1','age2'), data = lung)

table_forest <- SubMultiCoxCatX(Surv(time, status) ~ sex, var_subgroups = c("kk",'kk1','age2'), data = lung)

table_forest %>% tibble::tibble() -> table_forest.t
library(forestplot)
HR <- ifelse(!is.na(table_forest.t$HR), paste0(table_forest.t$HR, ' (', table_forest.t$Lower, '-', table_forest.t$Upper, ')'),NA)
tabletext <- cbind(c("Subgroup", table_forest.t$Variable),
                   c("With",  table_forest.t$With),
                   c("Without", table_forest.t$Without),
                   c("Hazard ratio (95% CI)",HR),
                   c("P value",table_forest.t$Pvalue),
                   c("P for interaction",  table_forest.t$`P for interaction`))

forestplot(labeltext = tabletext,
           mean = c(NA, as.numeric(table_forest.t$HR)),
           lower = c(NA, as.numeric(table_forest.t$Lower)),
           upper = c(NA, as.numeric(table_forest.t$Upper)),
           zero=1,
           boxsize=0.2,
           graph.pos= 4,#图放在第四列
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "13"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           xlab="Hazard ratio (95% CI)",
           # xticks=c(0,1,3) ,
           #----------------#字体
           is.summary=c(T, F,T,F,F,T,F,F,T,F,F,F),#T=粗体
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
           clip=c(0,10),
           #----------------#行间距、字间距/box形状
           ineheight=unit(6, 'mm'),
           line.margin=unit(8, 'mm'),
           colgap=unit(6, 'mm'),
           col=fpColors(zero = "#e22e2a",
                        box = '#048410',
                        lines = 'black'),
           fn.ci_norm="fpDrawCircleCI",
           title="Forestplot")
