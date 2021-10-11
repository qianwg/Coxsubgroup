#' @title Coxsubgroup: Sub-group analysis table for Cox model.
#' @description Sub-group analysis table for Cox model.
#' @param formula formula with survival analysis.
#' @param formulas formula list with survival analysis of different outcomes
#' @param var_subgroup 1 sub-group variable for analysis, Default: NULL
#' @param var_cov Variables for additional adjust, Default: NULL
#' @param data Data or svydesign in survey package.
#' @param decimal.hr Decimal for hazard ratio, Default: 2
#' @param decimal.percent Decimal for percent, Default: 1
#' @param decimal.pvalue Decimal for pvalue, Default: 3
#' @return Sub-group analysis table.
#' @details This result is used to make forestplot.
#' @examples
# library(survival);library(dplyr)
# lung %>%
#   mutate(status = as.integer(status == 1),
#          sex = factor(sex),
#          kk = factor(as.integer(pat.karno >= 70)),
#          kk1 = factor(as.integer(pat.karno >= 60))) -> lung
#
# lung %>%
#   mutate(age2 = cut(age, c(min(age), quantile(lung$age)[2], quantile(lung$age)[4], max(age)))) %>%
#   mutate(age2 = factor(age2)) ->
#   lung
#
# SubCoxConX(Surv(time, status) ~ age, var_subgroup = "kk", data = lung)
#
# SubMultiCoxConX(Surv(time, status) ~ age, var_subgroups = c("kk",'kk1','age2'), data = lung)
#
# SubCoxUniY(Surv(time, status) ~ age2, data = lung)
#
# formulas <- c(as.formula(Surv(time, status) ~ age2),
#               as.formula(Surv(time, status) ~ age2),
#               as.formula(Surv(time, status) ~ age2))
#
# SubCoxMulY(formulas, data = lung)
#
# SubCoxCatX(Surv(time, status) ~ sex, var_subgroup = 'kk', data = lung)
#
# SubMultiCoxCatX(Surv(time, status) ~ sex, var_subgroups = c("kk",'kk1','age2'), data = lung)

#' @seealso
#'  \code{\link[purrr]{safely}},\code{\link[purrr]{map}},\code{\link[purrr]{map2}}
#'  \code{\link[survival]{coxph}}
#'  \code{\link[survey]{svycoxph}}
#'  \code{\link[stats]{confint}}
#' @rdname Coxsubgroup
#' @importFrom purrr possibly map_dbl map map2
#' @importFrom dplyr select filter mutate bind_cols
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom survival coxph
#' @importFrom survey svycoxph regTermTest
#' @importFrom stats confint coefficients
#' @importFrom utils tail
#' @importFrom tableone CreateTableOne

library(survival)
library(dplyr)

#'@export
SubCoxConX <- function(formula, var_subgroup = NULL, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL

  possible_coxph <- purrr::possibly(survival::coxph, NA)
  possible_table <- purrr::possibly(table, NA)
  possible_confint <- purrr::possibly(stats::confint, NA)
  possible_rowone <- purrr::possibly(function(x){x[1, ]}, NA)
  possible_pv <- purrr::possibly(function(x){summary(x)[["coefficients"]][1, ] %>% tail(1)}, NA)

  if (!is.null(var_cov)){
    formula <- as.formula(paste0(deparse(formula), " + ", paste(var_cov, collapse = "+")))
  }

  if (is.null(var_subgroup)){

    model <- coxph(formula, data = data)
    HR <- round(exp(coef(model)), decimal.hr)[1]
    CI <- round(exp(confint(model)[1, ]), decimal.CI)
    pv <- round(summary(model)$coefficients[1, "Pr(>|z|)"], decimal_pvalue)


    tibble::tibble(Variable = 'Overall', HR = HR, Lower = CI[1], Upper = CI[2]) %>%
      mutate('Pvalue' = ifelse(pv >= 0.001, pv, "<0.001")) -> out.all
    return(out.all)
  }else{

    data %>% filter(!is.na(get(var_subgroup))) %>% split(.[[var_subgroup]]) %>% purrr::map(~possible_coxph(formula, data = ., x=T)) -> model
    # data %>% filter(!is.na(data[var_subgroup])) %>% split(data[[var_subgroup]]) %>% purrr::map(~possible_coxph(formula, data = ., x=T)) -> model

    data %>% filter(!is.na(get(var_subgroup))) %>% select(var_subgroup) %>% table %>% names -> label_val
    # data %>% filter(!is.na(data[var_subgroup])) %>% select(var_subgroup) %>% table %>% names -> label_val

    xlev <- survival::coxph(formula, data = data)$xlevels
    xlabel <- names(attr(model[[which(!is.na(model))[1]]]$x, 'contrast'))[1]

    model %>% purrr::map('coefficients', .default = NA) %>% purrr::map_dbl(1) %>% exp %>%  round(decimal.hr) -> HR
    model %>% purrr::map(possible_confint) %>% purrr::map(possible_rowone) %>% Reduce(rbind, .) %>% exp %>% round(decimal.CI) -> CI
    model %>% purrr::map(possible_pv) %>% purrr::map_dbl(~round(., decimal_pvalue)) -> pv

    tibble::tibble(Variable = paste(" ", label_val), HR = HR, Lower = CI[,1], Upper = CI[,2]) %>%
            mutate('Pvalue' = ifelse(pv >= 0.001, pv, "<0.001")) -> out

    table <- print(tableone::CreateTableOne(vars = var_subgroup, data = data), showAllLevels = T)
    table[-1,2] %>% tibble::tibble() -> prop

    out.final <- rbind(c(var_subgroup, rep(NA, ncol(out) - 1)), out) %>%
    cbind(rbind(c(NA),prop[1])) %>% rename(Prop = ".")
    return(out.final)
  }

}

#' @export
SubMultiCoxConX <- function(formula, var_subgroups = NULL, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL
  out.all <- SubCoxConX(formula, var_subgroup = NULL, var_cov = var_cov, data = data, decimal.hr = decimal.hr, decimal.CI=decimal.CI, decimal_pvalue= decimal_pvalue)

  table.all <- print(tableone::CreateTableOne(vars = var_subgroups, data = data), showAllLevels = T)
  prop.all <- table.all[1,2] %>% tibble::tibble()

  cbind(out.all,  prop.all[1]) %>% rename(Prop = ".") -> out.all2

  outlist <- purrr::map(var_subgroups, ~ SubCoxConX(formula, var_subgroup = ., var_cov = var_cov, data = data, decimal.hr = decimal.hr, decimal.CI=decimal.CI, decimal_pvalue= decimal_pvalue))

  return(rbind(out.all2, outlist %>% dplyr::bind_rows()))
}




#' @export
SubCoxUniY <- function(formula, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL


  possible_coxph <- purrr::possibly(survival::coxph, NA)
  possible_table <- purrr::possibly(table, NA)
  possible_confint <- purrr::possibly(stats::confint, NA)
  possible_rowone <- purrr::possibly(function(x){x[1, ]}, NA)
  possible_pv <- purrr::possibly(function(x){summary(x)[["coefficients"]][1, ] %>% tail(1)}, NA)

  if (!is.null(var_cov)){
    formula <- as.formula(paste0(deparse(formula), " + ", paste(var_cov, collapse = "+")))
  }

  model <- coxph(formula, data = data)
  xlev <- names(model$xlevels)
  HR <- round(exp(coef(model)), decimal.hr)
  CI <- round(exp(confint(model)), decimal.CI)
  pv <- round(summary(model)$coefficients[, "Pr(>|z|)"], decimal_pvalue)
  tibble::tibble('Variable' = xlev, "HR" = HR, Lower = CI[,1], Upper = CI[,2]) %>%
    mutate('Pvalue' = ifelse(pv >= 0.001, pv, "<0.001")) -> out1

  get_variable_name <- function(x){
    name <- paste0(xlev, '-', x)
    return(name)
  }

  total_name <- c(0:nrow(out1)+1)
  purrr::map(total_name,get_variable_name) -> names

  out.all <- rbind(c(NA,'1',1,1,NA), out1)
  out.all$Variable <- names %>% unlist()

  formula_status <- formula[[2]][[3]]
  out.all2 <- rbind(c(as.character(formula_status),NA,NA,NA,NA), out.all)
  return(out.all2)
}



#' @export
SubCoxMulY <- function(formulas = NULL, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL


  outlist <- purrr::map(formulas, ~ SubCoxUniY(formula = ., var_cov = var_cov, data = data, decimal.hr = decimal.hr, decimal.CI=decimal.CI, decimal_pvalue= decimal_pvalue))

  return(outlist %>% dplyr::bind_rows())
}



#' @export
SubCoxCatX <- function(formula, var_subgroup = NULL, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL

  possible_coxph <- purrr::possibly(survival::coxph, NA)
  possible_table <- purrr::possibly(table, NA)
  possible_confint <- purrr::possibly(stats::confint, NA)
  possible_rowone <- purrr::possibly(function(x){x[1, ]}, NA)
  possible_pv <- purrr::possibly(function(x){summary(x)[["coefficients"]][1, ] %>% tail(1)}, NA)

  if (!is.null(var_cov)){
    formula <- as.formula(paste0(deparse(formula), " + ", paste(var_cov, collapse = "+")))
  }

  if (is.null(var_subgroup)){
    model <- coxph(formula, data = data)
    HR <- round(exp(coef(model)), decimal.hr)[1]
    CI <- round(exp(confint(model)[1, ]), decimal.CI)
    pv <- round(summary(model)$coefficients[1, "Pr(>|z|)"], decimal_pvalue)


    tibble::tibble(Variable = 'Overall', HR = HR, Lower = CI[1], Upper = CI[2]) %>%
      mutate('Pvalue' = ifelse(pv >= 0.001, pv, "<0.001"), `P for interaction` = NA) -> out.all


    return(out.all)
  }else{
    data %>% filter(!is.na(get(var_subgroup))) %>%split(.[[var_subgroup]]) %>% purrr::map(~possible_coxph(formula, data = ., x=T)) -> model

    data %>% filter(!is.na(get(var_subgroup))) %>% select(var_subgroup) %>% table %>% names -> label_val

    xlev <- survival::coxph(formula, data = data)$xlevels
    xlabel <- names(attr(model[[which(!is.na(model))[1]]]$x, 'contrast'))[1]
    model.int <- possible_coxph(as.formula(gsub(xlabel, paste0(xlabel, "*", var_subgroup), deparse(formula))), data = data)
    pvs_int <- model.int %>% summary %>% coefficients
    pv_int <- round(pvs_int[nrow(pvs_int), ncol(pvs_int)], decimal_pvalue)

    model %>% purrr::map('coefficients', .default = NA) %>% purrr::map_dbl(1) %>% exp %>%  round(decimal.hr) -> HR
    model %>% purrr::map(possible_confint) %>% purrr::map(possible_rowone) %>% Reduce(rbind, .) %>% exp %>% round(decimal.CI) -> CI
    model %>% purrr::map(possible_pv) %>% purrr::map_dbl(~round(., decimal_pvalue)) -> pv

    tibble::tibble(Variable = paste(" ", label_val), HR = HR, Lower = CI[,1], Upper = CI[,2]) %>%
      mutate('Pvalue' = ifelse(pv >= 0.001, pv, "<0.001"), `P for interaction` = NA) -> out

    table <- jstable::CreateTableOne2(vars = var_subgroup, strata = xlabel, data = data)
    table[-1,2:3] %>% tibble::tibble() -> prop

    out.prop <- rbind(c(var_subgroup, rep(NA, ncol(out) - 2), ifelse(pv_int >= 0.001, pv_int, "<0.001")), out) %>%
      cbind(rbind(c(NA, NA), prop[[1]])) %>% rename(With = `1`, Without = `2`)
    # return(rbind(c(var_subgroup, rep(NA, ncol(out) - 2), ifelse(pv_int >= 0.001, pv_int, "<0.001")), out))
    return(out.prop)
  }

}

#' @export
SubMultiCoxCatX <- function(formula, var_subgroups = NULL, var_cov = NULL, data, decimal.hr = 3, decimal.CI = 3, decimal_pvalue = 3){
  . <- NULL
  out.all <- SubCoxCatX(formula, var_subgroup = NULL, var_cov = var_cov, data = data, decimal.hr = decimal.hr, decimal.CI=decimal.CI, decimal_pvalue= decimal_pvalue)

  table.all <- jstable::CreateTableOne2(vars = var_subgroups, strata = as.character(formula[[3]]), data = data)
  prop.all <- table.all[1,2:3] %>% tibble::tibble()

  cbind(out.all,  t(prop.all[[1]])) %>% rename(With = `1`, Without = `2`) -> out.all2

  outlist <- purrr::map(var_subgroups, ~ SubCoxCatX(formula, var_subgroup = ., var_cov = var_cov, data = data, decimal.hr = decimal.hr, decimal.CI=decimal.CI, decimal_pvalue= decimal_pvalue))

  return(rbind(out.all2, outlist %>% dplyr::bind_rows()))
}

