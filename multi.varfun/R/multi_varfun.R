#' Multilevel Variance Function Regression
#'
#'This package run multilevel variance function regression.
#' @param outcome Outcome variable.
#' @param indep List of independent variables (e.g., "age+sex+edu+gdp").
#' @param level2 ID-Variable for level 2 units (e.g., countryID).
#' @param level3 ID-Variable for level 3 units (e.g., countryID). Defaults to 1 denoting estimation of 2 level model.
#' @param rslope2 List of variables for which random slopes are estimated for level 2 (e.g., "edu+age"). Defaults to 1 denoting no random slope estimation.
#' @param rslope3 List of variables for which random slopes are estimated for level 3 (e.g., "edu+age"). Defaults to 1 denoting no random slope estimation.
#' @param data Name of dataset.
#' @keywords multilevel
#' @keywords variance function regression
#' @export
#' 
#' @examples 
#' ## Ethnic disparities in math competences: 
#' ethn_math <- multi_varfun(outcome="math", indep="migrant+edu_f+age+grade+sex", level2="cid", rslope="migrant", data=PISA)
#' summary(ethn_math$meanmodel)
#' summary(ethn_math$varmod)
#' 
#' @import lme4
#' @author Christoph Spörlein, \email{christoph.spoerlein@@uni-bamberg.de}
#' @references Western and Bloome 2009: Variance Function Regression for studying Inequality. Sociological Methodology 39:293-326.


multi_varfun <- function(outcome, indep=1, level2, level3=1, rslope2=1, rslope3=1, data){

  if (level3==1) {

  # keep relevant data
  form0 <- paste0(outcome,"~",indep,"+",level2)
  mod <- lm(form0, data=data)
  datanomiss <- mod$model

  form <- paste0(outcome,"~",indep,"+(1+",rslope2,"|",level2,")")
  form2 <- paste0("resmod~",indep,"+(1+",rslope2,"|",level2,")")

  meanmodel <- lmer(formula=form, data=datanomiss, REML=F)
  datanomiss$resmod <- residuals(meanmodel)^2
  varmod <- glmer(formula=form2, data=datanomiss, family=Gamma(link = "log"), nAGQ=0)
  datanomiss$fitmod <- fitted(varmod)

  datanomiss$LOGLIK <- -.5*(log(datanomiss$fitmod)+(datanomiss$resmod/datanomiss$fitmod))

  LLO <- sum(datanomiss$LOGLIK)

  DLL <- 1
  datanomiss$weight <- 1/datanomiss$fitmod


  #### loop

  while(DLL > .0001){
  model <- lmer(formula=form, weights=weight, data=datanomiss, REML=F)
  datanomiss$resmod <- residuals(model)^2
  varmod <- glmer(formula=form2, data=datanomiss, family=Gamma(link = "log"))
  datanomiss$fitmod <- fitted(varmod)
  datanomiss$weight <- 1/datanomiss$fitmod
  datanomiss$LOGLIK <- -.5*(log(datanomiss$fitmod)+(datanomiss$resmod/datanomiss$fitmod))
  LLN <- sum(datanomiss$LOGLIK)
  DLL <- LLN-LLO
  LLO <- LLN
  print(DLL)
  }

  results <- list("meanmodel"=meanmodel,"varmod"=varmod, "model"=model)
  return(results)

 } else {

  form0 <- paste0(outcome,"~",indep,"+",level2,"+",level3)
  mod <- lm(form0, data=data)
  datanomiss <- mod$model

  form <- paste0(outcome,"~",indep,"+(1+",rslope2,"|",level2,")+(1+",rslope3,"|",level3,")")
  form2 <- paste0("resmod~",indep,"+(1+",rslope2,"|",level2,")+(1+",rslope3,"|",level3,")")

  meanmodel <- lmer(formula=form, data=datanomiss, REML=F)
  datanomiss$resmod <- residuals(meanmodel)^2
  varmod <- glmer(formula=form2, data=datanomiss, family=Gamma(link = "log"), nAGQ=0)
  datanomiss$fitmod <- fitted(varmod)

  datanomiss$LOGLIK <- -.5*(log(datanomiss$fitmod)+(datanomiss$resmod/datanomiss$fitmod))

  LLO <- sum(datanomiss$LOGLIK)

  DLL <- 1
  datanomiss$weight <- 1/datanomiss$fitmod


  #### loop

  while(DLL > .0001){
  model <- lmer(formula=form, weights=weight, data=datanomiss, REML=F)
  datanomiss$resmod <- residuals(model)^2
  varmod <- glmer(formula=form2, data=datanomiss, family=Gamma(link = "log"))
  datanomiss$fitmod <- fitted(varmod)
  datanomiss$weight <- 1/datanomiss$fitmod
  datanomiss$LOGLIK <- -.5*(log(datanomiss$fitmod)+(datanomiss$resmod/datanomiss$fitmod))
  LLN <- sum(datanomiss$LOGLIK)
  DLL <- LLN-LLO
  LLO <- LLN
  print(DLL)
  }

  results <- list("meanmodel"=meanmodel,"varmod"=varmod, "model"=model)
  return(results)

}

}











