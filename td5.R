library(zoo)
library(tseries)
library(lmtest)
library(forecast)
library(ggplot2)
library(patchwork)

# Correction Q1 à 3 disponible https://github.com/AQLT/TP_Ensae2A/wiki

##########################
####### Question 1 #######
##########################
data <-read.csv("Data/data_tp5.csv",sep = ";")
data$dates[1] 
spread <- ts(data$spread, start = c(1986, 3), frequency = 12)
dspread <- diff(spread,1) # difference premiere


##########################
####### Question 2 #######
##########################
plot(cbind(spread,dspread))
# La série en niveau semble avoir une tendance déterministe (ou deux tendances)
# La série différenciée pourrait être stationnaire


##########################
####### Question 3 #######
##########################

# rmq : tester la présence d'une tendance par régression n'a pas de sens
# car si on a une racine unitaire on est dans le cas d'une spurious régression
y = cumsum(rnorm(n=100))
summary(lm(y ~ seq_along(y)))

library(urca)
library(fUnitRoots)
# Ici on teste la présence de racine unitaire
# adf invalid si on ne rajoute pas de variable explicative
# dans le doute c'est toujours mieux de rajouter la tendance et constante
adf <- adfTest(spread, type = "ct",lags = 0)

Qtests <- function(series, k = 24, fitdf=0) {
  t(sapply(1:k,function(l){
    if(l <= fitdf){
      b <- list(statistic = NA, p.value = NA)
    }else{
      b <- Box.test(series,"Ljung-Box",lag = l,
                    fitdf = fitdf
      )
    }
    data.frame(lag = l,
               b$statistic,
               b$p.value
    )
  }))
}
adfTest_valid <- function(series,kmax,type){
  k <- 0
  # On test ADF jusqu'à ce que les résidus ne soient pas autocorrélés
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,
                    24,
                    fitdf=length(adf@test$lm$coefficients))[,3]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")
    }else {
      cat("nope \n")
    }
    k <- k + 1
  }
  return(adf)
}

# On ne fait que le test à l'ordre 24
adfTest_valid2 <- function(series,kmax,type){
  k <- 0
  # On test ADF jusqu'à ce que les résidus ne soient pas autocorrélés
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,
                    24,
                    fitdf=length(adf@test$lm$coefficients))[,3]
    if (sum(pvals[24]<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")
    }else {
      cat("nope \n")
    }
    k <- k + 1
  }
  return(adf)
}

df1 <- adfTest_valid(spread,24,"ct")
df2 <- adfTest_valid2(spread,24,"ct")

# On teste ici le modèle
# ∆y_t = a + bt + γ y_t-1 +e_t
# tau3 correspond au test γ = 0
# phi2 correspond au test a = b = γ = 0
# phi3 correspond au test b = γ = 0
# voir https://new.mmf.lnu.edu.ua/wp-content/uploads/2018/03/enders_applied_econometric_time_series.pdf
summary(
  urca::ur.df(spread, type  = "trend",lags = 12,
              selectlags = "AIC")
) # 5 lags retenus par AIC

# Dans tous cas on ne rejette pas H0 (mais a priori pas de tendance linéaire)

# PP n'a pas besoin de beaucoup d'hypothèses parce que la stat calculée
# sur beta est construite de façon non paramétrique et qui corrige
# toute forme d'autocorrélation. Mais pas très bon quand peu de données
tseries::pp.test(spread) # on rejette pas hypothèse de présence de racine unitaire
tseries::kpss.test(spread, null = "Trend") # on rejette hypothèse de stationnarité

autoplot(dspread) / 
  (ggAcf(dspread) + labs(title = "ACF")  +
     ggPacf(dspread) + labs(title = "PACF") )

tseries::pp.test(dspread) # on rejette pas hypothèse de présence de racine unitaire
tseries::kpss.test(dspread, null = "Level")
df1 <- adfTest_valid(dspread,24,"c")
df2 <- adfTest_valid2(dspread,24,"c")
summary(
  urca::ur.df(dspread, type  = "none",lags = 12,
              selectlags = "AIC")
)


##########################
####### Question 4 #######
##########################

autoplot(dspread) / 
  (ggAcf(dspread) + labs(title = "ACF")  +
     ggPacf(dspread) + labs(title = "PACF") )

##########################
####### Question 5 #######
##########################

models_possibles <- expand.grid(p = c(0,1,2,3), d = 1, q = c(0, 1, 2, 3))

evaluation_model <- function(order, x = spread){
  # ici on utilise Arima plutôt que arima pour la fonction accuracy
  model <- forecast::Arima(x, order = order,include.constant = FALSE)
  qualite <- c(AIC(model), BIC(model), accuracy(model))
  names(qualite) <- c("AIC", "BIC", colnames(accuracy(model)))
  qualite
}

all_models <- data.frame(t(apply(models_possibles,1,evaluation_model)))
rownames (all_models) <- sprintf("ARIMA(%i,%i,%i)", models_possibles[,"p"],
                                 models_possibles[,"d"], models_possibles[,"q"])
all_models

rownames(all_models)[which.min(all_models$AIC)]
rownames(all_models)[which.min(all_models$BIC)]
rownames(all_models)[which.min(all_models$RMSE)]
apply(all_models,2,function(x) rownames(all_models)[which.min(x)])


##########################
####### Question 6 #######
##########################

# On peut retenir les deux modèles qui minimises AIC et BIC
# Rmq : le modèle retenu par auto.arima est un ARIMA(2,1,2) qui ne minimise pas l'AIC
auto.arima(spread,max.D = 0, max.P = 0, max.Q = 0)

arima310 <- arima(spread,c(3,1,0),include.mean=FALSE)
arima011 <- arima(spread,c(0,1,1),include.mean=FALSE)
lmtest::coeftest(
  arima310
) # modèle bien ajusté coef AR(3) significatif
lmtest::coeftest(
  arima011
) # modèle bien ajusté coef MA(1) significatif


##########################
####### Question 7 #######
##########################

Qtests(residuals(
  arima310
),
fitdf = 3)

Qtests(residuals(
  arima011
),
fitdf = 3)

##########################
####### Question 9 #######
##########################

(autoplot(residuals(arima310)) +
   geom_vline(xintercept = 2001+11/12,linetype="dashed",
              col = "red", alpha = 0.7) + 
   labs(title = "ARIMA(3,1,0)")) /
  (autoplot(residuals(arima011)) +
     geom_vline(xintercept = 2001+11/12,linetype="dashed",
                col = "red", alpha = 0.7)+ 
     labs(title = "ARIMA(0,1,1)"))
# Il y a vraisemblablement un point atypique

arima310_ind <- arima(spread,c(3,1,0),include.mean=FALSE,xreg = time(spread) == 2001+11/12)
arima011_ind <- arima(spread,c(0,1,1),include.mean=FALSE,xreg = time(spread) == 2001+11/12)

# L'indicatrice est significative
lmtest::coeftest(
  arima310_ind
) # modèle bien ajusté coef AR(3) significatif
lmtest::coeftest(
  arima011_ind
)

(autoplot(residuals(arima310_ind)) +
    geom_vline(xintercept = 2001+11/12,linetype="dashed",
               col = "red", alpha = 0.7) + 
    labs(title = "ARIMA(3,1,0)")) /
  (autoplot(residuals(arima011_ind)) +
     geom_vline(xintercept = 2001+11/12,linetype="dashed",
                col = "red", alpha = 0.7)+ 
     labs(title = "ARIMA(0,1,1)"))

##########################
####### Question 10 #######
##########################

# On peut faire un test de chow
# Voir https://cran.r-project.org/web/packages/strucchange/vignettes/strucchange-intro.pdf
# Voir également https://robjhyndman.com/hyndsight/structural-breaks/
# Code à vérifier
test_rupture <- function(order, date_break = 1995){
  glob_mod = arima(spread,
                   order = order, include.mean = FALSE)
  rss <- sum(residuals(glob_mod)^2)
  sigma2 <- glob_mod$sigma2
  k <- length(coef(glob_mod))
  n <- length(spread) - order[2]
  fit1 <- arima(window(spread, end = date_break),
                order = order, include.mean = FALSE)
  fit2 <- arima(window(spread, start = date_break+1/12),
                order = order, include.mean = FALSE)
  ess <- sum(c(residuals(fit1), residuals(fit2))^2)
  stats <- (rss - ess)/(n-2*k)
  stats <- stats/k
  c(stats = stats,
    qf = qf(0.05, df1 = k, df2 = n-2*k, lower.tail = FALSE),
    pval = 1- pf(stats, df1 = k, df2 = n-2*k)
  )
}

round(test_rupture(order = c(0,1,1),date_break = 2001+11/12), 3)
round(test_rupture(order = c(3,1,0),date_break = 2001+11/12), 3)

#############################
####### TP avec fable #######
#############################
library(fable)
library(lubridate)
library(dplyr)
library(feasts)
library(ggplot2)

fable_mod = lapply(sprintf("(value) ~ 0 + pdq(%i,%i,%i) + PDQ(0,0,0) ", models_possibles[,"p"],
                           models_possibles[,"d"], models_possibles[,"q"]), as.formula)
fable_mod = lapply(fable_mod, ARIMA)
names(fable_mod) <- rownames (all_models)
spread_ts = spread %>% 
  as_tsibble()
all_mod = do.call(function (...) model(spread_ts, ...), fable_mod)

all_mod %>%
  report() 

all_mod %>%
  report() %>% 
  mutate(across(AIC:BIC, ~ .model[which.min(.x)])) %>% 
  distinct(AIC, AICc, BIC)
