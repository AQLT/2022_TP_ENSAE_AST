library(zoo)
library(tseries)
library(lmtest)
library(forecast)
library(ggplot2)
library(urca)

##########################
####### Question 1 #######
##########################
xm_tot <-read.csv("Data/Donnees1.csv")
# C'est une série mensuelle, on prendre une date de départ au hasard
xm_tot <- ts(xm_tot, start = 2003, frequency = 12)
# Garder les dernières valeurs
xm <- xm_tot[1:(length(xm_tot)-4)]
xm <- ts(xm, start = start(xm_tot),
         frequency = frequency(xm_tot))

##########################
####### Question 3 #######
##########################
#Trace le graphique
plot(xm, type = "l")
# On remarque qu'il y a une saisonnalité assez nette

# abscisse pas forcément très lisible pour les objets ts()
acf(xm) 
# Notons différence avec un objet numeric
acf(c(xm))

# On peut égalementutiliser la fonction du package forecast
# qui commence l'ACF à 1
forecast::Acf(xm)

# Il y a également une fonction avec une sortie ggplot2
ggAcf(xm) + labs(title = "ACF") 

# Ce qui permet de combiner facilement plusieurs graphiques :
library(patchwork)
autoplot(xm) / (ggAcf(xm) + labs(title = "ACF")  +
                  ggPacf(xm) + labs(title = "PACF") ) 

desaison <- diff(xm, 12) # fonction aussi avec un objet numeric
# ou de manière équivalente
# desaison <- xm - lag(xm,-12)

##########################
####### Question 3 #######
##########################
# a priori non intégrée : pas de décroissance linéaire,
# coef non proche de 1
acf(desaison)
pacf(desaison)
autoplot(desaison) / 
  (ggAcf(desaison) + labs(title = "ACF")  +
     ggPacf(desaison) + labs(title = "PACF") )


##########################
####### Question 4 #######
##########################
# Il n'y a pas de tendance
plot(desaison, type = "l")

# Plusieurs tests de racine unitaire. Dans le cours
# ADF
# PP
# KPSS

# Test ADF
m <- fUnitRoots::adfTest(desaison, type = "c")
# ou 
# m <- tseries::adf.test(desaison)
m
# RMQ : modèle bien spécifié... On verra ça TD suivant
# res = residuals(m@test$lm)
# Box.test(res,12, fitdf=length(m@test$lm$coefficients))

# RMQ2 : pas très grave de mal spécifier ce modèle
m <- fUnitRoots::adfTest(desaison, type = "ct")
m
summary(m@test$lm)

# PP test
tseries::pp.test(desaison) # test utilisé dans corrigé
summary(urca::ur.pp(desaison, type="Z-tau"))
summary(urca::ur.pp(desaison, type="Z-tau",model = "trend"))

# KPSS
tseries::kpss.test(desaison)


##########################
####### Question 5 #######
##########################

# acf(desaison) # suggère MA(2)
# pacf(desaison) # suggère AR(3)
(ggAcf(desaison) + labs(title = "ACF")  + # Suggère MA(2)
   ggPacf(desaison) + labs(title = "PACF") ) # Suggère AR(3)
# Comme la série n'est pas intégrée les ordres p, d, q maximum sont :
# ARIMA(3,0,2)

# RMQ : par défaut il y a une constante, pas besoin de centrer les séries
model_maxi <- arima(desaison, order = c(3,0,2))
residus_maxi <- residuals(model_maxi)
# Parait bon :
ggAcf(residus_maxi) + ggPacf(residus_maxi)

# Pour l'application du test de Ljung-Box sur les résidus du modèle ARIMA
# il est suggéré de considéré p+q degrés de libertés
lbtest <- t(sapply(1:24,function(l){
  if(l <=  length(coef(model_maxi))){
    b <- list(statistic = NA, p.value = NA)
  }else{
    b <- Box.test(residus_maxi,"Ljung-Box",lag = l,
                  # il faut ajuster du degré de liberté, i.e. : nombre de coefficients estimés
                  fitdf = length(coef(model_maxi))
    )
  }
  
  data.frame(lag = l,
             b$statistic,
             b$p.value
  )
}))
lbtest # Il n'y a pas autocorrélation des résidus
# Remarque:portes n'utilise pas le bon nombre de degré de liberté (devrait être égal à lags-6)
# Le problème est qu'il n'utilise pas la constante
portes::LjungBox(model_maxi)
# Autre exemple, si l'on rajoute une variable supplémentaire (par exemple une indicatrice choisit au hasard)
# alors c'est encore le même nombre de degré de libertés :
# faire attention en utilisant ce package
portes::LjungBox(arima(desaison, order = c(3,0,2),
                       xreg = time(desaison)==2018))

# Utiliser plutôt
portes::LjungBox(residuals(model_maxi),order = 6)

# On pourrait également tester heteroscédasticité en appliquant le test sur le carré des résidus
lb2test <- t(sapply(1:24,function(l){
  if(l <=  length(coef(model_maxi))){
    b <- list(statistic = NA, p.value = NA)
  }else{
    b <- Box.test(residus_maxi^2,"Ljung-Box",lag = l,
                  fitdf = length(coef(model_maxi))
    )
  }
  
  data.frame(lag = l,
             b$statistic,
             b$p.value
  )
}))
lb2test
# Et normalité (suppose résidus indépendants et homoscédastiques)
tseries::jarque.bera.test(residus_maxi)

##########################
####### Question 6 #######
##########################
lmtest::coeftest(model_maxi) # aucun coefficient n'est significatif à part la constante
evaluation_model <- function(order, x = desaison, lags = 24,...){
  # ici on utilise Arima plutôt que arima pour la fonction accuracy
  model <- forecast::Arima(x, order = order,...)
  residus <- residuals(model)
  lbtest <- t(sapply(1:lags,function(l){
    if(l <=  length(coef(model))){
      b <- list(statistic = NA, p.value = NA)
    }else{
      b <- Box.test(residus,"Ljung-Box",lag = l,
                    fitdf = length(coef(model))
      )
    }
    data.frame(lag = l,
               b$statistic,
               b$p.value
    )
  }))
  # on ajoute un tryCatch pour éviter les erreurs
  ttest <- tryCatch(lmtest::coeftest(model), error = \(e) NULL)
  qualite <- c(AIC(model), BIC(model), accuracy(model))
  names(qualite) <- c("AIC", "BIC", colnames(accuracy(model)))
  list(model = model,
       ttest = ttest,
       lbtest = lbtest,
       qualite = qualite)
  
}

models_possibles <- expand.grid(p = c(0,1,2,3), d = 0, q = c(0, 1, 2))
models_evalues <- apply(models_possibles,1, evaluation_model)
names(models_evalues) <- sprintf("ARIMA(%i,%i,%i)", models_possibles[,"p"],
                                 models_possibles[,"d"], models_possibles[,"q"])
## Pour éviter de tout écrire à la main :
#cat(paste(sprintf("models_evalues$`%s`",names(models_evalues)),collapse = "\n"))


models_evalues$`ARIMA(0,0,0)`
# il n'y a pas indépendance des résidus
models_evalues$`ARIMA(1,0,0)`
# il n'y a pas indépendance des résidus jusqu'à 10  -> non retenu corrigé
models_evalues$`ARIMA(2,0,0)`
# il y a indépendance mais AR(2) non signif : modèle non ajusté
models_evalues$`ARIMA(3,0,0)`
# Tout est parfait
models_evalues$`ARIMA(0,0,1)`
# il n'y a pas indépendance des résidus (à 10 %) -> non retenu corrigé
models_evalues$`ARIMA(1,0,1)`
# Ma(1) non signif
models_evalues$`ARIMA(2,0,1)`
# Modèle qui parait plutôt bien
models_evalues$`ARIMA(3,0,1)`
# Aucun coef significatif
models_evalues$`ARIMA(0,0,2)`
# Modèle qui parait bien
models_evalues$`ARIMA(1,0,2)`
# AR(1) non signif : mal spécifié
models_evalues$`ARIMA(2,0,2)`
# AR(2) non signif : mal spécifié
models_evalues$`ARIMA(3,0,2)`
# AR(3) non signif : mal spécifié

nom_modeles_valides <- c("ARIMA(3,0,0)", "ARIMA(2,0,1)", "ARIMA(0,0,2)",
                         "ARIMA(1,0,0)", "ARIMA(0,0,1)" # deux modèles non retenus par le corrigé
)
models_valides <- models_evalues[nom_modeles_valides]
qualite_modeles_valides <- sapply(models_valides, function(x) x$qualite)
round(qualite_modeles_valides,4)
apply(qualite_modeles_valides,1,function(x) colnames(qualite_modeles_valides)[which.min(x)])
# En fonction du critère on ne choisit donc pas forcément le même modèle mais sur les critères 
# d'information "ARIMA(0,0,2)" meilleur

# NB : pas le même résultat avec auto.arima
auto.arima(desaison) 
auto.arima(xm) 
# Même modèle retenu avec TRAMO mais avec une transformation au log :
# TRAMO est utilisé pour faire de la désaisonnalisation
# RJDemetra::regarima_tramoseats(xm)

##########################
####### Question 7 #######
##########################
lapply(models_valides,function(x) forecast(x$model, h = 4))
prev <- lapply(models_valides,function(x){
  forecast(x$model, h = 4)$mean + tail(xm,12)[1:4] # on rajoute le lag12
} )
rmse <- sapply(prev,function(x)sqrt(mean((x- tail(xm_tot,4))^2)))
rmse

# Autre façon de faire en faisant directement un modèle ARIMA
models_possibles <- data.frame(p = c(3,2,0,1,0), d = 0, q = c(0,1,2,0,1))
models_evalues_2 <- apply(models_possibles,1, evaluation_model, x = xm, include.constant = TRUE,seasonal = c(0,1,0))
names(models_evalues_2) <- sprintf("ARIMA(%i,%i,%i)(0,1,0)", models_possibles[,"p"],
                                   models_possibles[,"d"], models_possibles[,"q"])
prev2 <- lapply(models_evalues_2,function(x) forecast(x$model, h = 4)$mean)
rmse2 <- sapply(prev2,function(x)sqrt(mean((x- tail(xm_tot,4))^2)))
rmse2


##########################
####### Question 8 #######
##########################

xm_tot <-read.csv("Data/Donnees2.csv")
# C'est une série mensuelle, on prendre une date de départ au hasard
xm_tot <- ts(xm_tot, start = 2003, frequency = 12)
# Garder les dernières valeurs
xm <- xm_tot[1:(length(xm_tot)-4)]
xm <- ts(xm, start = start(xm_tot),
         frequency = frequency(xm_tot))

# Il y a clairement une tendance linéaire
# La question est de savoir si la tendance est juste trend-stationnaire ou 
# s'il y a également une tendance stochastique
autoplot(xm) / 
  (ggAcf(xm) + labs(title = "ACF")  +
     ggPacf(xm) + labs(title = "PACF") )

autoplot(diff(xm,1)) / 
  (ggAcf(diff(xm,1)) + labs(title = "ACF")  +
     ggPacf(diff(xm,1)) + labs(title = "PACF") )

# Ici on teste modèle

# On teste ici le modèle
# ∆y_t = a + bt + γ y_t-1 +e_t
# tau3 correspond au test γ = 0
# phi2 correspond au test a = b = γ = 0
# phi3 correspond au test b = γ = 0
# voir https://new.mmf.lnu.edu.ua/wp-content/uploads/2018/03/enders_applied_econometric_time_series.pdf
summary(urca::ur.df(xm, type  = "trend",lags = 12,
                    selectlags = "AIC"))
# on rejette si stat < critical value
# Ici on ne rejette pas

# rmq : 
fUnitRoots::adfTest(xm, type = "ct",lags = 1)  # rejette
fUnitRoots::adfTest(xm, type = "ct",lags = 4)  # rejette pas : non stationnaire

tseries::pp.test(xm) # on rejette hypothèse de non stationnaire
tseries::kpss.test(xm, null = "Trend") # rejette hypothèse de stationnarité

# On va plutôt suivre PP (plus)

# À partir de R 4.1:
desaison = lm(xm ~ time(xm)) |> 
  residuals() |> 
  ts(start = start(xm), frequency = frequency(xm))
# Sinon utiliser ce code:
desaison = ts(residuals(
  lm(xm ~ time(xm))
),
start = start(xm),
frequency = frequency(xm)
)
tseries::pp.test(desaison) # on rejette hypothèse de non stationnaire
tseries::kpss.test(xm, null = "Level") ## rejette hypothèse de stationnarité
summary(urca::ur.df(xm, type  = "none",lags = 12,
                    selectlags = "AIC")) ## rejette pas hyp de non stationnarité

##########################
####### Question 8.5 #######
##########################

acf(desaison) # suggère MA(22)
pacf(desaison) # suggère AR(4)
autoplot(desaison) / 
  (ggAcf(desaison) + labs(title = "ACF")  + # MA(24)
     ggPacf(desaison) + labs(title = "PACF") ) # AR4
# Cela suggère plutôt un pur AR
# voir : 
y = ts(arima.sim(list(ar = c(0.8, -0.1334, 0.1, 0.2)), n = 200),
       start = 2000, frequency = 12)
autoplot(y) /
  (ggAcf(y) + labs(title = "ACF")  +
     ggPacf(y) + labs(title = "PACF") )

# Se restreindre à moins d'ordres MA que 24
models_possibles <- expand.grid(p = c(1,2,3,4), d = 0, q = c(0,1,2))
models_evalues <- apply(models_possibles,1, evaluation_model,include.mean = FALSE)
names(models_evalues) <- sprintf("ARIMA(%i,%i,%i)", models_possibles[,"p"],
                                 models_possibles[,"d"], models_possibles[,"q"])
## Pour éviter de tout écrire à la main :
# cat(paste(sprintf("models_evalues$`%s`",names(models_evalues)),collapse = "\n"))

models_evalues$`ARIMA(1,0,0)` # autocorrelation
models_evalues$`ARIMA(2,0,0)` # possible
models_evalues$`ARIMA(3,0,0)` # non signif AR(3)
models_evalues$`ARIMA(4,0,0)` # possible
models_evalues$`ARIMA(1,0,1)` # possible
models_evalues$`ARIMA(2,0,1)` # possible
models_evalues$`ARIMA(3,0,1)` # AR(3) non signif
models_evalues$`ARIMA(4,0,1)` # AR(4) non signif
models_evalues$`ARIMA(1,0,2)` # MA(2) non signif
models_evalues$`ARIMA(2,0,2)` # AR(2) et MA(2) non signif
models_evalues$`ARIMA(3,0,2)` # coefs non signifs
models_evalues$`ARIMA(4,0,2)` # coefs non signifs


nom_modeles_valides <- c("ARIMA(2,0,0)", "ARIMA(4,0,0)", "ARIMA(1,0,1)",
                         "ARIMA(2,0,1)"
)
models_valides <- models_evalues[nom_modeles_valides]
qualite_modeles_valides <- sapply(models_valides, function(x) x$qualite)
round(qualite_modeles_valides,4)
apply(qualite_modeles_valides,1,function(x) colnames(qualite_modeles_valides)[which.min(x)])

models_possibles <- data.frame(p = c(2,4,1,2), d = 0, q = c(0,0,1,1))
models_evalues_2 <- apply(models_possibles,1, evaluation_model,
                          x = xm,  xreg = time(xm))
names(models_evalues_2) <- sprintf("ARIMA(%i,%i,%i)", models_possibles[,"p"],
                                   models_possibles[,"d"], models_possibles[,"q"])
prev2 <- lapply(models_evalues_2,function(x) forecast(x$model, h = 4,xreg = tail(time(xm_tot),4))$mean)
# Rmq, coefs constante + tendance assez proche de lm(xm ~ time(xm))
rmse2 <- sapply(prev2,function(x)sqrt(mean((x- tail(xm_tot,4))^2)))
rmse2
qualite_modeles_valides2 = sapply(models_evalues_2, function(x) x$qualite)
round(qualite_modeles_valides2,4)
apply(qualite_modeles_valides2,1,function(x) colnames(qualite_modeles_valides2)[which.min(x)])

auto.arima(xm, max.D = 0, max.P = 0, max.Q = 0, 
           max.d = 0 # Sinon on retient un retard car algo utilise KPSS
)
