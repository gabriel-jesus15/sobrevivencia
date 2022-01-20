# Os dados sao referentes aos tempos de recidiva de pacientes que realizaram cirurgia 
# ULTDIAG: Tempos de recidiva ou censura,
# VCAN:Indicador de censura

library(survival)
library(dplyr)
library(ggfortify)
library(ranger)
install.packages('flexsurv')

#Carregamento dos dados
dados <- BD_AVR_PUL

dados$CATEATEND <- coalesce(dados$CATEATEND, 9)

shapiro.test(dados$ULTDIAG)

mKM_pul <- survfit(Surv(dados$ULTDIAG, dados$VCAN) ~ dados$CIRURGIA, se.fit = FALSE)
plot(mKM_pul, mark.time=T,xlab = "Tempo", ylab = "Survival function", cex.axis = 1.5,
     cex.lab = 1.5,lty=1:2)

dados2 <- na.omit(data.frame(tempos = dados$ULTDIAG, censuras = dados$VCAN, sexo = dados$SEXO, cate = dados$CATEATEND,
                    trat = dados$TRATAMENTO, diag = dados$DIAGPREV))

KM <- survfit(Surv(tempos, censuras) ~ as.factor(sexo) + as.factor(cate) 
        + as.factor(trat) + as.factor(diag), data = dados2, se.fit = FALSE)

KM <- survfit(Surv(tempos, censuras) ~ sexo + cate, data = dados2, se.fit = FALSE)

autoplot(KM)

### Curva de sobrevivencia para os dados Completos

# Kaplan-Meier Padrao

KM1 <- survfit(Surv(tempos, censuras) ~ 1, data = dados2)
plot(KM1)

### Curvas de sobrevivencia para as variaveis independentes

# Sexo
KM_sexo <- survfit(Surv(tempos, censuras) ~ sexo, data = dados2)
autoplot(KM_sexo)

survdiff(Surv(tempos, censuras) ~ sexo, data = dados2, rho= 0)

# Cate
KM_cate <- survfit(Surv(tempos, censuras) ~ cate, data = dados2)
autoplot(KM_cate)

survdiff(Surv(tempos, censuras) ~ cate, data = dados2, rho= 0)

# trat
KM_trat <- survfit(Surv(tempos, censuras) ~ trat, data = dados2)
plot(KM_trat)

survdiff(Surv(tempos, censuras) ~ trat, data = dados2, rho= 0)

# Diag
KM_diag <- survfit(Surv(tempos, censuras) ~ diag, data = dados2)
autoplot(KM_diag)

survdiff(Surv(tempos, censuras) ~ diag, data = dados2, rho= 0)

### Metodos Parametricos

## Escolha da Distribuicao

ajuste1 <- survreg(Surv(tempos, censuras) ~ 1, data = dados2, dist="exponential")
alpha1 <- exp(ajuste1$coefficients[1])
ajuste2 <- survreg(Surv(tempos, censuras) ~ 1, data = dados2, dist="weibull")
alpha2 <- exp(ajuste2$coefficients[1]);
gama <- 1/ajuste2$scale
ajuste3 <- survreg(Surv(tempos, censuras) ~ 1, data = dados2, dist="lognorm")

#S(T) VS ESTIMATIVAS

ekm <- survfit(Surv(tempos, censuras) ~ 1, data = dados2)
time <- ekm$time
st <- ekm$surv
ste <- exp(-time/alpha1)
stw <- exp(-(time/alpha2)^gama)
stln <- pnorm((-log(time)+ 4.778)/2.0347)
cbind(time,st,ste,stw,stln)

## Escolhendo a melhor distribuiÃ§Ã£o para o modelo de sobrevivencia pelo teste TRV

plot(KM1, conf.int=F, xlab="Tempos", ylab="S(t)", main= "Kaplan-Meier vs Exponencial vs Weibull vs LogNormal" )
lines(c(0,time),c(1,ste), lty=2, col='red')
lines(c(0,time),c(1,stw), lty=3, col='blue')
lines(c(0,time),c(1,stln), lty=4, col='green')
legend(5000,0.99,lty=c(1,2),c("Kaplan-Meier", "Exponencial", "Weibull", "LogNormal")
       ,bty="n",cex=0.8, col=c("black",'red','blue','green'))

### Ajuste do melhor modelo com todas as covariaveis escolhidas

# Geral
wei <- survreg(Surv(tempos, censuras)~ sexo + cate + trat + diag, data=dados2, dist="weibull")

# 1

wei_1 = survreg(Surv(tempos, censuras) ~ 1, data=dados2,dist="weibull")
summary(wei_1)

# Sexo

wei_sexo = survreg(Surv(tempos, censuras) ~ sexo, data=dados2, dist="weibull")
summary(wei_sexo)
# p valor < 0.1

# Cate

wei_cate = survreg(Surv(tempos, censuras) ~ cate, data=dados2, dist="weibull")
summary(wei_cate)

# Trat

wei_trat = survreg(Surv(tempos,censuras) ~ trat, data=dados2, dist="weibull")
summary(wei_trat)

# Diag

wei_diag = survreg(Surv(tempos, censuras) ~ diag, data=dados2, dist="weibull")
summary(wei_diag)

# Retirado sexo

wei_red <- survreg(Surv(tempos, censuras) ~ cate + trat + diag, data=dados2, dist="weibull")

### Testando para cada remoção

# - cate

wei_1 <- survreg(Surv(tempos, censuras) ~ trat + diag, data = dados2, dist="weibull")
anova(wei_red,wei_1)

# - trat

wei_2 <- survreg(Surv(tempos, censuras) ~ cate + diag, data = dados2, dist="weibull")
anova(wei_red,wei_2)

# - diag

wei_3 <- survreg(Surv(tempos, censuras) ~ cate + trat, data = dados2, dist="weibull")
anova(wei_red,wei_3)
# p valor < 0.1

### Comparando segunda remoção

wei_red2 <- survreg(Surv(tempos,censuras) ~ cate + trat, data = dados2, dist="weibull")

anova(wei,wei_red) # ok remover sexo

anova(wei,wei_red2) # problema ao remover diag 

### Comparando por interação

wei_prod1 <- survreg(Surv(tempos, censuras) ~ cate + trat + diag + (cate * trat), data = dados2, dist="weibull")
anova(wei_prod1 ,wei_red)
# p valor < alpha

wei_prod2 <- survreg(Surv(tempos, censuras) ~ cate + trat + diag + (cate * diag), data = dados2, dist="weibull")
anova(wei_prod2 ,wei_red)

wei_prod3 <- survreg(Surv(tempos, censuras) ~ cate + trat + diag + (diag * trat), data = dados2, dist="weibull")
anova(wei_prod3 ,wei_red)
# p valor < alpha

modelofinal <- survreg(Surv(tempos, censuras) ~ cate + trat + diag + (cate * trat) + (diag * trat), data = dados2, dist="weibull")
summary(modelofinal)

AIC(modelofinal)

modelofinal$coefficients

### Escolhendo melhor fit

fit1 <- coxph(Surv(tempos, censuras) ~ cate + trat + diag, data = dados2, method = 'breslow')
# diag alto

fit3 <- coxph(Surv(tempos, censuras) ~ cate + trat, data = dados2, method = 'breslow')
# sem problemas

fit4 <- coxph(Surv(tempos, censuras) ~ cate + trat + cate*trat, data = dados2, method = 'breslow')
# Algumas interações > q a = 0.1

fit_final <- fit3

AIC(fit_final)

summary(fit_final)

fit_final$loglik

zph <- cox.zph(fit_final)
zph 

# As variáveis rejeita a hipotese nula de proporcionalidade.

### Residuos 

plot(zph)

mresid <- resid(fit_final); 
plot(mresid)
