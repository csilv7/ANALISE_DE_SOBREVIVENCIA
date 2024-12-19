# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
# Análise de Sobrevivência - Modelagem de Dados com Censura
# Analista: Breno C R da Silva
# Data: 19/12/2024
# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------

# -----------------------
# [1] Ativação de Pacotes
# -----------------------

if (!require("survival")) install.packages("survival")
if (!require("ggplot2")) install.packages("ggplot2")

library(survival)
library(ggplot2)

# -----------------------
# [2] Simulação dos Dados
# -----------------------

# Definindo Semente
set.seed(123)

# Tamanho amostral
n <- 1000

# Parâmetros - Distribuição Weibull
gamma <- 2
alpha <- 1.5

# Geração de tempos de falha e censura
Tfalha <- rweibull(n, shape = gamma, scale = alpha)
Tcensu <- rexp(n, rate = 1)

# Determinação do tempo observado e indicador de censura
Tobservado <- pmin(Tfalha, Tcensu)
indCensura <- as.numeric(Tfalha <= Tcensu)

# Dados simulados
dados <- data.frame(Tempo = Tobservado, Censura = indCensura)

# Histograma
ggplot(data = dados, aes(x = Tempo)) +
  geom_histogram(bins = 15) +
  labs(title = "Histograma do Tempo de Sobrevivência - Dados Simulados",
       x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------
# [3] Estimação
# -------------

# -------------------------------
# [3.1] Estimação Não Paramátrica
# -------------------------------
# ---------------------------------
# [3.1.1] Estimador de Kaplan-Meier
# ---------------------------------

# Modelo de Kaplan-Meier
ekm <- survfit(Surv(Tempo, Censura) ~ 1, data = dados)

# Formatando como DataFrame
DataEKM <- data.frame(Time = ekm$time, Survival = ekm$surv, Type = "Kaplan-Meier")

# Visualização da curva
ggplot(data = DataEKM, aes(x = Time, y = Survival)) +
  geom_line(color = "blue", lwd = 1.5) +
  labs(title = "S(t) - Kaplan-Meier",
       x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# ---------------------------
# [3.2] Estimação Paramátrica
# ---------------------------
# --------------------------------
# [3.2.1] Distribuição Exponencial
# --------------------------------

# Função de Sobrevivência
Stexp <- function(t, alpha) exp(-alpha * t)

# EMV de α
emvExp <- sum(dados$Censura) / sum(dados$Tempo)

# EMV da Sobrevivência
EMVSurvExp <- Stexp(dados$Tempo, emvExp)

# Formatando como DataFrame
DataExp <- data.frame(Time = dados$Tempo, Survival = EMVSurvExp, Type = "Exponencial")

# Visualização da curva
ggplot(data = DataExp, aes(x = Time, y = Survival)) +
  geom_line(color = "red", lwd = 1.5) +
  labs(title = "S(t) - Exponencial",
       x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ----------------------------
# [3.2.1] Distribuição Weibull
# ----------------------------

# Função de Sobrevivência
StWeibull <- function(t, gamma, alpha) exp(-(t / alpha)^gamma)

# EMV de γ e α
# 1. Função Log-verossimilhança
logWeibull <- function(theta, dados){
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  t <- dados$Tempo   # Tempo de falha
  c <- dados$Censura # Variável indicadora
  
  logv <- (sum(c) * log(gamma)) - (gamma * log(alpha) * sum(c)) +
    (gamma - 1) * sum(c * log(t)) - sum((t / alpha)^gamma)
  return(-logv)
}

# 2. Otimizando
Theta0 <- c(1.5, 1)
estimate <- optim(par = Theta0, fn = logWeibull,
                  gr = NULL, method = "BFGS", hessian = TRUE, dados = dados)

estimate

# EMV da Sobrevivência
EMVSurvWeib <- StWeibull(dados$Tempo, estimate$par[1], estimate$par[2])

# Formatando como DataFrame
DataWeib <- data.frame(Time = dados$Tempo, Survival = EMVSurvWeib, Type = "Weibull")

# Visualização da curva
ggplot(data = DataWeib, aes(x = Time, y = Survival)) +
  geom_line(color = "red", lwd = 1.5) +
  labs(title = "S(t) - Weibull",
       x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --------------------
# [4] Análise Conjunta
# --------------------

# Unindo os dados para visualização
AllData <- rbind(DataEKM, DataExp, DataWeib)

# Gráfico com ggplot2
ggplot(AllData, aes(x = Time, y = Survival)) +
  geom_line(aes(color = Type, linetype = Type), lwd = 1.5) +
  labs(title = "S(t): Kaplan-Meier, Weibull e Exponencial",
       x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank()) +
  scale_color_manual(values = c("blue", "green", "red")) +
  scale_linetype_manual(values = c(1, 2, 3)) # c("solid", "dashed", "dotted")
