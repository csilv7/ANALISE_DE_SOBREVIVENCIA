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

library(survival)
library(ggplot2)

# -----------------------
# [2] Simulação dos Dados
# -----------------------

# Definindo Semente
set.seed(123)

# Parâmetros
n <- 1000       # Número total de observações
gamma <- 2      # Parâmetro shape da Weibull
alpha <- 1.5    # Parâmetro scale da Weibull
TaxaExp <- 1    # Taxa da distribuição exponencial
propCens <- 1/4 # Proporção desejada de censuras

# Vetores para armazenar os resultados
Tobservado <- numeric(n)
indCensura <- numeric(n)

# Contadores
nFalhas <- 0
nCensuras <- 0

# Loop para gerar os tempos
for (i in 1:n) {
  # Gerar um tempo de falha e um tempo de censura
  Tfalha <- rweibull(1, shape = gamma, scale = alpha)
  Tcensu <- rexp(1, rate = TaxaExp)
  
  # Verificar qual é o menor tempo
  if (Tfalha <= Tcensu) {
    Tobservado[i] <- Tfalha
    indCensura[i] <- 1 # Falha
    nFalhas <- nFalhas + 1
  } else {
    if (nCensuras < propCens * n) {
      Tobservado[i] <- Tcensu
      indCensura[i] <- 0 # Censura
      nCensuras <- nCensuras + 1
    } else {
      Tobservado[i] <- Tfalha
      indCensura[i] <- 1 # Falha
      nFalhas <- nFalhas + 1
    }
  }
}

# Verificar as proporções
cat("Proporção de falhas:", nFalhas / n, "\n")
cat("Proporção de censuras:", nCensuras / n, "\n")

# Dados simulados
dados <- data.frame(Tempo = Tobservado, Censura = indCensura)

# Histograma
ggplot(data = dados, aes(x = Tempo)) +
  geom_histogram(bins = 15, fill = "blue") +
  #xlim(c(0, max(dados$Tempo))) +
  labs(x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14)

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
DataEKM <- data.frame(Time = ekm$time, Survival = ekm$surv, Estimação = "Kaplan-Meier")


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
DataExp <- data.frame(Time = dados$Tempo, Survival = EMVSurvExp, 
                      Estimação = "Distrib Exponencial")

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
theta0 <- c(1.5, 1)
estimate <- optim(par = theta0, fn = logWeibull,
                  gr = NULL, method = "BFGS", 
                  hessian = TRUE, dados = dados)

# EMV da Sobrevivência
EMVSurvWeib <- StWeibull(dados$Tempo, estimate$par[1], estimate$par[2])

# Formatando como DataFrame
DataWeib <- data.frame(Time = dados$Tempo, Survival = EMVSurvWeib, 
                       Estimação = "Distrib Weibull")

# --------------------
# [4] Análise Conjunta
# --------------------

# Unindo os dados para visualização
AllData <- rbind(DataEKM, DataExp, DataWeib)

# Gráfico com ggplot2
ggplot(AllData, aes(x = Time, y = Survival)) +
  geom_line(aes(color = Estimação, lwd = Estimação)) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("blue", "green", "red")) +
  scale_linewidth_manual(values = c(1.5, 1.5, 1))


