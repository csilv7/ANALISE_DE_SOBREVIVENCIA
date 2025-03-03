# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
library(dplyr)
library(ggplot2)

set.seed(123)
n <- 1000 # Tamanho da amostra simulada

# ---------------------------
# [1] DISTRIBUIÇÃO EXPONENCIAL
# ---------------------------
# -------------
# [1.1] FUNÇÕES
# -------------
ftexp <- function(t, alpha) alpha * exp(-alpha * t)
Stexp <- function(t, alpha) exp(-alpha * t)
htexp <- function(t, alpha) rep(alpha, length(t))
Ltexp <- function(t, alpha) alpha * t

# ----------------------------------------
# [1.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
n <- 1000                  # Tamanho amostral
tempo <- rexp(n, rate = 1) # Simulando dados de uma exponencial
alphas <- c(1, 1.5, 2)     # Valores do parâmetro a serem avaliados

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(alphas, function(alpha) {
  data.frame(
    tempo = sort(tempo),
    ft = ftexp(sort(tempo), alpha),
    St = Stexp(sort(tempo), alpha),
    ht = htexp(sort(tempo), alpha),
    Lt = Ltexp(sort(tempo), alpha),
    alpha = factor(alpha)
  )
}))

# --------------------
# [1.3] FUNÇÃO GRÁFICA
# --------------------
PlotFunction <- function(dados, ft, label) {
  ggplot(data = dados, aes_string(x = "tempo", y = ft, color = "alpha")) +
    geom_line(size = 1.2) +
    labs(
      x = "Tempo",
      y = label, 
      color = expression(alpha)
    ) +
    scale_color_manual(
      values = c("red", "blue", "green"),
      labels = scales::parse_format()(levels(dados$alpha))
    ) +
    theme_minimal(base_size = 12)
}

# Plotando a função densidade de probabilidade
PlotFunction(dados, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
PlotFunction(dados, "St", "Função de Sobrevivência")

# Plotando a função de risco
PlotFunction(dados, "ht", "Função de Risco")

# Plotando a função de risco acumulado
PlotFunction(dados, "Lt", "Função de Risco Acumulado")

# ------------------------
# [2] DISTRIBUIÇÃO WEIBULL
# ------------------------
# -------------
# [2.1] FUNÇÕES
# -------------
ftexp <- function(t, gamma, alpha) {
  ft <- gamma*alpha^(-gamma)*t^(gamma-1)*exp(-(t/alpha)^gamma)
  return(ft)
}
Stexp <- function(t, gamma, alpha) {
  St <- exp(-(t/alpha)^gamma)
  return(St)
}
htexp <- function(t, gamma, alpha) {
  ht <- gamma*alpha^(-gamma)*t^(gamma-1)
  return(ht)
}
Ltexp <- function(t, gamma, alpha) {
  Lt <- (t/alpha)^gamma
  return(Lt)
}

# ----------------------------------------
# [2.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
n <- 1000                                  # Tamanho amostral
tempo <- rweibull(n, shape = 2, scale = 1) # Simulando dados de uma Weibull
alpha <- 1                                 # Fixo para simplificar
gammas <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # # Valores do parâmetro a serem avaliados

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(gammas, function(gamma) {
  data.frame(
    tempo = sort(tempo),
    ft = ftexp(sort(tempo), gamma, alpha),
    St = Stexp(sort(tempo), gamma, alpha),
    ht = htexp(sort(tempo), gamma, alpha),
    Lt = Ltexp(sort(tempo), gamma, alpha),
    gamma = factor(gamma)
  )
}))

# --------------------
# [2.3] FUNÇÃO GRÁFICA
# --------------------
PlotFunction <- function(dados, ft, label) {
  ggplot(data = dados, aes_string(x = "tempo", y = ft, color = "gamma")) +
    geom_line(size = 1.2) +
    labs(
      x = "Tempo",
      y = label, 
      color = expression(gamma)
    ) +
    scale_color_manual(
      values = c("red", "blue", "green", "purple", "orange", "brown"),
      labels = scales::parse_format()(levels(dados$gamma))
    ) +
    theme_minimal(base_size = 12)
}

# Plotando a função densidade de probabilidade
PlotFunction(dados, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
PlotFunction(dados, "St", "Função de Sobrevivência")

# Plotando a função de risco
PlotFunction(dados, "ht", "Função de Risco")

# Plotando a função de risco acumulado
PlotFunction(dados, "Lt", "Função de Risco Acumulado")

# ---------------------------
# [3] DISTRIBUIÇÃO LOG-NORMAL
# ---------------------------
# -------------
# [3.1] FUNÇÕES
# -------------
ftlnorm <- function(t, mu, sigma) {
  (1 / (t * sigma * sqrt(2 * pi))) * exp(-0.5 * ((log(t) - mu) / sigma)^2)
}

Stlnorm <- function(t, mu, sigma) {
  1 - pnorm((log(t) - mu) / sigma, lower.tail = TRUE)
}

htlnorm <- function(t, mu, sigma) {
  ftlnorm(t, mu, sigma) / Stlnorm(t, mu, sigma)
}

Ltlnorm <- function(t, mu, sigma) {
  -log(Stlnorm(t, mu, sigma))
}

# ----------------------------------------
# [3.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
n <- 1000                                  # Tamanho amostral
tempo <- rlnorm(n, meanlog = 0, sdlog = 1) # Simulando dados de uma Log-normal
mus <- c(0, 0.5, 1, 1.5, 2)                # Valores de mu
sigma <- 1                                 # Valor fixo de sigma

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(mus, function(mu) {
  data.frame(
    tempo = sort(tempo),
    ft = ftlnorm(sort(tempo), mu, sigma),
    St = Stlnorm(sort(tempo), mu, sigma),
    ht = htlnorm(sort(tempo), mu, sigma),
    Lt = Ltlnorm(sort(tempo), mu, sigma),
    mu = factor(mu)
  )
}))

# --------------------
# [3.3] FUNÇÃO GRÁFICA
# --------------------
PlotFunction <- function(dados, ft, label) {
  ggplot(data = dados, aes_string(x = "tempo", y = ft, color = "mu")) +
    geom_line(size = 1.2) +
    labs(
      x = "Tempo",
      y = label, 
      color = expression(mu)
    ) +
    scale_color_manual(
      values = c("red", "blue", "green", "purple", "orange"),
      labels = scales::parse_format()(levels(dados$mu))
    ) +
    theme_minimal(base_size = 12)
}

# Plotando a função densidade de probabilidade
PlotFunction(dados, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
PlotFunction(dados, "St", "Função de Sobrevivência")

# Plotando a função de risco
PlotFunction(dados, "ht", "Função de Risco")

# Plotando a função de risco acumulado
PlotFunction(dados, "Lt", "Função de Risco Acumulado")


################################################################################
################################################################################
################################################################################
################################################################################

# --------------------------
# [1] Configurações Iniciais
# --------------------------

# Semente
set.seed(123456789)

# Tamanho Amostral
n <- 1000

# ------------------------
# [2] Distribuição Weibull
# ------------------------
# -------------------------
# [2.1] Simulação dos Dados
# -------------------------

w_shape <- 2   # Parâmetro de Forma
w_scale <- 1.5 # Parâmetro de Escala

# Simulação
w_dados <- rweibull(n, shape = w_shape, scale = w_scale)

# -------------
# [3] Estimação
# -------------
# ----------------------------------------
# [3.1] Método Iterativo de Newton-Raphson
# ----------------------------------------
# ----------------------
# [3.1.1] Função `optim`
# ----------------------

# Definindo Função Log-verossimilhança
logVerossimil <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Função Log-verossimilhança
  flv <- log(gamma)*n - gamma*log(alpha)*n + (gamma - 1)*sum(log(t)) - sum((t/alpha)^gamma)
  
  # Retorna o valor simétrico
  return(-flv)
}

# Vetor de Parâmetros: Chute inicial
theta0 <- c(1, 1)

# Obtendo as estimativas
estimate <- optim(
  par = theta0,
  fn = logVerossimil,
  method = "BFGS",
  hessian = TRUE,
  times = w_dados
)

# Resultados
print(estimate)

# -----------------
# [3.1.1] Algoritmo
# -----------------
# -------------------------
# [3.1.1.1] Vetor Gradiente
# -------------------------
GRADIEN <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Derivadas Parciais
  DerivGamma <- n/gamma - log(alpha)*n + sum(log(t)) - sum(((t/alpha)^gamma)*log(t/alpha))
  DerivAlpha <- -(gamma/alpha)*n + gamma*(alpha^(-gamma-1))*sum(t^gamma)
  
  # Vetor Gradiente
  gradient <- c(DerivGamma, DerivAlpha)
  
  # Retornar
  return(gradient)
}

# -----------------------
# [2.2.2] Matriz Hessiana
# -----------------------
HESSIAN <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Derivadas de 2ª ordem
  D2Gamma <- - n/gamma^2 - sum(((t/alpha)^gamma)*(log(t/alpha)^2))
  D2Alpha <- - (gamma/alpha^2)*n - gamma*(gamma + 1)*(alpha^(-gamma-2))*sum(t^gamma)
  D2 <- -n*alpha + (alpha^(-gamma-1))*sum((t^gamma)*(gamma*log(t/alpha) + 1))
  
  # Matriz Hessiana
  H <- matrix(
    data = c(D2Gamma, D2, D2, D2Alpha),
    nrow = 2, ncol = 2
  )
  
  # Retornar
  return(H)
}

# --------------------------------------
# [3] Método Iterativo de Newton-Raphson
# --------------------------------------
theta0 <- c(1, 1)     # Chute Inicial
diff <- 1             # Diferença entre o passo atual e o passo anterior
error <- 10^(-8)      # Erro tolerável
id <- 1               # Contador da iteração

# Iteração
while (diff > error) {
  # Vetor Gradiente e Matriz Hessiana
  U <- GRADIEN(times = w_dados, theta = theta0) # Vetor Escore
  H <- HESSIAN(times = w_dados, theta = theta0) # Matriz Hessiana
  
  # Solução do sistema linear H %*% solution = U
  solution <- solve(H, U)
  
  # Atualização do Algoritmo
  theta1 <- theta0 - solution
  
  # Diferença entre os parâmetros
  diff <- max(abs(theta1 - theta0))
  
  # Imprimir resultados na tela
  cat("Iteração:", id, " -  Estimativa = (Forma:", theta1[1], ", Escala:", theta1[2], ") \n")
  
  # Controle do Algoritmo
  theta0 <- theta1
  id <- id + 1
}


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


