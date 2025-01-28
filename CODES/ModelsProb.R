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