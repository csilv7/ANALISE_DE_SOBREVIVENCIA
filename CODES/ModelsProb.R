# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
library(ggplot2)

set.seed(123)
n <- 1000 # Tamanho da amostra simulada

# ---------------------------
# [2] DISTRIBUIÇÃO EXPONENCIAL
# ---------------------------

# -------------
# [2.1] FUNÇÕES
# -------------
# As funções de sobrevivência, risco e risco acumulado são simplificadas
Stexp <- function(t, alpha) exp(-alpha * t)
htexp <- function(alpha) rep(alpha, length(t))
Ltexp <- function(t, alpha) alpha * t

# ----------------------------------------
# [2.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
tempo <- rexp(n, rate = 1) # Simulando dados de uma exponencial
alphas <- c(1, 1.5, 2)     # Valores de alpha a serem avaliados

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(alphas, function(alpha) {
  data.frame(
    tempo = tempo,
    St = Stexp(tempo, alpha),
    ht = htexp(alpha),
    Lt = Ltexp(tempo, alpha),
    alpha = factor(alpha)
  )
}))

# -------------
# [2.3] GRÁFICOS
# -------------

# Criando uma função para gerar gráficos
plot_func <- function(data, y_var, y_label, color_values, y_expression) {
  ggplot(data, aes(x = tempo, y = !!sym(y_var), color = alpha)) +
    geom_line(stat = "summary", fun = mean, size = 1) +
    labs(x = "Tempo", y = y_expression, color = expression(alpha)) +
    scale_color_manual(values = color_values,
                       labels = lapply(alphas, function(a) bquote(alpha == .(a)))) +
    theme_minimal()
}

# Função de Sobrevivência
plot_func(dados, "St", "S(t)", c("red", "blue", "green"), expression(S(t)))

# Função de Risco
plot_func(dados, "ht", expression(lambda(t)), c("red", "blue", "green"), expression(lambda(t)))

# Função de Risco Acumulado
plot_func(dados, "Lt", expression(Lambda(t)), c("red", "blue", "green"), expression(Lambda(t)))

# ------------------------
# [3] DISTRIBUIÇÃO WEIBULL
# ------------------------

# -------------
# [3.1] FUNÇÕES
# -------------

# Funções para Weibull
StWei <- function(t, alpha, gamma) exp(-(alpha * t)^gamma)
htWei <- function(t, alpha, gamma) gamma * (alpha^gamma) * t^(gamma - 1)
LtWei <- function(t, alpha, gamma) (alpha * t)^gamma

# ----------------------------------------
# [3.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------

# Simulando dados de uma Weibull
tempo <- rweibull(n, shape = 2, scale = 1)
alpha <- 1   # Fixo para simplificar
gammas <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0) # Valores de gamma

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(gammas, function(gamma) {
  data.frame(
    tempo = tempo,
    St = StWei(tempo, alpha, gamma),
    ht = htWei(tempo, alpha, gamma),
    Lt = LtWei(tempo, alpha, gamma),
    gamma = factor(gamma)
  )
}))

# -------------
# [3.3] GRÁFICOS
# -------------

# Função genérica para gráficos
plot_func <- function(data, y_var, y_label, color_values, y_expression) {
  ggplot(data, aes(x = tempo, y = !!sym(y_var), color = gamma)) +
    geom_line(stat = "summary", fun = mean, size = 1) +
    labs(x = "Tempo", y = y_expression, color = expression(gamma)) +
    scale_color_manual(values = color_values,
                       labels = lapply(gammas, function(g) bquote(gamma == .(g)))) +
    theme_minimal()
}

# Paleta de cores
color_values <- c("red", "blue", "green", "purple", "orange", "brown")

# Função de Sobrevivência
plot_func(dados, "St", expression(S(t)), color_values, expression(S(t)))

# Função de Risco
plot_func(dados, "ht", expression(lambda(t)), color_values, expression(lambda(t)))

# Função de Risco Acumulado
plot_func(dados, "Lt", expression(Lambda(t)), color_values, expression(Lambda(t)))

# --------------------------
# [4] DISTRIBUIÇÃO LOG-NORMAL
# --------------------------

library(ggplot2)
library(dplyr)

set.seed(123)
n <- 1000 # Tamanho da amostra simulada

# -------------
# [4.1] FUNÇÕES
# -------------

# Função densidade (f)
ftLogNormal <- function(t, mu, sigma) {
  (1 / (t * sigma * sqrt(2 * pi))) * exp(-0.5 * ((log(t) - mu) / sigma)^2)
}

# Função de Sobrevivência (S)
StLogNormal <- function(t, mu, sigma) {
  pnorm(-(log(t) - mu) / sigma, lower.tail = TRUE)
}

# Função de Risco (h)
htLogNormal <- function(t, mu, sigma) {
  ftLogNormal(t, mu, sigma) / StLogNormal(t, mu, sigma)
}

# Função de Risco Acumulado (Lambda)
LtLogNormal <- function(t, mu, sigma) {
  -log(StLogNormal(t, mu, sigma))
}

# ----------------------------------------
# [4.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------

# Simulando dados da distribuição log-normal
tempo <- rlnorm(n, meanlog = 0, sdlog = 1)
mus <- c(0, 0.5, 1) # Valores de mu
sigma <- 1         # Valor fixo de sigma

# Criando um Data Frame com valores das funções
dados <- do.call(rbind, lapply(mus, function(mu) {
  data.frame(
    tempo = tempo,
    ft = ftLogNormal(tempo, mu, sigma),
    St = StLogNormal(tempo, mu, sigma),
    ht = htLogNormal(tempo, mu, sigma),
    Lt = LtLogNormal(tempo, mu, sigma),
    mu = factor(mu)
  )
}))

# -------------
# [4.3] GRÁFICOS
# -------------

# Função genérica para gráficos
plot_func <- function(data, y_var, y_label, color_values, y_expression) {
  ggplot(data, aes(x = tempo, y = !!sym(y_var), color = mu)) +
    geom_line(stat = "summary", fun = mean, size = 1) +
    labs(x = "Tempo", y = y_expression, color = expression(mu)) +
    scale_color_manual(values = color_values,
                       labels = lapply(mus, function(m) bquote(mu == .(m)))) +
    theme_minimal()
}

# Paleta de cores
color_values <- c("red", "blue", "green")

# Função de Densidade
plot_func(dados, "ft", expression(f(t)), color_values, expression(f(t)))

# Função de Sobrevivência
plot_func(dados, "St", expression(S(t)), color_values, expression(S(t)))

# Função de Risco
plot_func(dados, "ht", expression(lambda(t)), color_values, expression(lambda(t)))

# Função de Risco Acumulado
plot_func(dados, "Lt", expression(Lambda(t)), color_values, expression(Lambda(t)))