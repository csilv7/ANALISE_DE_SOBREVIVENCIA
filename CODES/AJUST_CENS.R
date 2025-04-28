# -------------------
# ATIVAÇÃO DE PACOTES
# -------------------
library(ggplot2)
library(survival)
library(eha)

# ------------------
# [1] MODELO WEIBULL
# ------------------

# Função de sobrevivência
survival.weib <- function(times,shape.par,scale.par) 1 - pweibull(q=times, shape=shape.par, scale=scale.par)

# ---------------------------
# [1.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
shape.weib <- 2     # Parâmetro de forma
scale.weib <- 1.5   # Parâmetro de escala
rate.exp <- 1       # Parâmetro de taxa da exponencial (censura)

# Simulação dos dados
t <- rweibull(n, shape.weib, scale.weib) # Tempos de evento
c <- rexp(n, rate = rate.exp)            # Tempos censurados
times <- pmin(t, c)                      # Tempos observados
delta <- as.numeric(t <= c)              # Variável indicadora

# ----------------------------------------------
# [1.1.1] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.weib <- function(par, times, cens) {
  # Distição dos parâmetros
  gamma <- par[1] # Parâmetro de forma
  alpha <- par[2] # Parâmetro de escala
  
  # Função Log-verossimilhança
  ft <- dweibull(x = times, gamma, alpha)
  st <- 1 - pweibull(q = times, gamma, alpha)
  flv <- sum(cens * log(ft) + (1 - cens) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- c(1, 1)

# Otimização
ajust <- optim(
  par = init, fn = loglikelihood.weib, method = "BFGS", 
  hessian = TRUE, times = times, cens = delta
)

# Visualização
ajust

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Estimador de Kaplan-Meier
#ekm <- survfit(Surv(times, delta)~1)

# Organização dos dados
dados.weib <- data.frame(
  times = sort(times),
  st = survival.weib(sort(times), shape.weib, scale.weib),
  #st.ekm = ekm$surv,
  st.emv = survival.weib(sort(times), ajust$par[1], ajust$par[2])
)

# Plot da função de sobrevivência
ggplot(dados.weib, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  #geom_line(aes(y = st.ekm, color = "KM"), lwd = 1, lty = 2) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  #scale_color_manual(values = c("Verdadeiro" = "black", "KM" = "blue", "EMV" = "red")) +
  scale_color_manual(values = c("Verdadeiro" = "black", "EMV" = "red")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Sobrevida") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------
# [2] MODELO LOG-NORMAL
# ---------------------

# Função de sobrevivência
survival.lnorm <- function(times, loc.par, scale.par) 1 - pnorm(q=(log(times)-loc.par)/scale.par)

# ---------------------------
# [2.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
mu <- 1             # Parâmetro de locação
sigma <- 2          # Parâmetro de escala
rate.exp <- 1       # Parâmetro de taxa da exponencial (censura)

# Simulação dos dados
t <- rlnorm(n, meanlog = mu, sdlog = sigma) # Tempos de evento
c <- rexp(n, rate = rate.exp)               # Tempos censurados
times <- pmin(t, c)                         # Tempos observados
delta <- as.numeric(t <= c)                 # Variável indicadora

# ----------------------------------------------
# [2.1.1] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.lnorm <- function(par, times, cens) {
  # Distição dos parâmetros
  mu <- par[1] # Parâmetro de locação
  sigma <- par[2] # Parâmetro de escala
  
  # Função Log-verossimilhança
  ft <- dlnorm(x = times, mu, sigma)
  st <- 1 - pnorm(q=(log(times) - mu) / sigma)
  flv <- sum(cens * log(ft) + (1 - cens) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- c(0.5, 0.5)

# Otimização
ajust <- optim(
  par = init, fn = loglikelihood.lnorm, method = "BFGS", 
  hessian = TRUE, times = times, cens = delta
)

# Visualização
ajust

# ------------------------------------
# [2.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Estimador de Kaplan-Meier
#ekm <- survfit(Surv(times, delta)~1)

# Organização dos dados
dados.lnorm <- data.frame(
  times = sort(times),
  st = survival.lnorm(sort(times), mu, sigma),
  #st.ekm = ekm$surv,
  st.emv = survival.lnorm(sort(times), ajust$par[1], ajust$par[2])
)

# Plot da função de sobrevivência
ggplot(dados.lnorm, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  #geom_line(aes(y = st.ekm, color = "KM"), lwd = 1, lty = 2) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  #scale_color_manual(values = c("Verdadeiro" = "black", "KM" = "blue", "EMV" = "red")) +
  scale_color_manual(values = c("Verdadeiro" = "black", "EMV" = "red")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Sobrevida") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------------------
# [3] MODELO EXPONENCIAL POR PARTES
# ---------------------------------

# Função de sobrevivência
survival.pch <- function(times, cuts.points, levels.par) 1 - ppch(q=times, cuts=cuts.points, levels=levels.par)

# ---------------------------
# [3.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
rates <- c(0.5, 1, 1.5, 2) # Parâmetro de escala
breaks <- c(0.4, 1.2, 1.8) # Pontos de corte
rate.exp <- 1       # Parâmetro de taxa da exponencial (censura)

# Simulação dos dados
t <- rpch(n, cuts = breaks, levels = rates) # Tempos de evento
c <- rexp(n, rate = rate.exp)               # Tempos censurados
times <- pmin(t, c)                         # Tempos observados
delta <- as.numeric(t <= c)                 # Variável indicadora

# ----------------------------------------------
# [3.1.1] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.pch <- function(par, times, cens, cuts.points) {
  
  # Função Log-verossimilhança
  ft <- dpch(x = times, cuts = cuts.points, levels = par)
  st <- 1 - ppch(q = times, cuts = cuts.points, levels = par)
  flv <- sum(cens * log(ft) + (1 - cens) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- rep(1, length(rates))

# Otimização
ajust <- optim(par=init, fn = loglikelihood.pch, 
               gr = NULL, method = "BFGS", hessian = TRUE, 
               times=times, cens=delta, cuts.points=breaks)

# Visualização
ajust

# ------------------------------------
# [3.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Estimador de Kaplan-Meier
#ekm <- survfit(Surv(times, delta)~1)

# Organização dos dados
dados.pch <- data.frame(
  times = sort(times),
  st = survival.pch(sort(times), breaks, rates),
  #st.ekm = ekm$surv,
  st.emv = survival.pch(sort(times), breaks, ajust$par)
)

# Plot da função de sobrevivência
ggplot(dados.pch, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  #geom_line(aes(y = st.ekm, color = "KM"), lwd = 1, lty = 2) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  #scale_color_manual(values = c("Verdadeiro" = "black", "KM" = "blue", "EMV" = "red")) +
  scale_color_manual(values = c("Verdadeiro" = "black", "EMV" = "red")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Sobrevida") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------------------------------
# [4] MODELO EXPONENCIAL POR PARTES DE POTÊNCIA
# ---------------------------------------------

# Função de sobrevivência
survival.pchp <- function(times, cuts.points, levels.par, power.par) 1 - ppch(q=times, cuts=cuts.points, levels=levels.par)^power.par

# ---------------------------
# [4.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
rates <- c(0.5, 1, 1.5, 2) # Parâmetros de taxa
breaks <- c(0.4, 1.2, 1.8) # Pontos de corte
power <- 3/2               # Parâmetro de potência
rate.exp <- 1             # Parâmetro de taxa da exponencial (censura)

# Funções de Simulação
time <- function(t, cuts.points=cuts.points, rates.par=rates.par, power.par=power.par, u.unif) {
  surv <- 1 - ppch(q = t, cuts = cuts.points, levels = rates.par)^power.par
  return(surv - u.unif)
}

gen.pchp <- function(n, cuts.points, rates.par, power.par) {
  # Vetor para armazenar os tempos gerados
  pchp.times <- numeric(n)
  
  for (i in 1:n) {
    # Gerando um único valor de U para cada iteração
    u <- runif(1)
    
    # Encontrando a raiz para cada observação
    raiz <- uniroot(
      time, interval = c(0, 10000), cuts.points = cuts.points, 
      rates.par = rates.par, power.par = power.par, u.unif = u
    )
    
    pchp.times[i] <- raiz$root
  }
  
  return(pchp.times)
}

# Simulação dos dados
t <- gen.pchp(n, cuts.points=breaks, rates.par=rates, power.par=power) # Tempos de evento
c <- rexp(n, rate = rate.exp)               # Tempos censurados
times <- pmin(t, c)                         # Tempos observados
delta <- as.numeric(t <= c)                 # Variável indicadora

# ----------------------------------------------
# [4.1.1] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.pchp <- function(par, times, cens, cuts.points) {
  # Ajuste de Parâmetros
  n.par <- length(par)
  n.cuts <- length(cuts.points)
  rates.par <- par[1:(n.cuts + 1)]
  power.par <- par[n.par]
  
  # Ajutes de variáveis
  t <- times
  
  # Função Log-verossimilhança
  ft <- power.par*(ppch(q=t, cuts=cuts.points, levels=rates.par))^(power.par-1)*dpch(x=t, cuts=cuts.points, levels=rates.par)
  st <- 1 - ppch(q = t, cuts = cuts.points, levels = rates.par)^power.par
  flv <- sum(cens * log(ft) + (1 - cens) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- rep(1, length(rates) + 1)

# Otimização
ajust <- optim(par=init, fn = loglikelihood.pchp, 
               gr = NULL, method = "BFGS", hessian = TRUE, 
               times=times, cens=delta, cuts.points=breaks)

# Visualização
ajust

# ------------------------------------
# [4.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Estimador de Kaplan-Meier
#ekm <- survfit(Surv(times, delta)~1)

# Organização dos dados
dados.pchp <- data.frame(
  times = sort(times),
  st = survival.pchp(sort(times), breaks, rates, power),
  #st.ekm = ekm$surv,
  st.emv = survival.pchp(sort(times), breaks, ajust$par[1:(length(ajust$par)-1)], ajust$par[length(ajust$par)])
)

# Plot da função de sobrevivência
ggplot(dados.pchp, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  #geom_line(aes(y = st.ekm, color = "KM"), lwd = 1, lty = 2) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  #scale_color_manual(values = c("Verdadeiro" = "black", "KM" = "blue", "EMV" = "red")) +
  scale_color_manual(values = c("Verdadeiro" = "black", "EMV" = "red")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Sobrevida") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")