library(ggplot2)
library(survival)
library(eha)

# ----------------------
# [1] MODELO EXPONENCIAL
# ----------------------

# Função de sobrevivência
survival.exp <- function(times, rate.par) 1 - pexp(q = times, rate = rate.par)

# ---------------------------
# [1.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
rate <- 1.50        # Parâmetro de Taxa

# Simulação dos dados
times <- rexp(n, rate = rate)

# Estimador de Máxima Verossimilhança (EMV)
emv.exp <- n / sum(times)

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Organização dos dados
dados.exp <- data.frame(
  times = sort(times),
  st = survival.exp(sort(times), rate),
  st.emv = survival.exp(sort(times), emv.exp)
)

# Plot da função de sobrevivência
ggplot(dados.exp, aes(x = times)) +
  geom_line(aes(y = st, color = "rate"), lwd = 1) +
  geom_line(aes(y = st.emv, color = "emv.exp"), lwd = 1, lty = 4) +
  scale_color_manual(
    values = c("rate" = "black", "emv.exp" = "red"),
    labels = c(bquote(alpha == .(rate)), bquote(hat(alpha) == .(emv.exp)))
  ) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Parâmetro") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ------------------
# [2] MODELO WEIBULL
# ------------------

# Função de sobrevivência
survival.weib <- function(times,shape.par,scale.par) 1 - pweibull(q=times, shape=shape.par, scale=scale.par)

# ---------------------------
# [2.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
shape.weib <- 2     # Parâmetro de forma
scale.weib <- 1.5   # Parâmetro de escala

# Simulação dos dados
times <- rweibull(n, shape.weib, scale.weib)

# ----------------------------------------------
# [2.1.1] IMPLEMENTAÇÃO SEM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Vetor gradiente
GRADIEN <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Ajutes de variáveis
  t <- times
  
  # Derivadas Parciais
  DerivGamma <- n/gamma - log(alpha)*n + sum(log(t)) - sum(((t/alpha)^gamma)*log(t/alpha))
  DerivAlpha <- -(gamma/alpha)*n + gamma*(alpha^(-gamma-1))*sum(t^gamma)
  
  # Vetor Gradiente
  gradient <- c(DerivGamma, DerivAlpha)
  
  # Retornar
  return(gradient)
}

# Matriz Hessiana
HESSIAN <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Ajutes de variáveis
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

# Método Iterativo de Newton-Raphson
init <- c(1, 1)  # Chute Inicial
diff <- 1        # Diferença entre o passo atual e o passo anterior
error <- 10^(-8) # Erro tolerável
id <- 1          # Contador da iteração

# Iteração
while (diff > error) {
  # Vetor Gradiente e Matriz Hessiana
  U <- GRADIEN(times = times, theta = init) # Vetor Gradiente
  H <- HESSIAN(times = times, theta = init) # Matriz Hessiana
  
  # Solução do sistema linear H %*% solution = U
  solution <- solve(H, U)
  
  # Atualização do Algoritmo
  ajust <- init - solution
  
  # Diferença entre os parâmetros
  diff <- max(abs(ajust - init))
  
  # Imprimir resultados na tela
  #cat("Iteração:", id, " -  Estimativa = (Forma:", ajust[1], ", Escala:", ajust[2], ") \n")
  
  # Controle do Algoritmo
  init <- ajust
  id <- id + 1
}

# Impressão de resultados
cat("Número de Iterações Necessárias:", id, "\n")
cat("Estimativa para o parâmetro de forma:", ajust[1], "\n")
cat("Estimativa para o parâmetro de forma:", ajust[2], "\n")

# ----------------------------------------------
# [2.1.2] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.weib <- function(par, times) {
  # Distição dos parâmetros
  gamma <- par[1] # Parâmetro de forma
  alpha <- par[2] # Parâmetro de escala
  
  # Ajutes de variáveis
  t <- times
  c <- rep(1, length(t))
  
  # Função Log-verossimilhança
  ft <- dweibull(x = t, gamma, alpha)
  st <- 1 - pweibull(q = t, gamma, alpha)
  flv <- sum(c * log(ft) + (1 - c) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- c(1, 1)

# Otimização
ajust <- optim(
  par = init, fn = loglikelihood.weib,
  method = "BFGS", hessian = TRUE, times = times
)

# Visualização
ajust

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Organização dos dados
dados.weib <- data.frame(
  times = sort(times),
  st = survival.weib(sort(times), shape.weib, scale.weib),
  st.emv = survival.weib(sort(times), ajust$par[1], ajust$par[2])
)

# Plot da função de sobrevivência
ggplot(dados.weib, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  scale_color_manual(
    values = c("Verdadeiro" = "black", "EMV" = "red"),
    labels = c(
      bquote("Verdadeiro: " ~ gamma == .(shape.weib) ~ ";" ~ alpha == .(scale.weib)),
      bquote("EMV: " ~ hat(gamma) == .(ajust$par[1]) ~ ";" ~ hat(alpha) == .(ajust$par[2]))
    )
  ) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Parâmetro") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------
# [3] MODELO LOG-NORMAL
# ---------------------

# Função de sobrevivência
survival.lnorm <- function(times, loc.par, scale.par) 1 - pnorm(q=(log(times)-loc.par)/scale.par)

# ---------------------------
# [3.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
mu <- 0             # Parâmetro de locação
sigma <- 1          # Parâmetro de escala

# Simulação dos dados
times <- rlnorm(n, meanlog = mu, sdlog = sigma)

# ----------------------------------------------
# [2.1.2] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.lnorm <- function(par, times) {
  # Distição dos parâmetros
  mu <- par[1] # Parâmetro de locação
  sigma <- par[2] # Parâmetro de escala
  
  # Ajutes de variáveis
  t <- times
  c <- rep(1, length(t))
  
  # Função Log-verossimilhança
  ft <- dlnorm(x = t, mu, sigma)
  st <- 1 - pnorm(q=(log(t) - mu) / sigma)
  flv <- sum(c * log(ft) + (1 - c) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- c(1, 0.5)

# Otimização
ajust <- optim(
  par = init, fn = loglikelihood.lnorm,
  method = "BFGS", hessian = TRUE, times = times
)

# Visualização
ajust

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Organização dos dados
dados.lnorm <- data.frame(
  times = sort(times),
  st = survival.lnorm(sort(times), mu, sigma),
  st.emv = survival.lnorm(sort(times), ajust$par[1], ajust$par[2])
)

# Plot da função de sobrevivência
ggplot(dados.lnorm, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  scale_color_manual(
    values = c("Verdadeiro" = "black", "EMV" = "red"),
    labels = c(
      bquote("Verdadeiro: " ~ mu == .(mu) ~ ";" ~ sigma == .(sigma)),
      bquote("EMV: " ~ hat(mu) == .(ajust$par[1]) ~ ";" ~ hat(sigma) == .(ajust$par[2]))
    )
  ) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Parâmetro") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------------------
# [4] MODELO EXPONENCIAL POR PARTES
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

# Simulação dos dados
times <- rpch(n, cuts = breaks, levels = rates)

# ----------------------------------------------
# [2.1.2] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.pch <- function(par, times, cuts.points) {
  # Ajutes de variáveis
  t <- times
  c <- rep(1, length(t))
  
  # Função Log-verossimilhança
  ft <- dpch(x = times, cuts = cuts.points, levels = par)
  st <- 1 - ppch(q = times, cuts = cuts.points, levels = par)
  flv <- sum(c * log(ft) + (1 - c) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- rep(1, length(rates))

# Otimização
ajust <- optim(par=init, fn = loglikelihood.pch, 
               gr = NULL, method = "BFGS", hessian = TRUE, 
               times=times, cuts.points=breaks)

# Visualização
ajust

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Organização dos dados
dados.pch <- data.frame(
  times = sort(times),
  st = survival.pch(sort(times), breaks, rates),
  st.emv = survival.pch(sort(times), breaks, ajust$par)
)

# Plot da função de sobrevivência
ggplot(dados.pch, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  scale_color_manual(
    values = c("Verdadeiro" = "black", "EMV" = "red"),
  ) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Parâmetro") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---------------------------------------------
# [5] MODELO EXPONENCIAL POR PARTES DE POTÊNCIA
# ---------------------------------------------

# Função de sobrevivência
survival.pchp <- function(times, cuts.points, levels.par, power.par) 1 - ppch(q=times, cuts=cuts.points, levels=levels.par)^power.par

# ---------------------------
# [3.1] SIMULAÇÃO E ESTIMAÇÃO
# ---------------------------

# Parâmetros de simulação
set.seed(123456789) # Semente aleatória
n <- 1000           # Tamanho amostral
rates <- c(0.5, 1, 1.5, 2) # Parâmetros de taxa
breaks <- c(0.4, 1.2, 1.8) # Pontos de corte
power <- 3/2

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
times <- gen.pchp(n, cuts.points = breaks, rates.par = rates, power.par = power)

# ----------------------------------------------
# [2.1.2] IMPLEMENTAÇÃO COM FUNÇÃO DE OTIMIZAÇÃO
# ----------------------------------------------

# Definir a função log-verossimilhança
loglikelihood.pchp <- function(par, times, cuts.points) {
  # Ajuste de Parâmetros
  n.par <- length(par)
  n.cuts <- length(cuts.points)
  rates.par <- par[1:(n.cuts + 1)]
  power.par <- par[n.par]
  
  # Ajutes de variáveis
  t <- times
  c <- rep(1, length(t))
  
  # Função Log-verossimilhança
  ft <- power.par*(ppch(q=t, cuts=cuts.points, levels=rates.par))^(power.par-1)*dpch(x=t, cuts=cuts.points, levels=rates.par)
  st <- 1 - ppch(q = t, cuts = cuts.points, levels = rates.par)^power.par
  flv <- sum(c * log(ft) + (1 - c) * log(st))
  
  # Retorna o valor simétrico
  return(-flv)
}

# Chute Inicial
init <- rep(1, length(rates) + 1)

# Otimização
ajust <- optim(par=init, fn = loglikelihood.pchp, 
               gr = NULL, method = "BFGS", hessian = TRUE, 
               times=times, cuts.points=breaks)

# Visualização
ajust

# ------------------------------------
# [1.2] VISUALIZAÇÃO GRÁFICA DO AJUSTE
# ------------------------------------

# Organização dos dados
dados.pchp <- data.frame(
  times = sort(times),
  st = survival.pchp(sort(times), breaks, rates, power),
  st.emv = survival.pchp(sort(times), breaks, ajust$par[1:(length(ajust$par)-1)], ajust$par[length(ajust$par)])
)

# Plot da função de sobrevivência
ggplot(dados.pchp, aes(x = times)) +
  geom_line(aes(y = st, color = "Verdadeiro"), lwd = 1) +
  geom_line(aes(y = st.emv, color = "EMV"), lwd = 1, lty = 4) +
  scale_color_manual(
    values = c("Verdadeiro" = "black", "EMV" = "red"),
  ) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Parâmetro") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")