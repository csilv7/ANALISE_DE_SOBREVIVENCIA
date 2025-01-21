# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
# Análise de Sobrevivência - Modelos de Tempo de Vida Acele
# Analista: Breno C R da Silva
# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------

# -----------------------
# [1] Ativação de Pacotes
# -----------------------

library(ggplot2)

set.seed(123456789)

# ----------------------------
# [2] Distribuição Exponencial
# ----------------------------

myRexp <- function(n, taxa, coVariates, vecCoef) {
  # V.A. com Distribuição Uniforme(0, 1)
  U <- runif(n, 0, 1)
  
  # Combinação linear
  efeito <- rowSums(coVariates * vecCoef)
  
  # Tempo de Sobrevivência para TVA
  times <- - exp(efeito) * log(1 - U) / taxa
  return(times)
}

# ------------------------------
# [2.2] Simulação de Monte Carlo
# ------------------------------
# ---------------------------
# [2.2.1] Função de Simulação
# ---------------------------

myRexp <- function(n, rate, coVariates, vecCoef) {
  U <- runif(n, 0, 1) # Uniforme(0, 1)
  effect <- rowSums(coVariates * vecCoef) # Combinação linear
  times <- -exp(effect) * log(1 - U) / rate
  return(times)
}

# --------------------------
# [2.2.2] Simulação de Dados
# --------------------------

set.seed(123456789)
n <- 100
x1 <- rnorm(n, 0, 1)
x2 <- rbinom(n, 1, 0.5)
X <- cbind(1, x1, x2)
vecCoef <- c(1.5, 2/3, 2)
timesSurv <- myRexp(n, 1.5, X, vecCoef)

# Visualização
dados <- data.frame(Tempo = timesSurv)

ggplot(data = dados, aes(x = Tempo)) +
  geom_histogram(fill = "blue") +
  labs(title = "Histograma do Tempo de Sobrevivência - Simulação",
       x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---------------------------
# [3] Estimação de Parâmetros
# ---------------------------
# -------------------------
# [3.1] Log-Verossimilhança
# -------------------------

logVerossimil <- function(theta, coVariates, times) {
  n <- length(times)
  rate <- theta[1]
  beta <- theta[-1]
  
  effect <- rowSums(coVariates * beta)
  flv <- n * log(rate) - rate * sum(times / exp(effect))
  
  return(-flv) # Minimizar
}

# ---------------------------------
# [3.2] Algoritmo de Newton-Raphson
# ---------------------------------

theta0 <- c(1, 1, 0.5, 0.5) # Chute inicial
estimate <- optim(
  par = theta0,
  fn = logVerossimil,
  coVariates = X,
  times = timesSurv,
  method = "BFGS",
  hessian = TRUE
)

# Resultados
print(estimate)

# --------------
# [4] Comparação
# --------------

# Função de Sobrevivência
Stexp <- function(t, alpha, betas, X) exp(-alpha * (t/exp(rowSums(X * betas))))






# Formatando como DataFrame
DataExp <- data.frame(Time = timesSurv, Survival = Stexp(timesSurv, 1.5, vecCoef, X),
                      EMV_Survival = Stexp(timesSurv, estimate$par[1], estimate$par[-1], X),
                      Type = "Exponencial")

# Gráfico com ggplot2
ggplot(DataExp, aes(x = Time)) +
  geom_line(aes(y = EMV_Survival)) +
  geom_line(aes(y = Survival), lwd = 1.5) +
  labs(title = "S(t) Exponencial",
       x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_minimal(base_size = 14)
  
