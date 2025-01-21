# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Análise de Sobrevivência - Modelos de Tempo de Vida Acelerado
# Analista: Breno C R da Silva
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------

# -----------------------
# [1] Ativação de Pacotes
# -----------------------

library(ggplot2)

set.seed(123456789)

# ----------------------------
# [2] Distribuição Exponencial
# ----------------------------
# ------------------------------
# [2.2] Simulação de Monte Carlo
# ------------------------------
# ---------------------------
# [2.2.1] Função de Simulação
# ---------------------------

myRexp <- function(n, rate, coVariates, vecCoef) {
  # Uniforme(0, 1)
  U <- runif(n, 0, 1)
  
  # Combinação linear dos preditores lineares
  effect <- rowSums(coVariates * vecCoef)
  
  # Tempos de sobrevivência
  times <- - exp(effect) * log(1 - U) / rate
  
  return(times)
}

# --------------------------
# [2.2.2] Simulação de Dados
# --------------------------

set.seed(123456789)
n <- 1000                               # Tamanho Amostral
x1 <- rnorm(n, 0, 1)                    # Normal(0, 1)
x2 <- rbinom(n, 1, 0.5)                 # Bernoulli(0.5)
X <- cbind(1, x1, x2)                   # Vetor de variáveis explicativas
vecCoef <- c(1.5, 2/3, 2)               # Vetor de coeficientes betas
timesSurv <- myRexp(n, 1, X, vecCoef) # Simulando o Tempo de Sobrevivência

# Visualização
ggplot(data = data.frame(Tempo = timesSurv), aes(x = Tempo)) +
  geom_histogram(fill = "blue") +
  labs(x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14)

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
  
  return(-flv)
}

# ---------------------------------
# [3.2] Algoritmo de Newton-Raphson
# ---------------------------------

theta0 <- c(0.5, 1, 0.5, 0.5) # Chute inicial
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
Stexp <- function(t, alpha, betas, X) alpha * exp(- alpha * ( t / exp(rowSums(X * betas)) ) )


# Formatando como DataFrame
DataExp <- data.frame(Time = timesSurv, 
                      Survival = Stexp(timesSurv, 1.5, vecCoef, X),
                      EMV_Survival = Stexp(timesSurv, estimate$par[1], estimate$par[-1], X))

# Gráfico com ggplot2
ggplot(DataExp, aes(x = Time)) +
  geom_line(aes(y = Survival, color = "Verdadeira")) +
  geom_line(aes(y = EMV_Survival, color = "Estimada")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  scale_color_manual(name = "Curva", values = c("Verdadeira" = "blue", "Estimada" = "red")) +
  theme_minimal(base_size = 14)
