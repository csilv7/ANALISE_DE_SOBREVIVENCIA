# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Análise de Sobrevivência - Modelos de Tempo de Vida Acelerado
# Analista de Sobrevicência: Breno C R da Silva
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------

# -----------------------
# [1] Ativação de Pacotes
# -----------------------

library(ggplot2)
library(gridExtra)

library(survival)

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
  effect <- coVariates %*% vecCoef
  
  # Tempos de sobrevivência
  times <- - exp(effect) * log(1 - U) / rate
  
  return(times)
}

# --------------------------
# [2.2.2] Simulação de Dados
# --------------------------

set.seed(123456789)

# Tamanho da amostra e número de variáveis
n <- 1000
p <- 2

# Parâmetro da distribuição exponencial
taxa <- 1

# Matriz Design
x1 <- rnorm(n, 0, 1)
x2 <- rbinom(n, 1, 0.5)
X <- matrix(
  data = c(rep(1, n), x1, x2),
  nrow = n, ncol = p + 1
)

# Vetor de Coeficientes Betas
betas <- matrix(data = c(1.5, 2/3, 2), nrow = p + 1, ncol = 1)

# Simulando o Tempo de Sobrevivência
timesSurv <- myRexp(n, taxa, X, betas)

# ------------------------------
# [2.2.1] Visualizações Gráficas
# ------------------------------

# Histograma
g1 <- ggplot(data = data.frame(Tempo = timesSurv), aes(x = Tempo)) +
  geom_histogram(fill = "blue") +
  labs(title = "Histograma do Tempo de Sobrevida",
       x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
  )

# Função de Sobrevivência
cens <- rep(1, n)
ekm <- survfit(Surv(timesSurv, cens)~1)
ekm_data <- data.frame(
  time = ekm$time,
  survival = ekm$surv,
  lower = ekm$lower,
  upper = ekm$upper
)
g2 <- ggplot(ekm_data, aes(x = time, y = survival)) +
  geom_step(color = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Função de Sobrevida de Kaplan-Meier",
    x = "Tempo",
    y = "Probabilidade de Sobrevivência",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
  )

# Combine os gráficos em um subplot com duas colunas e uma linha
grid.arrange(g1, g2, ncol = 2, widths=c(1, 1.25))

# ---------------------------
# [3] Estimação de Parâmetros
# ---------------------------

# -----------------------------------
# [3.1] Ajuste pelo Pacote "Survival"
# -----------------------------------

fit1 <- survreg(formula = Surv(timesSurv, cens)~x1 + x2, 
                dist = "exponential")
summary(fit1)

# --------------------
# [3.1] Via otimização
# --------------------

# ---------------------------
# [3.1.1] Log-Verossimilhança
# ---------------------------

logVerossimil <- function(theta, coVariates, times) {
  # Comprimento do vetor de parâmetros
  nPar <- length(theta)
  
  # Parâmetro de Taxa e Preditores Lineares
  rate <- theta[1]
  betas <- theta[2:nPar]
  

  # Combinação linear dos preditores lineares
  effect <- coVariates %*% betas
  
  # Função de Log-verossimilhança
  flv <- n * log(rate) - rate * sum(times / exp(effect))
  
  return(-flv)
}

# --------------------
# [3.1.1] Função optim
# --------------------

theta0 <- rep(1.5, 4) # Chute inicial

# Aplicação do algoritmo
estimate <- optim(
  par = theta0,
  fn = logVerossimil,
  method = "BFGS",
  hessian = TRUE,
  coVariates = X,
  times = timesSurv
)

# Resultados
print(estimate)

# --------------
# [4] Comparação
# --------------

thetaEst <- estimate$par
nPar <- length(thetaEst)

# Função de Sobrevivência
Stexp <- function(t, alpha, betas, X) alpha * exp(- alpha * ( t / exp(X %*% betas) ) )


# Formatando como DataFrame
DataExp <- data.frame(Time = timesSurv, 
                      Survival = Stexp(timesSurv, taxa, betas, X),
                      EMV_Survival = Stexp(timesSurv, thetaEst[1], thetaEst[2:nPar], X))

# Gráfico com ggplot2
ggplot(DataExp, aes(x = Time)) +
  geom_line(aes(y = Survival, color = "Verdadeira")) +
  geom_line(aes(y = EMV_Survival, color = "Estimada")) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  scale_color_manual(name = "Curva", values = c("Verdadeira" = "blue", "Estimada" = "red")) +
  theme_minimal(base_size = 14)
