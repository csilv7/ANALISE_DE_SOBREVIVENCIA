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

n <- 1000

taxa <- 1

X <- matrix(
  data = c(rep(1, n), rnorm(n, 0, 1), rbinom(n, 1, 0.5)),
  nrow = n, ncol = 3
)

betas <- matrix(data = c(1.5, 2/3, 2), nrow = 3, ncol = 1)

timesSurv <- myRexp(n, taxa, X, betas) # Simulando o Tempo de Sobrevivência

# Visualização
ggplot(data = data.frame(Tempo = timesSurv), aes(x = Tempo)) +
  geom_histogram(fill = "blue") +
  labs(x = "Tempo", y = "Frequência") +
  theme_minimal(base_size = 14)



library(survival)

cens <- rep(1, n)

ekm <- survfit(Surv(timesSurv, cens)~1)

# Preparando os dados para o ggplot2
ekm_data <- data.frame(
  time = ekm$time,
  survival = ekm$surv,
  lower = ekm$lower,
  upper = ekm$upper
)

# Gráfico com ggplot2
ggplot(ekm_data, aes(x = time, y = survival)) +
  geom_step(color = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  labs(
    title = "Curva de Sobrevivência Kaplan-Meier",
    x = "Tempo",
    y = "Probabilidade de Sobrevivência",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )










