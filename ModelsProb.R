# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
library(ggplot2)

# -------------
# [2] Simulação
# -------------

# Lambda verdadeiro
lambda <- 2 / 3

# Definindo semente para reprodutibilidade
set.seed(123)

# Definindo o tamanho da amostra
n <- 1000

# Simulando
survival_times <- rexp(n, rate = lambda)

# -------------------------
# [3] Visualização do dados
# -------------------------
df <- data.frame(TIME = survival_times)

ggplot(data = df, aes(x = TIME)) +
  geom_histogram(bins = 20, fill = "gray", color = "black") +
  labs(x = "Dados Simulados", y = "Frequência") +
  theme_minimal()

