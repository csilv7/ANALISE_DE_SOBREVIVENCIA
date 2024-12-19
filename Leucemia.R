# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------

if (!require("survival")){
  install.packages("survival")
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
}

library(survival)
library(ggplot2)

# ---------------------------------
# [2] IMPORTAÇÃO E AJUSTE DOS DADOS
# ---------------------------------

# Caminho URL para os dados
url <- "https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"

# Leitura dos dados
dados <- read.table(url, header = TRUE)

# -------------
# [3] ESTIMAÇÃO
# -------------

# -------------------------------
# [3.1] ESTIMADOR DE KAPLAN-MEIER
# -------------------------------
ekm <- survfit(Surv(tempos, cens) ~ 1, data = dados)

# -----------------
# [4] VISUALIZAÇÃO
# -----------------

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
    caption = "Fonte: Dados de leucemia"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )
