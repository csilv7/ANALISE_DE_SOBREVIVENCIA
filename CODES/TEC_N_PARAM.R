# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
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

# ------------------------------------------------------------------------------

# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
library(survival)
library(ggplot2)
library(dplyr)

# -----------------------------
# [2] IMPORTAÇÃO E TRATAMENTO DOS DADOS
# -----------------------------

# Caminho URL para os dados
url <- "https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"

# Leitura dos dados
df <- read.table(url, header = TRUE)

# Decodificando a coluna r_6
df <- df %>%
  mutate(grupo = ifelse(r6 == 0, "Category Zero", "Category One"))

# ----------------------
# [3] TESTE DE LOGRANK
# ----------------------

# Criando o objeto de sobrevivência
dados_surv <- Surv(df$tempos, df$cens)

# Aplicando o Teste de Logrank
teste_logrank <- survdiff(dados_surv ~ grupo, data = df)
print(teste_logrank)

# -----------------------------
# [4] ESTIMADOR DE KAPLAN-MEIER
# -----------------------------

# Ajuste do modelo Kaplan-Meier por grupo
ekm <- survfit(dados_surv ~ grupo, data = df)

# -----------------------
# [5] GRÁFICO DAS CURVAS
# -----------------------

# Preparando os dados para o ggplot2
ekm_data <- data.frame(
  time = ekm$time,
  survival = ekm$surv,
  lower = ekm$lower,
  upper = ekm$upper,
  grupo = rep(levels(factor(df$grupo)), ekm$strata)
)

# Gráfico com ggplot2
ggplot(ekm_data, aes(x = time, y = survival, color = grupo, fill = grupo)) +
  geom_line(lwd = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x = "Tempo",
    y = "Probabilidade de Sobrevivência",
    color = "Grupo",
    fill = "Grupo",
    title = "Curvas de Sobrevivência de Kaplan-Meier por Grupo",
    caption = "Fonte: https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  )

