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

# Caminho url para os dados
url <- "https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"

# Leitura dos dados
dados <- read.table(url, header = TRUE)

# -------------
# [3] ESTIMAÇÃO
# -------------
# ---------------------------
# [3.1] ESTIMAÇÃO NÃO-PARAMÉTRICA
# ---------------------------

# ---------------------------------
# [3.1.1] ESTIMADOR DE KAPLAN-MEIER
# ---------------------------------

ekm <- survfit(Surv(tempos, cens) ~ 1, data = dados)
ekm