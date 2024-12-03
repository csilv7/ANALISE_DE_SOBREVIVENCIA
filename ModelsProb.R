library(ggplot2)

# Função de Sobrevivência (Exponencial)
St <- function(t, alpha) {
  exp(-alpha * t)
}

# Função de Risco (Constante - para simplificar)
ht <- function(alpha) {
  alpha
}

# Função de Risco Acumulado (Constante)
Lt <- function(t, alpha) {
  alpha * t
}

# Simulando dados de sobrevivência (exponencial)
set.seed(123)
n <- 1000
tempo <- rexp(n, rate = 1)  # 'rate' é o parâmetro da distribuição exponencial

# Calculando a função de sobrevivência para diferentes valores de alpha
fSt1 <- St(t = tempo, alpha = 1)
fSt2 <- St(t = tempo, alpha = 1.5)
fSt3 <- St(t = tempo, alpha = 2)

# Plotando as curvas de sobrevivência
ggplot() +
  geom_line(aes(x = tempo, y = fSt1, color = expression(paste(alpha, " = 1")))) +
  geom_line(aes(x = tempo, y = fSt2, color = "alpha = 1.5")) +
  geom_line(aes(x = tempo, y = fSt3, color = "alpha = 2")) +
  labs(x = "Tempo", y = "S(t)", color = "Valor de alpha") +
  theme_minimal()



library(ggplot2)

# Função de Sobrevivência (Exponencial)
St <- function(t, alpha) {
  exp(-alpha * t)
}

# Função de Risco (Constante - para simplificar)
ht <- function(alpha) {
  alpha
}

# Função de Risco Acumulado (Constante)
Lt <- function(t, alpha) {
  alpha * t
}

# Simulando dados de sobrevivência (exponencial)
set.seed(123)
n <- 1000
tempo <- rexp(n, rate = 1)

# Criando um data frame com valores de S(t) para diferentes alphas
dados <- data.frame(
  tempo = rep(tempo, 3),
  St = c(St(tempo, alpha = 1), St(tempo, alpha = 1.5), St(tempo, alpha = 2)),
  alpha = factor(rep(c(1, 1.5, 2), each = n))
)

# Plotando as curvas de sobrevivência com LaTeX nas legendas
ggplot(dados, aes(x = tempo, y = St, color = alpha)) +
  geom_line(stat = "summary", fun = mean, size = 1) +
  labs(x = "Tempo", y = "S(t)", 
    color = expression(alpha)) +
  scale_color_manual(
    values = c("red", "blue", "green"),
    labels = c(expression(alpha == 1), expression(alpha == 1.5), expression(alpha == 2))) +
  theme_minimal()

