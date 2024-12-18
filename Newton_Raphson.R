# ------------------------
# [2] Distribuição Weibull
# ------------------------

# Semente
set.seed(123)

# ---------------------------------------
# [2.1] Simulação e Visualização do Dados
# ---------------------------------------
# Tamanho da amostra
n <- 1000

# Parâmetros da distribuição
wShape <- 2
wScale <- 1.5

# Simulação
dadosWeibull <- rweibull(n, shape = wShape, scale = wScale)

# -------------
# [2.2] Funções
# -------------
# --------------------
# [2.2.1] Vetor Escore
# --------------------
score <- function(theta, dados) {
  gamma <- theta[1] # Parâmetro de Forma
  alpha <- theta[2] # Parâmetro de Escala
  
  n <- length(dados) # Tamanho da amostra
  t <- dados         # Dados observados
  
  # Derivadas Parciais
  # Em relação ao parâmetro de FORMA
  dGamma <- (n / gamma) - n * log(alpha) + sum(log(t)) - sum( (t/alpha)^gamma * log(t/alpha) )
  
  # Em relação ao parâmetro de ESCALA
  dAlpha <- - (n * gamma/alpha) + gamma*alpha^(- gamma - 1) * sum(t^gamma)
  
  # Retornar
  return(c(dGamma, dAlpha))
}

# -----------------------
# [2.2.2] Matriz Hessiana
# -----------------------
Hessian <- function(theta, dados) {
  gamma <- theta[1] # Parâmetro de Forma
  alpha <- theta[2] # Parâmetro de Escala
  
  n <- length(dados) # Tamanho da amostra
  t <- dados         # Dados observados
  
  dGamma2 <- - (n/gamma^2) - sum( (t/alpha)^gamma * log(t/alpha)^2 )
  dAlpha2 <- (n * gamma/alpha^2) - (gamma*(gamma + 1)) * (alpha^(- gamma - 2)) * sum(t^gamma)
  
  derivada <- - (n/alpha) - sum( (t^gamma / alpha^(gamma + 1)) * (gamma * log(t/alpha) + 1) )
  
  return(matrix(c(dGamma2, derivada, derivada, dAlpha2), nrow = 2, ncol = 2))
}

# --------------------------------------
# [3] Método Iterativo de Newton-Raphson
# --------------------------------------
# Definindo algumas grnadezas
theta0 <- c(1.5, 1)     # Chute Inicial
differ <- 1             # Diferença entre o passo atual e passo anterior
error <- 10^(-6)        # Erro tolerável
id <- 1                 # Contador (índice da iteração)

while(differ > error) {
  ESC <- score(theta = theta0, dados = dadosWeibull)   # Vetor Escore
  HES <- Hessian(theta = theta0, dados = dadosWeibull) # Matriz Hessiana
  
  soluc <- solve(HES, ESC)
  theta1 <- theta0 - soluc # Atulização do Algoritmo
  
  differ <- max(abs(theta1 - theta0))
  
  theta0 <- theta1
  
  cat("Iteração:", id, " ;  Estimativa = (Forma:", theta1[1], ", Escala:", theta1[2], ") \n")
  
  id <- id + 1
}