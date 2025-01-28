# --------------------------
# [1] Configurações Iniciais
# --------------------------

# Semente
set.seed(123456789)

# Tamanho Amostral
n <- 1000

# ------------------------
# [2] Distribuição Weibull
# ------------------------
# -------------------------
# [2.1] Simulação dos Dados
# -------------------------

w_shape <- 2   # Parâmetro de Forma
w_scale <- 1.5 # Parâmetro de Escala

# Simulação
w_dados <- rweibull(n, shape = w_shape, scale = w_scale)

# -------------
# [3] Estimação
# -------------
# ----------------------------------------
# [3.1] Método Iterativo de Newton-Raphson
# ----------------------------------------
# ----------------------
# [3.1.1] Função `optim`
# ----------------------

# Definindo Função Log-verossimilhança
logVerossimil <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Função Log-verossimilhança
  flv <- log(gamma)*n - gamma*log(alpha)*n + (gamma - 1)*sum(log(t)) - sum((t/alpha)^gamma)
  
  # Retorna o valor simétrico
  return(-flv)
}

# Vetor de Parâmetros: Chute inicial
theta0 <- c(1, 1)

# Obtendo as estimativas
estimate <- optim(
  par = theta0,
  fn = logVerossimil,
  method = "BFGS",
  hessian = TRUE,
  times = w_dados
)

# Resultados
print(estimate)

# -----------------
# [3.1.1] Algoritmo
# -----------------
# -------------------------
# [3.1.1.1] Vetor Gradiente
# -------------------------
U <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Derivadas Parciais
  DerivGamma <- n/gamma - log(alpha)*n + sum(log(t)) - sum((t/alpha)^gamma*log((t/alpha)^gamma))
  DerivAlpha <- -(gamma/alpha)*n + gamma*alpha^(-gamma-1)*sum(t^gamma)
  
  # Vetor Gradiente
  gradient <- c(DerivGamma, DerivAlpha)
    
  # Retornar
  return(gradient)
}

# -----------------------
# [2.2.2] Matriz Hessiana
# -----------------------
H <- function(times, theta) {
  # Número de observações
  n <- length(times)
  
  # Distição dos parâmetros
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  
  # Mudando o nome p/ facilitar a escrita
  t <- times
  
  # Derivadas de 2ª ordem
  D2Gamma <- n/gamma^2 - sum((t/alpha)^gamma*log((t/alpha)^gamma)^2)
  D2Alpha <- (gamma/alpha^2)*n - gamma*(gamma + 1)*alpha^(-gamma-2)*sum(t^gamma)
  D2 <- -n*alpha + sum((t^gamma/alpha^(gamma+1)))*(gamma*log(t/alpha)+1)
    
  # Matriz Hessiana
  H <- matrix(
    data = c(D2Gamma, D2, D2, D2Alpha),
    nrow = 2, ncol = 2
  )
  
  # Retornar
  return(H)
}

# --------------------------------------
# [3] Método Iterativo de Newton-Raphson
# --------------------------------------
theta0 <- c(1, 1) # Chute Inicial
diff <- 1           # Diferença entre o passo atual e passo anterior
error <- 10^(-8)    # Erro tolerável
id <- 1             # Contador da iteração

# Iteração
while(diff > error) {
  U <- U(times = w_dados, theta = theta0) # Vetor Escore
  H <- H(times = w_dados, theta = theta0) # Matriz Hessiana
  
  # Solução
  solution <- solve(H, U)
  
  # Atulização do Algoritmo
  theta1 <- theta0 - solution
  theta0 <- theta1
  
  # Diferença
  diff <- max(abs(theta1 - theta0))
  
  # Imprimir resultados na tela
  cat("Iteração:", id, " ;  Estimativa = (Forma:", theta1[1], ", Escala:", theta1[2], ") \n")
  
  # Controle do Algoritmo
  id <- id + 1
}