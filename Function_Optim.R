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

# Histograma
hist(x = dadosWeibull, xlab = "Dados Simulados", ylab = "Frequência",
     main = "Histograma dos Dados Simulados")

# --------------
# [3] Otimização
# --------------
# --------------------------------
# [3.1] Função Log-verossimilhança
# --------------------------------

logWeibull <- function(theta, dados){
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  n <- length(dados)
  t <- dados
  
  logverossimil <- (n * log(gamma)) - (gamma * log(alpha) * n) + (gamma - 1) * sum(log(t)) - sum((t/alpha)^gamma)
  return(-logverossimil)
}

# ------------------------------
# [3.2] Aplicando a função optim
# ------------------------------

theta0 <- c(1.5, 1) # Chute inicial
estimate <- optim(par = theta0, fn = logWeibull, gr = NULL , method = "BFGS" ,
                  hessian = TRUE, dados=dadosWeibull)
estimate

