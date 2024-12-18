# Função de Sobrevivência
StWeibull <- function(t, gamma, alpha) exp(-(t/alpha)^gamma)

n <- 1000

tempo <- c(rep(0, n))
deltc <- c(rep(0, n))

set.seed(12345)

for (id in (1:n)) {
  if ((table(deltc)[1] / n) <= 1/3) {
    tempo[id] <- rexp(1, rate = 1)
    deltc[id] <- 0
  } else {
    tempo[id] <- rweibull(1, shape = 2, scale = 1.5)
    deltc[id] <- 1
  }
}

hist(tempo)

# --------------------------
# Função Log-verossimilhança
# --------------------------

logWeibull <- function(theta, dados){
  gamma <- theta[1] # Parâmetro de forma
  alpha <- theta[2] # Parâmetro de escala
  n <- length(dados)
  t <- dados[1]
  c <- dados[2]
  
  logverossimil <- (sum(c) * log(gamma)) - (gamma * log(alpha) * sum(c)) + (gamma - 1) * sum(c * log(t)) - sum((t/alpha)^gamma)
  return(-logverossimil)
}

# ------------------------
# Aplicando a função optim
# ------------------------

theta0 <- c(1.5, 1) # Chute inicial
estimate <- optim(par = theta0, fn = logWeibull, gr = NULL , method = "BFGS" ,
                  hessian = TRUE, dados=c(tempo, deltc))
estimate
