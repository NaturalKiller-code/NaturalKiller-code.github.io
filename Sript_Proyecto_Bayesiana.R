# PROYECTO Regresión Bayesiana
# Emiliano-Ferro

# Subimos el dataset
rm(list=ls())

# Cargar los datos
Data <- read.csv("TFM.csv", header=TRUE, as.is = FALSE, sep = ",")

# Ver la estructura de los datos
str(Data)

# Hacer análisis de regresión lineal frecuentista
model <- lm(Data$Tasa.de.flujo.ma.ximo..l.min. ~ Data$Altura..cm., data = Data)

# Ver el resumen del modelo
summary(model)

# Gráfico de los puntos de los datos y la línea de regresión ajustada
plot(Data$Altura..cm., Data$Tasa.de.flujo.ma.ximo..l.min., 
     xlab = "Altura (cm)", ylab = "Tasa de flujo máximo (l/min)", 
     main = "Regresión Lineal Simple", 
     pch = 16, col = "blue")  # Puntos azules


# Agregar la línea de regresión
abline(model, col = "red", lwd = 2)  # Línea roja de regresión


plot(model) 

# Mostrar los residuales
rsdls <- residuals(model)

# Test de normalidad de los residuos (Shapiro-Wilk)
shapiro.test(rsdls)

# Intervalos de confianza para los coeficientes
confint(model)
  
########
  # Quitamos intercepto
  # Hacer análisis de regresión lineal frecuentista sin intercepto
  model2 <- lm(Data$Tasa.de.flujo.ma.ximo..l.min. ~ Data$Altura..cm. - 1, data = Data)
  
  # Ver el resumen del modelo
  summary(model2)
  
  # Gráfico de los puntos de los datos y la línea de regresión ajustada
  plot(Data$Altura..cm., Data$Tasa.de.flujo.ma.ximo..l.min., 
       xlab = "Altura (cm)", ylab = "Tasa de flujo máximo (l/min)", 
       main = "Regresión Lineal Simple (Sin Intercepto)", 
       pch = 16, col = "blue")  # Puntos azules
  
  # Agregar la línea de regresión ajustada
  abline(model2, col = "red", lwd = 2)  # Línea roja de regresión
  
  # Gráficos diagnósticos del modelo
  plot(model2)
  
  # Mostrar los residuales
  rsdls2 <- residuals(model2)
  
  # Test de normalidad de los residuos (Shapiro-Wilk)
  shapiro.test(rsdls2)
  
  # Intervalos de confianza para los coeficientes
  confint(model2)

  
########
# Regresión Líneal simple Bayesiano

# Importar datos
Data <- read.csv("TFM.csv", header = TRUE, as.is = FALSE, sep = ",")

# Modelo clásico con intercepto
model_clasico <- lm(Data$Tasa.de.flujo.ma.ximo..l.min. ~ Data$Altura..cm., data = Data)

# Resumen del modelo clásico
summary(model_clasico)

# Ajuste del modelo bayesiano
# Parámetros de la distribución previa
mu0 <- c(0, 0)  # Media previa para (β0, β1)
Sigma0 <- diag(c(1000, 1000))  # Matriz de varianza-covarianza previa (2x2)

# Matriz de diseño X (con intercepto)
X <- cbind(1, Data$Altura..cm.)  # Primera columna de unos para el intercepto
Y <- Data$Tasa.de.flujo.ma.ximo..l.min.
n <- length(Y)

# Cálculo de los parámetros posteriores
sigma_squared <- var(Y)  # Estimación de la varianza
Sigma1 <- solve(t(X) %*% X / sigma_squared + solve(Sigma0))  # Varianza posterior
mu1 <- Sigma1 %*% (t(X) %*% Y / sigma_squared + solve(Sigma0) %*% mu0)  # Media posterior

# Media a posteriori (estimadores MAP)
beta_MAP <- mu1

# Crear la gráfica comparativa
plot(Data$Altura..cm., Data$Tasa.de.flujo.ma.ximo..l.min., 
     xlab = "Altura (cm)", ylab = "Tasa de flujo máximo (l/min)", 
     main = "Comparación entre Modelos Clásico y Bayesiano",
     pch = 16, col = "blue")  # Puntos de datos en azul

# Agregar la línea ajustada del modelo clásico
abline(a = coef(model_clasico)[1], b = coef(model_clasico)[2], col = "red", lwd = 2, lty = 1)  # Línea roja sólida

# Agregar la línea ajustada del modelo bayesiano
abline(a = beta_MAP[1], b = beta_MAP[2], col = "green", lwd = 2, lty = 2)  # Línea verde discontinua

# Leyenda para identificar las líneas
legend("topright", legend = c("Modelo Clásico (Con intercepto)", "Modelo Bayesiano (MAP)"),
       col = c("red", "green"), lty = c(1, 2), lwd = 2)


########

# Distribuciones a posteriori marginales de los parámetros β0 y β1
# Generar muestras de la distribución posterior para β0 y β1
library(MASS)
set.seed(123)

# Generar 10,000 muestras de la distribución posterior conjunta de β0 y β1
samples <- mvrnorm(10000, mu = as.numeric(mu1), Sigma = as.matrix(Sigma1))

# Histograma de la distribución posterior marginal de β0
hist(samples[, 1], col = "blue", breaks = 50, probability = TRUE,
     xlab = expression(beta[0]), main = "Distribución posterior marginal de β0")

# Agregar una línea para la densidad de β0
lines(density(samples[, 1]), col = "red", lwd = 2)

# Histograma de la distribución posterior marginal de β1
hist(samples[, 2], col = "green", breaks = 50, probability = TRUE,
     xlab = expression(beta[1]), main = "Distribución posterior marginal de β1")

# Agregar una línea para la densidad de β1
lines(density(samples[, 2]), col = "red", lwd = 2)


# Distribución posterior conjunta de β0 y β1
# Cargar librería
library(MASS)

# Parámetros de la distribución posterior (asegúrate de que mu1 y Sigma1 estén calculados correctamente)
# mu1 debe ser un vector de longitud 2, y Sigma1 debe ser una matriz 2x2
set.seed(123)

# Generar muestras de la distribución posterior conjunta
samples <- mvrnorm(10000, mu = mu1, Sigma = Sigma1)

# Gráfico de dispersión de la distribución posterior conjunta
plot(samples[, 1], samples[, 2], col = rgb(0, 0, 1, 0.1), pch = 20, 
     xlab = expression(beta[0]), ylab = expression(beta[1]), 
     main = "Distribución posterior conjunta de β0 y β1")



# Cálculo de la correlación entre β0 y β1
correlation <- cor(samples[, 1], samples[, 2])
cat("Correlación entre β0 y β1:", correlation, "\n")

###########

# Algoritmo de Monte Carlo
# Cargar librerías necesarias
library(ggplot2)

# Definir los datos
altura <- c(174, 183, 176, 169, 183, 186, 178, 175, 172, 179, 171, 184, 200, 195, 176, 176, 190)
tasa_flujo <- c(733, 572, 500, 738, 616, 787, 866, 670, 550, 660, 575, 577, 783, 625, 470, 642, 856)

# Ajustar el modelo de regresión lineal clásico
modelo <- lm(tasa_flujo ~ altura)
predicciones <- predict(modelo, interval = "confidence")
predicciones_pred <- predict(modelo, interval = "prediction")

# Algoritmo de Monte Carlo (como el proporcionado anteriormente)
n_iterations <- 200000  # Total de iteraciones
burn_in <- 2000        # Burn-in
lag <- 5              # Thinning (lag)
samples <- matrix(NA, nrow = floor((n_iterations - burn_in) / lag), ncol = 2)  # Guardar muestras

# Valores iniciales
beta0_current <- coef(modelo)[1]  # β0 del modelo clásico
beta1_current <- coef(modelo)[2]  # β1 del modelo clásico
sigma <- summary(modelo)$sigma    # Desviación estándar del modelo clásico

# Algoritmo Metropolis
set.seed(42)
for (i in 1:n_iterations) {
  # Propuestas para β0 y β1
  beta0_star <- beta0_current + rnorm(1, mean = 0, sd = 0.1)
  beta1_star <- beta1_current + rnorm(1, mean = 0, sd = 0.1)
  
  # Likelihoods
  log_likelihood_current <- sum(dnorm(tasa_flujo - (beta0_current + beta1_current * altura), mean = 0, sd = sigma, log = TRUE))
  log_likelihood_star <- sum(dnorm(tasa_flujo - (beta0_star + beta1_star * altura), mean = 0, sd = sigma, log = TRUE))
  
  # Proporción de aceptación
  acceptance_ratio <- exp(log_likelihood_star - log_likelihood_current)
  if (runif(1) < acceptance_ratio) {
    beta0_current <- beta0_star
    beta1_current <- beta1_star
  }
  
  # Guardar muestras después de burn-in y lag
  if (i > burn_in && (i - burn_in) %% lag == 0) {
    idx <- (i - burn_in) / lag
    samples[idx, ] <- c(beta0_current, beta1_current)
  }
}

# Separar muestras de β0 y β1
beta0_samples <- samples[, 1]
beta1_samples <- samples[, 2]


# --- 1. Aproximación a la función de densidad para beta0 y beta1 ---
beta0_density <- density(beta0_samples)
beta1_density <- density(beta1_samples)

# --- 2. Crear histogramas con la densidad superpuesta ---
# Histograma con la densidad ajustada para beta0
ggplot(data.frame(beta0 = beta0_samples), aes(x = beta0)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7, color = "black") +
  geom_line(data = data.frame(x = beta0_density$x, y = beta0_density$y), aes(x = x, y = y), color = "red", lwd = 1) +
  ggtitle("Histograma y Densidad Posterior para β0") +
  xlab(expression(beta[0])) +
  ylab("Densidad") +
  theme_minimal()

# Histograma con la densidad ajustada para beta1
ggplot(data.frame(beta1 = beta1_samples), aes(x = beta1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "green", alpha = 0.7, color = "black") +
  geom_line(data = data.frame(x = beta1_density$x, y = beta1_density$y), aes(x = x, y = y), color = "red", lwd = 1) +
  ggtitle("Histograma y Densidad Posterior para β1") +
  xlab(expression(beta[1])) +
  ylab("Densidad") +
  theme_minimal()

# --- 3. Mostrar las estadísticas descriptivas ---
summary(beta0_samples)
summary(beta1_samples)

# --- 4. Superponer la curva de densidad sobre el histograma ---
par(mfrow = c(1, 2))  # Dividir la pantalla en 2 gráficos


# --- 5. Gráficas de Autocorrelación para β0 y β1 ---
par(mfrow = c(1, 2))  # Dividir la pantalla en 2 gráficos

# Autocorrelación para β0
acf(beta0_samples, main = "Autocorrelación de β0")

# Autocorrelación para β1
acf(beta1_samples, main = "Autocorrelación de β1")

# Histograma con la densidad ajustada para beta0
hist(beta0_samples, probability = TRUE, col = "blue", main = "Posterior de β0", xlab = expression(beta[0]))
lines(beta0_density, col = "red", lwd = 2)

# Histograma con la densidad ajustada para beta1
hist(beta1_samples, probability = TRUE, col = "green", main = "Posterior de β1", xlab = expression(beta[1]))
lines(beta1_density, col = "red", lwd = 2)

# --- 5. Graficar la cadena de MCMC para β0 y β1 ---
# Verificar la convergencia de las cadenas
plot(beta0_samples, type = "l", col = "blue", main = "Cadena de β0", ylab = "β0")
plot(beta1_samples, type = "l", col = "green", main = "Cadena de β1", ylab = "β1")

# Gráficas adicionales para el modelo clásico y Monte Carlo

# Predicciones usando Monte Carlo
predicciones_mc <- sapply(1:length(altura), function(i) {
  sapply(1:nrow(samples), function(j) samples[j, 1] + samples[j, 2] * altura[i])
})

# Bandas de confianza y predicción para Monte Carlo
lower_bound <- apply(predicciones_mc, 2, quantile, probs = 0.025)
upper_bound <- apply(predicciones_mc, 2, quantile, probs = 0.975)
predicciones_mc_mean <- apply(predicciones_mc, 2, mean)

# Predicciones de intervalos de predicción Monte Carlo
lower_pred <- apply(predicciones_mc + qnorm(0.025, 0, sigma), 2, quantile, probs = 0.025)
upper_pred <- apply(predicciones_mc + qnorm(0.975, 0, sigma), 2, quantile, probs = 0.975)

# Crear data frames para ggplot
data_clasica <- data.frame(altura, tasa_flujo, predicciones = predicciones[, "fit"], lower_conf = predicciones[, "lwr"], upper_conf = predicciones[, "upr"], lower_pred = predicciones_pred[, "lwr"], upper_pred = predicciones_pred[, "upr"])
data_mc <- data.frame(altura, tasa_flujo, predicciones_mc_mean, lower_bound, upper_bound, lower_pred, upper_pred)

# Gráfico de Regresión Lineal Clásica con Bandas de Predicción y Confianza
ggplot(data_clasica, aes(x = altura, y = tasa_flujo)) +
  geom_point(color = 'blue', size = 2) +
  geom_line(aes(y = predicciones), color = 'red', size = 1) +
  geom_ribbon(aes(ymin = lower_conf, ymax = upper_conf), fill = 'darkgray', alpha = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_pred, ymax = upper_pred), fill = 'gray', alpha = 0.3) +
  labs(title = "Regresión Lineal Clásica con Bandas de Confianza y Predicción",
       x = "Altura (cm)", y = "Tasa de Flujo Máximo (l/min)") +
  theme_minimal()

# Gráfico de Monte Carlo con Bandas de Predicción y Confianza
ggplot(data_mc, aes(x = altura, y = tasa_flujo)) +
  geom_point(color = 'blue', size = 2) +
  geom_line(aes(y = predicciones_mc_mean), color = 'green', size = 1) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = 'darkgray', alpha = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin= lower_pred, ymax = upper_pred), fill = 'gray', alpha = 0.3) + labs(title = "Monte Carlo: Bandas de Confianza y Predicción", x = "Altura (cm)", y = "Tasa de Flujo Máximo (l/min)") + theme_minimal()



# PREGUNTA 1
# Altura de interés
altura_interes <- 160

# Calcular las tasas de flujo máxima utilizando las muestras de MCMC
tasa_flujo_muestras <- beta0_samples + beta1_samples * altura_interes

# Calcular el intervalo de credibilidad al 95% para la media
credibilidad_95 <- quantile(tasa_flujo_muestras, probs = c(0.025, 0.975))
media_tasa_flujo <- mean(tasa_flujo_muestras)

# Imprimir resultados
cat("Tasa de flujo máxima promedio para una altura de 1.6 metros:", media_tasa_flujo, "\n")
cat("Intervalo de credibilidad al 95% para la media: [", credibilidad_95[1], ", ", credibilidad_95[2], "]\n")


# PREGUNTA 2

# Altura de interés (210 cm)
altura_interes2 <- 210

# Calcular las tasas de flujo máxima utilizando las muestras de MCMC para la nueva altura
tasa_flujo_muestras2 <- beta0_samples + beta1_samples * altura_interes2

# Calcular el intervalo de credibilidad al 95% para la media
credibilidad_95_2 <- quantile(tasa_flujo_muestras2, probs = c(0.025, 0.975))
media_tasa_flujo_2 <- mean(tasa_flujo_muestras2)

# Imprimir resultados
cat("Tasa de flujo máxima promedio para una altura de 2.1 metros:", media_tasa_flujo_2, "\n")
cat("Intervalo de credibilidad al 95% para la media: [", credibilidad_95_2[1], ", ", credibilidad_95_2[2], "]\n")



















# Algoritmo de Monte Carlo
# Cargar librerías necesarias
library(ggplot2)
library(coda)  # Para las gráficas de autocorrelación

# Definir los datos
altura <- c(174, 183, 176, 169, 183, 186, 178, 175, 172, 179, 171, 184, 200, 195, 176, 176, 190)
tasa_flujo <- c(733, 572, 500, 738, 616, 787, 866, 670, 550, 660, 575, 577, 783, 625, 470, 642, 856)

# Ajustar el modelo de regresión lineal clásico
modelo <- lm(tasa_flujo ~ altura)
predicciones <- predict(modelo, interval = "confidence")
predicciones_pred <- predict(modelo, interval = "prediction")

# Algoritmo de Monte Carlo
n_iterations <- 200000  # Total de iteraciones
burn_in <- 2000         # Burn-in
lag <- 5                # Thinning (lag)
samples <- matrix(NA, nrow = floor((n_iterations - burn_in) / lag), ncol = 2)  # Guardar muestras

# Valores iniciales
beta0_current <- coef(modelo)[1]  # β0 del modelo clásico
beta1_current <- coef(modelo)[2]  # β1 del modelo clásico
sigma <- summary(modelo)$sigma    # Desviación estándar del modelo clásico

# Algoritmo Metropolis
set.seed(42)
for (i in 1:n_iterations) {
  # Propuestas para β0 y β1
  beta0_star <- beta0_current + rnorm(1, mean = 0, sd = 0.1)
  beta1_star <- beta1_current + rnorm(1, mean = 0, sd = 0.1)
  
  # Likelihoods
  log_likelihood_current <- sum(dnorm(tasa_flujo - (beta0_current + beta1_current * altura), mean = 0, sd = sigma, log = TRUE))
  log_likelihood_star <- sum(dnorm(tasa_flujo - (beta0_star + beta1_star * altura), mean = 0, sd = sigma, log = TRUE))
  
  # Proporción de aceptación
  acceptance_ratio <- exp(log_likelihood_star - log_likelihood_current)
  if (runif(1) < acceptance_ratio) {
    beta0_current <- beta0_star
    beta1_current <- beta1_star
  }
  
  # Guardar muestras después de burn-in y lag
  if (i > burn_in && (i - burn_in) %% lag == 0) {
    idx <- (i - burn_in) / lag
    samples[idx, ] <- c(beta0_current, beta1_current)
  }
}

# Separar muestras de β0 y β1
beta0_samples <- samples[, 1]
beta1_samples <- samples[, 2]

# --- 1. Aproximación a la función de densidad para beta0 y beta1 ---
beta0_density <- density(beta0_samples)
beta1_density <- density(beta1_samples)

# --- 2. Crear histogramas con la densidad superpuesta ---
# Histograma con la densidad ajustada para beta0
ggplot(data.frame(beta0 = beta0_samples), aes(x = beta0)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7, color = "black") +
  geom_line(data = data.frame(x = beta0_density$x, y = beta0_density$y), aes(x = x, y = y), color = "red", lwd = 1) +
  ggtitle("Histograma y Densidad Posterior para β0") +
  xlab(expression(beta[0])) +
  ylab("Densidad") +
  theme_minimal()

# Histograma con la densidad ajustada para beta1
ggplot(data.frame(beta1 = beta1_samples), aes(x = beta1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "green", alpha = 0.7, color = "black") +
  geom_line(data = data.frame(x = beta1_density$x, y = beta1_density$y), aes(x = x, y = y), color = "red", lwd = 1) +
  ggtitle("Histograma y Densidad Posterior para β1") +
  xlab(expression(beta[1])) +
  ylab("Densidad") +
  theme_minimal()

# --- 3. Mostrar las estadísticas descriptivas ---
cat("Estadísticas descriptivas para β0:\n")
print(summary(beta0_samples))
cat("\nEstadísticas descriptivas para β1:\n")
print(summary(beta1_samples))

# --- 4. Graficar la cadena de MCMC para β0 y β1 ---
# Verificar la convergencia de las cadenas
par(mfrow = c(1, 2))  # Dividir la pantalla en 2 gráficos
plot(beta0_samples, type = "l", col = "blue", main = "Cadena de β0", ylab = "β0")
plot(beta1_samples, type = "l", col = "green", main = "Cadena de β1", ylab = "β1")

# --- 5. Gráficas de Autocorrelación para β0 y β1 ---
par(mfrow = c(1, 2))  # Dividir la pantalla en 2 gráficos

# Autocorrelación para β0
acf(beta0_samples, main = "Autocorrelación de β0")

# Autocorrelación para β1
acf(beta1_samples, main = "Autocorrelación de β1")

# Gráficas adicionales para el modelo clásico y Monte Carlo

# Predicciones usando Monte Carlo
predicciones_mc <- sapply(1:length(altura), function(i) {
  sapply(1:nrow(samples), function(j) samples[j, 1] + samples[j, 2] * altura[i])
})

# Bandas de confianza y predicción para Monte Carlo
lower_bound <- apply(predicciones_mc, 2, quantile, probs = 0.025)
upper_bound <- apply(predicciones_mc, 2, quantile, probs = 0.975)
predicciones_mc_mean <- apply(predicciones_mc, 2, mean)

# Predicciones de intervalos de predicción Monte Carlo
lower_pred <- apply(predicciones_mc + qnorm(0.025, 0, sigma), 2, quantile, probs = 0.025)
upper_pred <- apply(predicciones_mc + qnorm(0.975, 0, sigma), 2, quantile, probs = 0.975)

# Crear data frames para ggplot
data_clasica <- data.frame(altura, tasa_flujo, predicciones = predicciones[, "fit"], lower_conf = predicciones[, "lwr"], upper_conf = predicciones[, "upr"], lower_pred = predicciones_pred[, "lwr"], upper_pred = predicciones_pred[, "upr"])
data_mc <- data.frame(altura, tasa_flujo, predicciones_mc_mean, lower_bound, upper_bound, lower_pred, upper_pred)

# Gráfico de Regresión Lineal Clásica con Bandas de Predicción y Confianza
ggplot(data_clasica, aes(x = altura, y = tasa_flujo)) +
  geom_point(color = 'blue', size = 2) +
  geom_line(aes(y = predicciones), color = 'red', size = 1) +
  geom_ribbon(aes(ymin = lower_conf, ymax = upper_conf), fill = 'darkgray', alpha = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_pred, ymax = upper_pred), fill = 'gray', alpha = 0.3) +
  labs(title = "Regresión Lineal Clásica con Bandas de Confianza y Predicción",
       x = "Altura (cm)", y = "Tasa de Flujo Máximo (l/min)") +
  theme_minimal()

# Gráfico de Monte Carlo con Bandas de Pred