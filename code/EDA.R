library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(xtable)

montos <- read.csv("../data/montos_creditos.csv")

montos2 <- montos %>%
  gather(Tipo, monto, Industrial:Servicios) %>%
  filter(complete.cases(.))

resumen <- data.frame(Grupo = names(montos), 
                      Media = sapply(montos, mean, na.rm = T),
                      Varianza = sapply(montos, var, na.rm = T)) %>%
  cbind(t(sapply(montos, quantile, na.rm = T)))
names(resumen) <- c("Grupo", 
                    "Media", 
                    "Varianza", 
                    "Mínimo", 
                    "Cuantil 25", 
                    "Mediana", 
                    "Cuantil 75", 
                    "Máximo")
row.names(resumen) <- NULL

xtable(resumen)

ggplot(montos2) +
  geom_boxplot(aes(x = Tipo, y = monto)) +
  xlab("Grupo de riesgo") +
  ylab("Monto") +
  theme_few()

ggsave(filename = "../output/histograma_datos.jpg", plot = last_plot())
