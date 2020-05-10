library(RJDemetra)
x13_model <- x13(Y[,2])
plot(x13_model, type_chart = "sa-trend")
plot(x13_model$decomposition)

ts_model <- tramoseats(Y[,2])