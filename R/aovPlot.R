aov.residualplot <- function(formula, data, point.col = "blue", line.col = "grey",
                             line.lty = 2, ...) {
  library("lattice")
  exp.var <- all.vars(formula)[2]

  data.aov <- aov(formula, data = data)

  xyplot(data.aov$residuals ~ data[[exp.var]], data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = 0, col = line.col, lty = line.lty)
         },
         xlab = exp.var,
         ylab = "Residuals",
  )
}

aov.xyplot <- function(formula, data, point.col = "blue", line.col = "grey",
                       line.lty = 2, ...) {
    library("lattice")


  response.var <- all.vars(formula)[1]

  grand.mean <- mean(data[[response.var]], na.rm = TRUE)

  xyplot(formula, data = data,
         panel = function(x, y, col, ...) {
           panel.xyplot(x, y, col = point.col, ...)
           panel.abline(h = grand.mean, col = line.col, lty = line.lty)
         },
         ...
  )
}