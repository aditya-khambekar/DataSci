ak.aov = setRefClass(
    fields = list(aov = NULL),
    methods = list(
        residualsnormal = function(){
          qqnorm(aov$residuals)
          qqline(aov$residuals)
        },
        residualsd = function() {
          residual.df <- data.frame(
            Group = model.frame(aov)[[all.vars(aov$call$formula)[2]]],
            Residual = aov$residuals
          )
          result <- aggregate(Residual ~ Group, data = residual.df, FUN = sd)
          colnames(result) <- c("Group", "SD")

          return(result)
        },
        dataset.summary = function(){ # Note: Removed assignment arrow '<-'
          return(aov.dataset.summary(aov$call$formula, model.frame(aov)))
        }
        xyplot = function(point.col = "blue", line.col = "grey",
                             line.lty = 2, ...){
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
    )
)